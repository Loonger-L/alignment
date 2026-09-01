#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <zlib.h>
#include "mmpriv.h"

struct bam_ref_t { std::string name; int32_t len; };

static void put32(std::vector<uint8_t> &b, int32_t x)
{
	b.push_back(x); b.push_back(x >> 8); b.push_back(x >> 16); b.push_back(x >> 24);
}
static void put16(std::vector<uint8_t> &b, uint16_t x)
{
	b.push_back(x); b.push_back(x >> 8);
}
static int ref_id(const std::vector<bam_ref_t> &refs, const char *name)
{
	if (strcmp(name, "*") == 0) return -1;
	for (size_t i = 0; i < refs.size(); ++i)
		if (refs[i].name == name) return (int)i;
	return -1;
}
static int cigar_op(char c)
{
	const char *ops = "MIDNSHP=XB";
	const char *p = strchr(ops, c);
	return p? (int)(p - ops) : -1;
}
static int cigar_end(const char *cigar, int pos)
{
	int n = 0, end = pos;
	if (strcmp(cigar, "*") == 0) return end;
	for (const char *p = cigar; *p; ++p) {
		if (*p >= '0' && *p <= '9') n = n * 10 + *p - '0';
		else {
			if (strchr("MDN=X", *p)) end += n;
			n = 0;
		}
	}
	return end;
}
static uint16_t reg2bin(int beg, int end)
{
	--end;
	if (beg >> 14 == end >> 14) return 4681 + (beg >> 14);
	if (beg >> 17 == end >> 17) return 585 + (beg >> 17);
	if (beg >> 20 == end >> 20) return 73 + (beg >> 20);
	if (beg >> 23 == end >> 23) return 9 + (beg >> 23);
	if (beg >> 26 == end >> 26) return 1 + (beg >> 26);
	return 0;
}
static int bgzf_write(FILE *fp, const uint8_t *data, size_t len)
{
	while (len) {
		size_t n = len > 60000? 60000 : len;
		uint8_t compressed[65536];
		z_stream zs = {};
		if (deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) return -1;
		zs.next_in = (Bytef *)data; zs.avail_in = n;
		zs.next_out = compressed; zs.avail_out = sizeof(compressed);
		if (deflate(&zs, Z_FINISH) != Z_STREAM_END) { deflateEnd(&zs); return -1; }
		size_t clen = zs.total_out;
		deflateEnd(&zs);
		uint8_t hdr[] = {31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0,0,0};
		uint16_t bsize = (uint16_t)(clen + 25);
		hdr[16] = bsize; hdr[17] = bsize >> 8;
		if (fwrite(hdr, 1, sizeof(hdr), fp) != sizeof(hdr) ||
			fwrite(compressed, 1, clen, fp) != clen) return -1;
		uint32_t crc = crc32(0L, data, n);
		uint8_t tail[] = {(uint8_t)crc,(uint8_t)(crc>>8),(uint8_t)(crc>>16),(uint8_t)(crc>>24),
			(uint8_t)n,(uint8_t)(n>>8),(uint8_t)(n>>16),(uint8_t)(n>>24)};
		if (fwrite(tail, 1, sizeof(tail), fp) != sizeof(tail)) return -1;
		data += n; len -= n;
	}
	return 0;
}
static int bgzf_finish(FILE *fp)
{
	static const uint8_t eof[] = {31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0,27,0,3,0,0,0,0,0,0,0,0,0};
	return fwrite(eof, 1, sizeof(eof), fp) == sizeof(eof)? 0 : -1;
}
static int append_tag(std::vector<uint8_t> &b, const char *tag)
{
	if (strlen(tag) < 5 || tag[2] != ':' || tag[4] != ':') return -1;
	b.push_back(tag[0]); b.push_back(tag[1]); b.push_back(tag[3]);
	const char *v = tag + 5;
	switch (tag[3]) {
	case 'A': b.push_back(*v); break;
	case 'c': case 'C': b.push_back((uint8_t)strtol(v, 0, 10)); break;
	case 's': case 'S': put16(b, (uint16_t)strtol(v, 0, 10)); break;
	case 'i': case 'I': put32(b, (int32_t)strtol(v, 0, 10)); break;
	case 'f': { float f = (float)strtod(v, 0); uint32_t x; memcpy(&x, &f, 4); put32(b, x); break; }
	case 'Z': case 'H': b.insert(b.end(), v, v + strlen(v) + 1); break;
	default: return -1;
	}
	return 0;
}
static int append_record(std::vector<uint8_t> &out, char *line, const std::vector<bam_ref_t> &refs)
{
	char *f[4096], *p = line;
	int nf = 0;
	while (nf < 4096) {
		f[nf] = strsep(&p, "\t");
		if (f[nf] == 0) break;
		++nf;
	}
	if (nf < 11) return -1;
	int32_t pos = atoi(f[3]) - 1, l_seq = strcmp(f[9], "*")? (int32_t)strlen(f[9]) : 0;
	int32_t rid = ref_id(refs, f[2]), nrid = strcmp(f[6], "=") == 0? rid : ref_id(refs, f[6]);
	std::vector<uint8_t> b;
	put32(b, rid); put32(b, pos);
	uint32_t bin_mq_nl = ((uint32_t)(rid < 0? 0 : reg2bin(pos, cigar_end(f[5], pos))) << 16) |
		((uint32_t)atoi(f[4]) << 8) | (uint32_t)(strlen(f[0]) + 1);
	put32(b, bin_mq_nl);
	int nc = 0;
	if (strcmp(f[5], "*") != 0)
		for (const char *q = f[5]; *q; ++q) if (!(*q >= '0' && *q <= '9')) ++nc;
	put32(b, ((uint32_t)atoi(f[1]) << 16) | nc);
	put32(b, l_seq); put32(b, nrid); put32(b, atoi(f[7]) - 1); put32(b, atoi(f[8]));
	b.insert(b.end(), f[0], f[0] + strlen(f[0]) + 1);
	int n = 0;
	for (const char *q = f[5]; *q && *q != '*'; ++q)
		if (*q >= '0' && *q <= '9') n = n * 10 + *q - '0';
		else { int op = cigar_op(*q); if (op < 0) return -1; put32(b, n << 4 | op); n = 0; }
	for (int i = 0; i < l_seq; i += 2) {
		const char *base = " =ACMGRSVTWYHKDBN"; const char *a = strchr(base, f[9][i] & ~32);
		const char *d = i + 1 < l_seq? strchr(base, f[9][i+1] & ~32) : base;
		b.push_back(((a? a-base : 15) << 4) | (d? d-base : 0));
	}
	for (int i = 0; i < l_seq; ++i) b.push_back(strcmp(f[10], "*")? f[10][i] - 33 : 255);
	for (int i = 11; i < nf; ++i) if (append_tag(b, f[i]) < 0) return -1;
	put32(out, b.size()); out.insert(out.end(), b.begin(), b.end());
	return 0;
}
int mm_sam_to_bam(FILE *sam, const char *fn)
{
	char *line = 0; size_t cap = 0; ssize_t n;
	std::string hdr; std::vector<bam_ref_t> refs;
	while ((n = getline(&line, &cap, sam)) >= 0) {
		if (line[n-1] == '\n') line[--n] = 0;
		if (line[0] == '@') {
			hdr += line; hdr += '\n';
			if (strncmp(line, "@SQ\t", 4) == 0) {
				char *sn = strstr(line, "\tSN:"), *ln = strstr(line, "\tLN:");
				if (!sn || !ln) { free(line); return -1; }
				sn += 4; char *end = strchr(sn, '\t');
				refs.push_back({std::string(sn, end? end : sn + strlen(sn)), (int32_t)strtol(ln + 4, 0, 10)});
			}
		} else break;
	}
	free(line);
	if (fseek(sam, 0, SEEK_SET) != 0) return -1;
	std::vector<uint8_t> bam = {'B','A','M',1};
	put32(bam, hdr.size()); bam.insert(bam.end(), hdr.begin(), hdr.end()); put32(bam, refs.size());
	for (size_t i = 0; i < refs.size(); ++i) {
		put32(bam, refs[i].name.size() + 1); bam.insert(bam.end(), refs[i].name.begin(), refs[i].name.end());
		bam.push_back(0); put32(bam, refs[i].len);
	}
	FILE *fp = fopen(fn, "wb");
	if (!fp) return -1;
	int ret = 0;
	line = 0; cap = 0;
	while ((n = getline(&line, &cap, sam)) >= 0) {
		if (n > 0 && line[n-1] == '\n') line[n-1] = 0;
		if (line[0] != '@') {
			std::vector<uint8_t> record;
			if (append_record(record, line, refs) < 0) { ret = -1; break; }
			if (!bam.empty() && bam.size() + record.size() > 60000) {
				if (bgzf_write(fp, bam.data(), bam.size()) < 0) { ret = -1; break; }
				bam.clear();
			}
			bam.insert(bam.end(), record.begin(), record.end());
		}
	}
	free(line);
	if (ret == 0) ret = bgzf_write(fp, bam.data(), bam.size());
	if (ret == 0) ret = bgzf_finish(fp);
	if (fclose(fp) != 0) ret = -1;
	return ret;
}
