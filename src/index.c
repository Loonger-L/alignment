#include <stdlib.h>
#include <assert.h>
#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#else
#include <unistd.h>
#endif
#include <fcntl.h>
#include <stdio.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#define __STDC_LIMIT_MACROS
#include "kthread.h"
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kvec.h"
#include "khash.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

#include <algorithm>
#include <unordered_map>
#include <utility>
#include <vector>

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

KHASH_MAP_INIT_STR(str, uint32_t)

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

typedef struct mm_idx_bucket_s {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct {
	int32_t st, en, max; // max is not used for now
	int32_t score:30, strand:2;
} mm_idx_intv1_t;

typedef struct mm_idx_intv_s {
	int32_t n, m;
	mm_idx_intv1_t *a;
} mm_idx_intv_t;

mm_idx_t *mm_idx_init(int w, int k, int b, int flag)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	if (!(mm_dbg_flag & 1)) mi->km = km_init();
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	uint32_t i;
	if (mi == 0) return;
	if (mi->h) kh_destroy(str, (khash_t(str)*)mi->h);
	if (mi->B) {
		for (i = 0; i < 1U<<mi->b; ++i) {
			free(mi->B[i].p);
			free(mi->B[i].a.a);
			kh_destroy(idx, (idxhash_t*)mi->B[i].h);
		}
	}
	if (mi->I) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->I[i].a);
		free(mi->I);
	}
	if (!mi->km) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	} else km_destroy(mi->km);
	mm_kmer_weight_destroy(mi->weight_db);
	free(mi->B); free(mi->S); free(mi);
}

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

void mm_idx_stat(const mm_idx_t *mi)
{
	int n = 0, n1 = 0;
	uint32_t i;
	uint64_t sum = 0, len = 0;
	fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n", __func__, mi->k, mi->w, mi->flag&MM_I_HPC, mi->n_seq);
	if (mi->w < 25) //display a warning if user adjusted window size is too low
		fprintf(stderr, "[M::%s] warning: setting a lower window size (-w) can increase runtime and memory\n", __func__);
	for (i = 0; i < mi->n_seq; ++i)
		len += mi->seq[i].len;
	for (i = 0; i < 1U<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	for (i = 0; i < 1U<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		khint_t k;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k)) {
				sum += kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
				if (kh_key(h, k)&1) ++n1;
			}
	}
	fprintf(stderr, "[M::%s::%.3f*%.2f] distinct minimizers: %d (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf\n",
			__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), n, 100.0*n1/n, (double)sum / n, (double)len / sum);
}

int mm_idx_index_name(mm_idx_t *mi)
{
	khash_t(str) *h;
	uint32_t i;
	int has_dup = 0, absent;
	if (mi->h) return 0;
	h = kh_init(str);
	for (i = 0; i < mi->n_seq; ++i) {
		khint_t k;
		k = kh_put(str, h, mi->seq[i].name, &absent);
		if (absent) kh_val(h, k) = i;
		else has_dup = 1;
	}
	mi->h = h;
	if (has_dup && mm_verbose >= 2)
		fprintf(stderr, "[WARNING] some database sequences have identical sequence names\n");
	return has_dup;
}

int mm_idx_name2id(const mm_idx_t *mi, const char *name)
{
	khash_t(str) *h = (khash_t(str)*)mi->h;
	khint_t k;
	if (h == 0) return -2;
	k = kh_get(str, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t i, st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	if (en > mi->seq[rid].len) en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st;
	en1 = mi->seq[rid].offset + en;
	for (i = st1; i < en1; ++i)
		seq[i - st1] = mm_seq4_get(mi->S, i);
	return en - st;
}

int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

/*********************************
 * Sort and generate hash tables *
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int n, n_keys;
	size_t j, start_a, start_p;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>8>>mi->b<<1, &absent);
			assert(absent && j == start_a + n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == (int32_t)start_p);

	// deallocate and clear b->a
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}

/******************
 * Generate index *
 ******************/

#include <string.h>
#include <zlib.h>
#include "bseq.h"

/*
 * A single distinct-frequency group.
 *
 * All canonical k-mers with identical reference occurrence count share
 * one group and therefore share the same score interval.
 */
struct mm_kmer_weight_group_t {
	uint64_t count;  // raw number of occurrences in the entire reference
	double lower;    // inclusive lower endpoint of [lower, lower + width)
	double width;    // interval width
};

/*
 * First prototype implementation.
 *
 * This is intentionally an exact hash table. It is not the final compact
 * representation; its purpose is to validate the graded weighting algorithm.
 */
struct mm_kmer_weight_db_s {
	std::unordered_map<uint64_t, uint32_t> group_id_by_kmer;
	std::vector<mm_kmer_weight_group_t> groups;
	uint64_t max_count;
};

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t *seq;
	mm128_v a;
} step_t;

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->B[a[i].x>>8&mask].a;
		kv_push(mm128_t, 0, *p, a[i]);
	}
}

/* Forward declarations for helpers defined later in this file. */
static inline uint64_t mm_canonical_kmer(uint64_t forward, uint64_t reverse);
static struct mm_kmer_weight_db_s *mm_kmer_weight_build_from_reference(const char *fn, int k, int is_hpc);
static struct mm_kmer_weight_db_s *mm_kmer_weight_build_from_sequences(int n, const char **seq, int k, int is_hpc);

static void mm_count_sequence_kmers(
	const char *str,
	int len,
	int k,
	int is_hpc,
	std::unordered_map<uint64_t, uint64_t> &counts)
{
	uint64_t shift1 = 2 * (k - 1);
	uint64_t mask = (1ULL << (2 * k)) - 1;
	uint64_t kmer[2] = {0, 0};

	/*
	 * Same purpose as tiny_queue_t in sketch.c:
	 * record homopolymer run lengths while -H is active.
	 *
	 * Winnowmap requires k <= 28, therefore 32 slots are sufficient.
	 */
	int span_queue[32];
	int queue_front = 0;
	int queue_count = 0;

	int l = 0;
	int kmer_span = 0;

	assert(str != 0);
	assert(len >= 0);
	assert(k > 0 && k <= 28);

	for (int i = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];

		/*
		 * Match mm_sketch():
		 * an ambiguous base breaks the current k-mer run.
		 */
		if (c >= 4) {
			l = 0;
			kmer_span = 0;
			queue_front = 0;
			queue_count = 0;
			kmer[0] = 0;
			kmer[1] = 0;
			continue;
		}

		if (is_hpc) {
			int skip_len = 1;

			/*
			 * Match mm_sketch() HPC behavior:
			 * collapse a homopolymer run to one base for k-mer construction,
			 * but preserve its original length in kmer_span.
			 */
			if (i + 1 < len &&
			    seq_nt4_table[(uint8_t)str[i + 1]] == c) {
				for (skip_len = 2; i + skip_len < len; ++skip_len) {
					if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
						break;
				}

				/*
				 * After this update, the outer loop increments i once more,
				 * so all bases in the run are skipped.
				 */
				i += skip_len - 1;
			}

			span_queue[(queue_front + queue_count) & 0x1f] = skip_len;
			++queue_count;
			kmer_span += skip_len;

			if (queue_count > k) {
				kmer_span -= span_queue[queue_front];
				queue_front = (queue_front + 1) & 0x1f;
				--queue_count;
			}
		} else {
			kmer_span = l + 1 < k ? l + 1 : k;
		}

		/* Generate forward and reverse-complement encoded k-mers. */
		kmer[0] = (kmer[0] << 2 | c) & mask;
		kmer[1] = (kmer[1] >> 2) | (3ULL ^ c) << shift1;

		/*
		 * Match mm_sketch(): symmetric k-mers are ignored because their
		 * strand is ambiguous.
		 *
		 * Do not increment l here: this intentionally matches the current
		 * implementation of mm_sketch().
		 */
		if (kmer[0] == kmer[1])
			continue;

		++l;

		/*
		 * Match mm_sketch() eligibility conditions for a real candidate.
		 */
		if (l >= k && kmer_span < 256) {
			uint64_t canonical = mm_canonical_kmer(kmer[0], kmer[1]);
			++counts[canonical];
		}
	}
}



static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch
		if (s->seq) {
			uint32_t old_m, m;
			assert((uint64_t)p->mi->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->mi->seq
			old_m = p->mi->n_seq, m = p->mi->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m)
				p->mi->seq = (mm_idx_seq_t*)krealloc(p->mi->km, p->mi->seq, m * sizeof(mm_idx_seq_t));
			// make room for p->mi->S
			if (!(p->mi->flag & MM_I_NO_SEQ)) {
				uint64_t sum_len, old_max_len, max_len;
				for (i = 0, sum_len = 0; i < s->n_seq; ++i) sum_len += s->seq[i].l_seq;
				old_max_len = (p->sum_len + 7) / 8;
				max_len = (p->sum_len + sum_len + 7) / 8;
				kroundup64(old_max_len); kroundup64(max_len);
				if (old_max_len != max_len) {
					p->mi->S = (uint32_t*)realloc(p->mi->S, max_len * 4);
					memset(&p->mi->S[old_max_len], 0, 4 * (max_len - old_max_len));
				}
			}
			// populate p->mi->seq
			for (i = 0; i < s->n_seq; ++i) {
				mm_idx_seq_t *seq = &p->mi->seq[p->mi->n_seq];
				uint32_t j;
				if (!(p->mi->flag & MM_I_NO_NAME)) {
					seq->name = (char*)kmalloc(p->mi->km, strlen(s->seq[i].name) + 1);
					strcpy(seq->name, s->seq[i].name);
				} else seq->name = 0;
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				// copy the sequence
				if (!(p->mi->flag & MM_I_NO_SEQ)) {
					for (j = 0; j < seq->len; ++j) { // TODO: this is not the fastest way, but let's first see if speed matters here
						uint64_t o = p->sum_len + j;
						int c = seq_nt4_table[(uint8_t)s->seq[i].seq[j]];
						mm_seq4_set(p->mi->S, o, c);
					}
				}
				// update p->sum_len and p->mi->n_seq
				p->sum_len += seq->len;
				s->seq[i].rid = p->mi->n_seq++;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			mm_bseq1_t *t = &s->seq[i];
			if (t->l_seq > 0)
				mm_sketch(0, t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, p->mi->flag&MM_I_HPC, &s->a, p->mi);
			else if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING] the length database sequence '%s' is 0\n", t->name);
			free(t->seq); free(t->name);
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		mm_idx_add(p->mi, s->a.n, s->a.a);
		kfree(0, s->a.a); free(s);
	}
    return 0;
}

uint64_t encodeKmer(const std::string &str)
{
	uint64_t kmer[2] = {0,0};
	int k = str.length();
	uint64_t shift1 = 2 * (k - 1);

	for (int i = 0; i < k; ++i) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		kmer[0] = kmer[0] << 2 | c;  
		kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1;
	}

	return kmer[0] < kmer[1]? kmer[0] : kmer[1];
}

static inline uint64_t mm_canonical_kmer(uint64_t forward, uint64_t reverse)
{
	return forward < reverse ? forward : reverse;
}

mm_idx_t *mm_idx_gen(
	const char *fn,
	int w,
	int k,
	int b,
	int flag,
	int mini_batch_size,
	int n_threads,
	uint64_t batch_size)
{
	if (fn == 0)
		return 0;

	/*
	 * Pass 1:
	 * Count all valid canonical k-mers over the complete reference.
	 */
	mm_kmer_weight_db_s *weight_db =
		mm_kmer_weight_build_from_reference(fn, k, flag & MM_I_HPC);

	if (weight_db == 0)
		return 0;

	/*
	 * Pass 2:
	 * Re-open the reference and perform the existing minimizer indexing pass.
	 */
	mm_bseq_file_t *fp = mm_bseq_open(fn);

	if (fp == 0) {
		fprintf(stderr, "[ERROR] failed to reopen reference '%s' for indexing\n", fn);
		mm_kmer_weight_destroy(weight_db);
		return 0;
	}

	pipeline_t pl;
	memset(&pl, 0, sizeof(pipeline_t));

	pl.mini_batch_size =
		(uint64_t)mini_batch_size < batch_size ?
		mini_batch_size : batch_size;

	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.mi = mm_idx_init(w, k, b, flag);
	pl.mi->weight_db = weight_db;

	kt_pipeline(n_threads < 3 ? n_threads : 3,
	            worker_pipeline,
	            &pl,
	            3);

	mm_bseq_close(fp);

	if (mm_verbose >= 3) {
		fprintf(stderr,
		        "[M::%s::%.3f*%.2f] collected weighted minimizers\n",
		        __func__,
		        realtime() - mm_realtime0,
		        cputime() / (realtime() - mm_realtime0));
	}

	mm_idx_post(pl.mi, n_threads);

	if (mm_verbose >= 3) {
		fprintf(stderr,
		        "[M::%s::%.3f*%.2f] sorted minimizers\n",
		        __func__,
		        realtime() - mm_realtime0,
		        cputime() / (realtime() - mm_realtime0));
	}

	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int flag, int n_threads)
{
	return mm_idx_gen(fn, w, k, 14, flag, 1 << 18,
	                  n_threads, UINT64_MAX);
}

mm_idx_t *mm_idx_str(
	int w,
	int k,
	int is_hpc,
	int bucket_bits,
	int n,
	const char **seq,
	const char **name)
{
	uint64_t sum_len = 0;
	mm128_v a = {0, 0, 0};
	mm_idx_t *mi;
	khash_t(str) *h;
	int i;
	int flag = 0;

	if (n <= 0 || seq == 0)
		return 0;

	if (k <= 0 || k > 28) {
		fprintf(stderr,
		        "[ERROR] k-mer size must be between 1 and 28; got %d\n",
		        k);
		return 0;
	}

	if (is_hpc)
		flag |= MM_I_HPC;

	if (name == 0)
		flag |= MM_I_NO_NAME;

	if (bucket_bits < 0)
		bucket_bits = 14;

	/*
	 * Pass 0:
	 * Compute total reference length before allocating the packed sequence
	 * storage mi->S. This is retained from the original implementation.
	 */
	for (i = 0; i < n; ++i) {
		if (seq[i] == 0) {
			fprintf(stderr,
			        "[ERROR] null sequence pointer at index %d\n",
			        i);
			return 0;
		}

		sum_len += strlen(seq[i]);
	}

	/*
	 * Pass 1:
	 *
	 * Count canonical k-mers across the COMPLETE set of input sequences,
	 * group identical counts, calculate each group's [lower, upper)
	 * interval, and build:
	 *
	 *     canonical k-mer -> group_id
	 *     group_id -> { count, lower, width }
	 *
	 * This pass must finish before mm_sketch() is called. Otherwise, a
	 * k-mer occurring in later sequences would not have its final global
	 * count when an earlier sequence is indexed.
	 */
	mm_kmer_weight_db_s *weight_db =
		mm_kmer_weight_build_from_sequences(n, seq, k, is_hpc);

	if (weight_db == 0)
		return 0;

	/*
	 * Create the standard minimizer index and attach the global weight DB.
	 * From this point, mm_sketch() can query every reference k-mer's
	 * group metadata through mi->weight_db.
	 */
	mi = mm_idx_init(w, k, bucket_bits, flag);

	if (mi == 0) {
		mm_kmer_weight_destroy(weight_db);
		return 0;
	}

	mi->weight_db = weight_db;

	mi->n_seq = n;
	mi->seq = (mm_idx_seq_t*)kcalloc(
		mi->km, n, sizeof(mm_idx_seq_t));

	/*
	 * This allocation is unchanged from the original implementation.
	 * Four bits are stored for each reference base.
	 */
	mi->S = (uint32_t*)calloc((sum_len + 7) / 8, 4);

	if (mi->seq == 0 || mi->S == 0) {
		fprintf(stderr,
		        "[ERROR] failed to allocate memory for in-memory reference index\n");
		mm_idx_destroy(mi);
		return 0;
	}

	mi->h = h = kh_init(str);

	/*
	 * Pass 2:
	 *
	 * Copy sequence metadata and packed bases into the index, then call
	 * mm_sketch(). The new applyWeight() will use mi->weight_db to obtain
	 * the k-mer frequency group and calculate its graded score.
	 */
	for (i = 0, sum_len = 0; i < n; ++i) {
		const char *s = seq[i];
		mm_idx_seq_t *p = &mi->seq[i];
		uint32_t j;

		if (!(mi->flag & MM_I_NO_NAME)) {
			int absent;

			/*
			 * The original implementation assumes name[i] is valid when
			 * name != NULL. Preserve that behavior but provide a clearer
			 * error for a malformed input array.
			 */
			if (name[i] == 0) {
				fprintf(stderr,
				        "[ERROR] null sequence name at index %d\n",
				        i);
				mm_idx_destroy(mi);
				free(a.a);
				return 0;
			}

			p->name = (char*)kmalloc(mi->km, strlen(name[i]) + 1);
			strcpy(p->name, name[i]);

			kh_put(str, h, p->name, &absent);
			assert(absent);
		} else {
			p->name = 0;
		}

		p->offset = sum_len;
		p->len = (uint32_t)strlen(s);

		/*
		 * Copy the reference into Winnowmap's packed 4-bit sequence store.
		 */
		for (j = 0; j < p->len; ++j) {
			int c = seq_nt4_table[(uint8_t)s[j]];
			uint64_t offset = sum_len + j;

			mm_seq4_set(mi->S, offset, c);
		}

		sum_len += p->len;

		if (p->len > 0) {
			a.n = 0;

			/*
			 * This is the existing minimizer-index construction call.
			 * The difference is that `mi->weight_db` already exists,
			 * so applyWeight() now uses:
			 *
			 * k-mer -> group_id -> {lower, width} -> graded score.
			 */
			mm_sketch(
				0,
				s,
				p->len,
				w,
				k,
				(uint32_t)i,
				is_hpc,
				&a,
				mi);

			mm_idx_add(mi, (int)a.n, a.a);
		}
	}

	free(a.a);

	/*
	 * Existing final index construction:
	 * sort minimizers and build lookup buckets.
	 */
	mm_idx_post(mi, 1);

	return mi;
}
/*************
 * index I/O *
 *************/

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t sum_len = 0;
	uint32_t x[5], i;

	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 5, fp);
	for (i = 0; i < mi->n_seq; ++i) {
		if (mi->seq[i].name) {
			uint8_t l = strlen(mi->seq[i].name);
			fwrite(&l, 1, 1, fp);
			fwrite(mi->seq[i].name, 1, l, fp);
		} else {
			uint8_t l = 0;
			fwrite(&l, 1, 1, fp);
		}
		fwrite(&mi->seq[i].len, 4, 1, fp);
		sum_len += mi->seq[i].len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ))
		fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
	fflush(fp);
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	char magic[4];
	uint32_t x[5], i;
	uint64_t sum_len = 0;
	mm_idx_t *mi;

	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 5, fp) != 5) return 0;
	mi = mm_idx_init(x[0], x[1], x[2], x[4]);
	mi->n_seq = x[3];
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		mm_idx_seq_t *s = &mi->seq[i];
		fread(&l, 1, 1, fp);
		if (l) {
			s->name = (char*)kmalloc(mi->km, l + 1);
			fread(s->name, 1, l, fp);
			s->name[l] = 0;
		}
		fread(&s->len, 4, 1, fp);
		s->offset = sum_len;
		sum_len += s->len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, fp);
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, fp);
		fread(&size, 4, 1, fp);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
			fread(x, 8, 2, fp);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ)) {
		mi->S = (uint32_t*)malloc((sum_len + 7) / 8 * 4);
		fread(mi->S, 4, (sum_len + 7) / 8, fp);
	}
	return mi;
}

int64_t mm_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	int64_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
#ifdef WIN32
	if ((off_end = _lseeki64(fd, 0, SEEK_END)) >= 4) {
		_lseeki64(fd, 0, SEEK_SET);
#else
	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
#endif // WIN32
		ret = read(fd, magic, 4);
		if (ret == 4 && strncmp(magic, MM_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
{
	int64_t is_idx;
	mm_idx_reader_t *r;
	is_idx = mm_idx_is_idx(fn);
	if (is_idx < 0) return 0; // failed to open the index
	r = (mm_idx_reader_t*)calloc(1, sizeof(mm_idx_reader_t));
    r->fn = strdup(fn);
	r->is_idx = is_idx;
	if (opt) r->opt = *opt;
	else mm_idxopt_init(&r->opt);
	if (r->is_idx) {
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
	} else r->fp.seq = mm_bseq_open(fn);
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

void mm_idx_reader_close(mm_idx_reader_t *r)
{
	if (r == 0)
		return;

	if (r->is_idx)
		fclose(r->fp.idx);
	else if (r->fp.seq)
		mm_bseq_close(r->fp.seq);

	if (r->fp_out)
		fclose(r->fp_out);

	free(r->fn);
	free(r);
}

mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
{
	mm_idx_t *mi;
	if (r->is_idx) {
		mi = mm_idx_load(r->fp.idx);
		if (mi && mm_verbose >= 2 && (mi->k != r->opt.k || mi->w != r->opt.w || (mi->flag&MM_I_HPC) != (r->opt.flag&MM_I_HPC)))
			fprintf(stderr, "[WARNING]\033[1;31m Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\033[0m\n");
	} else {
	if (r->done)
		return 0;

	mi = mm_idx_gen(r->fn,
	                r->opt.w,
	                r->opt.k,
	                r->opt.bucket_bits,
	                r->opt.flag,
	                r->opt.mini_batch_size,
	                n_threads,
	                r->opt.batch_size);

	r->done = 1;
    }
	if (mi) {
		if (r->fp_out) mm_idx_dump(r->fp_out, mi);
		mi->index = r->n_parts++;
	}
	return mi;
}

int mm_idx_reader_eof(const mm_idx_reader_t *r) // TODO: in extremely rare cases, mm_bseq_eof() might not work
{
	return r->is_idx? (feof(r->fp.idx) || ftell(r->fp.idx) == r->idx_size) : mm_bseq_eof(r->fp.seq);
}

#include <ctype.h>
#include <zlib.h>
#include "ksort.h"
#include "kseq.h"
KSTREAM_DECLARE(gzFile, gzread)

#define sort_key_bed(a) ((a).st)
KRADIX_SORT_INIT(bed, mm_idx_intv1_t, sort_key_bed, 4)

mm_idx_intv_t *mm_idx_read_bed(const mm_idx_t *mi, const char *fn, int read_junc)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	mm_idx_intv_t *I;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	I = (mm_idx_intv_t*)calloc(mi->n_seq, sizeof(*I));
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		mm_idx_intv_t *r;
		mm_idx_intv1_t t = {-1,-1,-1,-1,0};
		char *p, *q, *bl, *bs;
		int32_t i, id = -1, n_blk = 0;
		for (p = q = str.s, i = 0;; ++p) {
			if (*p == 0 || *p == '\t') {
				int32_t c = *p;
				*p = 0;
				if (i == 0) { // chr
					id = mm_idx_name2id(mi, q);
					if (id < 0) break; // unknown name; TODO: throw a warning
				} else if (i == 1) { // start
					t.st = atol(q); // TODO: watch out integer overflow!
					if (t.st < 0) break;
				} else if (i == 2) { // end
					t.en = atol(q);
					if (t.en < 0) break;
				} else if (i == 4) { // BED score
					t.score = atol(q);
				} else if (i == 5) { // strand
					t.strand = *q == '+'? 1 : *q == '-'? -1 : 0;
				} else if (i == 9) {
					if (!isdigit(*q)) break;
					n_blk = atol(q);
				} else if (i == 10) {
					bl = q;
				} else if (i == 11) {
					bs = q;
					break;
				}
				if (c == 0) break;
				++i, q = p + 1;
			}
		}
		if (id < 0 || t.st < 0 || t.st >= t.en) continue;
		r = &I[id];
		if (i >= 11 && read_junc) { // BED12
			int32_t st, sz, en;
			st = strtol(bs, &bs, 10); ++bs;
			sz = strtol(bl, &bl, 10); ++bl;
			en = t.st + st + sz;
			for (i = 1; i < n_blk; ++i) {
				mm_idx_intv1_t s = t;
				if (r->n == r->m) {
					r->m = r->m? r->m + (r->m>>1) : 16;
					r->a = (mm_idx_intv1_t*)realloc(r->a, sizeof(*r->a) * r->m);
				}
				st = strtol(bs, &bs, 10); ++bs;
				sz = strtol(bl, &bl, 10); ++bl;
				s.st = en, s.en = t.st + st;
				en = t.st + st + sz;
				if (s.en > s.st) r->a[r->n++] = s;
			}
		} else {
			if (r->n == r->m) {
				r->m = r->m? r->m + (r->m>>1) : 16;
				r->a = (mm_idx_intv1_t*)realloc(r->a, sizeof(*r->a) * r->m);
			}
			r->a[r->n++] = t;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return I;
}

int mm_idx_bed_read(mm_idx_t *mi, const char *fn, int read_junc)
{
	int32_t i;
	if (mi->h == 0) mm_idx_index_name(mi);
	mi->I = mm_idx_read_bed(mi, fn, read_junc);
	if (mi->I == 0) return -1;
	for (i = 0; i < mi->n_seq; ++i) // TODO: eliminate redundant intervals
		radix_sort_bed(mi->I[i].a, mi->I[i].a + mi->I[i].n);
	return 0;
}

int mm_idx_bed_junc(const mm_idx_t *mi, int32_t ctg, int32_t st, int32_t en, uint8_t *s)
{
	int32_t i, left, right;
	mm_idx_intv_t *r;
	memset(s, 0, en - st);
	if (mi->I == 0 || ctg < 0 || ctg >= mi->n_seq) return -1;
	r = &mi->I[ctg];
	left = 0, right = r->n;
	while (right > left) {
		int32_t mid = left + ((right - left) >> 1);
		if (r->a[mid].st >= st) right = mid;
		else left = mid + 1;
	}
	for (i = left; i < r->n; ++i) {
		if (st <= r->a[i].st && en >= r->a[i].en && r->a[i].strand != 0) {
			if (r->a[i].strand > 0) {
				s[r->a[i].st - st] |= 1, s[r->a[i].en - 1 - st] |= 2;
			} else {
				s[r->a[i].st - st] |= 8, s[r->a[i].en - 1 - st] |= 4;
			}
		}
	}
	return left;
}

static std::unordered_map<uint64_t, uint64_t>
mm_count_reference_kmers(const char *fn, int k, int is_hpc)
{
	std::unordered_map<uint64_t, uint64_t> counts;

	mm_bseq_file_t *fp = mm_bseq_open(fn);
	if (fp == 0) {
		fprintf(stderr, "[ERROR] failed to open reference '%s' for k-mer counting\n", fn);
		return counts;
	}

	for (;;) {
		int n_seq = 0;
		mm_bseq1_t *seq = mm_bseq_read(fp, 1 << 18, 0, &n_seq);

		if (seq == 0)
			break;

		for (int i = 0; i < n_seq; ++i) {
			if (seq[i].l_seq > 0)
				mm_count_sequence_kmers(seq[i].seq, seq[i].l_seq, k,
				                        is_hpc, counts);

			free(seq[i].seq);
			free(seq[i].name);
		}

		free(seq);
	}

	mm_bseq_close(fp);
	return counts;
}    
    
static mm_kmer_weight_db_s *
mm_kmer_weight_build_from_counts(
	std::unordered_map<uint64_t, uint64_t> &counts)
{
	if (counts.empty()) {
		fprintf(stderr,
		        "[ERROR] no valid canonical k-mers were found in the reference\n");
		return 0;
	}

	/*
	 * Convert the mutable counting hash table to a sortable compact vector.
	 *
	 * pair.first  = canonical k-mer
	 * pair.second = raw occurrence count
	 */
	std::vector<std::pair<uint64_t, uint64_t> > entries;
	entries.reserve(counts.size());

	uint64_t max_count = 0;

	for (std::unordered_map<uint64_t, uint64_t>::const_iterator it =
		     counts.begin();
	     it != counts.end();
	     ++it) {
		entries.push_back(std::make_pair(it->first, it->second));

		if (it->second > max_count)
			max_count = it->second;
	}

	assert(max_count > 0);

	/*
	 * Deterministic ordering:
	 *
	 * 1. raw count ascending;
	 * 2. canonical k-mer ascending for stable behavior/debugging.
	 */
	std::sort(entries.begin(), entries.end(),
		[](const std::pair<uint64_t, uint64_t> &a,
		   const std::pair<uint64_t, uint64_t> &b) {
			if (a.second != b.second)
				return a.second < b.second;

			return a.first < b.first;
		});

	mm_kmer_weight_db_s *db = new mm_kmer_weight_db_s;
	db->max_count = max_count;
	db->group_id_by_kmer.reserve(entries.size());

	/*
	 * `entries` now holds everything needed from `counts`.
	 * Release the count hash table before building the permanent dictionary.
	 */
	counts.clear();
	counts.rehash(0);

	size_t begin = 0;

	while (begin < entries.size()) {
		uint64_t current_count = entries[begin].second;
		size_t end = begin + 1;

		/* Locate the end of this equal-count group. */
		while (end < entries.size() &&
		       entries[end].second == current_count)
			++end;

		if (db->groups.size() >= UINT32_MAX) {
			fprintf(stderr,
			        "[ERROR] too many distinct k-mer count groups "
			        "for uint32_t group IDs\n");
			delete db;
			return 0;
		}

		uint32_t group_id = (uint32_t)db->groups.size();

		/*
		 * r_current = current_count / max_count
		 */
		double r_current =
			(double)current_count / (double)max_count;

		double lower;

		if (group_id == 0) {
			/*
			 * First frequency group:
			 *
			 * interval = [0, r0)
			 */
			lower = 0.0;
		} else {
			/*
			 * Subsequent frequency groups:
			 *
			 * interval = [r_previous * r_current, r_current)
			 */
			uint64_t previous_count =
				db->groups[group_id - 1].count;

			double r_previous =
				(double)previous_count / (double)max_count;

			lower = r_previous * r_current;
		}

		mm_kmer_weight_group_t group;
		group.count = current_count;
		group.lower = lower;
		group.width = r_current - lower;

		/*
		 * Sanity checks for the required half-open interval [0, 1).
		 */
		assert(group.lower >= 0.0);
		assert(group.width > 0.0);
		assert(group.lower + group.width <= 1.0);

		db->groups.push_back(group);

		/*
		 * Every k-mer in this count group receives the same group ID.
		 * No k-mer stores its own floating-point interval.
		 */
		for (size_t i = begin; i < end; ++i) {
			db->group_id_by_kmer.emplace(entries[i].first, group_id);
		}

		begin = end;
	}

	fprintf(stderr,
	        "[M::%s] built graded k-mer weight database: "
	        "distinct-kmers=%zu, frequency-groups=%zu, max-count=%" PRIu64 "\n",
	        __func__,
	        db->group_id_by_kmer.size(),
	        db->groups.size(),
	        db->max_count);

	if (mm_verbose >= 3) {
		size_t preview =
			db->groups.size() < 10 ? db->groups.size() : 10;

		for (size_t i = 0; i < preview; ++i) {
			const mm_kmer_weight_group_t &group = db->groups[i];

			fprintf(stderr,
			        "[M::%s] group=%zu count=%" PRIu64
			        " interval=[%.12g, %.12g)\n",
			        __func__,
			        i,
			        group.count,
			        group.lower,
			        group.lower + group.width);
		}
	}

	return db;
}
    
static mm_kmer_weight_db_s *
mm_kmer_weight_build_from_sequences(
	int n,
	const char **seq,
	int k,
	int is_hpc)
{
	std::unordered_map<uint64_t, uint64_t> counts;

	if (n <= 0 || seq == 0) {
		fprintf(stderr,
		        "[ERROR] invalid sequence array for k-mer weight construction\n");
		return 0;
	}

	fprintf(stderr,
	        "[M::%s::%.3f*%.2f] counting canonical k-mers from "
	        "%d in-memory reference sequence(s)\n",
	        __func__,
	        realtime() - mm_realtime0,
	        cputime() / (realtime() - mm_realtime0),
	        n);

	for (int i = 0; i < n; ++i) {
		if (seq[i] == 0) {
			fprintf(stderr,
			        "[ERROR] null sequence pointer at index %d\n",
			        i);
			return 0;
		}

		int len = (int)strlen(seq[i]);

		if (len > 0)
			mm_count_sequence_kmers(seq[i], len, k, is_hpc, counts);
	}

	fprintf(stderr,
	        "[M::%s::%.3f*%.2f] counted %zu distinct canonical k-mers\n",
	        __func__,
	        realtime() - mm_realtime0,
	        cputime() / (realtime() - mm_realtime0),
	        counts.size());

	return mm_kmer_weight_build_from_counts(counts);
}
    
static mm_kmer_weight_db_s *
mm_kmer_weight_build_from_reference(
	const char *fn,
	int k,
	int is_hpc)
{
	fprintf(stderr,
	        "[M::%s::%.3f*%.2f] counting reference k-mers\n",
	        __func__,
	        realtime() - mm_realtime0,
	        cputime() / (realtime() - mm_realtime0));

	std::unordered_map<uint64_t, uint64_t> counts =
		mm_count_reference_kmers(fn, k, is_hpc);

	fprintf(stderr,
	        "[M::%s::%.3f*%.2f] counted %zu distinct canonical k-mers\n",
	        __func__,
	        realtime() - mm_realtime0,
	        cputime() / (realtime() - mm_realtime0),
	        counts.size());

	return mm_kmer_weight_build_from_counts(counts);
}
    
int mm_kmer_weight_get(
	const mm_idx_t *mi,
	uint64_t kmer,
	double *lower,
	double *width)
{
	if (mi == 0 || mi->weight_db == 0 ||
	    lower == 0 || width == 0)
		return 0;

	const mm_kmer_weight_db_s *db = mi->weight_db;

	std::unordered_map<uint64_t, uint32_t>::const_iterator it =
		db->group_id_by_kmer.find(kmer);

	if (it == db->group_id_by_kmer.end())
		return 0;

	uint32_t group_id = it->second;

	assert(group_id < db->groups.size());

	const mm_kmer_weight_group_t &group =
		db->groups[group_id];

	*lower = group.lower;
	*width = group.width;

	return 1;
}

void mm_kmer_weight_destroy(mm_kmer_weight_db_s *db)
{
	delete db;
}