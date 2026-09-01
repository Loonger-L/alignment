# Guide: Checking, Reviewing, and Merging Coding-Agent Changes in `Loonger-L/alignment`

This guide details how to check task status, locate Pull Requests (PRs), review changes, and merge completed coding-agent branches into the `main` branch of `Loonger-L/alignment`.

---

## 1. Checking Task Status & Locating the PR

When the Copilot Coding Agent completes a code modification task, it pushes the changes to a dedicated branch and creates or updates a Pull Request.

### Locating via GitHub Web UI
1. Navigate to the repository on GitHub: [`https://github.com/Loonger-L/alignment`](https://github.com/Loonger-L/alignment).
2. Click on the **Pull requests** tab.
3. Locate the PR created by the coding agent (e.g., `copilot-swe-agent[bot]`).
   - **Target PR**: [PR #1: Fix build: forward declarations in index.c, remove meryl, fix C++11 hex float](https://github.com/Loonger-L/alignment/pull/1)
   - **Branch**: `copilot/fix-kmer-weighting-changes`
   - **Base Branch**: `main`

### Locating via Command Line (Git / GitHub CLI)
```bash
# List remote branches created by agents
git fetch origin
git branch -r | grep copilot

# Or using GitHub CLI
gh pr list
```

---

## 2. Reviewing the Code Changes

Before merging, review the modifications to ensure they align with requirements and compile cleanly.

### Summary of Completed Changes in PR #1 (`copilot/fix-kmer-weighting-changes`)
1. **`src/index.c` — Forward Declarations Fix**
   - Added forward declarations for `mm_canonical_kmer`, `mm_kmer_weight_build_from_reference`, and `mm_kmer_weight_build_from_sequences` before their first usage in `mm_count_sequence_kmers`, `mm_idx_gen`, and `mm_idx_str`.

2. **`Makefile` — Removed Meryl Build Dependency**
   - Removed external meryl compilation commands (`+$(MAKE) -C ext/meryl/src`) and `BUILDSTACKTRACE=0` exports because exact k-mer counting is now performed internally from reference sequences.

3. **`src/sketch.c` — C++11 Compatibility Fix**
   - Replaced hexadecimal float literal `0x1.0p-64` (C++17 feature) with `ldexp(1.0, -64)` to adhere strictly to `-std=c++11`.

### Local Inspection & Build Testing
To review and build the branch locally:
```bash
# 1. Fetch and checkout the agent branch
git fetch origin copilot/fix-kmer-weighting-changes:copilot/fix-kmer-weighting-changes
git checkout copilot/fix-kmer-weighting-changes

# 2. View differences against main
git diff main..copilot/fix-kmer-weighting-changes

# 3. Clean and compile
make clean && make -j8

# 4. Verify binary execution
./bin/winnowmap --help
```

---

## 3. Merging the PR into `main`

Once satisfied with the review and test build, merge the changes into `main`.

### Method A: GitHub Web UI (Recommended)
1. Open [PR #1](https://github.com/Loonger-L/alignment/pull/1) on GitHub.
2. If the PR status is **Draft**, click the **"Ready for review"** button at the bottom of the PR page.
3. Click **"Merge pull request"** (or choose "Squash and merge" / "Rebase and merge").
4. Click **"Confirm merge"**.
5. (Optional) Click **"Delete branch"** to clean up `copilot/fix-kmer-weighting-changes`.

### Method B: Git Command Line
```bash
# 1. Checkout and update main
git checkout main
git pull origin main

# 2. Merge the agent's branch
git merge origin/copilot/fix-kmer-weighting-changes -m "Merge coding-agent fixes for k-mer weighting build"

# 3. Push to main branch on GitHub
git push origin main
```

### Method C: GitHub CLI (`gh`)
```bash
# Mark PR as ready and merge
gh pr ready 1
gh pr merge 1 --merge
```

---

## Summary Checklist for Merging Agent Work
- [x] Agent completed task & created PR (#1)
- [x] Code builds cleanly with `make -j8`
- [x] PR converted from Draft to Ready for Review
- [x] PR merged into `main` branch
