# Blocked LD computation in `make-bld`

This document describes how the original `make-bld` LD computation worked, why it
was slow on dense SNP panels, and how it was restructured into the current
blocked (BLAS-3) implementation while preserving output compatibility.

Relevant code: `src/bfile.cpp`, function `ld_report()` (and helpers `getMaxNum()`,
`initX()`, `makex_xVec_subset()`).

---

## 1. What `make-bld` produces

`smr --bfile <prefix> --make-bld --r|--r2 --ld-wind <Kb> --out <out>` computes,
for every SNP `i` in the (possibly chromosome-restricted) panel, the LD
correlation (or squared correlation) between `i` and every following SNP `j`
within `--ld-wind` kilobases, and writes them to `<out>.bld` (binary) plus
`<out>.esi` (text SNP information).

`.bld` layout (unchanged by this work):

| field                 | type     | content                                                               |
|-----------------------|----------|-----------------------------------------------------------------------|
| 16 × int              | header   | indicator (`r`/`r2`), individual count, SNP count, window Kb, padding |
| uint64                | `valnum` | total number of LD values                                             |
| uint64 × (`snpNum+1`) | `cols`   | cumulative value counts per SNP (`cols[i+1] = cols[i] + ldnum_i`)     |
| float × `valnum`      | values   | LD r (or r²) per SNP pair, ordered by increasing `j`                  |

`getMaxNum()` precomputes `cols` and returns `m`, the smallest power of two
greater than the maximum per-SNP partner count. `m` determines the width of the
genotype buffer.

## 2. The original algorithm

The original code kept a **rolling genotype window** `X` (individuals × `m`,
float). For SNP `i`, its genotype lives in column `i & (m-1)`; the columns hold
SNPs `[i, i+m)` at any time. Per SNP `i`:

```cpp
VectorXf ldv = X.col(start) / (X.rows() - 1);   // start = i & (m-1)
ldv = X.transpose() * ldv;                       // gemv over the WHOLE X
```

then the in-window band `[st..ed]` was written to the file, and the next
genotype (SNP `i+m`) was decoded into the freed column (`makex_xVec_subset`).
LD values are therefore Pearson-style dot products of *centered/standardized*
genotype columns, produced by `initX()`/`makex_xVec_subset()`.

### Why it is slow at production scale

On the UK Biobank chr1 test panel (**600,843 SNPs × 10,000 individuals**, densest
4 Mb window = 14,847 SNPs → `m = 16384`):

* `X` = 10,000 × 16,384 floats = **655 MB**.
* The gemv `X.transpose() * ldv` **re-reads all 655 MB once per SNP**:
  600K × 655 MB ≈ **~400 TB of memory traffic** (hours at tens of GB/s).
* The flops (~60-200 TF) are secondary; the kernel is a BLAS-2 gemv, which is
  memory-bound and runs at a fraction of peak BLAS-3 throughput.
* Measured (`--ld-wind 100`, `m = 512`): **10 min 26 s** for the full panel.
  At the production default `--ld-wind 4000` (`m = 16384`) the same run is
  infeasible in practice (many hours).

The unit costs that remain even at small `m`: per-SNP genotype decode
(`makex_xVec_subset`, bit-unpacking through `vector<bool>` proxies), one-time
`calcu_mu`, and the `.bld` writes.

## 3. The blocked algorithm (this work)

### 3.0 What `X` contains, and the role of `X.rows() - 1` (worked example)

The rolling buffer `X` holds **standardized genotype columns**, produced by
`initX()` (initial `m` columns) and `makex_xVec_subset()` (rolling updates):

1. `calcu_mu` computes per-SNP mean dosage `μ_k` over **non-missing** individuals.
2. Per individual `i` and SNP `k`: dosage `g_i ∈ {0,1,2}` (allele count, flipped
   to the reference allele: `g_i → 2 - g_i` when `allele1 != ref_A`), then
   `z_i = g_i - μ_k`; **missing genotypes become `z_i = 0`** (mean imputation).
3. Each column is scaled by `1/sd_k`, where `sd_k = sqrt(Σ z_i² / (n-1))` is the
   empirical standard deviation (`n = X.rows()` = kept individuals).

So for SNPs `a` and `b`, the Pearson correlation is

```
r = Σ_i x_i,a · x_i,b / (n - 1)
```

because the columns already have `Σ x²/(n-1) = 1`. The original loop divides
**one** column by `(X.rows() - 1)` *before* the dot product, which produces `r`
directly; squaring elementwise produces `r²`.

**Example** (2 SNPs, 5 individuals, one missing genotype for SNP A):

```
individual:        1    2    3    4    5
SNP A dosage g:    0    1    2    1   NA
SNP B dosage g:    1    1    2    0    2

μ_A = 4/4 = 1.0            μ_B = 6/5 = 1.2
z_A = [-1, 0, +1, 0, 0]    z_B = [-0.2, -0.2, +0.8, -1.2, +0.8]   (NA -> 0)
Σz_A² = 2.0  -> sd_A = sqrt(2.0/4) = 0.7071
Σz_B² = 2.8  -> sd_B = sqrt(2.8/4) = 0.8367

x_A = [-1.4142, 0, +1.4142, 0, 0]
x_B = [-0.2390, -0.2390, +0.9565, -1.4348, +0.9565]

Σ x_A·x_B = (-1.4142)(-0.2390) + 0 + (1.4142)(0.9565) + 0 + 0 = 0.6906
r = 0.6906 / (5 - 1) = 0.1726
```

Note the convention: the denominator is `n-1` over **all** kept individuals,
and `μ`, `sd` include mean-imputed zeros — i.e. the GCTA/PLINK-style LD, *not*
the textbook pairwise-complete correlation (which would give `r = 0.5` here
using only the 4 complete pairs). The `.bld` values follow the former
convention consistently.

The computation is restructured to **one BLAS-3 gemm per block of SNPs** instead
of one BLAS-2 gemv per SNP. `blk = min(512, m)`; SNPs are processed in blocks
`[i0, i0+b)`.

### Correctness constraints

Two subtleties make a naive blocked version *wrong* (an initial version of this
patch had exactly this bug and produced LD values off by ~0.15 on average):

1. **Rolling updates overwrite the block's own columns.** In the original loop
   the update for SNP `i` (decoding SNP `i+m` into column `i & (m-1)`) happens
   *after* SNP `i`'s LD computation. If all `b` updates run *before* the block's
   gemm, the buffer no longer holds the block's own genotypes.
2. **The partner-column space collides.** A partner SNP `q` maps to column
   `q & (m-1)`. After pre-updates, that column holds SNP `q+m`, not `q` —
   so partners *within* the block would read stale genotypes.

### The scheme

For each block `[i0, i0+b)` (`start0 = i0 & (m-1)`):

1. **Snapshot**: `Xb = X.middleCols(start0, b) / (X.rows() - 1)` — a copy of the
   block's own genotypes, pre-scaled exactly like the original `ldv`.
2. **C_main**: `C = X.transpose() * Xb` (before any updates). Valid for partners
   `q ∈ [i0, i0+m)`.
3. **Rolling updates** (unchanged order): decode SNPs `[i0+m, i0+b+m)` into the
   rolling columns.
4. **C_extra**: `C_extra = X.middleCols(start0, b).transpose() * Xb` — the
   newly decoded columns, valid for partners `q ∈ [i0+m, i0+m+b)`.
5. **Assemble each SNP's `ldv` into a scratch vector** so it is *exactly* what
   the original loop produced:
   `scratch = C.col(t)` for all `m` entries, then patch the tail
   `scratch[start0 .. start0+t) = C_extra.col(t)[0..t)`.

   The cur-space partition is exact: forward in-block partners map to
   `[start0+t, m)` (C_main), tail partners `q ≥ i0+m` map to
   `[start0, start0+t)` (C_extra); the two ranges never overlap, and the tail
   never exceeds `t` entries because a SNP's band has at most `m-1` partners.
6. The **band scan and segment writes are byte-for-byte the original logic**
   (same `cur = (start+j) & (m-1)` mapping, same mid-band flush at `ed == m-1`,
   same partner order). `vc == valnum` accounting is unchanged.

For `--r2`, squaring is applied to `C` and `C_extra` once per block, identical to
the original per-element square.

`bitmod == false` (when `m` is capped at the SNP count) works unchanged: no
rolling updates happen, `C_extra` is empty, all partners come from `C_main`.

### Why it is fast

* The rolling `X` is now read **once per block of 512 SNPs** instead of once per
  SNP: traffic drops from ~400 TB to ~0.5-1 TB (full panel, `m = 16384`).
* The kernel is BLAS-3 `sgemm` (MKL, multithreaded) instead of memory-bound
  BLAS-2 `sgemv`: the same ~200 TF of math executes at hundreds of GFLOP/s.
* Genotype decode order and count are unchanged (still once per SNP), and
  `getMaxNum()`/header/cols generation are untouched.

## 4. Output compatibility

| Item                    | Result                                                                                                        |
|-------------------------|---------------------------------------------------------------------------------------------------------------|
| `.esi`                  | byte-identical                                                                                                |
| `.bld` header (16 ints) | byte-identical                                                                                                |
| `valnum`, `cols`        | byte-identical                                                                                                |
| LD values               | **drift ≤ ~2e-05 absolute** (max), mean ~1e-06, p99 ~1e-05                                                    |
| Downstream `.smr`       | identical SNP sets (`nsnp_HEIDI` unchanged on all 1523 probes); `p_HEIDI` shifts at 5th-6th significant digit |

The drift comes from summation-order differences between `sgemm` (blocked,
possibly FMA) and `sgemv` (sequential), and from pre-scaling `Xb` by `1/(n-1)`
before the product. It is at the level of float rounding error and is
scientifically irrelevant (typical HEIDI thresholds are 0.01-0.05), but it does
mean `.bld` files are **not bit-identical** to the old implementation.

## 5. Measured performance

Full UK Biobank chr1 panel (600,843 SNPs, 10,000 individuals):

| Command                         | old (gemv)                                | new (blocked)   |
|---------------------------------|-------------------------------------------|-----------------|
| `--make-bld --r --ld-wind 100`  | **10 min 26 s**                           | **2 min 59 s**  |
| `--make-bld --r --ld-wind 4000` | infeasible (hours; killed at 30 min, <5%) | **23 min 19 s** |

Notes:

* At small `m` the shared costs dominate (per-SNP decode, one-time `calcu_mu`,
  `.bld` writes — 26 GB at `--ld-wind 4000`), so the ratio is ~3.5×.
  At the production default (`--ld-wind 4000`, `m = 16384`) where the gemv
  traffic dominates, the speedup is in the **~50-100×** range.
* Remaining cost breakdown at `--ld-wind 4000` (new): gemm ~16 min (BLAS-3,
  197 TF), genotype decode ~2 min, scratch assembly ~1 min, file I/O ~2 min.

## 6. Future work (not done)

* **OpenMP over blocks**: blocks are dependent through the rolling update, but
  the gemm itself can be parallelized per block; the decode phase could also be
  parallelized (independent columns). MKL already multithreads the gemm.
* **`--chr` full-genome runs**: the same algorithm applies per chromosome;
  runtime scales with panel size and window density, not chromosome count.
* **Float → double accumulation** in the gemm would reduce the ~1e-5 absolute
  error to ~1e-8 at ~2× the flops; currently deemed unnecessary.
