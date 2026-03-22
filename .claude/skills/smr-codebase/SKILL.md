---
name: smr-codebase
description: >
  Teaches AI coding agents how to work effectively with the SMR
  (Summary-data-based Mendelian Randomization) C++ codebase. Covers
  architecture, data formats, build system, statistical methods, and
  common pitfalls.
triggers:
  - SMR
  - smr codebase
  - mendelian randomization
  - eQTL
  - BESD
  - HEIDI
  - GWAS summary
---

# SMR Codebase Skill

## Project Overview

**SMR** (Summary-data-based Mendelian Randomization) is a C++17 command-line tool that tests for pleiotropic associations between molecular traits (gene expression, DNA methylation, protein abundance) and complex disease traits using GWAS and xQTL summary-level data.

The core scientific question: *Is a SNP's effect on a phenotype mediated through gene expression?*

- **Methods**: SMR test + HEIDI (Heterogeneity in Dependent Instruments) test
- **Version**: 1.4.0
- **Lab**: Yang Lab, Westlake University
- **License**: MIT
- **Reference**: Zhu et al. (2016), Nature Genetics

---

## Tech Stack and Key Dependencies

| Component | Detail |
|-----------|--------|
| Language | C++17 (required) |
| Linear algebra | Eigen 3.4.0 (auto-fetched via CMake FetchContent) |
| Compression | zlib 1.3.1 (auto-fetched via CMake FetchContent) |
| Parallelism | OpenMP (required; controls `--thread-num`) |
| Build (primary) | CMake ≥ 3.16 |
| Build (legacy) | GNU Make + GCC 13.2.0 (Homebrew on macOS) |

---

## Directory Structure

```
SMR/
├── src/                    # All .cpp implementation files
│   ├── SMR.cpp             # CLI entry point, option parsing, dispatch (~1,328 lines)
│   ├── SMR_data.cpp        # Core data I/O, filtering, allele checks (~7,530 lines)
│   ├── SMR_data_p1.cpp     # BESD file creation, lookup, effect-size estimation (~4,528 lines)
│   ├── SMR_data_p2.cpp     # Dense BESD creation, data manipulation (~4,180 lines)
│   ├── SMR_data_p3.cpp     # Additional data processing (~3,961 lines)
│   ├── SMR_plot.cpp        # Plotting / visualization (~2,481 lines)
│   ├── CommFunc.cpp        # Common utilities: file I/O, vector ops (~487 lines)
│   ├── StatFunc.cpp        # P-value calculations (~797 lines)
│   ├── StrFunc.cpp         # String parsing (~311 lines)
│   ├── bfile.cpp           # PLINK bed/bim/fam binary format support (~1,580 lines)
│   └── dcdflib.cpp         # CDF library: beta/gamma/chi2/normal/F/t (~9,205 lines)
│
├── include/                # All .hpp/.h header files (mirrors src/)
│   ├── SMR.hpp             # CLI option flags and global state
│   ├── SMR_data.hpp        # Core data structures (gwasData, eqtlInfo, bInfo, SMRWK, SMRRLT)
│   ├── config.h            # Version string: SMR_VERSION "1.4.0"
│   └── ...                 # One header per src file
│
├── cmake/
│   └── hpc-toolchain.cmake # HPC cluster toolchain settings
├── scripts/hpc/build.sh    # HPC build script
├── CMakeLists.txt          # Primary build definition (155 lines)
├── Makefile                # Legacy build (46 lines)
└── docs/changlog.txt       # Changelog
```

---

## Build Commands

### CMake (preferred)
```bash
mkdir build && cd build
cmake ..                          # Release build
cmake -DCMAKE_BUILD_TYPE=DEBUG .. # Debug build
cmake -DBUILD_STATIC=ON ..        # Static binary
make -j$(nproc)
```
CMake auto-downloads Eigen and zlib — no manual dependency install needed.

### Make (legacy, macOS/Linux)
```bash
make smr                          # Shared binary
make smr_static                   # Static binary
make clean
make install                      # Copies to /usr/local/bin
```
Override dependency paths:
```bash
make EIGEN_PATH=/custom/eigen ZLIB_INCLUDE=/custom/zlib/include ZLIB_LIB=/custom/zlib/lib
```

### HPC Cluster
```bash
bash scripts/hpc/build.sh
```

---

## Core Data Structures (in `include/SMR_data.hpp`)

All major structs use **Struct-of-Arrays (SoA)** layout — separate vectors per field. Do not convert to Array-of-Structs; this pattern exists for cache performance and SIMD.

| Struct | Purpose |
|--------|---------|
| `gwasData` | GWAS summary stats: SNP names, bp positions, alleles, freq, beta, SE, p-values |
| `eqtlInfo` | eQTL data: probes, SNPs, effect sizes; supports sparse and dense BESD formats |
| `bInfo` | Binary genotype data from PLINK bed/bim/fam; members prefixed with `_` (e.g., `_chr`, `_snp_name`) |
| `ldInfo` | Linkage disequilibrium information |
| `SMRWK` | Working buffer for one SMR analysis unit: effect sizes, freqs, Z-scores |
| `SMRRLT` | Result record: probe, top SNP, effect sizes, SMR p-value, HEIDI p-value |

---

## Key Namespaces and Coding Conventions

- `SMRDATA` — primary namespace for data structures and analysis functions
- `StatFunc` — statistical test functions
- `CommFunc` — common utilities (file ops, vector math)
- `StrFunc` — string parsing

**Naming rules**:
- Struct/class member variables: underscore prefix (`_chr`, `_bp`, `_freq`)
- Functions: camelCase (`read_bimfile()`, `allele_check()`, `heidi_test_new()`)
- Global variables: lowercase with underscores (`thread_num`, `mute`, `MAX_NUM_LD`)

**Platform guards** appear throughout — always respect them:
```cpp
#ifdef _WIN64   // Windows 64-bit
#ifdef __APPLE__ // macOS
#ifdef __linux__
```

---

## BESD File Format (Critical — Read Before Touching I/O)

BESD (Binary Effect Size Data) is SMR's native eQTL binary format. Two variants:

| Format | Constant | Description |
|--------|----------|-------------|
| Sparse | `SPARSE_FILE_TYPE_3` | Stores only significant SNP–probe pairs; memory-efficient |
| Dense | `DENSE_FILE_TYPE_3` | Full matrix; faster random access |

**Files per dataset**: `.besd` (binary data) + `.esi` (SNP info) + `.epi` (probe info)

Key functions:
- `read_besd()` / `read_besd_sparse()` — data loading
- `make_besd()` — BESD creation from text eQTL summary
- `query_besd()` — region-based lookup

**Do not change the binary layout constants** (`SPARSE_FILE_TYPE_3`, `DENSE_FILE_TYPE_3`). Existing user data files depend on them.

---

## Statistical Methods

### SMR Test
Tests if GWAS-top SNP effect is consistent with being mediated through expression:
- Combines GWAS β and eQTL β via Wald ratio
- Produces SMR chi-square test statistic and p-value

### HEIDI Test
Distinguishes pleiotropy from linkage by testing heterogeneity across multiple SNPs:
- `heidi_test()` — standard HEIDI
- `heidi_test_new()` — updated method
- `heidi_test_ref_new()` — reference-panel-based
- Controlled by `--heidi-off` (skip entirely)

### LD Calculation
- `cor_calc()` — pairwise LD calculation
- `ld_calc_o2m()` — one-to-many LD
- `get_square_ldpruning_idxes()` — LD pruning
- `MAX_NUM_LD` — maximum number of LD SNPs (tunable via `--max-num-ld`)

### CDF / P-value Functions (`dcdflib.cpp`, `StatFunc.cpp`)
Do not reimplement or replace these — the 9,205-line `dcdflib.cpp` is a proven numerical library. All p-value calculations must go through `StatFunc` wrappers.

---

## Data Flow: GWAS + eQTL → SMR Result

```
CLI (SMR.cpp)
  └─ parse options & dispatch
        ├─ read_gwas_summary()      → gwasData
        ├─ read_beqtl_summary()     → eqtlInfo
        ├─ read_bimfile() / bInfo   → reference LD panel (optional)
        │
        ├─ allele_check()           → harmonize alleles between datasets
        ├─ filter_by_maf/chr/snp    → prune SNPs
        │
        ├─ smr_test()               → fills SMRWK working buffer
        │     └─ heidi_test_*()     → heterogeneity p-value
        │
        └─ output_smr_result()      → .smr text file (SMRRLT records)
```

---

## Common Pitfalls

1. **SoA layout assumption** — All vector accesses across `gwasData` / `bInfo` / `eqtlInfo` are indexed in sync. Adding a new field means adding a parallel vector *and* keeping all resize/filter operations in sync.

2. **BESD format constants** — Never change `SPARSE_FILE_TYPE_3` / `DENSE_FILE_TYPE_3`. Changing these breaks backward compatibility with user data.

3. **Platform-specific `ftell` overflow** — Windows 64-bit has a known `ftell` overflow for large BESD files; the fix uses `_ftelli64`. Don't simplify this away.

4. **Thread safety with OpenMP** — Global variables like `thread_num` and `mute` affect OpenMP loops. Be cautious modifying shared state in parallel regions.

5. **Allele strand ambiguity** — A/T and C/G SNPs require special handling in `allele_check()`. This logic is subtle — read it fully before touching.

6. **`dcdflib.cpp` is not to be modified** — It is a third-party numerical library. Route new statistics through `StatFunc` wrappers.

7. **No test suite exists** — Changes must be validated manually with known input/output datasets. Be conservative with refactors.

---

## Major CLI Options Reference

```
Input:
  --bfile            PLINK binary genotype (LD reference)
  --gwas-summary     GWAS summary statistics file
  --beqtl-summary    Binary eQTL summary (BESD format)
  --eqtl-summary     Text eQTL summary file

Analysis:
  --smr              Run SMR analysis
  --heidi-off        Skip HEIDI test
  --trans            Trans-eQTL mode
  --plot             Generate SMR locus plots
  --make-besd        Create BESD from text eQTL file

Filtering:
  --maf              Minor allele frequency cutoff
  --ld-upper-limit   LD r² upper threshold (HEIDI)
  --ld-lower-limit   LD r² lower threshold (HEIDI)
  --max-num-ld       Max LD SNPs for HEIDI (MAX_NUM_LD)
  --extract-snp / --exclude-snp
  --extract-probe / --exclude-probe
  --chr / --snp-chr / --probe-chr

Output:
  --out              Output file prefix
  --thread-num       OpenMP thread count
```

---

## Adding New Features — Checklist

- [ ] Add CLI flag in `include/SMR.hpp` (option struct)
- [ ] Parse flag in `src/SMR.cpp` (option parsing block)
- [ ] Implement logic in the appropriate `SMR_data*.cpp` module
- [ ] Declare in matching header under `include/`
- [ ] If new data field: add parallel vector to the relevant struct AND update all filter/resize call sites
- [ ] If new output column: update both the writer and the result struct `SMRRLT`
- [ ] Test with real GWAS + BESD data; no automated tests exist

---

## Resources

- Official documentation: https://yanglab.westlake.edu.cn/software/smr/
- BESD format spec: documented in SMR website tutorials
- Changelog: `docs/changlog.txt`
