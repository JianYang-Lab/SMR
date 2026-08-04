# SMR

The SMR software tool was originally developed to implement the SMR & HEIDI methods to test for pleiotropic association between the expression level of a gene and a complex trait of interest using summary-level data from GWAS and expression quantitative trait loci (eQTL) studies ([Zhu et al. 2016 Nature Genetics](https://www.nature.com/articles/ng.3538/)). The SMR & HEIDI methodology can be interpreted as an analysis to test if the effect size of a SNP on the phenotype is mediated by gene expression. This tool can therefore be used to prioritize genes underlying GWAS hits for follow-up functional studies. The methods are applicable to all kinds of molecular QTL (xQTL) data, including DNA methylation QTL (mQTL) and protein abundance QTL (pQTL).


## Requirement

* Eigen

* libz

## Installation

### Install requirments
Dependencies (Eigen, zlib, spdlog, fmt, mimalloc) are fetched automatically by CMake — no manual installation needed. A C++17 compiler, CMake ≥ 3.16, and OpenMP are required.

### Build

#### using make

Simply, in smr directory (delegates to `scripts/local/build.sh`):

```
make
```

#### using cmake

In SMR directory,

```shell
mkdir build
cd build
cmake ..
make
```

The default build type is Release. To build with debug information,

```shell
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

Useful CMake options:

- `-DSMR_BLAS_BACKEND=mkl|openblas|none` (default `mkl`): BLAS backend for Eigen. `mkl` uses Intel oneMKL (`-DMKL_DIR=...`/`-DMKL_ROOT=...`); `openblas` uses OpenBLAS for non-Intel systems (AMD etc.) — a system install is used when found (hint with `-DOPENBLAS_ROOT=...`), otherwise OpenBLAS v0.3.30 is fetched and built from source (statically linked, `DYNAMIC_ARCH=ON` for runtime CPU dispatch); `none` uses Eigen's native backend.
- `-DBUILD_WITH_MKL=ON/OFF` (default ON): accelerate Eigen with Intel oneMKL; falls back to Eigen's native backend if MKL is not found. Use `-DMKL_DIR=...` or `-DMKL_ROOT=...` to point at your oneMKL installation.
- `-DBUILD_WITH_MIMALLOC=ON/OFF` (default ON): link mimalloc for faster memory allocation.
- `-DSMR_VERBOSE_CONFIGURE=ON`: print extended CMake configuration diagnostics.

Helper build scripts: `scripts/local/build.sh [-f|--fresh] [-g|--generate]` (Linux/macOS) and `scripts/hpc/build.sh` (HPC clusters). The scripts accept the `SMR_BLAS_BACKEND` environment variable (default `mkl`), e.g. `SMR_BLAS_BACKEND=openblas ./scripts/local/build.sh -f -g`.


## Usage
visit https://yanglab.westlake.edu.cn/software/smr/ for software's document.

## Using SMR with AI Tools

The SMR release package includes a built-in MCP (Model Context Protocol) server, enabling you to run SMR analyses through AI tools such as **Claude**, **Codex**, and **OpenCode** via natural language conversation.

### Quick Start

1. Download and extract the SMR release package.
2. Configure the MCP server for your AI tool (see `CLAUDE_CODE_USAGE.md` in the release package for details).
3. Start chatting with the AI — for example: *"Run an SMR analysis with bfile data/1kg_eur, gwas-summary data/bmi.ma, beqtl-summary data/eqtl.besd, output to results/smr_out"*.

### Available MCP Tools

| Tool | Function |
|------|----------|
| `run_smr_analysis` | Run the main SMR test (Wald ratio + HEIDI) |
| `smr_multi` | Run set-based (multi-SNP) SMR analysis |
| `make_besd` | Create a sparse BESD file from a text eQTL summary |
| `make_besd_dense` | Create a dense BESD file from a text eQTL summary |
| `query_besd` | Query a BESD file for significant SNP-probe associations |
| `show_sample_size` | Display the sample size stored in a BESD file |
| `recode_besd` | Convert a BESD file to COJO/SMR text format |
| `plot_smr` | Generate SMR locus plots |
| `make_bld` | Create a BLD (binary LD) file from a PLINK bfile |
| `update_freq` | Update allele frequencies in a BESD file |
| `meta_analysis` | Run meta-analysis of multiple eQTL studies (MeCS) |
| `combine_besd` | Combine multiple BESD files into one |
| `count_cis` | Count cis-eQTL in a BESD file |
| `count_trans` | Count trans-eQTL in a BESD file |
| `update_epi_esi` | Update EPI/ESI annotation files in a BESD file |
| `get_version` | Get SMR version information |
| `run_raw_command` | Pass arbitrary SMR command-line arguments (escape hatch) |

For full configuration instructions for Claude Desktop, Claude Code, OpenCode, and Codex, see **`CLAUDE_CODE_USAGE.md`** in the release package.
