# SMR

The SMR software tool was originally developed to implement the SMR & HEIDI methods to test for pleiotropic association between the expression level of a gene and a complex trait of interest using summary-level data from GWAS and expression quantitative trait loci (eQTL) studies ([Zhu et al. 2016 Nature Genetics](https://www.nature.com/articles/ng.3538/)). The SMR & HEIDI methodology can be interpreted as an analysis to test if the effect size of a SNP on the phenotype is mediated by gene expression. This tool can therefore be used to prioritize genes underlying GWAS hits for follow-up functional studies. The methods are applicable to all kinds of molecular QTL (xQTL) data, including DNA methylation QTL (mQTL) and protein abundance QTL (pQTL).


## Requirement

* Eigen

* libz

## Installation

### Install requirments
Install Eigen and libz. If you want build smr staticlly, you need static library of those them. And the static version of omp is needed too.

### Build

#### using make

Simply, in smr directoy,

```
make
```

If you want compile static version, type:

```
make smr_static
```

If your Eigen or libz is not installed in standard head file and library searching path, you can
specific them as following:

```
make EIGEN_PATH="where the Eigen's head file located" ZLIB_INCLUDE="path of zlib head file" ZLIB_LIB="path of zlib library"
```

To build smr by debug mode, using:

```
make DEBUG=ON
```

#### using cmake

In SMR directory,

```shell
mkdir build
cd build
cmake ..
make
```

To turn on debug,

```shell
cmake -DCMAKE_BUILD_TYPE=DEBUG ..
```

To build static version of SMR,

```shell
cmake -DBUILD_STATIC=ON ..
```

If you need specify Eigne path you can use `-Deigen_path="path_to_Eigen"` when run `cmake ..`. To specify zlib head and library files, using `-Dzlib_path=path_to_your_zlib`. If your zlib head file and library located in different directory, you can use `-Dzlib_include_path=where_zlib_head_files` and `-Dzlib_lib_path=where_zlib_lib` to specify them.


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
