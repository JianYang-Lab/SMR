#!/usr/bin/env python3
"""
SMR MCP Server

Exposes SMR (Summary-data-based Mendelian Randomization) CLI functionality
as MCP (Model Context Protocol) tools, enabling conversational invocation
through Claude or other LLM clients.

Usage:
    python3 smr_mcp_server.py

Environment variables:
    SMR_BIN          Path to the smr binary (default: auto-detect)
    SMR_MCP_TIMEOUT  Command timeout in seconds (default: 3600)
    SMR_MCP_CWD      Working directory for smr execution (default: inherited)
"""

import json
import os
import shutil
import subprocess
import sys
import traceback

# ─── Server Configuration ────────────────────────────────────────────────────

SERVER_NAME = "smr-mcp-server"
SERVER_VERSION = "1.0.0"
PROTOCOL_VERSION = "2025-06-18"

DEFAULT_TIMEOUT = int(os.environ.get("SMR_MCP_TIMEOUT", "3600"))
DEFAULT_CWD = os.environ.get("SMR_MCP_CWD")

# Maximum output size before truncation (in characters)
MAX_OUTPUT = 50000

# ─── SMR Binary Discovery ────────────────────────────────────────────────────


def find_smr_binary() -> str:
    """Locate the smr binary by checking env var, PATH, and common locations.

    Search order:
      1. SMR_BIN environment variable
      2. PATH lookup
      3. Relative to this script — covers both source-tree and release-package layouts:

         Source tree:
           SMR/scripts/mcp/smr_mcp_server.py
           SMR/build/smr                  ← candidate

         Release package:
           smr-1.4.2-linux-x86_64/mcp/smr_mcp_server.py
           smr-1.4.2-linux-x86_64/smr     ← candidate (one level up)
    """
    # 1. SMR_BIN environment variable
    env_bin = os.environ.get("SMR_BIN")
    if env_bin:
        if os.path.isfile(env_bin):
            return os.path.abspath(env_bin)
        return env_bin  # Return as-is; execution will give a clear error

    # 2. PATH lookup
    path_bin = shutil.which("smr")
    if path_bin:
        return path_bin

    # 3. Common locations relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    candidates = [
        # Release package: smr sits one level above mcp/
        os.path.join(script_dir, "..", "smr"),
        # Source tree: smr built in build/ two levels above mcp/
        os.path.join(script_dir, "..", "..", "build", "smr"),
        os.path.join(script_dir, "..", "..", "smr"),
        # CWD fallbacks
        "./build/smr",
        "./smr",
    ]

    for candidate in candidates:
        if os.path.isfile(candidate):
            return os.path.abspath(candidate)

    return "smr"  # Fallback — execution will fail with a helpful message


SMR_BIN = find_smr_binary()


# ─── SMR Command Execution ───────────────────────────────────────────────────


def execute_smr(args: list, timeout: int = DEFAULT_TIMEOUT) -> dict:
    """
    Execute the smr binary with given arguments.

    Returns a dict with:
        - success: bool
        - command: str (the full command executed)
        - stdout: str
        - stderr: str
        - returncode: int
        - output_files: list of str (files that exist after execution)
    """
    cmd = [SMR_BIN] + args
    result = {
        "success": False,
        "command": " ".join(cmd),
        "stdout": "",
        "stderr": "",
        "returncode": -1,
        "output_files": [],
    }

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=DEFAULT_CWD,
        )
        result["stdout"] = proc.stdout or ""
        result["stderr"] = proc.stderr or ""
        result["returncode"] = proc.returncode
        result["success"] = proc.returncode == 0

    except subprocess.TimeoutExpired as e:
        result["stderr"] = f"Command timed out after {timeout} seconds."
        if e.stdout:
            result["stdout"] = e.stdout.decode() if isinstance(e.stdout, bytes) else e.stdout
        if e.stderr:
            result["stderr"] += "\n" + (
                e.stderr.decode() if isinstance(e.stderr, bytes) else e.stderr
            )
    except FileNotFoundError:
        result["stderr"] = (
            f"ERROR: SMR binary not found at '{SMR_BIN}'.\n\n"
            "To fix this:\n"
            "  1. Build SMR first:  make build   (or:  mkdir build && cd build && cmake .. && make)\n"
            "  2. Or set the SMR_BIN environment variable to the full path of the smr binary.\n"
        )
    except Exception as e:
        result["stderr"] = f"Unexpected error executing SMR: {e}"

    # Detect output files created by --out
    out_prefix = None
    for i, arg in enumerate(args):
        if arg == "--out" and i + 1 < len(args):
            out_prefix = args[i + 1]
            break

    if out_prefix:
        for ext in [".smr", ".plot", ".bld", ".besd", ".esi", ".epi", ".log"]:
            candidate = out_prefix + ext
            if os.path.isfile(candidate):
                result["output_files"].append(candidate)

    return result


def format_result(result: dict) -> str:
    """Format an SMR execution result as a human-readable string."""
    parts = []

    parts.append(f"Command: {result['command']}")
    parts.append(f"Exit code: {result['returncode']}")

    if result["success"]:
        parts.append("Status: SUCCESS")
    else:
        parts.append("Status: FAILED")

    # Combine stdout and stderr, truncating if too long
    combined = ""
    if result["stdout"]:
        combined += result["stdout"]
    if result["stderr"]:
        if combined:
            combined += "\n"
        combined += result["stderr"]

    if len(combined) > MAX_OUTPUT:
        truncated = combined[:MAX_OUTPUT]
        remaining = len(combined) - MAX_OUTPUT
        combined = f"{truncated}\n\n... [output truncated, {remaining} more characters] ..."

    if combined:
        parts.append(f"\n--- Output ---\n{combined}")

    if result["output_files"]:
        parts.append("\n--- Output Files ---")
        for f in result["output_files"]:
            size = os.path.getsize(f) if os.path.isfile(f) else 0
            parts.append(f"  {f}  ({size} bytes)")

    return "\n".join(parts)


# ─── Tool Implementations ────────────────────────────────────────────────────
# Each function corresponds to one MCP tool.  The function receives the
# tool's arguments as keyword arguments and returns a string result.


def tool_run_smr_analysis(
    bfile: str,
    gwas_summary: str,
    beqtl_summary: str,
    out: str,
    maf: float = 0.01,
    thread_num: int = 1,
    heidi_off: bool = False,
    ld_upper_limit: float = 0.9,
    ld_lower_limit: float = 0.05,
    max_num_ld: int = 500,
    trans: bool = False,
    diff_freq: float = 0.2,
    diff_freq_prop: float = 0.05,
) -> str:
    """Run the main SMR analysis.

    Tests whether a SNP's effect on a phenotype is mediated through gene
    expression, using GWAS summary statistics and an eQTL reference panel.
    Output is written to ``<out>.smr``.
    """
    args = [
        "--bfile", bfile,
        "--gwas-summary", gwas_summary,
        "--beqtl-summary", beqtl_summary,
        "--out", out,
        "--maf", str(maf),
        "--thread-num", str(thread_num),
        "--ld-upper-limit", str(ld_upper_limit),
        "--ld-lower-limit", str(ld_lower_limit),
        "--max_num_ld", str(max_num_ld),
        "--diff-freq", str(diff_freq),
        "--diff-freq-prop", str(diff_freq_prop),
        "--smr",
    ]
    if heidi_off:
        args.append("--heidi-off")
    if trans:
        args.append("--trans")
    return format_result(execute_smr(args))


def tool_make_besd(
    eqtl_summary: str,
    out: str,
    cis_wind: int = 2000,
    trans_wind: int = 1000,
    peqtl_trans: float = 5e-8,
    peqtl_cis: float = 5e-8,
    peqtl_other: float = 1e-5,
    thread_num: int = 1,
) -> str:
    """Create a sparse BESD file from a text eQTL summary.

    The BESD (Binary Effect Size Data) format stores eQTL summary
    statistics efficiently.  Only significant SNP-probe pairs are kept
    in the sparse variant.
    """
    args = [
        "--eqtl-summary", eqtl_summary,
        "--out", out,
        "--cis-wind", str(cis_wind),
        "--trans-wind", str(trans_wind),
        "--peqtl-trans", str(peqtl_trans),
        "--peqtl-cis", str(peqtl_cis),
        "--peqtl-other", str(peqtl_other),
        "--thread-num", str(thread_num),
        "--make-besd",
    ]
    return format_result(execute_smr(args))


def tool_make_besd_dense(
    eqtl_summary: str,
    out: str,
    cis_wind: int = 2000,
    trans_wind: int = 1000,
    peqtl_trans: float = 5e-8,
    peqtl_cis: float = 5e-8,
    peqtl_other: float = 1e-5,
    thread_num: int = 1,
) -> str:
    """Create a dense BESD file from a text eQTL summary.

    The dense variant stores the full SNP-probe effect-size matrix,
    trading memory for faster random access.
    """
    args = [
        "--eqtl-summary", eqtl_summary,
        "--out", out,
        "--cis-wind", str(cis_wind),
        "--trans-wind", str(trans_wind),
        "--peqtl-trans", str(peqtl_trans),
        "--peqtl-cis", str(peqtl_cis),
        "--peqtl-other", str(peqtl_other),
        "--thread-num", str(thread_num),
        "--make-besd-dense",
    ]
    return format_result(execute_smr(args))


def tool_plot_smr(
    bfile: str,
    gwas_summary: str,
    beqtl_summary: str,
    out: str,
    maf: float = 0.01,
    thread_num: int = 1,
    heidi_off: bool = False,
    ld_upper_limit: float = 0.9,
    ld_lower_limit: float = 0.05,
    max_num_ld: int = 500,
) -> str:
    """Generate SMR locus plots.

    Creates publication-quality plots showing the SMR association,
    eQTL effect sizes, and HEIDI test results for each significant locus.
    """
    args = [
        "--bfile", bfile,
        "--gwas-summary", gwas_summary,
        "--beqtl-summary", beqtl_summary,
        "--out", out,
        "--maf", str(maf),
        "--thread-num", str(thread_num),
        "--ld-upper-limit", str(ld_upper_limit),
        "--ld-lower-limit", str(ld_lower_limit),
        "--max-num-ld", str(max_num_ld),
        "--plot",
    ]
    if heidi_off:
        args.append("--heidi-off")
    return format_result(execute_smr(args))


def tool_query_besd(
    beqtl_summary: str,
    out: str,
    p_value: float = 5e-8,
    maf: float = 0.01,
    thread_num: int = 1,
) -> str:
    """Query a BESD file for significant SNP-probe associations.

    Returns all SNP-probe pairs with eQTL p-value below the specified
    threshold.  Useful for exploring the contents of a BESD dataset.
    """
    args = [
        "--beqtl-summary", beqtl_summary,
        "--out", out,
        "--maf", str(maf),
        "--thread-num", str(thread_num),
        "--query", str(p_value),
    ]
    return format_result(execute_smr(args))


def tool_show_sample_size(
    beqtl_summary: str,
) -> str:
    """Display the sample size stored in a BESD file."""
    args = [
        "--beqtl-summary", beqtl_summary,
        "--show-n",
    ]
    return format_result(execute_smr(args))


def tool_recode_besd(
    beqtl_summary: str,
    out: str,
    maf: float = 0.01,
    thread_num: int = 1,
) -> str:
    """Convert a BESD file to COJO / SMR text format.

    The output text file can be used with GCTA-COJO or re-imported
    as a text eQTL summary.
    """
    args = [
        "--beqtl-summary", beqtl_summary,
        "--out", out,
        "--maf", str(maf),
        "--thread-num", str(thread_num),
        "--recode",
    ]
    return format_result(execute_smr(args))


def tool_smr_multi(
    bfile: str,
    gwas_summary: str,
    beqtl_summary: str,
    out: str,
    maf: float = 0.01,
    thread_num: int = 1,
    heidi_off: bool = False,
    ld_upper_limit: float = 0.9,
    ld_lower_limit: float = 0.05,
    max_num_ld: int = 500,
) -> str:
    """Run set-based (multi-SNP) SMR analysis.

    Performs SMR using multiple SNPs as instruments simultaneously,
    providing increased statistical power when LD structure is complex.
    """
    args = [
        "--bfile", bfile,
        "--gwas-summary", gwas_summary,
        "--beqtl-summary", beqtl_summary,
        "--out", out,
        "--maf", str(maf),
        "--thread-num", str(thread_num),
        "--ld-upper-limit", str(ld_upper_limit),
        "--ld-lower-limit", str(ld_lower_limit),
        "--max_num_ld", str(max_num_ld),
        "--smr-multi",
    ]
    if heidi_off:
        args.append("--heidi-off")
    return format_result(execute_smr(args))


def tool_make_bld(
    bfile: str,
    out: str,
    maf: float = 0.01,
    thread_num: int = 1,
) -> str:
    """Create a BLD (binary LD) file from a PLINK bfile.

    The BLD format stores LD information in a compact binary form,
    used as a reference LD panel for SMR/HEIDI analyses.
    """
    args = [
        "--bfile", bfile,
        "--out", out,
        "--maf", str(maf),
        "--thread-num", str(thread_num),
        "--make-bld",
    ]
    return format_result(execute_smr(args))


def tool_update_freq(
    beqtl_summary: str,
    out: str,
    freq_file: str,
    thread_num: int = 1,
) -> str:
    """Update allele frequencies in a BESD file.

    Replaces the stored allele frequencies with those from an external
    frequency file (e.g., from a reference panel).
    """
    args = [
        "--beqtl-summary", beqtl_summary,
        "--out", out,
        "--thread-num", str(thread_num),
        "--update-freq", freq_file,
    ]
    return format_result(execute_smr(args))


def tool_meta_analysis(
    besd_flist: str,
    out: str,
    method: str = "mecs",
    pmecs: float = 0.01,
    nmecs: int = 100,
    thread_num: int = 1,
) -> str:
    """Run meta-analysis of multiple eQTL studies.

    Combines effect sizes across multiple BESD files using either
    standard meta-analysis or MeCS (Mendelian randomization with
    Correlated SNPs) method.
    """
    args = [
        "--besd-flist", besd_flist,
        "--out", out,
        "--thread-num", str(thread_num),
    ]
    if method == "mecs":
        args.extend(["--mecs", "--pmecs", str(pmecs), "--nmecs", str(nmecs)])
    else:
        args.append("--meta")
    return format_result(execute_smr(args))


def tool_combine_besd(
    besd_flist: str,
    out: str,
    thread_num: int = 1,
) -> str:
    """Combine multiple BESD files into one.

    Merges effect-size data from multiple BESD files listed in the
    flist file into a single output BESD file.
    """
    args = [
        "--besd-flist", besd_flist,
        "--out", out,
        "--thread-num", str(thread_num),
    ]
    return format_result(execute_smr(args))


def tool_count_cis(
    beqtl_summary: str,
    out: str,
    p_threshold: float = 5e-8,
    cis_window: int = 2000,
    thread_num: int = 1,
) -> str:
    """Count cis-eQTL in a BESD file.

    Outputs the top SNP in each cis-region, where cis is defined
    as within cis_window Kb of the probe.
    """
    args = [
        "--beqtl-summary", beqtl_summary,
        "--out", out,
        "--cis-wind", str(cis_window),
        "--peqtl-cis", str(p_threshold),
        "--thread-num", str(thread_num),
        "--descriptive-cis",
    ]
    return format_result(execute_smr(args))


def tool_count_trans(
    beqtl_summary: str,
    out: str,
    p_threshold: float = 5e-8,
    cis_window: int = 2000,
    trans_window: int = 1000,
    thread_num: int = 1,
) -> str:
    """Count trans-eQTL in a BESD file.

    Outputs the top SNP in each trans-region, where trans is defined
    as outside the cis_window Kb region around the probe.
    """
    args = [
        "--beqtl-summary", beqtl_summary,
        "--out", out,
        "--cis-wind", str(cis_window),
        "--trans-wind", str(trans_window),
        "--peqtl-trans", str(p_threshold),
        "--thread-num", str(thread_num),
        "--descriptive-trans",
    ]
    return format_result(execute_smr(args))


def tool_update_epi_esi(
    beqtl_summary: str,
    out: str,
    epi_file: str = "",
    esi_file: str = "",
    thread_num: int = 1,
) -> str:
    """Update EPI (probe info) and/or ESI (SNP info) files in a BESD file.

    Replaces the stored probe or SNP annotation with external files.
    At least one of epi_file or esi_file must be provided.
    """
    args = [
        "--beqtl-summary", beqtl_summary,
        "--out", out,
        "--thread-num", str(thread_num),
    ]
    if epi_file:
        args.extend(["--update-epi", epi_file])
    if esi_file:
        args.extend(["--update-esi", esi_file])
    return format_result(execute_smr(args))


def tool_get_version() -> str:
    """Get SMR version information."""
    result = execute_smr([], timeout=30)
    # SMR exits with code 1 when run with no arguments (prints usage and exits),
    # but the version banner is still in stdout. Accept exit codes 0 and 1.
    if result["returncode"] in (0, 1):
        return result["stdout"]
    return format_result(result)


def tool_run_raw_command(args: list) -> str:
    """Run SMR with arbitrary command-line arguments.

    This is an escape hatch for advanced usage not covered by the
    typed tool wrappers.  Pass the arguments exactly as they would
    appear on the command line (without the ``smr`` prefix).

    Example::

        run_raw_command(["--bfile", "ref", "--gwas-summary", "g.ma",
                         "--beqtl-summary", "e.besd", "--out", "o",
                         "--smr", "--heidi-off", "--trans"])
    """
    return format_result(execute_smr(args))


# ─── Tool Registry ───────────────────────────────────────────────────────────

TOOLS = [
    {
        "name": "run_smr_analysis",
        "description": (
            "Run SMR (Summary-data-based Mendelian Randomization) analysis. "
            "Tests whether a SNP's effect on a phenotype is mediated through "
            "gene expression using GWAS and eQTL summary-level data. "
            "Produces an .smr output file with SMR test statistics and HEIDI p-values.\n\n"
            "Example: run_smr_analysis(bfile='data/ld_ref', "
            "gwas_summary='data/gwas.ma', beqtl_summary='data/eqtl.besd', "
            "out='results/smr_output')"
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "bfile": {
                    "type": "string",
                    "description": "PLINK binary genotype file prefix for LD reference "
                    "(e.g., 'data/1kg_eur'). The tool will look for .bed, .bim, .fam files.",
                },
                "gwas_summary": {
                    "type": "string",
                    "description": "GWAS summary statistics file in MaCH/Minimac format "
                    "(columns: SNP A1 A2 freq beta se p N).",
                },
                "beqtl_summary": {
                    "type": "string",
                    "description": "Binary eQTL summary file (BESD format).",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix. Results will be written to <out>.smr.",
                },
                "maf": {
                    "type": "number",
                    "description": "Minor allele frequency cutoff (0-0.5).",
                    "default": 0.01,
                    "minimum": 0,
                    "maximum": 0.5,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
                "heidi_off": {
                    "type": "boolean",
                    "description": "Skip the HEIDI heterogeneity test.",
                    "default": False,
                },
                "ld_upper_limit": {
                    "type": "number",
                    "description": "LD r-squared upper threshold for HEIDI.",
                    "default": 0.9,
                    "minimum": 0,
                    "maximum": 1,
                },
                "ld_lower_limit": {
                    "type": "number",
                    "description": "LD r-squared lower threshold for HEIDI.",
                    "default": 0.05,
                    "minimum": 0,
                    "maximum": 1,
                },
                "max_num_ld": {
                    "type": "integer",
                    "description": "Maximum number of LD SNPs for HEIDI (500-10000).",
                    "default": 500,
                    "minimum": 500,
                    "maximum": 10000,
                },
                "trans": {
                    "type": "boolean",
                    "description": "Run in trans-eQTL mode.",
                    "default": False,
                },
                "diff_freq": {
                    "type": "number",
                    "description": "Allele frequency difference threshold for QC.",
                    "default": 0.2,
                    "minimum": 0,
                    "maximum": 1,
                },
                "diff_freq_prop": {
                    "type": "number",
                    "description": "Proportion of SNPs allowed to exceed diff_freq.",
                    "default": 0.05,
                    "minimum": 0,
                    "maximum": 1,
                },
            },
            "required": ["bfile", "gwas_summary", "beqtl_summary", "out"],
        },
    },
    {
        "name": "make_besd",
        "description": (
            "Create a sparse BESD file from a text eQTL summary. "
            "BESD (Binary Effect Size Data) is SMR's native binary format "
            "for storing eQTL summary statistics efficiently.\n\n"
            "Example: make_besd(eqtl_summary='data/eqtl.txt', out='data/eqtl.besd')"
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "eqtl_summary": {
                    "type": "string",
                    "description": "Text eQTL summary file path.",
                },
                "out": {
                    "type": "string",
                    "description": "Output BESD file prefix.",
                },
                "cis_wind": {
                    "type": "integer",
                    "description": "Cis window size in Kb.",
                    "default": 2000,
                    "minimum": 1,
                },
                "trans_wind": {
                    "type": "integer",
                    "description": "Trans window size in Kb.",
                    "default": 1000,
                    "minimum": 0,
                },
                "peqtl_trans": {
                    "type": "number",
                    "description": "P-value threshold for trans-eQTL.",
                    "default": 5e-8,
                    "minimum": 0,
                    "maximum": 1,
                },
                "peqtl_cis": {
                    "type": "number",
                    "description": "P-value threshold for cis-eQTL.",
                    "default": 5e-8,
                    "minimum": 0,
                    "maximum": 1,
                },
                "peqtl_other": {
                    "type": "number",
                    "description": "P-value threshold for other eQTL.",
                    "default": 1e-5,
                    "minimum": 0,
                    "maximum": 1,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["eqtl_summary", "out"],
        },
    },
    {
        "name": "make_besd_dense",
        "description": (
            "Create a dense BESD file from a text eQTL summary. "
            "The dense variant stores the full SNP-probe effect-size matrix, "
            "trading memory for faster random access."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "eqtl_summary": {
                    "type": "string",
                    "description": "Text eQTL summary file path.",
                },
                "out": {
                    "type": "string",
                    "description": "Output BESD file prefix.",
                },
                "cis_wind": {
                    "type": "integer",
                    "description": "Cis window size in Kb.",
                    "default": 2000,
                    "minimum": 1,
                },
                "trans_wind": {
                    "type": "integer",
                    "description": "Trans window size in Kb.",
                    "default": 1000,
                    "minimum": 0,
                },
                "peqtl_trans": {
                    "type": "number",
                    "description": "P-value threshold for trans-eQTL.",
                    "default": 5e-8,
                    "minimum": 0,
                    "maximum": 1,
                },
                "peqtl_cis": {
                    "type": "number",
                    "description": "P-value threshold for cis-eQTL.",
                    "default": 5e-8,
                    "minimum": 0,
                    "maximum": 1,
                },
                "peqtl_other": {
                    "type": "number",
                    "description": "P-value threshold for other eQTL.",
                    "default": 1e-5,
                    "minimum": 0,
                    "maximum": 1,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["eqtl_summary", "out"],
        },
    },
    {
        "name": "plot_smr",
        "description": (
            "Generate SMR locus plots. Creates publication-quality plots "
            "showing the SMR association, eQTL effect sizes, and HEIDI test "
            "results for each significant locus.\n\n"
            "Example: plot_smr(bfile='data/ld_ref', gwas_summary='data/gwas.ma', "
            "beqtl_summary='data/eqtl.besd', out='results/plot')"
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "bfile": {
                    "type": "string",
                    "description": "PLINK binary genotype file prefix for LD reference.",
                },
                "gwas_summary": {
                    "type": "string",
                    "description": "GWAS summary statistics file.",
                },
                "beqtl_summary": {
                    "type": "string",
                    "description": "Binary eQTL summary file (BESD format).",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix for plots.",
                },
                "maf": {
                    "type": "number",
                    "description": "Minor allele frequency cutoff.",
                    "default": 0.01,
                    "minimum": 0,
                    "maximum": 0.5,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
                "heidi_off": {
                    "type": "boolean",
                    "description": "Skip HEIDI test.",
                    "default": False,
                },
                "ld_upper_limit": {
                    "type": "number",
                    "description": "LD r-squared upper threshold.",
                    "default": 0.9,
                    "minimum": 0,
                    "maximum": 1,
                },
                "ld_lower_limit": {
                    "type": "number",
                    "description": "LD r-squared lower threshold.",
                    "default": 0.05,
                    "minimum": 0,
                    "maximum": 1,
                },
                "max_num_ld": {
                    "type": "integer",
                    "description": "Maximum LD SNPs for HEIDI.",
                    "default": 500,
                    "minimum": 500,
                    "maximum": 10000,
                },
            },
            "required": ["bfile", "gwas_summary", "beqtl_summary", "out"],
        },
    },
    {
        "name": "query_besd",
        "description": (
            "Query a BESD file for significant SNP-probe associations. "
            "Returns all SNP-probe pairs with eQTL p-value below the "
            "specified threshold."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "beqtl_summary": {
                    "type": "string",
                    "description": "Binary eQTL summary file (BESD format).",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix.",
                },
                "p_value": {
                    "type": "number",
                    "description": "P-value threshold for query (default: 5e-8).",
                    "default": 5e-8,
                    "minimum": 0,
                    "maximum": 1,
                },
                "maf": {
                    "type": "number",
                    "description": "Minor allele frequency cutoff.",
                    "default": 0.01,
                    "minimum": 0,
                    "maximum": 0.5,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["beqtl_summary", "out"],
        },
    },
    {
        "name": "show_sample_size",
        "description": (
            "Display the sample size stored in a BESD file. "
            "Useful for verifying data integrity before analysis."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "beqtl_summary": {
                    "type": "string",
                    "description": "Binary eQTL summary file (BESD format).",
                },
            },
            "required": ["beqtl_summary"],
        },
    },
    {
        "name": "recode_besd",
        "description": (
            "Convert a BESD file to COJO/SMR text format. "
            "The output text file can be used with GCTA-COJO or "
            "re-imported as a text eQTL summary."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "beqtl_summary": {
                    "type": "string",
                    "description": "Binary eQTL summary file (BESD format).",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix.",
                },
                "maf": {
                    "type": "number",
                    "description": "Minor allele frequency cutoff.",
                    "default": 0.01,
                    "minimum": 0,
                    "maximum": 0.5,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["beqtl_summary", "out"],
        },
    },
    {
        "name": "smr_multi",
        "description": (
            "Run set-based (multi-SNP) SMR analysis. "
            "Performs SMR using multiple SNPs as instruments simultaneously, "
            "providing increased statistical power when LD structure is complex."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "bfile": {
                    "type": "string",
                    "description": "PLINK binary genotype file prefix for LD reference.",
                },
                "gwas_summary": {
                    "type": "string",
                    "description": "GWAS summary statistics file.",
                },
                "beqtl_summary": {
                    "type": "string",
                    "description": "Binary eQTL summary file (BESD format).",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix.",
                },
                "maf": {
                    "type": "number",
                    "description": "Minor allele frequency cutoff.",
                    "default": 0.01,
                    "minimum": 0,
                    "maximum": 0.5,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
                "heidi_off": {
                    "type": "boolean",
                    "description": "Skip HEIDI test.",
                    "default": False,
                },
                "ld_upper_limit": {
                    "type": "number",
                    "description": "LD r-squared upper threshold.",
                    "default": 0.9,
                    "minimum": 0,
                    "maximum": 1,
                },
                "ld_lower_limit": {
                    "type": "number",
                    "description": "LD r-squared lower threshold.",
                    "default": 0.05,
                    "minimum": 0,
                    "maximum": 1,
                },
                "max_num_ld": {
                    "type": "integer",
                    "description": "Maximum LD SNPs for HEIDI.",
                    "default": 500,
                    "minimum": 500,
                    "maximum": 10000,
                },
            },
            "required": ["bfile", "gwas_summary", "beqtl_summary", "out"],
        },
    },
    {
        "name": "make_bld",
        "description": (
            "Create a BLD (binary LD) file from a PLINK bfile. "
            "The BLD format stores LD information in a compact binary form, "
            "used as a reference LD panel for SMR/HEIDI analyses."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "bfile": {
                    "type": "string",
                    "description": "PLINK binary genotype file prefix.",
                },
                "out": {
                    "type": "string",
                    "description": "Output BLD file prefix.",
                },
                "maf": {
                    "type": "number",
                    "description": "Minor allele frequency cutoff.",
                    "default": 0.01,
                    "minimum": 0,
                    "maximum": 0.5,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["bfile", "out"],
        },
    },
    {
        "name": "update_freq",
        "description": (
            "Update allele frequencies in a BESD file. "
            "Replaces the stored allele frequencies with those from an "
            "external frequency file."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "beqtl_summary": {
                    "type": "string",
                    "description": "Binary eQTL summary file (BESD format).",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix.",
                },
                "freq_file": {
                    "type": "string",
                    "description": "External allele frequency file.",
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["beqtl_summary", "out", "freq_file"],
        },
    },
    {
        "name": "meta_analysis",
        "description": (
            "Run meta-analysis of multiple eQTL studies. "
            "Combines effect sizes across multiple BESD files listed in a "
            "flist file. Supports standard meta-analysis or MeCS method."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "besd_flist": {
                    "type": "string",
                    "description": "File listing BESD file prefixes to combine.",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix.",
                },
                "method": {
                    "type": "string",
                    "enum": ["mecs", "meta"],
                    "description": "Meta-analysis method: 'mecs' (MeCS) or 'meta' (standard). "
                    "Note: 'meta' flag is not in SMR's FLAGS_VALID_CK validation list "
                    "and may be rejected. Use 'mecs' for reliability.",
                    "default": "mecs",
                },
                "pmecs": {
                    "type": "number",
                    "description": "P-value threshold for MeCS (only used with method='mecs').",
                    "default": 0.01,
                    "minimum": 0,
                    "maximum": 1,
                },
                "nmecs": {
                    "type": "integer",
                    "description": "Number of common SNPs for MeCS correlation.",
                    "default": 100,
                    "minimum": 1,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["besd_flist", "out"],
        },
    },
    {
        "name": "combine_besd",
        "description": (
            "Combine multiple BESD files into one. "
            "Merges effect-size data from multiple BESD files listed in "
            "a flist file into a single output BESD file."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "besd_flist": {
                    "type": "string",
                    "description": "File listing BESD file prefixes to combine.",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix.",
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["besd_flist", "out"],
        },
    },
    {
        "name": "count_cis",
        "description": (
            "Count cis-eQTL in a BESD file. "
            "Outputs the top SNP in each cis-region, where cis is defined "
            "as within cis_window Kb of the probe."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "beqtl_summary": {
                    "type": "string",
                    "description": "Binary eQTL summary file (BESD format).",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix.",
                },
                "p_threshold": {
                    "type": "number",
                    "description": "P-value threshold for cis-eQTL.",
                    "default": 5e-8,
                    "minimum": 0,
                    "maximum": 1,
                },
                "cis_window": {
                    "type": "integer",
                    "description": "Cis window size in Kb.",
                    "default": 2000,
                    "minimum": 1,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["beqtl_summary", "out"],
        },
    },
    {
        "name": "count_trans",
        "description": (
            "Count trans-eQTL in a BESD file. "
            "Outputs the top SNP in each trans-region, where trans is "
            "defined as outside the cis_window Kb region around the probe."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "beqtl_summary": {
                    "type": "string",
                    "description": "Binary eQTL summary file (BESD format).",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix.",
                },
                "p_threshold": {
                    "type": "number",
                    "description": "P-value threshold for trans-eQTL.",
                    "default": 5e-8,
                    "minimum": 0,
                    "maximum": 1,
                },
                "cis_window": {
                    "type": "integer",
                    "description": "Cis window size in Kb (defines the cis boundary).",
                    "default": 2000,
                    "minimum": 1,
                },
                "trans_window": {
                    "type": "integer",
                    "description": "Trans window size in Kb.",
                    "default": 1000,
                    "minimum": 0,
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["beqtl_summary", "out"],
        },
    },
    {
        "name": "update_epi_esi",
        "description": (
            "Update EPI (probe info) and/or ESI (SNP info) files in a BESD file. "
            "Replaces the stored probe or SNP annotation with external files. "
            "At least one of epi_file or esi_file must be provided."
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "beqtl_summary": {
                    "type": "string",
                    "description": "Binary eQTL summary file (BESD format).",
                },
                "out": {
                    "type": "string",
                    "description": "Output file prefix.",
                },
                "epi_file": {
                    "type": "string",
                    "description": "External EPI (probe annotation) file.",
                    "default": "",
                },
                "esi_file": {
                    "type": "string",
                    "description": "External ESI (SNP annotation) file.",
                    "default": "",
                },
                "thread_num": {
                    "type": "integer",
                    "description": "Number of OpenMP threads.",
                    "default": 1,
                    "minimum": 1,
                },
            },
            "required": ["beqtl_summary", "out"],
        },
    },
    {
        "name": "get_version",
        "description": "Get SMR version information and build details.",
        "inputSchema": {"type": "object", "properties": {}},
    },
    {
        "name": "run_raw_command",
        "description": (
            "Run SMR with arbitrary command-line arguments. "
            "This is an escape hatch for advanced usage not covered by the "
            "typed tool wrappers. Pass the arguments exactly as they would "
            "appear on the command line (without the 'smr' prefix).\n\n"
            "Example: run_raw_command(['--bfile', 'ref', '--gwas-summary', 'g.ma', "
            "'--beqtl-summary', 'e.besd', '--out', 'o', '--smr', '--heidi-off', '--trans'])"
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "args": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of command-line arguments (without the 'smr' prefix).",
                },
            },
            "required": ["args"],
        },
    },
]


# ─── Tool Dispatch ───────────────────────────────────────────────────────────

TOOL_FUNCTIONS = {
    "run_smr_analysis": tool_run_smr_analysis,
    "make_besd": tool_make_besd,
    "make_besd_dense": tool_make_besd_dense,
    "plot_smr": tool_plot_smr,
    "query_besd": tool_query_besd,
    "show_sample_size": tool_show_sample_size,
    "recode_besd": tool_recode_besd,
    "smr_multi": tool_smr_multi,
    "make_bld": tool_make_bld,
    "update_freq": tool_update_freq,
    "meta_analysis": tool_meta_analysis,
    "combine_besd": tool_combine_besd,
    "count_cis": tool_count_cis,
    "count_trans": tool_count_trans,
    "update_epi_esi": tool_update_epi_esi,
    "get_version": tool_get_version,
    "run_raw_command": tool_run_raw_command,
}


# ─── MCP Protocol Implementation ─────────────────────────────────────────────


def send_message(msg: dict) -> None:
    """Write a JSON-RPC message to stdout."""
    sys.stdout.write(json.dumps(msg))
    sys.stdout.write("\n")
    sys.stdout.flush()


def handle_initialize(req_id, params):
    """Handle the initialize request."""
    send_message(
        {
            "jsonrpc": "2.0",
            "id": req_id,
            "result": {
                "protocolVersion": PROTOCOL_VERSION,
                "capabilities": {"tools": {"listChanged": False}},
                "serverInfo": {
                    "name": SERVER_NAME,
                    "version": SERVER_VERSION,
                },
            },
        }
    )


def handle_tools_list(req_id, params):
    """Handle the tools/list request."""
    send_message(
        {
            "jsonrpc": "2.0",
            "id": req_id,
            "result": {"tools": TOOLS},
        }
    )


def handle_tools_call(req_id, params):
    """Handle the tools/call request."""
    tool_name = params.get("name")
    arguments = params.get("arguments", {})

    if tool_name not in TOOL_FUNCTIONS:
        send_message(
            {
                "jsonrpc": "2.0",
                "id": req_id,
                "error": {
                    "code": -32601,
                    "message": f"Unknown tool: {tool_name}",
                },
            }
        )
        return

    func = TOOL_FUNCTIONS[tool_name]

    try:
        result_text = func(**arguments)
        send_message(
            {
                "jsonrpc": "2.0",
                "id": req_id,
                "result": {
                    "content": [{"type": "text", "text": result_text}],
                    "isError": "FAILED" in result_text,
                },
            }
        )
    except TypeError as e:
        # Argument validation error
        error_text = f"Error calling tool '{tool_name}': {e}"
        send_message(
            {
                "jsonrpc": "2.0",
                "id": req_id,
                "result": {
                    "content": [{"type": "text", "text": error_text}],
                    "isError": True,
                },
            }
        )
    except Exception as e:
        error_text = (
            f"Error calling tool '{tool_name}': {e}\n\n"
            f"Traceback:\n{traceback.format_exc()}"
        )
        send_message(
            {
                "jsonrpc": "2.0",
                "id": req_id,
                "result": {
                    "content": [{"type": "text", "text": error_text}],
                    "isError": True,
                },
            }
        )


def handle_request(msg: dict) -> None:
    """Route an incoming JSON-RPC message to the appropriate handler."""
    method = msg.get("method")
    req_id = msg.get("id")
    params = msg.get("params", {})

    if method == "initialize":
        handle_initialize(req_id, params)
    elif method == "notifications/initialized":
        # No response needed for notifications
        pass
    elif method == "tools/list":
        handle_tools_list(req_id, params)
    elif method == "tools/call":
        handle_tools_call(req_id, params)
    elif method == "ping":
        send_message(
            {"jsonrpc": "2.0", "id": req_id, "result": {}}
        )
    else:
        if req_id is not None:
            send_message(
                {
                    "jsonrpc": "2.0",
                    "id": req_id,
                    "error": {
                        "code": -32601,
                        "message": f"Method not found: {method}",
                    },
                }
            )


def main() -> None:
    """Main MCP server loop — reads JSON-RPC messages from stdin."""
    # Log server start to stderr (stdout is reserved for MCP protocol)
    print(f"[{SERVER_NAME}] v{SERVER_VERSION} starting...", file=sys.stderr)
    print(f"[{SERVER_NAME}] SMR binary: {SMR_BIN}", file=sys.stderr)

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue

        try:
            msg = json.loads(line)
            handle_request(msg)
        except json.JSONDecodeError as e:
            # Send error response for malformed JSON
            send_message(
                {
                    "jsonrpc": "2.0",
                    "id": None,
                    "error": {
                        "code": -32700,
                        "message": f"Parse error: {e}",
                    },
                }
            )
        except Exception as e:
            # Catch-all to prevent the server from crashing
            print(f"[{SERVER_NAME}] Unexpected error: {e}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)

    print(f"[{SERVER_NAME}] shutting down.", file=sys.stderr)


if __name__ == "__main__":
    main()
