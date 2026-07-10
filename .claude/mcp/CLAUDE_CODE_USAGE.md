# Using SMR with AI Tools (Claude, Codex, OpenCode)

The SMR release package includes a built-in MCP (Model Context Protocol) server that lets you invoke SMR functionality through AI tools such as Claude, Codex, and OpenCode via natural language conversation.

## Prerequisites

- Python 3.8+ (the MCP server depends only on the Python standard library; no packages need to be installed)
- The SMR binary executable (already included in this release package)
## Release Package Structure

The release package is a single directory (folder name may vary) containing the SMR binary and the MCP server script:

```
<SMR_DIR>/
├── smr                          # Pre-compiled SMR binary
└── mcp/
    └── smr_mcp_server.py        # MCP server script
```

In the configuration examples below, replace `<SMR_DIR>` with the actual absolute path to the release directory on your system (e.g., `/home/user/smr-1.4.2`).

> **Note:** The MCP server auto-detects the `smr` binary relative to its own location (`../smr`), so the folder name does not matter. You only need the absolute path for the AI client configuration.

## Available Tools

The MCP server exposes 17 tools. The AI automatically selects the appropriate tool based on your conversational intent:

| Tool | Function | SMR Flags |
|------|----------|-----------|
| `run_smr_analysis` | Run the main SMR test (Wald ratio + HEIDI) | `--smr` |
| `smr_multi` | Run set-based (multi-SNP) SMR analysis | `--smr-multi` |
| `make_besd` | Create a sparse BESD file from a text eQTL summary | `--make-besd` |
| `make_besd_dense` | Create a dense BESD file from a text eQTL summary | `--make-besd-dense` |
| `query_besd` | Query a BESD file for significant SNP-probe associations | `--query` |
| `show_sample_size` | Display the sample size stored in a BESD file | `--show-n` |
| `recode_besd` | Convert a BESD file to COJO/SMR text format | `--recode` |
| `plot_smr` | Generate SMR locus plots | `--plot` |
| `make_bld` | Create a BLD (binary LD) file from a PLINK bfile | `--make-bld` |
| `update_freq` | Update allele frequencies in a BESD file | `--update-freq` |
| `meta_analysis` | Run meta-analysis of multiple eQTL studies (MeCS) | `--mecs` |
| `combine_besd` | Combine multiple BESD files into one | `--besd-flist` |
| `count_cis` | Count cis-eQTL in a BESD file | `--descriptive-cis` |
| `count_trans` | Count trans-eQTL in a BESD file | `--descriptive-trans` |
| `update_epi_esi` | Update EPI/ESI annotation files in a BESD file | `--update-epi` / `--update-esi` |
| `get_version` | Get SMR version information | (run with no arguments) |
| `run_raw_command` | Pass arbitrary SMR command-line arguments directly | (escape hatch) |

## Configuration

### 1. Claude Desktop (macOS)

Edit `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "smr": {
      "command": "python3",
      "args": ["<SMR_DIR>/mcp/smr_mcp_server.py"],
      "env": {
        "SMR_BIN": "<SMR_DIR>/smr"
      }
    }
  }
}
```

Restart Claude Desktop after saving.

### 2. Claude Code

Create `.mcp.json` in your project root:

```json
{
  "mcpServers": {
    "smr": {
      "command": "python3",
      "args": ["<SMR_DIR>/mcp/smr_mcp_server.py"],
      "env": {
        "SMR_BIN": "<SMR_DIR>/smr"
      }
    }
  }
}
```

Or add the same entry to the user-level config at `~/.claude.json`.

### 3. OpenCode

Create or edit `opencode.json` in your project root:

```json
{
  "$schema": "https://opencode.ai/config.json",
  "mcp": {
    "smr": {
      "type": "local",
      "command": ["python3", "<SMR_DIR>/mcp/smr_mcp_server.py"],
      "enabled": true,
      "env": {
        "SMR_BIN": "<SMR_DIR>/smr"
      }
    }
  }
}
```

Restart OpenCode after saving.

### 4. Codex (OpenAI)

Add the following to `~/.codex/config.toml`:

```toml
[mcp_servers.smr]
command = "python3"
args = ["<SMR_DIR>/mcp/smr_mcp_server.py"]
env = { SMR_BIN = "<SMR_DIR>/smr" }
```

Restart Codex after saving.

> **Note:** The `SMR_BIN` environment variable is optional. If omitted, the MCP server will auto-detect the `smr` binary at `<SMR_DIR>/smr` (one directory above the script). You only need `SMR_BIN` if you want to point to a binary in a non-standard location.

## Usage Examples

Once configured, simply converse with the AI in natural language:

**Run an SMR analysis:**
> "Run an SMR analysis with bfile data/1kg_eur, gwas-summary data/bmi.ma, beqtl-summary data/eqtl.besd, output to results/smr_out"

The AI will automatically call the `run_smr_analysis` tool, generating and executing a command like:
```
smr --bfile data/1kg_eur --gwas-summary data/bmi.ma \
    --beqtl-summary data/eqtl.besd --out results/smr_out \
    --maf 0.01 --thread-num 1 --ld-upper-limit 0.9 \
    --ld-lower-limit 0.05 --max_num_ld 500 \
    --diff-freq 0.2 --diff-freq-prop 0.05 --smr
```

**Create a BESD file:**
> "Convert data/eqtl.txt to BESD format, output to data/eqtl.besd"

The AI will call the `make_besd` tool.

**Query BESD contents:**
> "Show me the sample size in data/eqtl.besd"

The AI will call the `show_sample_size` tool.

**Advanced usage (escape hatch):**
> "Use run_raw_command to run: smr --bfile ref --gwas-summary g.ma --beqtl-summary e.besd --out o --smr --heidi-off --trans"

When the preset tools don't cover your needs, you can pass arbitrary SMR command-line arguments via `run_raw_command`.

## Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `SMR_BIN` | Full path to the SMR binary. Optional — auto-detected if not set. | `../smr` relative to the script |
| `SMR_MCP_TIMEOUT` | Timeout for a single command execution (seconds) | 3600 |
| `SMR_MCP_CWD` | Working directory for SMR execution | Inherits current directory |

If `SMR_BIN` is not set, the MCP server searches for the binary in the following order:

1. The `smr` binary inside the release package (one directory above the script — auto-detected based on the script's own location, so the folder name does not matter)
2. The `build/smr` path in a source tree
3. The `smr` found in the system PATH

## Troubleshooting

### "SMR binary not found"

The SMR binary could not be located. Set the `SMR_BIN` environment variable to the full path of the smr binary.

### "Permission denied"

The smr binary does not have execute permission. Run:
```bash
chmod +x smr
```

### "error while loading shared libraries"

The system is missing required dependency libraries. Ensure that `libgomp` (OpenMP runtime) and `libz` are installed.

### MCP server fails to start

Check your Python version:
```bash
python3 --version
```
Python 3.8 or later is required.

## Further Information

- SMR official documentation: https://yanglab.westlake.edu.cn/software/smr/
- MCP protocol specification: https://modelcontextprotocol.io/
