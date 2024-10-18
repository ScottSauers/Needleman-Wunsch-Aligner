**Project structure:**
- `src/`
  - `aligner.rs`: The main aligner program for the Needleman-Wunsch algorithm.
  - `analysis.rs`: A script that automates pre-set alignments and analyses.
  - `utils.rs`: Functions that have to be used by both the analysis and aligner file.
- `Cargo.toml`: Rust project configuration file with dependencies.

**Download the project:**

```bash
git clone https://github.com/ScottSauers/Needleman-Wunsch-Aligner.git
cd aligner
```

**Build and run:**

```bash
cargo build --release
```

Example run:

```bash
./target/release/aligner -q <a_query_file> -r <reference_file> -o <output_file> -g <gap_penalty_amount> -p <mismatch_penalty> -m <match_score> -t <sequence_type>
```

**Params:**

- `-q, --query`: Query sequence file in FASTA format.
- `-r, --reference`: Reference sequence file in FASTA format.
- `-o, --output`: Output alignment file.
- `-g, --gap`: Gap penalty (negative int).
- `-p, --mismatch`: Mismatch penalty (negative int).
- `-m, --match`: Match score (positive int).
- `-t, --type`: Sequence type (nucleotide or aminoacid).
- `-u, --unpenalized`: Unpenalized start and end gaps. Omit for penalized.


**Analysis**

The analysis script does some alignments automatically.

Run it like this:

```bash
cargo run --release --bin analysis
```
