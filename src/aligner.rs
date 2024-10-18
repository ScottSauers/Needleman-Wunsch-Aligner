use clap::{Arg, Command};
use std::error::Error;
use std::fs;
use std::io::{Write};
use std::path::Path;

mod utils;
use crate::utils::{read_fasta_sequence};

// Holds alignment result
struct AlignmentResult {
    alignment_score: i32,
    align1: String,
    align2: String,
    alignment_visualization: String,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("Sequence Aligner")
        .arg(
            Arg::new("query")
                .short('q')
                .long("query")
                .value_name("FILE")
                .help("Query sequence file in FASTA format")
                .required(true)
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            Arg::new("reference")
                .short('r')
                .long("reference")
                .value_name("FILE")
                .help("Reference sequence file in FASTA format")
                .required(true)
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("Output alignment file")
                .required(true)
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            Arg::new("gap_penalty")
                .short('g')
                .long("gap")
                .value_name("INT")
                .help("Gap penalty (negative integer)")
                .required(true)
                .allow_hyphen_values(true)
                .value_parser(clap::value_parser!(i32)),
        )
        .arg(
            Arg::new("mismatch_penalty")
                .short('p')
                .long("mismatch")
                .value_name("INT")
                .help("Mismatch penalty (negative integer)")
                .required(true)
                .allow_hyphen_values(true)
                .value_parser(clap::value_parser!(i32)),
        )
        .arg(
            Arg::new("match_score")
                .short('m')
                .long("match")
                .value_name("INT")
                .help("Match score (positive integer)")
                .required(true)
                .value_parser(clap::value_parser!(i32)),
        )
        .arg(
            Arg::new("unpenalized_end_gaps")
                .short('u')
                .long("unpenalized")
                .help("Unpenalized start and end gaps")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("sequence_type")
                .short('t')
                .long("type")
                .value_name("TYPE")
                .help("Sequence type: 'nucleotide' or 'aminoacid'")
                .required(true)
                .value_parser(clap::value_parser!(String)),
        )
        .get_matches();

    let query_file = matches.get_one::<String>("query").unwrap();
    let reference_file = matches.get_one::<String>("reference").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();
    let gap_penalty: i32 = *matches.get_one::<i32>("gap_penalty").unwrap();
    let mismatch_penalty: i32 = *matches.get_one::<i32>("mismatch_penalty").unwrap();
    let match_score: i32 = *matches.get_one::<i32>("match_score").unwrap();
    let unpenalized_end_gaps = matches.get_flag("unpenalized_end_gaps");
    let sequence_type_input = matches.get_one::<String>("sequence_type").unwrap();
    let sequence_type = sequence_type_input.to_lowercase();
    println!("Sequence Type: {}", sequence_type);
    println!("Unpenalized End Gaps: {}", unpenalized_end_gaps);

    check_and_download_file(query_file)?;
    check_and_download_file(reference_file)?;

    if sequence_type == "nucleotide" {
        let (query_header, query_sequence) = read_fasta_sequence(query_file)?;
        let (reference_header, reference_sequence) = read_fasta_sequence(reference_file)?;
    
        // No translation
        let alignment = needleman_wunsch(
            &reference_sequence,
            &query_sequence,
            match_score,
            mismatch_penalty,
            gap_penalty,
            unpenalized_end_gaps,
        );
    
        write_alignment_output(output_file, &alignment, &reference_header, &query_header)?;
    } else if sequence_type == "aminoacid" {
        let (query_header, query_aa_sequence) = read_fasta_sequence(query_file)?;
        let (reference_header, reference_aa_sequence) = read_fasta_sequence(reference_file)?;
    
        let alignment = needleman_wunsch(
            &reference_aa_sequence,
            &query_aa_sequence,
            match_score,
            mismatch_penalty,
            gap_penalty,
            unpenalized_end_gaps,
        );
        
        write_alignment_output(output_file, &alignment, &reference_header, &query_header)?;
    } else {
        return Err("Invalid sequence type: specify 'nucleotide' or 'aminoacid'.".into());
    }


    Ok(())
}

// Download if file nonexistent
fn check_and_download_file(file_path: &str) -> Result<(), Box<dyn Error>> {
    if Path::new(file_path).exists() {
        println!("File '{}' already exists.", file_path);
    } else {
        println!("File '{}' not found. Downloading...", file_path);
        let base_url = "https://raw.githubusercontent.com/ScottSauers/Needleman-Wunsch-Aligner/main/";
        let url = format!("{}{}", base_url, file_path);

        let response = reqwest::blocking::get(&url)?;

        if response.status().is_success() {
            let content = response.text()?;
            fs::write(file_path, content)?;
            println!("File '{}' downloaded successfully.", file_path);
        } else {
            return Err(format!(
                "Failed to download file '{}'. HTTP Status: {}",
                file_path,
                response.status()
            )
            .into());
        }
    }
    Ok(())
}


// Needleman-Wunsch algorithm, global and semi-global alignment
fn needleman_wunsch(
    seq1: &str,
    seq2: &str,
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
    unpenalized_end_gaps: bool,
) -> AlignmentResult {
    let m = seq1.len();
    let n = seq2.len();

    // Init score and traceback matrices
    // The first row and column are initialized based on whether end gaps are penalized
    let mut score_matrix = vec![vec![0; n + 1]; m + 1];
    let mut trace_matrix = vec![vec![' '; n + 1]; m + 1];

    // Init first row and column
    for i in 0..=m {
        score_matrix[i][0] = if unpenalized_end_gaps { 0 } else { i as i32 * gap_penalty };
        trace_matrix[i][0] = 'U'; // Up
    }
    for j in 0..=n {
        score_matrix[0][j] = if unpenalized_end_gaps { 0 } else { j as i32 * gap_penalty };
        trace_matrix[0][j] = 'L'; // Left
    }
    trace_matrix[0][0] = '0';

    // Fill score and traceback matrices
    for i in 1..=m {
        let c1 = seq1.chars().nth(i - 1).unwrap();
        for j in 1..=n {
            let c2 = seq2.chars().nth(j - 1).unwrap();
            let match_mismatch = if c1 == c2 { match_score } else { mismatch_penalty };
            let diag_score = score_matrix[i - 1][j - 1] + match_mismatch;
            let up_score = score_matrix[i - 1][j]
                + if unpenalized_end_gaps && (i == m) { 0 } else { gap_penalty };
            let left_score = score_matrix[i][j - 1]
                + if unpenalized_end_gaps && (j == n) { 0 } else { gap_penalty };
            let max_score = diag_score.max(up_score).max(left_score);
            score_matrix[i][j] = max_score;
            if max_score == diag_score {
                trace_matrix[i][j] = 'D'; // Diagonal
            } else if max_score == up_score {
                trace_matrix[i][j] = 'U'; // Up
            } else {
                trace_matrix[i][j] = 'L'; // Left
            }
        }
    }

    // Determine the traceback starting point
    let (max_score, (start_i, start_j)) = if unpenalized_end_gaps {
        // Semi-global alignment to find the maximum score in the last row and column
        find_max_in_last_row_and_column(&score_matrix, m, n)
    } else {
        // Global alignment using the bottom-right corner
        (score_matrix[m][n], (m, n))
    };

    // Traceback to get alignment
    let mut align1 = String::new();
    let mut align2 = String::new();
    let mut alignment_visualization = String::new();
    let mut i = start_i;
    let mut j = start_j;

    while i > 0 || j > 0 {
        if i > 0 && j > 0 && trace_matrix[i][j] == 'D' {
            let c1 = seq1.chars().nth(i - 1).unwrap();
            let c2 = seq2.chars().nth(j - 1).unwrap();
            align1.push(c1);
            align2.push(c2);
            if c1 == c2 {
                alignment_visualization.push('|');
            } else {
                alignment_visualization.push('x');
            }
            i -= 1;
            j -= 1;
        } else if i > 0 && trace_matrix[i][j] == 'U' {
            let c1 = seq1.chars().nth(i - 1).unwrap();
            align1.push(c1);
            align2.push('_');
            alignment_visualization.push(' ');
            i -= 1;
        } else if j > 0 && trace_matrix[i][j] == 'L' {
            let c2 = seq2.chars().nth(j - 1).unwrap();
            align1.push('_');
            align2.push(c2);
            alignment_visualization.push(' ');
            j -= 1;
        } else {
            break;
        }
    }

    // Reverse the alignments for correct ordering
    align1 = align1.chars().rev().collect();
    align2 = align2.chars().rev().collect();
    alignment_visualization = alignment_visualization.chars().rev().collect();

    // Final alignment score
    let alignment_score = max_score;

    AlignmentResult {
        alignment_score,
        align1,
        align2,
        alignment_visualization,
    }
}

// Helper function to find the max score in the last row and last column
fn find_max_in_last_row_and_column(
    score_matrix: &Vec<Vec<i32>>,
    m: usize,
    n: usize,
) -> (i32, (usize, usize)) {
    let mut max_score = score_matrix[m][n];
    let mut max_pos = (m, n);

    // Check last row
    for j in 0..=n {
        if score_matrix[m][j] > max_score {
            max_score = score_matrix[m][j];
            max_pos = (m, j);
        }
    }

    // Check last column
    for i in 0..=m {
        if score_matrix[i][n] > max_score {
            max_score = score_matrix[i][n];
            max_pos = (i, n);
        }
    }

    (max_score, max_pos)
}

// Write alignment output to file
fn write_alignment_output(
    output_file: &str,
    alignment: &AlignmentResult,
    reference_header: &str,
    query_header: &str,
) -> Result<(), Box<dyn Error>> {
    let mut file = fs::File::create(output_file)?;
    writeln!(file, "{}", alignment.alignment_score)?;
    writeln!(file, "{}", reference_header)?;
    writeln!(file, "{}", alignment.align1.replace(' ', "_"))?;
    writeln!(file, "{}", alignment.alignment_visualization)?;
    writeln!(file, "{}", alignment.align2.replace(' ', "_"))?;
    writeln!(file, "{}", query_header)?;
    Ok(())
}
