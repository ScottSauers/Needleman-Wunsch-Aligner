use statrs::distribution::{ContinuousCDF, Normal};
use std::error::Error;
use std::fs;
use std::io::{self, BufRead};
use std::path::Path;
use std::process::Command;

mod utils;
use crate::utils::{read_fasta_sequence, translate_sequence, save_sequence_to_file};

const SEQUENCE_TYPE_NUCLEOTIDE: &str = "nucleotide";
const SEQUENCE_TYPE_AMINOACID: &str = "aminoacid";

fn main() -> Result<(), Box<dyn Error>> {
    let aligner_path = "./target/release/aligner";

    if !Path::new(aligner_path).exists() {
        println!(
            "Error: Aligner binary not found at '{}'. Build it using 'cargo build --release'.",
            aligner_path
        );
        return Ok(());
    }

    println!("Running alignment with penalties for start/end gaps. Query is pfizer_mrna.fna. Reference is sars_spike_protein.fna. Check output in question2_output.txt. Gap penalty is -2, mismatch penalty is -1, and match score is 1.");

    let alignment2_output = "question2_output.txt";
    let status = Command::new(&aligner_path)
        .args(&[
            "--query=pfizer_mrna.fna",
            "--reference=sars_spike_protein.fna",
            "--output=question2_output.txt",
            "--gap=-2",
            "--mismatch=-1",
            "--match=1",
            "--type", SEQUENCE_TYPE_NUCLEOTIDE,
        ])
        .status()?;
    if status.success() {
        println!("Alignment with penalties for start/end gaps completed successfully.");
    } else {
        println!("Warning: alignment with penalties for start/end gaps failed.");
    }

    println!("Running alignment with free start/end gaps. Query is pfizer_mrna.fna. Reference is sars_spike_protein.fna. Check output at question3_output.txt. Gap penalty is -2, mismatch penalty is -1, and match score is 1. Flag unpenalized is active.");

    let alignment3_output = "question3_output.txt";
    let status = Command::new(&aligner_path)
        .args(&[
            "--query=pfizer_mrna.fna",
            "--reference=sars_spike_protein.fna",
            "--output=question3_output.txt",
            "--gap=-2",
            "--mismatch=-1",
            "--match=1",
            "--unpenalized",
            "--type", SEQUENCE_TYPE_NUCLEOTIDE, 
        ])
        .status()?;
    if status.success() {
        println!("Alignment with free start/end gaps completed successfully.");
    } else {
        println!("Warning: alignment with free start/end gaps did not complete successfully.");
    }

    println!("Getting alignment scores...");
    let score2 = read_alignment_score(alignment2_output)?;
    let score3 = read_alignment_score(alignment3_output)?;
    println!(
        "Alignment with penalties for start/end gaps score: {}",
        score2
    );
    println!("Alignment with free start/end gaps score: {}", score3);

    println!("Analyzing Alignment with penalties for start/end gaps...");
    let (matches2, mismatches2, gaps2) = parse_alignment(alignment2_output)?;
    let total_mismatches2 = mismatches2 + gaps2;
    println!("Alignment with penalties for start/end gaps (With penalties for start/end gaps):");
    println!("Matches: {}", matches2);
    println!("Mismatches: {}", mismatches2);
    println!("Gaps: {}", gaps2);
    println!("Total Mismatches (including gaps): {}", total_mismatches2);

    println!("Finding alignment with free start/end gaps...");
    let (matches3, mismatches3, gaps3) = parse_alignment(alignment3_output)?;
    let total_mismatches3 = mismatches3 + gaps3;
    println!("Alignment with free start/end gaps:");
    println!("Matches: {}", matches3);
    println!("Mismatches: {}", mismatches3);
    println!("Gaps: {}", gaps3);
    println!("Total Mismatches (including gaps): {}", total_mismatches3);

    println!("Reading nucleotide sequences.");
    let (ref_header, ref_sequence) = read_fasta_sequence("sars_spike_protein.fna")?;
    let (query_header, query_sequence) = read_fasta_sequence("pfizer_mrna.fna")?;

    println!("Translating sequences to amino acids.");
    let ref_aa_sequence = translate_sequence(&ref_sequence)?;
    let query_aa_sequence = translate_sequence(&query_sequence)?;

    save_sequence_to_file("sars_spike_protein.aa", &ref_header, &ref_aa_sequence)?;
    save_sequence_to_file("pfizer_mrna.aa", &query_header, &query_aa_sequence)?;

    println!(
        "Running amino acid alignment... Gap penalty -2, mismatch penalty -1, and match score 1."
    );

    let alignment_aa_output = "question8_output.txt";

    let status = Command::new(&aligner_path)
        .args(&[
            "--query=pfizer_mrna.aa",
            "--reference=sars_spike_protein.aa",
            "--output=question8_output.txt",
            "--gap=-2",
            "--mismatch=-1",
            "--match=1",
            "--type", SEQUENCE_TYPE_AMINOACID,
        ])
        .status()?;
    if status.success() {
        println!("Amino Acid Alignment completed successfully.");
    } else {
        println!("Warning: amino acid alignment failed.");
    }

    if Path::new(alignment_aa_output).exists() {
        let (matches_aa, mismatches_aa, gaps_aa) = parse_alignment(alignment_aa_output)?;
        let total_mismatches_aa = mismatches_aa + gaps_aa;
        println!("Amino acid alignment:");
        println!("Matches: {}", matches_aa);
        println!("Mismatches: {}", mismatches_aa);
        println!("Gaps: {}", gaps_aa);
        println!("Total mismatches including gaps: {}", total_mismatches_aa);

        println!("Differences between a.a. sequences:");
        let differences_aa = extract_differences_aa(alignment_aa_output)?;
        if differences_aa.is_empty() {
            println!("No differences found");
        } else {
            for diff in differences_aa {
                println!("{}", diff);
            }
        }
    } else {
        println!("No amino acid file '{}' found.", alignment_aa_output);
    }

    let (gc_real, gc_real_total, gc_real_count) = calculate_gc_content("sars_spike_protein.fna")?;
    let (gc_pfizer, gc_pfizer_total, gc_pfizer_count) = calculate_gc_content("pfizer_mrna.fna")?;
    println!(
        "sars_spike_protein.fna GC Content: {:.2}% (GC Count: {}, Total: {})",
        gc_real, gc_real_count, gc_real_total
    );
    println!(
        "pfizer_mrna.fna GC Content: {:.2}% (GC Count: {}, Total: {})",
        gc_pfizer, gc_pfizer_count, gc_pfizer_total
    );

    println!("Performing Z-Test on GC Content...");

    let z_score = perform_z_test(
        gc_real_count,
        gc_real_total,
        gc_pfizer_count,
        gc_pfizer_total,
    )?;
    println!("Z-Score: {:.4}", z_score);
    let normal_dist = Normal::new(0.0, 1.0)?;
    let p_value = 2.0 * (1.0 - normal_dist.cdf(z_score.abs())).max(0.0);
    println!(
        "P-Value: {}",
        if p_value == 0.0 {
            "< 1e-10".to_string()
        } else {
            format!("{:.4}", p_value)
        }
    );
    if p_value < 0.05 {
        println!("Result: Significant difference in GC content.");
    } else {
        println!("Result: No significant difference in GC content.");
    }

    Ok(())
}

fn read_alignment_score(file_path: &str) -> Result<i32, Box<dyn Error>> {
    let file = fs::File::open(file_path)?;
    let mut reader = io::BufReader::new(file);
    let mut score_line = String::new();
    reader.read_line(&mut score_line)?;
    let score: i32 = score_line.trim().parse()?;
    Ok(score)
}

fn parse_alignment(file_path: &str) -> Result<(usize, usize, usize), Box<dyn Error>> {
    let file = fs::File::open(file_path)?;
    let reader = io::BufReader::new(file);
    let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;
    if lines.len() < 5 {
        println!("File '{}' needs to have at least six lines.", file_path);
        return Ok((0, 0, 0));
    }
    let visual_line = &lines[3];
    let matches = visual_line.matches('|').count();
    let mismatches = visual_line.matches('x').count();
    let gaps = visual_line.matches(' ').count();
    Ok((matches, mismatches, gaps))
}

fn extract_differences_aa(file_path: &str) -> Result<Vec<String>, Box<dyn Error>> {
    let file = fs::File::open(file_path)?;
    let reader = io::BufReader::new(file);
    let mut lines = reader.lines();
    lines.next();

    let _header1 = lines.next().unwrap_or(Ok(String::new()))?;
    let seq1 = lines.next().unwrap_or(Ok(String::new()))?;
    let visual = lines.next().unwrap_or(Ok(String::new()))?;
    let seq2 = lines.next().unwrap_or(Ok(String::new()))?;
    let _header2 = lines.next().unwrap_or(Ok(String::new()))?;

    let mut differences = Vec::new();
    for (i, ((c1, c2), v)) in seq1
        .chars()
        .zip(seq2.chars())
        .zip(visual.chars())
        .enumerate()
    {
        if v != '|' {
            differences.push(format!("Position {}: {} vs {}", i + 1, c1, c2));
        }
    }
    Ok(differences)
}

fn calculate_gc_content(file_path: &str) -> Result<(f64, f64, f64), Box<dyn Error>> {
    let file = fs::File::open(file_path)?;
    let reader = io::BufReader::new(file);
    let mut gc_count = 0.0;
    let mut total = 0.0;
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 && line.starts_with('>') {
            continue;
        }
        for c in line.chars() {
            match c.to_ascii_uppercase() {
                'G' | 'C' => {
                    gc_count += 1.0;
                    total += 1.0;
                }
                'A' | 'T' | 'U' => {
                    total += 1.0;
                }
                _ => (),
            }
        }
    }
    if total == 0.0 {
        println!("Failed. No valid nucleotides found in {}.", file_path);
        return Ok((0.0, 0.0, 0.0));
    }
    let gc_content = (gc_count / total) * 100.0;
    Ok((gc_content, total, gc_count))
}

fn perform_z_test(
    gc_count1: f64,
    total1: f64,
    gc_count2: f64,
    total2: f64,
) -> Result<f64, Box<dyn Error>> {
    let p1 = gc_count1 / total1;
    let p2 = gc_count2 / total2;
    let p = (gc_count1 + gc_count2) / (total1 + total2);
    let se = ((p * (1.0 - p)) * (1.0 / total1 + 1.0 / total2)).sqrt();
    if se == 0.0 {
        return Err("Standard error is zero, cannot perform Z-Test.".into());
    }
    let z = (p1 - p2) / se;

    println!("GC Count 1: {}", gc_count1);
    println!("Total 1: {}", total1);
    println!("GC Count 2: {}", gc_count2);
    println!("Total 2: {}", total2);

    Ok(z)
}
