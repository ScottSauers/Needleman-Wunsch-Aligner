use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::io::{self, BufRead, Write};

// Translate DNA/RNA sequence to amino acid sequence
pub fn translate_sequence(dna_sequence: &str) -> Result<String, Box<dyn Error>> {
    let codon_table = get_codon_table();

    let mut aa_sequence = String::new();

    // Convert to uppercase for fun (already should be upppercase) and replace 'U' with 'T' to handle RNA sequences in a way that aligns with the existing codon table
    let dna_sequence = dna_sequence.to_uppercase().replace('U', "T");

    let start_index = dna_sequence.find("ATG");

    // Check if a start codon 'ATG' is found
    if start_index.is_none() {
        return Err(format!(
            "Start codon 'ATG' not found in the provided sequence: '{}'. Length: {}",
            dna_sequence,
            dna_sequence.len()
        )
        .into());
    }

    let mut i = start_index.unwrap();

    while i + 3 <= dna_sequence.len() {
        let codon = &dna_sequence[i..i + 3];
        let amino_acid = codon_table.get(codon).unwrap_or(&"X");
        if *amino_acid == "*" {
            break;
        }
        aa_sequence.push_str(amino_acid);
        i += 3;
    }

    Ok(aa_sequence)
}

// Read a sequence from a FASTA file
pub fn read_fasta_sequence(file_path: &str) -> Result<(String, String), Box<dyn Error>> {
    let file = fs::File::open(file_path)?;
    let reader = io::BufReader::new(file);

    let mut header = String::new();
    let mut sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            header = line;
        } else {
            sequence.push_str(&line.trim());
        }
    }

    Ok((header, sequence))
}

pub fn save_sequence_to_file(
    file_path: &str,
    header: &str,
    sequence: &str,
) -> Result<(), Box<dyn Error>> {
    let mut file = fs::File::create(file_path)?;
    writeln!(file, "{}", header)?;
    writeln!(file, "{}", sequence)?;
    Ok(())
}

fn get_codon_table() -> HashMap<&'static str, &'static str> {
    let mut codon_table = HashMap::new();
    codon_table.insert("TTT", "F");
    codon_table.insert("TTC", "F");
    codon_table.insert("TTA", "L");
    codon_table.insert("TTG", "L");
    codon_table.insert("CTT", "L");
    codon_table.insert("CTC", "L");
    codon_table.insert("CTA", "L");
    codon_table.insert("CTG", "L");
    codon_table.insert("ATT", "I");
    codon_table.insert("ATC", "I");
    codon_table.insert("ATA", "I");
    codon_table.insert("ATG", "M"); // Start
    codon_table.insert("GTT", "V");
    codon_table.insert("GTC", "V");
    codon_table.insert("GTA", "V");
    codon_table.insert("GTG", "V");
    codon_table.insert("TCT", "S");
    codon_table.insert("TCC", "S");
    codon_table.insert("TCA", "S");
    codon_table.insert("TCG", "S");
    codon_table.insert("CCT", "P");
    codon_table.insert("CCC", "P");
    codon_table.insert("CCA", "P");
    codon_table.insert("CCG", "P");
    codon_table.insert("ACT", "T");
    codon_table.insert("ACC", "T");
    codon_table.insert("ACA", "T");
    codon_table.insert("ACG", "T");
    codon_table.insert("GCT", "A");
    codon_table.insert("GCC", "A");
    codon_table.insert("GCA", "A");
    codon_table.insert("GCG", "A");
    codon_table.insert("TAT", "Y");
    codon_table.insert("TAC", "Y");
    codon_table.insert("TAA", "*");
    codon_table.insert("TAG", "*");
    codon_table.insert("CAT", "H");
    codon_table.insert("CAC", "H");
    codon_table.insert("CAA", "Q");
    codon_table.insert("CAG", "Q");
    codon_table.insert("AAT", "N");
    codon_table.insert("AAC", "N");
    codon_table.insert("AAA", "K");
    codon_table.insert("AAG", "K");
    codon_table.insert("GAT", "D");
    codon_table.insert("GAC", "D");
    codon_table.insert("GAA", "E");
    codon_table.insert("GAG", "E");
    codon_table.insert("TGT", "C");
    codon_table.insert("TGC", "C");
    codon_table.insert("TGA", "*");
    codon_table.insert("TGG", "W");
    codon_table.insert("CGT", "R");
    codon_table.insert("CGC", "R");
    codon_table.insert("CGA", "R");
    codon_table.insert("CGG", "R");
    codon_table.insert("AGT", "S");
    codon_table.insert("AGC", "S");
    codon_table.insert("AGA", "R");
    codon_table.insert("AGG", "R");
    codon_table.insert("GGT", "G");
    codon_table.insert("GGC", "G");
    codon_table.insert("GGA", "G");
    codon_table.insert("GGG", "G");
    codon_table
}
