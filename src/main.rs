use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;

const AMINO_ACIDS: &str = "ACDEFGHIKLMNPQRSTVWY";

fn read_fasta(file_path: &str) -> io::Result<Vec<String>> {
    let mut sequences = Vec::new();
    let mut current_seq = String::new();

    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_seq.is_empty() {
                sequences.push(current_seq.clone());
                current_seq.clear();
            }
        } else {
            current_seq.push_str(&line);
        }
    }

    if !current_seq.is_empty() {
        sequences.push(current_seq);
    }

    Ok(sequences)
}

fn calculate_pwm(sequences: &[String]) -> Vec<HashMap<char, f64>> {
    if sequences.is_empty() {
        return Vec::new();
    }

    let seq_length = sequences[0].len();
    let num_sequences = sequences.len();
    let mut pwm = vec![HashMap::new(); seq_length];

    for seq in sequences {
        assert_eq!(seq.len(), seq_length, "Sequences must have the same length");

        for (i, residue) in seq.chars().enumerate() {
            let count = pwm[i].entry(residue).or_insert(0.0);
            *count += 1.0;
        }
    }

    for position in &mut pwm {
        for aa in AMINO_ACIDS.chars() {
            let count = position.entry(aa).or_insert(0.0);
            *count /= num_sequences as f64;
        }
    }

    pwm
}

fn write_pwm_to_tsv(pwm: &[HashMap<char, f64>], output_path: &str) -> io::Result<()> {
    let mut file = File::create(output_path)?;

    // Write header
    let header: Vec<String> = AMINO_ACIDS.chars().map(|aa| aa.to_string()).collect();
    writeln!(file, "Position\t{}", header.join("\t"))?;

    // Write each position's data
    for (i, position) in pwm.iter().enumerate() {
        let row: Vec<String> = AMINO_ACIDS
            .chars()
            .map(|aa| format!("{:.3}", position.get(&aa).unwrap_or(&0.0)))
            .collect();
        writeln!(file, "{}\t{}", i + 1, row.join("\t"))?;
    }

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <fasta_file_path> <output_tsv_path>", args[0]);
        std::process::exit(1);
    }

    let fasta_file = &args[1];
    let output_file = &args[2];

    if !Path::new(fasta_file).exists() {
        eprintln!("Error: File '{}' does not exist.", fasta_file);
        std::process::exit(1);
    }

    match read_fasta(fasta_file) {
        Ok(sequences) => {
            let pwm = calculate_pwm(&sequences);
            if let Err(e) = write_pwm_to_tsv(&pwm, output_file) {
                eprintln!("Error writing to TSV file: {}", e);
            } else {
                println!("PWM written to '{}'", output_file);
            }
        }
        Err(e) => eprintln!("Error reading FASTA file: {}", e),
    }
}
