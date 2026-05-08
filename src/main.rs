//! A command-line tool for sketching genomic FASTQ files and assigning
//! incoming samples to reference index partitions.
//!
//! Sketching is performed with MinHash via [SourMash]. Samples can be
//! assigned using round-robin rotation, similarity scoring, or an
//! asymmetrical redistribution strategy.
//!
//! [SourMash]: https://sourmash.readthedocs.io

mod sketch;

use clap::{Parser, Subcommand};
use sketch::{
    compare_sketches, make_initial_sketch, merge_sketches, read_sketch, read_sketches_from_dir,
    select_most_similar_sketch, sketch_dir_files, write_sketch,
    run_round_robin, run_similarity, write_assignments, write_results,
    validate_fastq_dir, run_asymmetrical_assignment,
};
use std::fs;
use sourmash::signature::SigsTrait;

/// Genomic sketch assignment tool.
#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Sketch all FASTQ files in a directory and report basic statistics.
    ///
    /// Merges all sketches into a single `merged.sig` file written to the
    /// current directory, then prints the hash count and similarity between
    /// the first two sketches.
    SketchFiles {
        /// Directory containing input FASTQ files.
        #[arg(long, default_value_t = String::from("fastq_files"), short = 'd')]
        fastq_dir: String,

        /// Scaled factor for MinHash sketching. Lower values retain more hashes.
        #[arg(long, default_value_t = 1000, short = 's')]
        scaled: u32,

        /// k-mer length used for hashing.
        #[arg(long, default_value_t = 21, short = 'k')]
        ksize: u32,
    },

    /// Build a reference sketch index from a directory of FASTQ files.
    ///
    /// Partitions files into `num_index` groups and writes one sketch per
    /// partition to `sig_dir`.
    BuildIndex {
        /// Directory containing input FASTQ files.
        #[arg(long, default_value_t = String::from("fastq_files"), short = 'd')]
        fastq_dir: String,

        /// Output directory for the generated index sketch files.
        #[arg(long, default_value_t = String::from("initial_index"), short = 'o')]
        sig_dir: String,

        /// Number of index partitions to create.
        #[arg(long, default_value_t = 5, short = 'n')]
        num_index: u32,

        /// Scaled factor for MinHash sketching. Lower values retain more hashes.
        #[arg(long, default_value_t = 1000, short = 's')]
        scaled: u32,

        /// k-mer length used for hashing.
        #[arg(long, default_value_t = 21, short = 'k')]
        ksize: u32,
    },

    /// Find the reference index sketch most similar to a single FASTQ file.
    ///
    /// Sketches `fastq_file_path` and compares it against all sketches in
    /// `sig_dir`, printing the name and score of the closest match.
    ///
    /// The `scaled` and `ksize` values must match those used when building
    /// the index.
    FindMostSimilarIndex {
        /// Path to the query FASTQ file.
        fastq_file_path: String,

        /// Directory containing reference index sketch files.
        #[arg(long, default_value_t = String::from("initial_index"), short = 'd')]
        sig_dir: String,

        /// Scaled factor used when the index was built.
        #[arg(long, default_value_t = 1000, short = 's')]
        scaled: u32,

        /// k-mer length used when the index was built.
        #[arg(long, default_value_t = 21, short = 'k')]
        ksize: u32,
    },

    /// Assign incoming FASTQ files to index partitions using round-robin rotation.
    ///
    /// Files are cycled evenly across partitions regardless of sequence content.
    /// Pass `--make-sketch` to sketch incoming files before assignment.
    /// Writes assignments to `<output>/round_robin_assignments.csv`.
    RunRoundRobin {
        /// Directory containing incoming FASTQ files to assign.
        #[arg(long, default_value_t = String::from("incoming_fastq"), short = 'd')]
        incoming_dir: String,

        /// Directory containing the reference sketch index.
        #[arg(long, default_value_t = String::from("initial_index"), short = 'e')]
        sig_dir: String,

        /// Output directory for updated sketches and the assignment CSV.
        #[arg(long, default_value_t = String::from("results/round_robin_sketches"), short = 'o')]
        output: String,

        /// Scaled factor for MinHash sketching.
        #[arg(long, default_value_t = 1000, short = 's')]
        scaled: u32,

        /// k-mer length used for hashing.
        #[arg(long, default_value_t = 21, short = 'k')]
        ksize: u32,

        /// Sketch incoming files before running the assignment.
        #[arg(long, default_value_t = false, short = 'm')]
        make_sketch: bool,
    },

    /// Assign incoming FASTQ files to index partitions by MinHash similarity.
    ///
    /// Each file is sketched and compared against all reference sketches; it
    /// is assigned to the closest match. Writes assignments to
    /// `<output>/similarity_assignments.csv`.
    RunSimilarity {
        /// Directory containing incoming FASTQ files to assign.
        #[arg(long, default_value_t = String::from("incoming_fastq"), short = 'd')]
        incoming_dir: String,

        /// Directory containing the reference sketch index.
        #[arg(long, default_value_t = String::from("initial_index"), short = 'e')]
        sig_dir: String,

        /// Output directory for updated sketches and the assignment CSV.
        #[arg(long, default_value_t = String::from("results/sim_final_sketches"), short = 'o')]
        output: String,

        /// Scaled factor for MinHash sketching.
        #[arg(long, default_value_t = 1000, short = 's')]
        scaled: u32,

        /// k-mer length used for hashing.
        #[arg(long, default_value_t = 21, short = 'k')]
        ksize: u32,
    },

    /// Redistribute samples from an existing assignment file asymmetrically.
    ///
    /// Reads a prior assignment CSV and reassigns samples without re-sketching.
    /// Writes results to `<output>/asymmetrical_assignments.csv`.
    RunAsymmetrical {
        /// Path to an existing assignment CSV from a previous run.
        #[arg(long, default_value_t = String::from("existing_assignment_file"), short = 'a')]
        existing_assignment_file: String,

        /// Output directory for the new assignment CSV.
        #[arg(long, default_value_t = String::from("results/asymmetrical_final_sketches"), short = 'o')]
        output: String,
    },

    /// Check a directory of FASTQ files for sketching compatibility.
    ///
    /// Prints the paths of any files that cannot be sketched at the given
    /// k-mer length.
    ValidateFastqDir {
        /// Directory of FASTQ files to validate.
        #[arg(long, default_value_t = String::from("fastq_files"), short = 'd')]
        fastq_dir: String,

        /// k-mer length to validate against.
        #[arg(long, default_value_t = 21, short = 'k')]
        ksize: u32,
    },

    /// Run both round-robin and similarity assignments and compare the results.
    ///
    /// Executes [`RunRoundRobin`] and [`RunSimilarity`] on the same input,
    /// writes individual assignment CSVs, and produces a side-by-side
    /// `comparison.csv` in `output_dir`.
    ///
    /// The index is read from disk twice to avoid cloning large sketch data.
    RunExperiment {
        /// Directory containing incoming FASTQ files.
        #[arg(long, default_value_t = String::from("incoming_fastq"), short = 'd')]
        incoming_dir: String,

        /// Directory containing the reference sketch index.
        #[arg(long, default_value_t = String::from("initial_index"), short = 'e')]
        sig_dir: String,

        /// Root output directory. Sub-directories are created automatically.
        #[arg(long, default_value_t = String::from("results"), short = 'o')]
        output_dir: String,

        /// Scaled factor for MinHash sketching.
        #[arg(long, default_value_t = 1000, short = 's')]
        scaled: u32,

        /// k-mer length used for hashing.
        #[arg(long, default_value_t = 21, short = 'k')]
        ksize: u32,

        /// Number of reference index partitions.
        #[arg(long, default_value_t = 5, short = 'n')]
        num_index: u32,

        /// Sketch incoming files before running assignments.
        #[arg(long, default_value_t = false, short = 'm')]
        make_sketch: bool,
    },
}

fn main() {
    let args = Args::parse();
    match args.command {
        Command::SketchFiles { fastq_dir, scaled, ksize } => {
            println!("Reading from {a}\nSketching with SourMash!", a = fastq_dir);
            let sketches = sketch_dir_files(&fastq_dir, scaled, ksize);
            let merged = merge_sketches(&sketches, scaled, ksize);
            let filename = "merged.sig";
            println!("Merged sketch contains {} hashes", merged.size());
            write_sketch(filename, &merged);
            let read_merged = read_sketch(filename);
            println!(
                "Read the merged sketch result contains {} hashes",
                read_merged.size()
            );
            let res = compare_sketches(&sketches[0], &sketches[1]);
            println!("similarity {}", res);
        }
        Command::BuildIndex { fastq_dir, sig_dir, num_index, scaled, ksize } => {
            println!("Building index from {fastq_dir} saving to {sig_dir}");
            make_initial_sketch(&fastq_dir, num_index, scaled, ksize, &sig_dir);
        }
        Command::FindMostSimilarIndex { fastq_file_path, sig_dir, scaled, ksize } => {
            let sketches = read_sketches_from_dir(&sig_dir);
            let most_similar_sketch =
                select_most_similar_sketch(&sketches, &fastq_file_path, scaled, ksize);
            println!(
                "Most similar sketch {}, {}",
                most_similar_sketch.0, most_similar_sketch.1
            );
        }
        Command::ValidateFastqDir { fastq_dir, ksize } => {
            println!("Checking {fastq_dir} for invalid fastq files");
            println!("The following files are invalid at k-mer length {}:", ksize);
            validate_fastq_dir(&fastq_dir, ksize);
        }
        Command::RunRoundRobin { incoming_dir, sig_dir, output, scaled, ksize, make_sketch } => {
            println!("Running round robin from {incoming_dir}");
            let sketches = read_sketches_from_dir(&sig_dir);
            let assignments =
                run_round_robin(&incoming_dir, make_sketch, sketches, scaled, ksize, &output);
            write_assignments(&format!("{output}/round_robin_assignments.csv"), &assignments);
            println!("Done. Assignments written to {output}");
        }
        Command::RunSimilarity { incoming_dir, sig_dir, output, scaled, ksize } => {
            println!("Running similarity assignment from {incoming_dir}");
            let sketches = read_sketches_from_dir(&sig_dir);
            let assignments =
                run_similarity(&incoming_dir, sketches, scaled, ksize, &output);
            write_assignments(&format!("{output}/similarity_assignments.csv"), &assignments);
            println!("Done. Assignments written to {output}");
        }
        Command::RunAsymmetrical { existing_assignment_file, output } => {
            let assignments =
                run_asymmetrical_assignment(&existing_assignment_file, &output);
            write_assignments(&format!("{output}/asymmetrical_assignments.csv"), &assignments);
        }
        Command::RunExperiment {
            incoming_dir, sig_dir, output_dir, scaled, ksize, num_index, make_sketch,
        } => {
            println!("Running full experiment from {incoming_dir}");

            // Read the index from disk twice rather than cloning — sketch data can be large.
            let rr_sketches = read_sketches_from_dir(&sig_dir);
            let sim_sketches = read_sketches_from_dir(&sig_dir);

            let rr_assignments = run_round_robin(
                &incoming_dir,
                make_sketch,
                rr_sketches,
                scaled,
                ksize,
                &format!("{output_dir}/rr_final_sketches"),
            );
            let sim_assignments = run_similarity(
                &incoming_dir,
                sim_sketches,
                scaled,
                ksize,
                &format!("{output_dir}/sim_final_sketches"),
            );

            fs::create_dir_all(&output_dir).unwrap();
            write_assignments(
                &format!("{output_dir}/round_robin_assignments.csv"),
                &rr_assignments,
            );
            write_assignments(
                &format!("{output_dir}/similarity_assignments.csv"),
                &sim_assignments,
            );
            write_results(
                &format!("{output_dir}/comparison.csv"),
                &rr_assignments,
                &sim_assignments,
                num_index as usize,
            );
            println!("Done. Results written to {output_dir}");
        }
    }
}
