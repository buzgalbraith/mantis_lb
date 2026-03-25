
mod sketch;
use sketch::{read_sketch, write_sketch, sketch_dir_files, merge_sketches, compare_sketches};
use clap::Parser;
use sourmash::signature::SigsTrait;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Directory with FastQ files 
    #[arg(short='d', long)]
    fastq_dir: String,
}



fn main() {
    let args = Args::parse();
    println!("Reading from {a}\nSketching with SourMash!", a=args.fastq_dir);
    let sketches = sketch_dir_files(&args.fastq_dir);
    let merged = merge_sketches(&sketches);
    let filename = "merged.sig";
    println!("Merged sketch contains {} hashes", merged.size());
    write_sketch(filename, &merged);
    let read_merged = read_sketch(filename);
    println!("Read the merged sketch result contains {} hashes", read_merged.size());
    let res = compare_sketches(&sketches[0], &sketches[1]); 
    println!("similarity {}", res);
    }
