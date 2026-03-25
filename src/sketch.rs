use needletail::{Sequence, parse_fastx_file};
use rand::random;
use sourmash::encodings::HashFunctions;
use sourmash::prelude::ToWriter;
use sourmash::{signature::SigsTrait, sketch::minhash::KmerMinHash};
use std::fs;

// wrapper for sketching an entire fastQ file
pub fn sketch_file(path: &str, scaled: u32, ksize: u32, seed: Option<u64>) -> KmerMinHash {
    let seed = seed.unwrap_or(random());
    let mut mh: KmerMinHash = KmerMinHash::new(
        scaled, // scaled size
        ksize,  // k-mer size
        HashFunctions::Murmur64Dna,
        seed,
        false, // track abundance
        0,     // if 0 use scaled
    );

    let mut reader = parse_fastx_file(path).expect("valid path/file");

    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let norm_seq = seqrec.normalize(false);
        mh.add_sequence(&norm_seq, true).unwrap();
    }
    println!("Sketch  contains {} hashes", mh.size());
    mh
}

// Sketch all fastQ files in a directory
pub fn sketch_dir_files(
    fastq_dir: &str,
    scaled: u32,
    ksize: u32,
    seed: Option<u64>,
) -> Vec<KmerMinHash> {
    let paths = fs::read_dir(fastq_dir).unwrap();
    let mut sketches: Vec<KmerMinHash> = Vec::new();
    for path in paths {
        let path = path.unwrap().path();
        let ext = path.extension().and_then(|e| e.to_str());
        if ext == Some("fastq") || ext == Some("fastq.gz") {
            sketches.push(sketch_file(
                path.to_str().expect("missing_path"),
                scaled,
                ksize,
                seed,
            ));
        }
    }
    sketches
}

// merge a vector of sketches into a single representation
pub fn merge_sketches(
    sketches: &Vec<KmerMinHash>,
    scaled: u32,
    ksize: u32,
    seed: u64,
) -> KmerMinHash {
    // seed must be provided and match the merge
    let mut merged: KmerMinHash =
        KmerMinHash::new(scaled, ksize, HashFunctions::Murmur64Dna, seed, false, 0);
    for sketch in sketches {
        merged.merge(sketch).expect("error");
    }
    merged
}
// super loose wrapper for getting sketch similarly, only here for the sake of completeness
pub fn compare_sketches(sketch_a: &KmerMinHash, sketch_b: &KmerMinHash) -> f64 {
    sketch_a.similarity(sketch_b, false, false).expect("error")
}

// write a sketch to a sig file
pub fn write_sketch(path: &str, sketch: &KmerMinHash) {
    let file = fs::File::create(path).unwrap(); // create, not open
    let mut writer = std::io::BufWriter::new(file);
    sketch.to_writer(&mut writer).expect("error");
}
// read a sketch from a sig file
pub fn read_sketch(path: &str) -> KmerMinHash {
    let file = fs::File::open(path).unwrap();
    let reader = std::io::BufReader::new(file);
    KmerMinHash::from_reader(reader).expect("missing")
}

// make n initial sketches, just using a round-robin for load balancing
pub fn make_initial_sketch(
    fastq_dir: &str,
    n: u32,
    scaled: u32,
    ksize: u32,
    sig_dir: &str,
) -> Vec<(u64, KmerMinHash)> {
    let mut sketches: Vec<(u64, KmerMinHash)> = Vec::new();
    for _ in 0..n {
        let seed: u64 = random();
        sketches.push((
            seed,
            KmerMinHash::new(scaled, ksize, HashFunctions::Murmur64Dna, seed, false, 0),
        ));
    }
    let paths = fs::read_dir(fastq_dir).unwrap();
    for (i, path) in paths.enumerate() {
        let idx: usize = i % n as usize;
        let path = path.unwrap().path();
        let ext = path.extension().and_then(|e| e.to_str());
        if ext == Some("fastq") || ext == Some("fastq.gz") {
            let file_sketch = sketch_file(
                path.to_str().expect("missing_path"),
                scaled,
                ksize,
                Some(sketches[idx].0),
            );
            sketches[idx].1.merge(&file_sketch).unwrap();
        }
    }
    fs::create_dir_all(sig_dir).expect("could not create sig dir");
    // write out results
    for i in 0..n {
        write_sketch(
            format!("{sig_dir}/cluster_sketch_{}.sig", i).as_str(),
            &sketches[i as usize].1,
        );
    }
    sketches
}
