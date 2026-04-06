use needletail::{Sequence, parse_fastx_file};
use rand::random;
use sourmash::encodings::HashFunctions;
use sourmash::prelude::ToWriter;
use sourmash::{signature::SigsTrait, sketch::minhash::KmerMinHash};
use std::fs;
use std::collections::HashMap;

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

pub fn read_sketches_from_dir(sketches_dir: &str) -> Vec<KmerMinHash> {
    let paths = fs::read_dir(sketches_dir).unwrap();
    let mut sketches: Vec<KmerMinHash> = Vec::new();
    for path in paths {
        sketches.push(read_sketch(path.unwrap().path().to_str().expect("error")));
    }
    sketches
}

// make n initial sketches, just using a round-robin for load balancing
pub fn make_initial_sketch(
    fastq_dir: &str,
    n: u32,
    scaled: u32,
    ksize: u32,
    sig_dir: &str,
) -> Vec<KmerMinHash> {
    let mut sketches: Vec<KmerMinHash> = Vec::new();
    for _ in 0..n {
        let seed: u64 = random();
        sketches.push(
            KmerMinHash::new(scaled, ksize, HashFunctions::Murmur64Dna, seed, false, 0)
        );
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
                Some(sketches[idx].seed()),
            );
            sketches[idx].merge(&file_sketch).unwrap();
        }
    }
    fs::create_dir_all(sig_dir).expect("could not create sig dir");
    // write out results
    for (i, sketch) in sketches.iter().enumerate(){
        write_sketch(
            format!("{sig_dir}/cluster_sketch_{}.sig", i).as_str(),
            sketch,
        );
    }
    sketches
}

pub fn write_sketches_to_dir(sketches: &Vec<KmerMinHash>, dir: &str) {
    fs::create_dir_all(dir).expect("could not create dir");
    for (i, sketch) in sketches.iter().enumerate() {
        write_sketch(
            &format!("{dir}/cluster_sketch_{}.sig", i),
            sketch,
        );
    }
}

pub fn run_round_robin(
    incoming_dir: &str,
    mut cluster_sketches: Vec<KmerMinHash>,
    scaled: u32,
    ksize: u32,
    final_sig_dir: &str,
) -> HashMap<String, usize> {
    let n = cluster_sketches.len();
    let mut assignments: HashMap<String, usize> = HashMap::new();

    let mut paths: Vec<_> = fs::read_dir(incoming_dir).unwrap().filter_map(|p| {
        let path = p.unwrap().path();
        let ext = path.extension().and_then(|e| e.to_str()).map(str::to_owned);
        if ext.as_deref() == Some("fastq") || ext.as_deref() == Some("fastq.gz"){
            Some(path)
        } else{
            None
        }
    }).collect();

    paths.sort();

    for(i, path) in paths.iter().enumerate() {
        let idx = i % n;
        let seed = cluster_sketches[idx].seed();
        let sketch = sketch_file(path.to_str().unwrap(), scaled, ksize, Some(seed));
        cluster_sketches[idx].merge(&sketch).unwrap();
        let filename = path.file_name().unwrap().to_str().unwrap().to_string();
        assignments.insert(filename, idx);
    }
    write_sketches_to_dir(&cluster_sketches, final_sig_dir);
    assignments
}



pub fn select_most_similar_sketch(
    sketches: &Vec<KmerMinHash>,
    fastq_file_path: &str,
    scaled: u32,
    ksize: u32,
) -> (usize, f64, KmerMinHash) {
    // initialize as empty
    let mut most_similar: (usize, f64, KmerMinHash) = (
        0,
        0.00,
        KmerMinHash::new(0, 0, HashFunctions::Murmur64Dna, 0, false, 0),
    );
    for (i, sketch) in sketches.iter().enumerate() {
        // use seed of current sketch //
        let new_sketch = sketch_file(fastq_file_path, scaled, ksize, Some(sketch.seed()));
        let cur_sim = new_sketch.similarity(sketch, false, false).expect("error");
        if cur_sim > most_similar.1 {
            most_similar = (i, cur_sim, new_sketch);
        }
        println!("Cluster {i} has sim {cur_sim}")
    }
    most_similar
}

pub fn run_similarity(
    incoming_dir: &str,
    mut cluster_sketches: Vec<KmerMinHash>,
    scaled: u32,
    ksize: u32,
    final_sig_dir: &str,
) -> HashMap<String, usize> {
    let mut assignments: HashMap<String, usize> = HashMap::new();

    let paths: Vec<_> = fs::read_dir(incoming_dir).unwrap().filter_map(|p|{
        let path = p.unwrap().path();
        let ext = path.extension().and_then(|e| e.to_str()).map(str::to_owned);
        if ext.as_deref() == Some("fastq") || ext.as_deref() == Some("fastq.gz"){
            Some(path)
        }else{
            None
        }
    }).collect();

    for path in paths.iter(){
        let(best_idx, _, sketch) = select_most_similar_sketch(
            &cluster_sketches, 
            path.to_str().unwrap(),
            scaled,
            ksize
        );
        cluster_sketches[best_idx].merge(&sketch).unwrap();
        let filename = path.file_name().unwrap().to_str().unwrap().to_string();
        assignments.insert(filename, best_idx);
    }
    write_sketches_to_dir(&cluster_sketches, final_sig_dir);
    assignments
}

pub fn write_assignments(path: &str, assignments: &HashMap<String, usize>) {
    use std::io::Write;

    let mut entries: Vec<_> = assignments.iter().collect();
    entries.sort_by(|a, b| a.0.cmp(b.0));

    let mut file = fs::File::create(path).unwrap();
    writeln!(file, "filename, cluster").unwrap();
    for(filename, cluster) in &entries {
        writeln!(file, "{},{}", filename, cluster).unwrap();
    }
}

pub fn write_results(
    path: &str,
    round_robin: &HashMap<String, usize>,
    similarity: &HashMap<String, usize>,
    n: usize,
){
    use std::io::Write;

    let mut rr_counts = vec![0usize; n];
    for cluster in round_robin.values(){
        rr_counts[*cluster] += 1;
    }

    let mut sim_counts = vec![0usize; n];
    for cluster in similarity.values() {
        sim_counts[*cluster] += 1;
    }

    let mut file = fs::File::create(path).unwrap();
    writeln!(file, "strategy, cluster, file_count").unwrap();
    for (i, count) in rr_counts.iter().enumerate(){
        writeln!(file, "round_robin,{},{}", i, count).unwrap();
    }
    for (i, count) in sim_counts.iter().enumerate() {
        writeln!(file, "similarity,{},{}", i, count).unwrap();
    }
}

/*
pub fn load_ballance_new_fastq_files(
    fastq_dir: &str,
    n: u32,
    scaled: u32,
    ksize: u32,
    sig_dir: &str,
) {
    let mut sketches = read_sketches_from_dir(&sig_dir);
    let paths = fs::read_dir(fastq_dir).unwrap();
    for path in paths {
        let most_similar = select_most_similar_sketch(
            &sketches,
            path.unwrap().path().to_str().expect("error"),
            scaled,
            ksize,
        );
        sketches[most_similar.0]
            .merge(&most_similar.2)
            .expect("error")
    }
    for i in 0..n {
        write_sketch(
            format!("{sig_dir}/cluster_sketch_{}.sig", i).as_str(),
            &sketches[i as usize],
        );
    }
}
*/

