use needletail::{Sequence, parse_fastx_file};
use sourmash::encodings::HashFunctions;
use sourmash::prelude::ToWriter;
use sourmash::{signature::SigsTrait, sketch::minhash::KmerMinHash};
use std::fs;
use std::io::BufRead;
use std::path::PathBuf;
use std::collections::HashMap;
use rayon::prelude::*;
use rand::prelude::*;
use rand::distr::weighted::WeightedIndex;
/// Seed used for all MinHash sketches.
///
/// A fixed seed ensures reproducible hashes across runs.
const TESTING_SEED: u64 = 101;

/// Sketch a single FASTQ file into a [`KmerMinHash`].
///
/// Reads every record in `path`, normalises the sequence, and adds it to a
/// new MinHash sketch parameterised by `scaled` and `ksize`.
pub fn sketch_file(path: &str, scaled: u32, ksize: u32) -> KmerMinHash {
    let mut mh: KmerMinHash = KmerMinHash::new(
        scaled, // scaled size
        ksize,  // k-mer size
        HashFunctions::Murmur64Dna,
        TESTING_SEED,
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

/// Sketch every `.fastq` or `.fastq.gz` file in `fastq_dir`.
///
/// Returns one [`KmerMinHash`] per file in directory-iteration order.
/// Files with other extensions are silently skipped.
pub fn sketch_dir_files(
    fastq_dir: &str,
    scaled: u32,
    ksize: u32,
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
            ));
        }
    }
    sketches
}

/// Return `true` if `path` contains at least one valid k-mer of length `ksize`.
///
/// A k-mer is considered valid when it is fully unambiguous — i.e. it contains
/// no `N` or `n` bases. Files where every read is shorter than `ksize` or
/// entirely ambiguous return `false`.
pub fn check_valid_fastq(path: &str, ksize: u32) -> bool {
    let mut reader = parse_fastx_file(path).expect("valid path/file");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let norm_seq = seqrec.normalize(false);
        let has_valid_kmer = norm_seq.len() >= ksize as usize
            && norm_seq
                .windows(ksize as usize)
                .any(|kmer| !kmer.contains(&b'N') && !kmer.contains(&b'n'));
        if has_valid_kmer {
            return true; // at least one valid line
        }
    }
    false // no valid lines
}

/// Print the path of every invalid FASTQ file in `fastq_dir`.
///
/// A file is invalid if [`check_valid_fastq`] returns `false` for it at the
/// given `ksize`. Files with other extensions are silently skipped.
pub fn validate_fastq_dir(
    fastq_dir: &str,
    ksize: u32,
) {
    let paths = fs::read_dir(fastq_dir).unwrap();
    for path in paths {
        let path = path.unwrap().path();
        let ext = path.extension().and_then(|e| e.to_str());
        if ext == Some("fastq") || ext == Some("fastq.gz") {
            let is_valid = check_valid_fastq(
                path.to_str().expect("missing_path"),
                ksize,
            );
            // print out any paths that are found to be invalid
            if !is_valid {
                println!("{}", path.to_str().expect("missing_path"));
            }
        }
    }
}

/// Merge a slice of sketches into a single [`KmerMinHash`].
///
/// All sketches must have been created with the same `scaled` and `ksize`
/// values, which are also used to initialise the merged result.
pub fn merge_sketches(
    sketches: &Vec<KmerMinHash>,
    scaled: u32,
    ksize: u32,
) -> KmerMinHash {
    let mut merged: KmerMinHash =
        KmerMinHash::new(scaled, ksize, HashFunctions::Murmur64Dna, TESTING_SEED, false, 0);
    for sketch in sketches {
        merged.merge(sketch).expect("error");
    }
    merged
}

/// Return the Jaccard similarity between two sketches.
pub fn compare_sketches(sketch_a: &KmerMinHash, sketch_b: &KmerMinHash) -> f64 {
    sketch_a.similarity(sketch_b, false, false).expect("error")
}

/// Serialise a sketch to a `.sig` file at `path`.
pub fn write_sketch(path: &str, sketch: &KmerMinHash) {
    let file = fs::File::create(path).unwrap(); // create, not open
    let mut writer = std::io::BufWriter::new(file);
    sketch.to_writer(&mut writer).expect("error");
}

/// Deserialise a sketch from a `.sig` file at `path`.
pub fn read_sketch(path: &str) -> KmerMinHash {
    let file = fs::File::open(path).unwrap();
    let reader = std::io::BufReader::new(file);
    KmerMinHash::from_reader(reader).expect("missing")
}

/// Read every `.sig` file in `sketches_dir` and return the sketches.
pub fn read_sketches_from_dir(sketches_dir: &str) -> Vec<KmerMinHash> {
    let paths = fs::read_dir(sketches_dir).unwrap();
    let mut sketches: Vec<KmerMinHash> = Vec::new();
    for path in paths {
        sketches.push(read_sketch(path.unwrap().path().to_str().expect("error")));
    }
    sketches
}
// sketch files in parallel 
pub fn parallel_sketch_files(k: usize, 
    fastq_list: &str, 
    scaled: u32,
    ksize: u32,
    )->Vec<KmerMinHash>{
    // k number of threads
    // let files: Vec<PathBuf> = fs::read_dir(fastq_dir)
    //     .unwrap()
    //     .filter_map(|e| e.ok())
    //     .map(|e| e.path())
    //     .filter(|p| p.is_file())
    //     .collect();
    let files: Vec<PathBuf> = std::io::BufReader::new(std::fs::File::open(fastq_list).unwrap())
        .lines()
        .filter_map(|l| l.ok())
        .map(|l| PathBuf::from(l.trim()))
        .filter(|p| p.is_file())
        .collect();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(k)
        .build()           // returns Result<ThreadPool>, not global
        .unwrap();

    pool.install(|| {
        files.par_iter().map(|path| sketch_file(path.to_str().expect("missing"), scaled, ksize)).collect()
    })
}

pub fn parallel_sketch_files_with_names(k: usize, 
    fastq_list: &str, 
    scaled: u32,
    ksize: u32,
    )->Vec<(KmerMinHash, PathBuf)>{
    // k number of threads
    let files: Vec<PathBuf> = std::io::BufReader::new(std::fs::File::open(fastq_list).unwrap())
        .lines()
        .filter_map(|l| l.ok())
        .map(|l| PathBuf::from(l.trim()))
        .filter(|p| p.is_file())
        .collect();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(k)
        .build()           // returns Result<ThreadPool>, not global
        .unwrap();
        pool.install(|| {
            files
                .par_iter()
                .map(|path| {
                    let sketch = sketch_file(path.to_str().expect("missing"), scaled, ksize);
                    (sketch, path.clone())
                })
                .collect::<Vec<(KmerMinHash, PathBuf)>>()
        })
}


/// Build `n` initial reference sketches from the FASTQ files in `fastq_dir`.
///
/// Files are assigned to partitions via round-robin, then each partition's
/// sketch is written to `sig_dir` as `cluster_sketch_<i>.sig`. The sketches
/// are also returned for immediate use.
pub fn sketch_initial_index(
    fastq_list: &str,
    n: u32,
    scaled: u32,
    ksize: u32,
    sig_dir: &str,
    num_threads:usize  ,
) -> Vec<KmerMinHash> {
    let mut worker_sketches: Vec<KmerMinHash> = Vec::new();
    for _ in 0..n {
        worker_sketches.push(
            KmerMinHash::new(scaled, ksize, HashFunctions::Murmur64Dna, TESTING_SEED, false, 0)
        );
    }
    let sequence_sketches = parallel_sketch_files(num_threads, fastq_list, scaled, ksize);
    for (i, sketch) in sequence_sketches.into_iter().enumerate() {
        let idx: usize = i % n as usize;
        worker_sketches[idx].merge(&sketch).unwrap();
    }
    fs::create_dir_all(sig_dir).expect("could not create sig dir");
    // write out results
    for (i, sketch) in worker_sketches.iter().enumerate(){
        write_sketch(
            format!("{sig_dir}/cluster_sketch_{}.sig", i).as_str(),
            sketch,
        );
    }
    worker_sketches
}

/// Write a slice of sketches to `dir` as `cluster_sketch_<i>.sig` files.
///
/// Creates `dir` if it does not already exist.
pub fn write_sketches_to_dir(sketches: &Vec<KmerMinHash>, dir: &str) {
    println!("{dir}");
    fs::create_dir_all(dir).expect("could not create dir");
    for (i, sketch) in sketches.iter().enumerate() {
        write_sketch(
            &format!("{dir}/cluster_sketch_{}.sig", i),
            sketch,
        );
    }
}

/// Assign incoming FASTQ files to cluster sketches using round-robin rotation.
///
/// Files in `incoming_dir` are sorted and then cycled across `cluster_sketches`
/// by index. If `make_sketch` is `true`, each file is sketched and merged into
/// its assigned cluster before the updated sketches are written to
/// `final_sig_dir`. If `make_sketch` is `false`, only the directory is created.
///
/// Returns a map of filename → cluster index.
pub fn run_round_robin(
    fastq_list: &str,
    make_sketch: bool,
    mut cluster_sketches: Vec<KmerMinHash>,
    scaled: u32,
    ksize: u32,
    final_sig_dir: &str,
) -> HashMap<String, usize> {
    let n = cluster_sketches.len();
    let mut assignments: HashMap<String, usize> = HashMap::new();

    // let mut paths: Vec<_> = fs::read_dir(incoming_dir).unwrap().filter_map(|p| {
    //     let path = p.unwrap().path();
    //     let ext = path.extension().and_then(|e| e.to_str()).map(str::to_owned);
    //     if ext.as_deref() == Some("fastq") || ext.as_deref() == Some("fastq.gz"){
    //         Some(path)
    //     } else{
    //         None
    //     }
    // }).collect();
    let mut files: Vec<PathBuf> = std::io::BufReader::new(std::fs::File::open(fastq_list).unwrap())
        .lines()
        .filter_map(|l| l.ok())
        .map(|l| PathBuf::from(l.trim()))
        .filter(|p| p.is_file())
        .collect();

    files.sort();

    for(i, path) in files.iter().enumerate() {
        let idx = i % n;
        let filename = path.file_name().unwrap().to_str().unwrap().to_string();
        if make_sketch{
            let sketch = sketch_file(path.to_str().unwrap(), scaled, ksize);
            cluster_sketches[idx].merge(&sketch).unwrap();
        }
        assignments.insert(filename, idx);
    }
    if make_sketch{
        write_sketches_to_dir(&cluster_sketches, final_sig_dir);
    }
    else{
        fs::create_dir_all(final_sig_dir).expect("could not create dir");
    }
    assignments
}

/// Find the cluster sketch most similar to a single FASTQ file.
///
/// Sketches `fastq_file_path` once, then compares it against every sketch in
/// `sketches`. Returns a tuple of `(best_index, best_similarity, query_sketch)`.
/// Similarity scores for each cluster are printed to stdout.
pub fn select_most_similar_sketch(
    sketches: &Vec<KmerMinHash>,
    new_sketch: KmerMinHash
) -> (usize, f64, KmerMinHash) {
    // initialize as empty
    let mut most_similar: (usize, f64, KmerMinHash) = (
        0,
        0.00,
        new_sketch,
    );
    for (i, sketch) in sketches.iter().enumerate() {
        let cur_sim = most_similar.2.similarity(sketch, false, false).expect("error");
        if cur_sim > most_similar.1 {
            most_similar.0 = i;
            most_similar.1 = cur_sim;
        }
        println!("Cluster {i} has sim {cur_sim}")
    }
    most_similar
}

/// Assign incoming FASTQ files to cluster sketches by MinHash similarity.
///
/// Each file in `incoming_dir` is sketched and assigned to the most similar
/// cluster via [`select_most_similar_sketch`]. The winning cluster's sketch is
/// then updated by merging the new file's sketch into it, so later assignments
/// reflect the growing clusters. Updated sketches are written to `final_sig_dir`.
///
/// Returns a map of filename → cluster index.
pub fn run_similarity(
    fastq_list: &str,
    mut cluster_sketches: Vec<KmerMinHash>,
    scaled: u32,
    ksize: u32,
    num_threads:usize,
    final_sig_dir: &str,
) -> HashMap<String, usize> {
    let mut assignments: HashMap<String, usize> = HashMap::new();
    let new_sketches = parallel_sketch_files_with_names(num_threads, fastq_list, scaled, ksize);
    for (new_sketch, path) in new_sketches.iter(){
        // let new_sketch = read_sketch(path.to_str().expect("Missing"));
        let(best_idx, _, sketch) = select_most_similar_sketch(
            &cluster_sketches,
            new_sketch.clone()
            
        );
        cluster_sketches[best_idx].merge(&sketch).unwrap();
        let base_name = path.file_prefix().unwrap().to_str().unwrap().to_string();
        assignments.insert(format!("{base_name}.fastq"), best_idx);
    }
    write_sketches_to_dir(&cluster_sketches, final_sig_dir);
    assignments
}

/// Write an assignment map to a CSV file at `path`.
///
/// Entries are sorted by filename before writing. The output format is:
/// ```text
/// filename, cluster
/// sample_a.fastq,0
/// sample_b.fastq,2
/// ```
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

/// Write a comparison of round-robin and similarity assignment counts to a CSV.
///
/// For each strategy, counts how many files were assigned to each cluster and
/// writes one row per cluster. The output format is:
/// ```text
/// strategy, cluster, file_count
/// round_robin,0,12
/// similarity,0,9
/// ```
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

/// Reassign files from an existing assignment CSV using weighted random sampling.
///
/// Reads `existing_assignment_file`, counts how many files were originally
/// assigned to each cluster, then randomly reassigns all files using those
/// counts as sampling weights. Files are shuffled before assignment to avoid
/// ordering bias. The output directory `dir` is created if it does not exist.
///
/// Returns a map of filename → new cluster index.
pub fn run_weighted_random_assignment(
    existing_assignment_file: &str,
    dir: &str,
) -> HashMap<String, usize> {
    fs::create_dir_all(dir).expect("could not create dir");
    let file = std::fs::File::open(existing_assignment_file);
    let reader = std::io::BufReader::new(file.expect("here"));
    let mut files: Vec<String> = Vec::new();
    let mut original_assignments: Vec<u128> = Vec::new();
    for (i, line) in reader.lines().enumerate() {
            // skip the first line
            if i == 0 {
                continue;
            }
            else{
                let a = line.expect("here");
                let line_val = a.split(",");
                let lines: Vec<&str> = line_val.collect();
                if let Some(f_name) = lines.first() {
                    files.push(f_name.to_string());
                }
                if let Some(assignment) = lines.last() {
                    let value: u128 = assignment.trim().parse().expect("not a valid number");
                    original_assignments.push(value);
                }
            }
    }
    let mut m: HashMap<u128, usize> = HashMap::new();
    let mut choices: Vec<u128> = Vec::new();
    let mut weights: Vec<usize> = Vec::new();
    let mut assignments: HashMap<String, usize> = HashMap::new();
    // count original assignments per cluster to use as sampling weights
    for x in original_assignments {
        *m.entry(x).or_default() += 1;
    }
    for x in m{
        choices.push(x.0);
        weights.push(x.1)
    }

    // sample a cluster for each file using the weighted distribution
    let dist = WeightedIndex::new(&weights).unwrap();
    let mut rng = rand::rng();
    files.shuffle(& mut rng);
    for path in files.iter() {
        path.to_string();
        let idx = choices[dist.sample(&mut rng)] as usize;
        assignments.insert(path.to_string(), idx);
    }
    return assignments;
}
