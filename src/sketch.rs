use sourmash::prelude::ToWriter;
use sourmash::{signature::SigsTrait, sketch::minhash::KmerMinHash};
use sourmash::encodings::HashFunctions;
use needletail::{Sequence, parse_fastx_file};
use std::fs;

pub fn sketch_file(path: &str) -> KmerMinHash {
    let mut mh: KmerMinHash = KmerMinHash::new(
        1000,    // scaled size
        20,   // k-mer size
        HashFunctions::Murmur64Dna,
        42,   
        false, // track abundance
        0, // if 0 use scaled
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

pub fn sketch_dir_files(fastq_dir:&str)->Vec<KmerMinHash>{
    let paths = fs::read_dir(fastq_dir).unwrap();
    let mut sketches: Vec<KmerMinHash> = Vec::new();
    for path in paths {
        sketches.push( sketch_file(path.unwrap().path().to_str().expect("missing_path")));
    }
    return sketches;
}

pub fn merge_sketches(sketches:&Vec<KmerMinHash>)->KmerMinHash{
    let mut merged: KmerMinHash = KmerMinHash::new(
        1000,    // num=0 → use scaled mode
        20,   // ksize
        HashFunctions::Murmur64Dna,
        42,   // seed
        false, // track abundance
        0, // scaled
    );
    for sketch in sketches{
        merged.merge(sketch).expect("error");
    }
    return  merged;
}

pub fn compare_sketches(sketch_a:&KmerMinHash, sketch_b:&KmerMinHash)->f64{
    return sketch_a.similarity(sketch_b, false, false).expect("error");
}


pub fn write_sketch(path:&str, sketch:&KmerMinHash){
    let file = fs::File::create(path).unwrap();  // create, not open
    let mut writer = std::io::BufWriter::new(file);
    sketch.to_writer(& mut writer).expect("error");
    
}
pub fn read_sketch(path:&str)->KmerMinHash{
    let file = fs::File::open(path).unwrap();
    let reader = std::io::BufReader::new(file);
    let mh = KmerMinHash::from_reader(reader).expect("missing");
    return  mh;
}
