use needletail::parse_fastx_file;
use clap::Parser;
use rayon::prelude::*;
use core::fmt;
use std::ffi::{CStr, CString, c_char};
use std::path::PathBuf;
use std::fmt::{write};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    filepath: PathBuf,
    #[arg(short, long)]
    cores: usize,
    #[arg(short, long)]
    verbose: bool,
}

struct SolvedSeq {
    ss:String,
    mfe:f32,
}

impl SolvedSeq {
    fn new(ss:String, mfe:f32) -> Self {
        Self {ss, mfe}
    }
}

impl fmt::Display for SolvedSeq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "(db: {}, mfe:{})", self.ss, self.mfe)
    }
}

unsafe extern "C" {
    unsafe fn vrna_fold(seq: *const c_char, structure: *mut c_char) -> f32;
}

fn rfold(seqs: &Vec<String>) -> Vec<SolvedSeq> {
    seqs.par_iter()
        .map(|seq| {
            let processed_string = CString::new(seq.as_str())
            .expect("Error converting Rust String to CString");
            let mut ss_buffer = vec![0u8; seq.as_bytes().len() + 1];
            // RNAFold Thermo Calculations
            let mfe = unsafe {
                vrna_fold(processed_string.as_ptr(), ss_buffer.as_mut_ptr() as *mut c_char)
            };
            let c_stucture: &CStr = unsafe {
                CStr::from_ptr(ss_buffer.as_ptr() as *const c_char)
            };

            // Load Solution
            let ss_string: String = c_stucture.to_string_lossy().into_owned();
            println!("Folded a seq");
            SolvedSeq::new(ss_string, mfe)
        })
        .collect()
}

fn main() {
    let args = Args::parse();

    // Configure the ThreadPool
    rayon::ThreadPoolBuilder::new()
    .num_threads(args.cores)
    .build_global()
    .expect("ThreadPool Initialization Error");

    if args.verbose {
        println!("Starting RNAFold with {} cores", args.cores)
    }

    let path: &PathBuf = &args.filepath;
    let mut reader = parse_fastx_file(path)
    .expect("NeedleTail Reader Recieved Ivalid Filepath");

    let mut sequences: Vec<String> = Vec::new();

    while let Some(record) = reader.next() {        
        // Processing raw bytes to CString format
        let record= record
        .expect("Recieved Invalid Sequence in fasta");
        let stringrecord: String = String::from_utf8(record.seq().into_owned())
        .expect("Invalid UTF-8 bytes recieved");
        sequences.push(stringrecord);
    }

    let data = rfold(&sequences);

    for seq in data{
        print!("{} \n", seq);
    }

}
