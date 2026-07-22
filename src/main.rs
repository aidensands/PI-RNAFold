use needletail::parse_fastx_file;
use clap::Parser;
use rayon::prelude::*;
use core::fmt;
use std::ffi::{CStr, CString, c_char};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    filepath: PathBuf,
    #[arg(short, long)]
    cores: usize,
    #[arg(short, long)]
    verbose: bool,
    #[arg(short, long)]
    output: Option<PathBuf>
}

struct SolvedSeq {
    id:String,
    ss:String,
    mfe:f32,
}

impl SolvedSeq {
    fn new(id:String, ss:String, mfe:f32) -> Self {
        Self {id, ss, mfe}
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

fn rfold(seqs: &Vec<(String, String)>, verbose: bool) -> Vec<SolvedSeq> {
    seqs.par_iter()
        .map(|seq| {
            let processed_string = CString::new(seq.0.as_str())
            .expect("Error converting Rust String to CString");
            let mut ss_buffer = vec![0u8; seq.0.as_bytes().len() + 1];
            // RNAFold Thermo Calculations
            let mfe = unsafe {
                vrna_fold(processed_string.as_ptr(), ss_buffer.as_mut_ptr() as *mut c_char)
            };
            let c_stucture: &CStr = unsafe {
                CStr::from_ptr(ss_buffer.as_ptr() as *const c_char)
            };

            // Load Solution
            let ss_string: String = c_stucture.to_string_lossy().into_owned();

            if verbose {
                println!("Thread Folded Seq: {}", seq.1);
            }

            SolvedSeq::new(seq.1.clone(), ss_string, mfe)
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

    let mut sequences: Vec<(String, String)> = Vec::new();

    while let Some(record) = reader.next() {        
        // Processing raw bytes to CString format
        let record= record
        .expect("Recieved Invalid Sequence in fasta");
        let stringrecord: String = String::from_utf8(record.seq().into_owned())
        .expect("Invalid UTF-8 bytes recieved");
        let seqname = str::from_utf8(record.id())
        .expect("Issue Loading Sequence ID");
        sequences.push((stringrecord, seqname.to_owned()));
    }

    let data = rfold(&sequences, args.verbose);

    if let Some(output_path) = args.output {
        let mut writer = csv::Writer::from_path(output_path)
        .expect("Failed to create csv output writer");

        writer.write_record(&["id", "structure", "mfe"])
        .expect("Failed to write header record");

        for seq in data {
            writer.write_record(&[seq.id, seq.ss, seq.mfe.to_string()])
            .expect("Failed to write data record");
        }

        writer.flush().expect("Failed to flush csv writer");

    }

}

