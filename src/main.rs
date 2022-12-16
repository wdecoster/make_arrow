use arrow::datatypes::{DataType, Field, Schema};
use clap::Parser;
use log::info;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::{bam, bam::Read, htslib};
use std::fs::File;
use std::path::PathBuf;
use std::sync::Arc;
use unzip_n::unzip_n;

use arrow::{
    self,
    array::{Float64Array, UInt64Array, UInt8Array},
    ipc::writer::FileWriter,
    record_batch::RecordBatch,
};

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[command(author, version, about="Tool to extract metrics from cram or bam to an arrow file", long_about = None)]
struct Cli {
    /// cram or bam file to check
    #[arg(value_parser)]
    input: String,

    /// Number of parallel decompression threads to use
    #[arg(short, long, value_parser, default_value_t = 4)]
    threads: usize,

    /// Output file name
    #[arg(short, long, value_parser, default_value_t = String::from("read_metrics.arrow"))]
    output: String,
}

fn main() {
    env_logger::init();
    let args = Cli::parse();
    is_file(&args.input).unwrap_or_else(|_| panic!("Input file {} is invalid", args.input));
    info!("Collected arguments");
    extract(&args.input, args.output, args.threads)
}

fn is_file(pathname: &str) -> Result<(), String> {
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}

// -qualities
// -aligned qualities

pub fn extract(bam_path: &String, output_path: String, threads: usize) {
    let mut bam = bam::Reader::from_path(&bam_path).expect("Error opening BAM.\n");
    bam.set_threads(threads)
        .expect("Failure setting decompression threads");
    unzip_n!(4);
    let (lengths, aligned_lengths, identities, mapqs): (Vec<u64>, Vec<u64>, Vec<f64>, Vec<u8>) =
        bam.rc_records()
            .map(|r| r.expect("Failure parsing Bam file"))
            .filter(|read| read.flags() & (htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY) as u16 == 0)
            .map(|read| {
                (
                    read.seq_len() as u64,
                    (read.reference_end() - read.reference_start()) as u64,
                    gap_compressed_identity(&read) * 100.0,
                    read.mapq(),
                )
            })
            .unzip_n_vec();
    save_as_arrow(output_path, lengths, aligned_lengths, identities, mapqs);
}

pub fn save_as_arrow(
    filename: String,
    lengths: Vec<u64>,
    aligned_lengths: Vec<u64>,
    identities: Vec<f64>,
    mapqs: Vec<u8>,
) {
    let identities_array = Arc::new(Float64Array::from(identities)) as _;
    let lengths_array = Arc::new(UInt64Array::from(lengths)) as _;
    let aligned_lengths_array = Arc::new(UInt64Array::from(aligned_lengths)) as _;
    let mapqs_array = Arc::new(UInt8Array::from(mapqs)) as _;
    let batch = RecordBatch::try_from_iter([
        ("identities", identities_array),
        ("lengths", lengths_array),
        ("aligned_lengths", aligned_lengths_array),
        ("mapQ", mapqs_array),
    ])
    .unwrap();

    let schema = Schema::new(vec![
        Field::new("identities", DataType::Float64, false),
        Field::new("lengths", DataType::UInt64, false),
        Field::new("aligned_lengths", DataType::UInt64, false),
        Field::new("mapQ", DataType::UInt8, false),
    ]);
    let buffer = File::create(filename).expect("create arrow file error");

    let mut writer = FileWriter::try_new(buffer, &schema).expect("create arrow file writer error");

    writer.write(&batch).expect("write arrow batch error");
    writer.finish().expect("finish write arrow error");
}

/// Calculates the gap-compressed identity
/// based on https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
/// recent minimap2 version have that as the de tag
/// if that is not present it is calculated from CIGAR and NM
fn gap_compressed_identity(record: &std::rc::Rc<rust_htslib::bam::Record>) -> f64 {
    match get_de_tag(record) {
        Some(v) => v as f64,
        None => {
            let mut matches = 0;
            let mut gap_size = 0;
            let mut gap_count = 0;
            for entry in record.cigar().iter() {
                match entry {
                    Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                        matches += *len;
                    }
                    Cigar::Del(len) | Cigar::Ins(len) => {
                        gap_size += *len;
                        gap_count += 1;
                    }
                    _ => (),
                }
            }
            1.0 - ((get_nm_tag(record) - gap_size + gap_count) as f64
                / (matches + gap_count) as f64)
        }
    }
}

fn get_nm_tag(record: &bam::Record) -> u32 {
    match record.aux(b"NM") {
        Ok(value) => match value {
            Aux::U8(v) => u32::from(v),
            Aux::U16(v) => u32::from(v),
            Aux::U32(v) => v,
            _ => panic!("Unexpected type of Aux {:?}", value),
        },
        Err(_e) => panic!("Unexpected result while trying to access the NM tag"),
    }
}

/// Get the de:f tag from minimap2, which is the gap compressed sequence divergence
/// Which is converted into percent identity with 100 * (1 - de)
/// This tag can be absent if the aligner version is not quite recent
fn get_de_tag(record: &bam::Record) -> Option<f32> {
    match record.aux(b"de") {
        Ok(value) => match value {
            Aux::Float(v) => Some(100.0 * (1.0 - v)),
            _ => panic!("Unexpected type of Aux {:?}", value),
        },
        Err(_e) => None,
    }
}

#[cfg(test)]
#[ctor::ctor]
fn init() {
    env_logger::init();
}

#[test]
fn verify_app() {
    use clap::CommandFactory;
    Cli::command().debug_assert()
}

#[test]
fn test_extract() {
    extract(
        &"test-data/small-test-phased.bam".to_string(),
        "test.arrow".to_string(),
        8,
    )
}
