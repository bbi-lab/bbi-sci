#![allow(unused_parens)]

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

//use std::error::Error;
use std::path::Path;
use std::fs::File;
// use std::io::Read;
use std::io::BufWriter;
use std::io::Write;
use std::str;
use std::collections::{HashMap, HashSet};

use process_hashes::barcode_utils;
extern crate clap;
use clap::{Arg, Command};
use seq_io::fastq;
use seq_io::fastq::Record;
use sprs;
use sprs::CsMat;
use detect_compression;


// url: https://www.reddit.com/r/rust/comments/jv3q3e/how_to_select_between_reading_from_a_file_and/

fn process_fastq_file<R: std::io::Read>(hash_edit_distance: usize, hash_whitelist: &Vec<HashMap<String, String>>, cells: &mut HashSet<String>, hash_counts: &mut Vec<u64>, hashdict: &mut HashMap<String, HashMap<String, HashSet<String>>>, fastq_reader: &mut seq_io::fastq::Reader<R>, num_hash: &mut u64) -> Result<(), Box<dyn std::error::Error>> {

  /*
  ** Loop through input reads.
  */
  let mut header: &str;
  let mut seq: &str;
  let mut hashbc: &str;
  let mut hashval: &String = &"none".to_string();
  let mut polya: &str;
  let mut header_toks: Vec<&str>;
  let mut len_header_toks: usize;
  let mut cell_barc: String;
  let mut umi: String;
  let mut is_hash: bool;
  let mut i_edit_distance: usize = 0;

  while let Some(result) = fastq_reader.next() {

    let record = result.expect("bad status reading fastq file");
    header = std::str::from_utf8(record.head()).expect("bad status getting fastq read header");
    seq = std::str::from_utf8(record.seq()).expect("bad status getting fastq read sequence");
 
    polya = &seq[11..15];

    is_hash = false;
    if(polya == "AAAA") {
      hashbc = &seq[0..10];
      for i in (0..hash_edit_distance+1) {
        if let Some(hashval_tmp) = hash_whitelist[i].get(hashbc) {
          is_hash = true;
          hashval = hashval_tmp;
          i_edit_distance = i;
          break
        }
      }
    }

    if(!is_hash) {
      continue;
    }

    header_toks = header.split("|").collect();
    len_header_toks = header_toks.len();
//    cell_barc = header_toks[2..len_header_toks-1].join("_").to_string();
    cell_barc = header_toks[2..len_header_toks-1].join("_");
    umi = header_toks[len_header_toks-1].to_string();


    /*
    ** The following line if for checking the read parsing. Keep it
    ** commented out unless necessary for testing.
    */
    // println!("hash_dict: {}|{}|{}", hashval, cell_barc, umi);


    /*
    ** Record distinct cells.
    */
    cells.insert(cell_barc.clone());

    /*
    ** Count distinct hash reads.
    */
    let mut umi_new: usize = 0;
    let _result = update_nested_maps(hashval, cell_barc, umi, hashdict, &mut umi_new);
    if(umi_new != 0) {
      hash_counts[i_edit_distance] += 1;
    }

    *num_hash += 1;
  }

  Ok(())
}


#[inline(always)]
fn update_nested_maps(hash_key1: &String, hash_key2: String, set_key1: String, map1: &mut HashMap<String, HashMap<String, HashSet<String>>>, umi_new: &mut usize) -> Result<(), Box<dyn std::error::Error>> {
  /*
  ** We expect that map1 was initialized with hash_key1 so we don't need
  ** to use map1.entry(key1.to_string()).or_insert(HashMap::new()); instead,
  ** we return an error if the hash_key1 isn't in map1.
  let map2 = map1.entry(key1.to_string()).or_insert_with(|| HashMap::new()).expect("unable to insert key/value pair");
  */

  let map2 = map1.get_mut(hash_key1).expect("hash barcode not in hash_sheet file.");
  let set1 = map2.entry(hash_key2).or_insert_with(|| HashSet::new());
  if(set1.insert(set_key1) == true) {
    *umi_new += 1;
  }

  Ok(())
}


fn make_sparse_matrix(cells: &HashSet<String>, hash_lookup: &HashMap<String, String>, map: &mut HashMap<String, HashMap<String, HashSet<String>>>) -> Result<(Vec<String>, Vec<String>, CsMat<usize>), Box<dyn std::error::Error>> {
//   -> Result<(Vec<String>, Vec<String>, sprs::TriMatBase<Vec<usize>, Vec<usize>>), Box<dyn std::error::Error>> {
 

  /*
  ** Note: it's possible to use the indexing operator
  **       with a HashMap get() operation, as long as
  **       one can be certain that the key exists in
  **       the HashMap.
  */
  let mut num_hash: usize = 0;
  let mut cell_counter: HashSet<String> = HashSet::new();
  for key1 in map.keys() {
    if(!map[key1].is_empty()) {
      num_hash += 1;
      for cell in map[key1].keys() {
        cell_counter.insert(cell.to_string());
      }
    }
  }

  let num_cell = cell_counter.len();
  if(num_cell != cells.len()) {
    panic!("Error: inconsistent cell counts! {} {}", num_cell, cells.len());
  }

  /*
  ** Columns are cells.
  ** Rows are hash reads.
  */
  let mut data: sprs::TriMatBase<Vec<usize>, Vec<usize>> = sprs::TriMat::with_capacity((num_hash, num_cell), 100);

  /*
  ** Hash reads.
  */
  let mut hash_seqs: Vec<String> = Vec::with_capacity(num_hash);
  let mut row_names: Vec<String> = Vec::with_capacity(num_hash);
  for hash_seq in map.keys() {
    if(!map.get(hash_seq).unwrap().is_empty()) {
      hash_seqs.push(hash_seq.to_string());
      row_names.push(hash_lookup.get(hash_seq).unwrap().to_string());
    }
  }

  /*
  ** Cells.
  */
  let mut col_names: Vec<String> = Vec::with_capacity(num_cell);
  for col_name in cells {
    col_names.push(col_name.to_string());
  }

  /*
  ** Load triplet matrix.
  */
  for (irow, hash_seq) in hash_seqs.iter().enumerate() {
    for (icol, col_name) in col_names.iter().enumerate() {
      let cell_map = map.get(hash_seq as &str).unwrap();

      /*
      ** There are not UMIs for every hash and cell combination so
      ** cell_map.get() is expect to return None at times.
      */
      let result = cell_map.get(col_name as &str);
      if let Some(set) = result {
        data.add_triplet(irow, icol, set.len());
      }
    }
  }

  /*
  ** Make CSC matrix.
  */
  let mat_sparse: CsMat<usize> = data.to_csc();

  Ok((row_names, col_names, mat_sparse))
}


/*
** Write sparse matrix to Matrix Market file.
*/
fn write_matrix(key: &String, row_names: Vec<String>, col_names: Vec<String>, smat: &sprs::CsMat<usize> ) -> Result<(), std::io::Error> {

  let mut file_name: String;
  let mut path: &Path;
  let mut file: File;
  let mut writer: BufWriter<File>;

  /*
  ** Write Matrix Market file.
  */
  file_name = format!("{}.hashumis.mtx", *key);
  path = Path::new(&file_name);
  let _result = sprs::io::write_matrix_market(path, smat).expect(&format!("error while writing file {}", file_name));

  /*
  ** Write row name file.
  */
  file_name = format!("{}.hashumis_hashes.txt", *key);
  path = Path::new(&file_name);
  file = File::create(path).expect(&format!("unable to open file {}", file_name));
  writer = BufWriter::new(file);
  let _ = writer.write_all(row_names.join("\n").as_bytes()).expect(&format!("error while writing file {}", file_name));
  let _ = writer.flush();
  std::mem::drop(writer);

  /*
  ** Write column name file.
  */
  file_name = format!("{}.hashumis_cells.txt", *key);
  path = Path::new(&file_name);
  file = File::create(path).expect(&format!("unable to open file {}", file_name));
  writer = BufWriter::new(file);
  let _ = writer.write_all(col_names.join("\n").as_bytes()).expect(&format!("error while writing file {}", file_name));
  let _ = writer.flush();
  std::mem::drop(writer);

  Ok(())
}


#[allow(dead_code)]
fn dump_nested_maps(map: &mut HashMap<String, HashMap<String, HashSet<String>>>) {
  for key1 in map.keys() {
    for key2 in map.get(key1).unwrap().keys() {
      for key3 in map.get(key1).unwrap().get(key2).unwrap().into_iter() {
        println!("table: key1: {}  key2: {}  key3: {}", key1, key2, key3);
      }
    }
  }
}


/// Find RNA-seq hash reads in a fastq file and write counts as a Matrix Market file.
///
/// Arguments:
///- -s hash sheet file in TSV format where the columns are *barcode_name*  *barcode_sequence*
///- -f fastq file names separated by spaces or '-' for stdin. The files may be gzip compressed.
///- -d hash edit distance allows for n discrepancies. n=1 by default.
///- -k output file root name
///
/// Return:
///
/// Nothing.
///
/// Reading fastq files rather than from stdin appears to be substantially faster.
///
fn main() {
  /*
  ** Get command line arguments.
  */
  let clarg = Command::new("process_hashes")
        .version("0.1.0")
        .about("Finds hash sequence reads in fastq file.")
        .arg(Arg::new("hash_sheet")   //required=true, no default
                  .required(true)
                  .short('s')
                  .long("hash_sheet")
                  .help("Path to hash sample sheet."))
        .arg(Arg::new("fastq")   // required=true
                  .required(true)
                  .num_args(1..)
                  .value_delimiter(' ')
                  .value_terminator("--")
                  .short('f')
                  .long("fastq")
                  .help("Input fastq filenames separated by spaces or '-' for stdin."))
        .arg(Arg::new("hash_edit_distance")  // required=false, default=1
                  .required(false)
                  .default_value("1")
                  .short('d')
                  .long("hash_edit_distance")
                  .value_parser(clap::value_parser!(usize))
                  .help("Allowed edit distance for hashes."))
        .arg(Arg::new("key")   //required=true
                  .required(true)
                  .short('k')
                  .long("key")
                  .help("Key for file name for output file."))
        .get_matches();

  let hash_sheet: String = clarg.get_one::<String>("hash_sheet").unwrap().to_string();
  let fastq_filenames: Vec<_> = clarg.get_many::<String>("fastq").unwrap().collect();
  let hash_edit_distance: usize = *clarg.get_one::<usize>("hash_edit_distance").unwrap();
  let key: String = clarg.get_one::<String>("key").unwrap().to_string();

  /*
  ** Set up hash whitelist map.
  */
  let hash_lookup = barcode_utils::read_barcode_file(&hash_sheet).unwrap();
  let barcodes: Vec<String> = hash_lookup.keys().cloned().collect();
  let hash_whitelist: Vec<HashMap<String, String>> = barcode_utils::construct_mismatch_to_whitelist_map(barcodes, 1, true).unwrap();

  /*
  ** Initialize counters.
  */

  /*
  ** Hash reads per cell counter.
  */
  let mut hashdict: HashMap<String, HashMap<String, HashSet<String>>> = HashMap::with_capacity(256);
  for hashseq in hash_lookup.keys() {
    hashdict.insert(hashseq.to_string(), HashMap::new());
  }

  /*
  ** Cell counter.
  ** Each element is a distinct cell barcode. Use cells.get(&seq)
  ** to get a reference to the barcode sequence. Here we use
  ** barcode sequence references as keys in maps indexed by 'cell'.
  */
  let mut cells: HashSet<String> = HashSet::with_capacity(100000);

  /*
  ** Total hashes assigned counter.
  */
  let mut hash_counts: Vec<u64> = vec![0; hash_edit_distance+1];

  /*
  ** Total hash counter used for diagnostics.
  */
  let mut num_hash: u64 = 0;

  /*
  ** Set up file handle and reader from either stdin or specified file.
  */
  //t reader: Box<dyn std::io::Read> = 
  if(fastq_filenames[0] == "-") {
    let reader = Box::new(std::io::stdin());
    let mut fastq_reader = fastq::Reader::new(reader);
    let _ = process_fastq_file(hash_edit_distance, &hash_whitelist, &mut cells, &mut hash_counts, &mut hashdict, &mut fastq_reader, &mut num_hash);
  }
  else {
    for fastq_filename in fastq_filenames {
      let reader = Box::new(detect_compression::DetectReader::open(fastq_filename).unwrap());
      let mut fastq_reader = fastq::Reader::new(reader);
      let _ = process_fastq_file(hash_edit_distance, &hash_whitelist, &mut cells, &mut hash_counts, &mut hashdict, &mut fastq_reader, &mut num_hash);
    }
  };

  // dump_nested_maps(&mut hashdict);

  /*
  ** Make a sparse, CSC, matrix from hashdict.
  */
  let (row_names, col_names, smat) = make_sparse_matrix(&cells, &hash_lookup, &mut hashdict).unwrap();

  /*
  ** Write matrix.
  */
  let _ = write_matrix(&key, row_names, col_names, &smat);


  let file_name = format!("{}_hash.log", key);
  let path = Path::new(&file_name);
  let file = File::create(path).expect(&format!("unable to open file {}", file_name));
  let mut writer = BufWriter::new(file);
  for i_dist in (0..hash_edit_distance+1) {
    writeln!(writer, "Hash UMIs detected with {} correction: {}", i_dist, hash_counts[i_dist]).expect("error writing log file");
  }
  let _ = writer.flush();

  /*
  ** Let OS free memory.
  */
  std::process::exit(0);
}

