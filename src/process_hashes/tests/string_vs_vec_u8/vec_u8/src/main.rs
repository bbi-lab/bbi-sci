#![allow(unused_parens)]

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::env;
use std::path::Path;
use std::fs::File;
use std::io::{Read, Write, BufWriter};
use std::str;
//use std::collections::{HashMap, HashSet};
use std::collections::hash_map::Entry;

use rustc_hash::FxHashMap;
use vec_u8::barcode_utils;
use seq_io::fastq::{Reader, Record};
use detect_compression;
use itertools::Itertools;


#[inline(always)]
fn update_nested_maps(hash_key1: &[u8], hash_key2: Vec<u8>, hash_key3: &[u8], map1: &mut FxHashMap<Vec<u8>, FxHashMap<Vec<u8>, FxHashMap<Vec<u8>, u64>>>) -> Result<(), Box<dyn std::error::Error>> {
  /*
  ** We expect that map1 was initialized with hash_key1 so we don't need
  ** to use map1.entry(key1.to_string()).or_insert(FxHashMap::new()); instead,
  ** we return an error if the hash_key1 isn't in map1.
  let map2 = map1.entry(key1.to_string()).or_insert_with(|| FxHashMap::default()).expect("unable to insert key/value pair");
  */

  /*
  ** map1 key is hash barcode (hash_key1).
  */
  let map2 = map1.get_mut(hash_key1).expect("hash barcode not in hash_sheet file.");

  /*
  ** map2 key is cell barcode (hash_key2).
  */
  let map3 = map2.entry(hash_key2).or_insert_with(|| FxHashMap::default());

  /*
  ** map3 key is UMI barcode (hash_key3).
  */
  let umi_count = map3.entry(hash_key3.to_vec()).or_insert(0);
  *umi_count += 1;

  Ok(())
}


fn main() {
  let hashseq_filename: String = String::from("/home/brent/work/data_sets/zebrafish/fastq_hashed/RNA3-072-a/timecourse_hash2.txt");

  let mut header: &[u8];
  let mut seq: &[u8];
  let mut hashbc: &[u8];
  let mut hashval: &[u8] = &[];
  let mut polya: &[u8];
  let mut header_toks: Vec<&[u8]>;
  let mut len_header_toks: usize;
  let mut cell_barc: Vec<u8>;
  let mut umi: &[u8];
  let hash_edit_distance: usize = 1;
  let mut i_edit_distance;
  let mut is_hash: bool;

  let hash_lookup = barcode_utils::read_barcode_file(&hashseq_filename).unwrap();
  let barcodes: Vec<String> = hash_lookup.keys().cloned().collect();
  let hash_whitelist_string: Vec<FxHashMap<String, String>> = barcode_utils::construct_mismatch_to_whitelist_map(barcodes, 1, true).unwrap();
  let hash_whitelist: Vec<FxHashMap<Vec<u8>, Vec<u8>>> = barcode_utils::mismatch_to_whitelist_map_as_u8(hash_whitelist_string).expect("oops");

  let mut hashdict: FxHashMap<Vec<u8>, FxHashMap<Vec<u8>, FxHashMap<Vec<u8>, u64>>> = FxHashMap::default();
  for hashseq in hash_lookup.keys() {
    hashdict.insert(hashseq.clone().into_bytes(), FxHashMap::default());
  }

//std::process::exit(0);

  for (i, arg) in env::args().enumerate() {
    println!("i: {}  arg: {}", i, arg);
    if(i == 0) {
      continue;
    }
    if(i > 1) {
      break;
    }
    let fastq_filename: String = String::from(arg);
    let reader = Box::new(detect_compression::DetectReader::open(fastq_filename).unwrap());
    let mut fastq_reader = Reader::new(reader);

    while let Some(result) = fastq_reader.next() {
     let record = result.expect("bad status reading fastq file");
      header = record.head();
      seq = record.seq();

      /*
      ** Is this a hash read?
      */
      polya = &seq[11..15];
      is_hash = false;
      if(polya == b"AAAA") {
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

      /*
      ** This is a hash read!
      */
      header_toks = header.split(|i| *i == b'|').collect();
      len_header_toks = header_toks.len();
      /*
      ** The join method allocates heap memory to store cell_barc
      ** so I pass cell_barc to update_nested_maps as a Vec
      ** rather than a reference, which would need to be converted
      ** back to a Vec for the call to entry(). (update_nested_maps
      ** should be inlined anyway.)
      */
      cell_barc = header_toks[2..len_header_toks-1].join(&b'_');
      umi = header_toks[len_header_toks-1];

      let _result = update_nested_maps(hashval, cell_barc, umi, &mut hashdict);
    }
  }
  std::process::exit(0);
}
