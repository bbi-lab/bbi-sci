/*
** Turn off warnings about unused parentheses.
*/
#![allow(unused_parens)]


//! Reproduce Andrew Hill's barcode_utils.py in Rust, a bit at a time.
pub mod barcode_utils {

  use std::io;
  use std::error::Error;
  use std::collections::HashMap;
  use itertools::Itertools;
  use csv::{ReaderBuilder, Trim};
  use serde::Deserialize;
  use regex::Regex;

  /// Read a barcode TSV file that has the format
  ///
  /// barcode name\tbarcode sequence
  ///
  /// Arguments:
  ///- file_path: a reference to the input file path
  ///
  /// Return:
  ///
  /// A hash map: the key is the barcode sequence and the value is
  /// the barcode name.
  ///
  /// Notes:
  ///
  ///- the barcode sequence names are edit to replace [ )(=_/] with [.].
  ///
  /// Based on URL: https://stackoverflow.com/questions/78639668/fast-reading-from-a-tsv-file-in-rust
  #[derive(Deserialize)]
  struct Record {
    hash_name: String,
    hash_barcode: String,
  }

  pub fn read_barcode_file(file_path: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
  let fp = std::fs::File::open(file_path)?;
  let buf_reader = io::BufReader::new(fp);
  let mut tsv_reader = ReaderBuilder::new()
                         .has_headers(false)
                         .trim(Trim::Fields)
                         .delimiter(b'\t')
                         .comment(Some(b'#'))
                         .from_reader(buf_reader);

  let mut hash_map: HashMap<String, String> = HashMap::new();

  /*
  ** Build regex for replacing characters in names.
  */
  let re = Regex::new(r"[ )(=_/-]").unwrap();

  /*
  ** Read the file and store the names and sequences.
  */
  for result in tsv_reader.deserialize() {
     let record: Record = result?;
     hash_map.entry(record.hash_barcode.to_owned()).or_insert_with(|| re.replace_all(&(record.hash_name), ".").to_string());
  }

  Ok(hash_map)
}


  /// Generate a vector of mismatched sequences to a given sequence. Must contain only ACGT.
  /// This is based heavily on a Biostars answer.
  /// Curiously, it returns the input sequence when num_mismatches is zero; otherwise, it does not.
  ///
  /// Arguments:
  ///- sequence: a reference to the input string.
  ///- num_mismatches: the number of mismatches in each output sequence.
  ///- allow_n: include Ns in the output sequences.
  ///
  /// Return:
  ///
  /// A Result enum containing a vector of mismatched sequence Strings.
  ///
  pub fn generate_mismatches(sequence: &str, num_mismatches: usize, allow_n: bool) -> Result<Vec<String>, Box<dyn Error>> {

    let sequence_upper: String = sequence.to_uppercase();

    if(num_mismatches == 0) {
      return(Ok(vec!(sequence_upper)));
    }

    let mut letters = String::from("ACGT");
    if(allow_n == true) {
      letters.push('N');
    }

    /*
    ** Return a Result of vector of strings.
    */
    let mut mismatches: Vec<String> = Vec::new();

    let seq_len = sequence_upper.len();

    /*
    ** Make iterators of vectors of indices where mismatches
    ** may show up in the input sequence, and iterator through
    ** these vectors.
    */
    let combination_iter = (0..seq_len).combinations(num_mismatches);
    for locs in combination_iter {
      /*
      ** sequence_vector begins with a vector of vectors where each base
      ** in the sequence is an inner vector. That is, for sequence ATGCTA,
      ** each loop iteration gives the save vector of vectors, at this point.
      **   sequence_vector: [['A'], ['T'], ['G'], ['C'], ['T'], ['A']]
      **   ...
      */
      let mut sequence_vector: Vec<Vec<char>> = sequence_upper.chars().map(|c| vec!(c)).collect();

      /*
      ** Here inner vectors are replaced with the three bases that are not
      ** in the original sequence. For example,
      **   [['A'], ['T'], ['G'], ['C'], ['T'], ['A']]
      ** becomes
      **   [['A'], ['A','C','G'], ['G'], ['C'], ['T'], ['A']]
      ** at mismatch index 1, which is the second inner
      ** vector.
      */
      for loc in locs {
        let orig_char = sequence_upper.chars().nth(loc as usize).unwrap();
        sequence_vector[loc] = letters.chars().filter(|c| *c != orig_char).collect();
      }

      /*
      ** And the cartesian product of the set of inner vectors expands
      ** to the desired mismatch sequences when the product vectors are
      ** converted to strings.
      */
      for vector_set in sequence_vector.iter().multi_cartesian_product() {
        let primer_string: String = vector_set.into_iter().collect();
        mismatches.push(primer_string);
      }
    }

    Ok(mismatches)
  }


  /// Construct a precomputed set  of all mismatches within a specified
  /// edit distance and the barcode whitelist.
  ///
  /// Arguments:
  /// * whitelist: xxx
  ///
  /// Return:
  ///
  /// A result enum containing a vector of maps of mismatched sequences to
  /// their whitelist sequences.
  ///
  pub fn construct_mismatch_to_whitelist_map(whitelist: Vec<String>, edit_distance: usize, allow_n: bool) -> Result<Vec<HashMap<String, String>>, Box<dyn Error>> {

    /*
    ** Set whitelist sequences to upper-case.
    */
    let whitelist_upper = whitelist.iter().map(|s| s.to_uppercase()).collect::<Vec<String>>();

    /*
    ** mismatch_to_whitelist_map is a vector of hash maps where the vector
    ** index is the number of substitutions in the mismatch sequence, and
    ** the hash maps are keyed by the sequences with substitutions and the
    ** values are the original whitelist sequences (no mismatches).
    */
    let mut mismatch_to_whitelist_map: Vec<HashMap<String, String>> = Vec::with_capacity(edit_distance+1);
    for _i in (0..edit_distance+1) {
      mismatch_to_whitelist_map.push(HashMap::new());
    }

    /*
    ** Set the zero mismatch sequence maps where the key and
    ** value are the same.
    */
    mismatch_to_whitelist_map[0] = whitelist_upper.iter().map(|s| (s.to_owned(), s.to_owned())).collect::<HashMap<String, String>>();


    /*
    ** Track  conflicts where mismatches map to different sequences.
    */
    let mut conflicting_mismatches: Vec<String> = Vec::new();

    /*
    ** Doesn't really matter as correction function will never see it,
    ** but exclude any perfect matches to actual seqs by mismatches.
    */
    conflicting_mismatches.extend(whitelist_upper.clone());


    for mismatch_count in (1..edit_distance+1) {

      for sequence in &whitelist_upper {
        /*
        ** Generate all possible mismatches in range.
        */
        let mismatches = generate_mismatches(&sequence, mismatch_count, allow_n)?;

        for mismatch in mismatches.iter() {
          if(mismatch_to_whitelist_map[mismatch_count].contains_key::<str>(&mismatch) == true) {
            conflicting_mismatches.push(mismatch.to_string());
          }
          mismatch_to_whitelist_map[mismatch_count].insert(mismatch.to_string(), sequence.clone());
        }
      }

      /*
      ** Go back and remove any conflicting mismatches.
      */
      for mismatch in conflicting_mismatches.clone().into_iter().unique() {
        if(mismatch_to_whitelist_map[mismatch_count].contains_key::<str>(&mismatch) == true) {
          mismatch_to_whitelist_map[mismatch_count].remove(&mismatch);
        }
      }
    }

    Ok(mismatch_to_whitelist_map)
  }

}

