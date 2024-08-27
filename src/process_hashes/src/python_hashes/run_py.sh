#!/bin/bash

root_dir="/home/brent/work/data_sets/zebrafish/fastq_hashed/RNA3-072-a"

#zcat $root_dir/demux_out/SeahubZ01.fq.part1-L001.fastq.gz | cargo run -- -s $root_dir/timecourse_hash2.txt -f - -k foop

# usage: process_hashes.py [-h] --hash_sheet HASH_SHEET --fastq [FASTQ] [--hash_edit_distance HASH_EDIT_DISTANCE] --key KEY

process_hashes.py --hash_sheet $root_dir/timecourse_hash2.txt --fastq <(zcat  $root_dir/demux_out/SeahubZ01.fq.part1-L001.fastq.gz) --key foo

