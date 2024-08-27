#!/bin/bash

root_dir="/home/brent/work/data_sets/zebrafish/fastq_hashed/RNA3-072-a"

/home/brent/src/rust/process_hashes/target/release/process_hashes -s $root_dir/timecourse_hash2.txt -f $root_dir/demux_out/SeahubZ01.fq.part1-L001.fastq.gz -k foop

#zcat $root_dir/demux_out/SeahubZ01.fq.part1-L001.fastq.gz | /home/brent/src/rust/process_hashes/target/release/process_hashes -s $root_dir/timecourse_hash2.txt -f - -k foop

#/home/brent/src/rust/process_hashes/target/release/process_hashes -s $root_dir/timecourse_hash2.txt -f SeahubZ01.fq.part1-L001.fastq -k foop

