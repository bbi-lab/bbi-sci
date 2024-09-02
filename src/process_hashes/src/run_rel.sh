#!/bin/bash

root_dir="/home/brent/work/data_sets/zebrafish/fastq_hashed/RNA3-072-a"

#/home/brent/git/bbi-sci/src/process_hashes/target/release/process_hashes -s $root_dir/timecourse_hash2.txt -f $root_dir/demux_out/SeahubZ01.fq.part1-L001.fastq.gz -k foop -n SeahubZ01

/home/brent/git/bbi-sci/src/process_hashes/target/release/process_hashes -s $root_dir/timecourse_hash2.txt -f $root_dir/demux_out/*.fastq.gz -k foop -n SeahubZ01

#zcat $root_dir/demux_out/SeahubZ01.fq.part1-L001.fastq.gz | /home/brent/git/bbi-sci/src/process_hashes/target/release/process_hashes -s $root_dir/timecourse_hash2.txt -f - -k foop -n SeahubZ01

#/home/brent/git/bbi-sci/src/process_hashes/target/release/process_hashes -s $root_dir/timecourse_hash2.txt -f SeahubZ01.fq.part1-L001.fastq -k foop -n SeahubZ01

