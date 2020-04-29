
// Parse input parameters
params.help = false
params.samples = false
params.star_file = "$baseDir/bin/star_file.txt"
params.gene_file = "$baseDir/bin/gene_file.txt"
params.umi_cutoff = 100
params.rt_barcode_file="default"
params.max_cores = 16
params.hash_list = false
params.max_wells_per_sample = 20

//print usage
if (params.help) {
    log.info ''
    log.info 'BBI sci-RNA-seq Pipeline'
    log.info '--------------------------------'
    log.info ''
    log.info 'For reproducibility, please specify all parameters to a config file'
    log.info 'by specifying -c CONFIG_FILE.config.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run bbi-sci -c CONFIG_FILE'
    log.info ''
    log.info 'Help: '
    log.info '    --help                                     Show this message and exit.'
    log.info ''
    log.info 'Required parameters (specify in your config file):'
    log.info '    params.output_dir = OUTPUT DIRECTORY       Output directory.'
    log.info '    params.sample_sheet = SAMPLE_SHEET_PATH    Sample sheet of the format described in the README.'
    log.info '    params.demux_out = DEMUX OUTPUT DIR        Path to the demux_out folder from the bbi-dmux run.'
    log.info '    params.level = 3                           2 or 3 level sci?'
    log.info ''
    log.info 'Optional parameters (specify in your config file):'
    log.info '    params.rt_barcode_file = "default"         The path to a custom RT barcode file. If "default", default BBI barcodes will be used.'
    log.info '    params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.'
    log.info '    process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.'
    log.info '    process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted. '
    log.info '    params.samples = [sample1, sample2]        Add to only run certain samples from trimming on. Default is to run all.'
    log.info '    params.star_file = PATH/TO/FILE            File with the genome to star maps, similar to the one included with the package.'
    log.info '    params.gene_file = PATH/TO/FILE            File with the genome to gene model maps, similar to the one included with the package.'
    log.info '    params.umi_cutoff = 100                    The umi cutoff to be called a cell in matrix output.'
    log.info '    params.hash_list = false                   Path to a tab-delimited file with at least two columns, first the hash name and second the hash barcode sequence. Default is false to indicate no hashing.'
    log.info '    params.max_wells_per_sample = 20           The maximum number of wells per sample - if a sample is in more wells, the fastqs will be split then reassembled for efficiency.'
    log.info ''
    log.info 'Issues? Contact hpliner@uw.edu'
    exit 1
}

// check required options
if (!params.output_dir || !params.sample_sheet || !params.level || !params.demux_out) {
    exit 1, "Must include config file using -c CONFIG_FILE.config that includes output_dir, sample_sheet, level and demux_out"
}


/*************

Process: check_sample_sheet
 Inputs:
    params.sample_sheet
    params.star_file
    params.level
    params.max_wells_per_sample
    params.rt_barcode_file

 Outputs:
    good_sample_sheet - corrected csv sample sheet
    logfile - running log
    log_piece1 - piece of log to be concatenated for full log

 Summary:
    Check and process sample sheet - check_sample_sheet.py
    Start log

 Downstream:
    gather_info
    trim_fastqs
    combine_logs

*************/

process check_sample_sheet {
    cache 'lenient'
    module 'modules:java/latest:modules-init:modules-gs:python/3.6.4'


    output:
        file "*.csv" into good_sample_sheet
        file '*.log' into log_check_sample
        file 'start.txt' into log_piece1

    """
    printf "BBI bbi-sci Pipeline Log\n\n" > start.log
    printf "Run started at: \$(date)\n\n" >> start.log

    printf "***** BEGIN PIPELINE *****: \n\n" >> start.log
    printf "** Start process 'check_sample_sheet' at: \$(date)\n\n" >> start.log
    printf "    Process versions: 
        \$(python --version)\n\n" >> start.log
    printf "    Process command: 
        check_sample_sheet.py --sample_sheet $params.sample_sheet 
            --star_file $params.star_file --level $params.level\n\n" >> start.log

    check_sample_sheet.py --sample_sheet $params.sample_sheet --star_file $params.star_file --level $params.level --rt_barcode_file $params.rt_barcode_file --max_wells_per_samp $params.max_wells_per_sample

    printf "** End process 'check_sample_sheet' at: \$(date)\n\n" >> start.log
    cp start.log start.txt
    """

}

// Generate a sample list with fixed naming
samp_file = file(params.sample_sheet)
def samp_list = []
for (line in samp_file.readLines()) {
  samp_list.add(line.split(",")[1])
}
samp_list = samp_list.collect{"$it".replaceAll(/\s/, ".").replaceAll(/_/, ".").replaceAll(/-/, ".").replaceAll(/\\//, ".")}
samp_list.removeElement("Sample.ID")
if (params.samples != false) {
  samp_list.intersect(params.samples.collect{"$it".replaceAll(/\s/, ".").replaceAll(/_/, ".").replaceAll(/-/, ".").replaceAll(/\\//, ".")})
}


/*************

Process: trim_fastqs
 Inputs:
    input_fastq - all fastq files from params.demux_out folder
    logfile - running log
 
 Outputs:
    trim_out - output folder from trimming - stops here
    key - sample id
    name - file id (including lane and split info)
    trimmed_fastq - trimmed, gzipped fastq
    logfile - running log
    log_piece2 - piece of log to be concatenated for full log
    input_fastq - all fastq files from params.demux_out folder

 Summary:
    Trim fastqs - trim_galore
    Continue log
    
 Downstream:
    gather_info
    process_hashes

 Notes:
    Only moves forward if sample is in samp_list - where params.samples comes in

*************/

process trim_fastqs {
    cache 'lenient'
    memory '12G'
    module 'java/latest:modules:modules-init:modules-gs:python/2.7.3:cutadapt/1.8.3:trim_galore/0.4.1'

    input:
        file input_fastq from Channel.fromPath("${params.demux_out}/*.fastq.gz")
        file logfile from log_check_sample

    output:
        file "trim_out" into trim_output
        set val(key), val(name), file("trim_out/*.fq.gz"), file('trim.log'), file('*trim.txt') into trimmed_fastqs
        set val(key), file(input_fastq) into fastqs_out

    when:
	!((input_fastq.name.split(/-L[0-9]{3}/)[0].split(/\.fq.part/)[0]) in "Undetermined") && ((input_fastq.name.split(/-L[0-9]{3}/)[0].split(/\.fq.part/)[0]) in samp_list)

    script:
    name = input_fastq.baseName - ~/.fastq/
    key = input_fastq.baseName.split(/-L[0-9]{3}/)[0].split(/\.fq.part/)[0]

    """
    cat ${logfile} > trim.log
    printf "** Start process 'trim_fastqs' for $input_fastq at: \$(date)\n\n" > piece.log
    printf "    Process versions: 
        " >> piece.log
    python --version &>> piece.log
    printf "        trim_galore \$(trim_galore -v | grep version | awk '{\$1=\$1;print}')
        cutadapt version \$(cutadapt --version)\n\n" >> piece.log

    printf "    Process command: 
        trim_galore $input_fastq -a AAAAAAAA --three_prime_clip_R1 1 
            --gzip -o ./trim_out/\n
    Process output:\n" >> piece.log
    mkdir trim_out
    trim_galore $input_fastq \
        -a AAAAAAAA \
        --three_prime_clip_R1 1 \
        --gzip \
        -o ./trim_out/
    cat trim_out/*trimming_report.txt | sed '/Overview of/,+RUN d' >> piece.log
    printf "** End process 'trim_fastqs' at: \$(date)\n\n" >> piece.log
    cp piece.log ${name}_trim.txt
    cat piece.log >> trim.log
    """
}


/*************

Process: gather_info
 Inputs:
    input_fastq - all fastq files from params.demux_out folder
    logfile - running log
 
 Outputs:
    trim_out - output folder from trimming - stops here
    key - sample id
    name - file id (including lane and split info)
    trimmed_fastq - trimmed, gzipped fastq
    logfile - running log
    log_piece2 - piece of log to be concatenated for full log
    input_fastq - all fastq files from params.demux_out folder

 Summary:
    Trim fastqs - trim_galore
    Continue log
    
 Downstream:
    gather_info
    process_hashes

 Notes:
    Only moves forward if sample is in samp_list - where params.samples comes in

*************/

process gather_info {
    cache 'lenient'
    memory '1G'
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'

    input:
        file good_sample_sheet
        set val(key), val(name), file(trimmed_fastq), file(logfile), file(log_piece2) from trimmed_fastqs

    output:
        set val(key), val(name), env(star_path), env(star_mem), file(trimmed_fastq), file(logfile), file(log_piece2) into align_prepped
        set val(key), env(gtf_path) into gtf_info
        set val(key), env(gtf_path) into gtf_info2

    """
    spec=`awk 'BEGIN {FS=",";OFS=","};{sub(" ", ".", \$2);sub("/", ".", \$2);sub("-", ".", \$2);sub("_", ".", \$2);split(\$2,a,"_fq_part");print(\$1, a[1], \$3)}' $good_sample_sheet | awk 'BEGIN {FS=","}; \$2=="$key" {print \$3}' | uniq`
    star_mem=`awk -v var="\$spec" '\$1==var {print \$3}' $params.star_file | uniq`
    star_path=`awk -v var="\$spec" '\$1==var {print \$2}' $params.star_file | uniq`
    gtf_path=`awk -v var="\$spec" '\$1==var {print \$2}' $params.gene_file | uniq`    

    """
}

fastqs_out
    .groupTuple()
    .set { for_hash }

save_hash_cell = {params.output_dir + "/" + it - ~/.hashumis_cells.txt/ + "/" + it}
save_hash_hash = {params.output_dir + "/" + it - ~/.hashumis_hashes.txt/ + "/" + it}
save_hash_mtx = {params.output_dir + "/" + it - ~/.hashumis.mtx/ + "/" + it}


process process_hashes {
    cache 'lenient'
    memory '12G'
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_cell, pattern: "*hashumis_cells.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_hash, pattern: "*hashumis_hashes.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_mtx, pattern: "*.mtx", mode: 'copy'

    input:
        set key, file(fqs) from for_hash

    output:
        file("*hash.log") into hash_logs
        set file("*mtx"), file("*hashumis_cells.txt"), file("*hashumis_hashes.txt") into hash_mats

    when:
        params.hash_list != false

    """
     process_hashes.py --hash_sheet $params.hash_list \
         --fastq <(zcat $fqs) --key $key 

    """
}

if (params.max_cores < 8) {
    cores_align = params.max_cores
} else {
    cores_align = 8
}

process align_reads {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:STAR/2.5.2b'
    memory { mem.toInteger()/cores_align + " GB" }
    penv 'serial'
    cpus cores_align

    input:
        set val(key), val(name), val(star), val(mem), file(input_file), file(logfile), file(log_piece2) from align_prepped

    output:
        file "align_out" into align_output
        set val(key), val(name), file("align_out/*Aligned.out.bam"), file('align.log'), file(log_piece2), file("*align.txt") into aligned_bams

    """
    cat ${logfile} > align.log
    printf "** Start process 'align_reads' for $input_file at: \$(date)\n\n" > piece.log
    printf "    Process versions: 
        \$(STAR --version)\n\n" >> piece.log

    printf "    Process command: 
        STAR --runThreadN $cores_align --genomeDir $star 
            --readFilesIn $input_file --readFilesCommand zcat --outFileNamePrefix ./align_out/${name} 
            --outSAMtype BAM Unsorted --outSAMmultNmax 1 --outSAMstrandField intronMotif\n
    Process output:\n" >> piece.log

    mkdir align_out 
    STAR \
        --runThreadN $cores_align \
        --genomeDir $star \
        --readFilesIn $input_file \
        --readFilesCommand zcat \
        --outFileNamePrefix ./align_out/${name}  \
        --outSAMtype BAM Unsorted \
        --outSAMmultNmax 1 \
        --outSAMstrandField intronMotif

    cat align_out/*Log.final.out >> piece.log

    printf "\n** End process 'align_reads' at: \$(date)\n\n" >> piece.log

    cp piece.log ${name}_align.txt
    cat piece.log >> align.log

    """

}

if (params.max_cores < 10) {
    cores_sf = params.max_cores
} else {
    cores_sf = 10
}

process sort_and_filter {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:samtools/1.4'
    memory '1 GB'
    penv 'serial'
    cpus cores_sf

    input:
        set val(key), val(name), file(aligned_bam), file(logfile), file(log_piece2), file(log_piece3) from aligned_bams
    
    output:
        set val(key), file("*.bam") into sorted_bams
        file "*_sf.log" into bam_logs
        set val(key), file(log_piece2), file(log_piece3), file("*_sf.txt") into log_pieces

    """
    printf "** Start process 'sort_and_filter' for $aligned_bam at: \$(date)\n\n" > ${name}_piece.log
    printf "    Process versions: 
        \$(samtools --version | tr '\n' ' ')\n\n" >> ${name}_piece.log
    printf "    Process command: 
        samtools view -bh -q 30 -F 4 '$aligned_bam' 
            | samtools sort -@ $cores_sf - > '${name}.bam'\n\n" >> ${name}_piece.log

    samtools view -bh -q 30 -F 4 "$aligned_bam" \
        | samtools sort -@ $cores_sf - \
        > "${name}.bam"

    printf "    Process stats:
        sort_and_filter starting reads: \$(samtools view -c $aligned_bam)
        sort_and_filter ending reads  : \$(samtools view -c ${name}.bam)\n\n" >> ${name}_piece.log
    printf "** End process 'sort_and_filter' at: \$(date)\n\n" >> ${name}_piece.log

    cp ${name}_piece.log ${name}_sf.txt
    cat ${logfile} > ${name}_sf.log
    cat ${name}_piece.log >> ${name}_sf.log

    """
}

log_pieces
    .groupTuple()
    .set { logs_to_combine }

process combine_logs {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs'
    memory '1 GB'

    input:
        file log_piece1
        set val(key), file(log2), file(log3), file(log4) from logs_to_combine

    output:
        set val(key), file("*_pre.log") into log_premerge

    """
    cat $log_piece1 $log2 $log3 $log4 > ${key}_pre.log

    """
}

sorted_bams
    .groupTuple()
    .set { bams_to_merge }
    
log_premerge.join(bams_to_merge).set{for_merge_bams}

save_bam = {params.output_dir + "/" + it - ~/.bam/ + "/" + it}

process merge_bams {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:samtools/1.4'
    publishDir = [path: "${params.output_dir}/", saveAs: save_bam, pattern: "*.bam", mode: 'copy' ]

    input:
        set key, file(logfile), file(bam_set) from for_merge_bams

    output:
        set key, file("*.bam"), file("merge_bams.log") into sample_bams
        set key, file("*.read_count.txt") into read_count

    """
    cat ${logfile} > merge_bams.log
    printf "** Start process 'merge_bams' at: \$(date)\n\n" >> merge_bams.log
    printf "    Process versions: 
        \$(samtools --version | tr '\n' ' ')\n\n" >> merge_bams.log
    printf "    Process command: 
        samtools merge ${key}.bam $bam_set\n\n" >> merge_bams.log

    samtools merge ${key}.bam $bam_set

    samtools view "${key}.bam" \
    | cut -d '|' -f 2 \
    | datamash -g 1 count 1 \
    | sort -k1,1 -S 5G \
    | datamash -g 1 sum 2 \
    > "${key}.read_count.txt"

    printf "** End process 'merge_bams' at: \$(date)\n\n" >> merge_bams.log
    """
}

sample_bams.join(gtf_info).set{assign_prepped}

process split_bam {
    cache 'lenient'
    memory '1 GB'
    module 'java/latest:modules:modules-init:modules-gs:samtools/1.4:bamtools/2.2.3:bedtools/2.26.0:python/3.6.4'

    input:
        set key, file(input_bam), file(logfile), val(gtf_path) from assign_prepped

    output:
        set key, file("split_bams/*.bam"), val(gtf_path) into split_bams mode flatten
        set key, file("remove_dups.log") into bam_and_log
        file input_bam into output 
        
    """
    cat ${logfile} > remove_dups.log
    printf "** Start processes 'remove duplicates, assign_genes, umi_rollup' at: \$(date)\n\n" >> remove_dups.log
    printf "    Process versions:
        \$(bedtools --version)
        \$(samtools --version | tr '\n' ' ')
        \$(bamtools --version | grep bamtools)
        \$(python --version)\n\n" >> remove_dups.log

    echo '    Process command:
        mkdir split_bams
        bamtools split -in $input_bam -reference -stub split_bams/split
        
        rmdup.py --bam in_bam --output_bam out.bam

        bedtools bamtobed -i out.bam -split \
                | sort -k1,1 -k2,2n -k3,3n -S 5G \
                > "in_bam.bed"

        bedtools map \
            -a in_bam.bed \
            -b exon_index \
            -nonamecheck -s -f 0.95 -c 7 -o distinct -delim "|" \
        | bedtools map \
            -a - -b gene_index \
            -nonamecheck -s -f 0.95 -c 4 -o distinct -delim "|" \
        | sort -k4,4 -k2,2n -k3,3n -S 5G\
        | datamash \
            -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8 \
        | assign-reads-to-genes.py gene_index \
        | awk \$3 == "exonic" || \$3 == "intronic" {{
                split(\$1, arr, "|")
                printf "%s_%s_%s\t%s\t%s\\n", arr[3], arr[4], arr[5], \$2, \$3
        }} \
        | sort -k2,2 -k1,1 -S 5G > in_bam.txt

        cat in_bed > key.bed
        sort -m -k2,2 -k1,1 in_assign | datamash -g 1,2 count 2 \
        | gzip > key.gz
        ' >> remove_dups.log 


    printf "    Process stats:
        remove_dups starting reads: \$(samtools view -c $input_bam)" >> remove_dups.log

    mkdir split_bams
    bamtools split -in $input_bam -reference -stub split_bams/split
    cd split_bams
    ls | grep "_[0-9A-Za-z\\.]\\{3,\\}.bam\$" | samtools merge split.REFnonstand.bam -b -
    ls | grep "_[0-9A-Za-z\\.]\\{3,\\}.bam\$" | xargs -d"\\n" rm
    mv split.REFnonstand.bam split.REF_nonstand.bam
    """

}

/**
Assign genes:
1. First use bedtools map to map the dedupped bed file to all exons with options:
-s forced strandedness
-f 0.95 95% of read must overlap exon
-c 7 map the name of the gene
-o distinct concatenate list of gene names
-delim "|" custom delimiter
-nonamecheck Don't error if there are different naming conventions for the chromosomes
2. Use bedtools map to map output to gene index
-s forced strandedness
-f 0.95 95% of read must overlap exon
-c 4 map the name of the cell name
-o distinct concatenate list of gene names
-delim "|" custom delimiter
-nonamecheck Don't error if there are different naming conventions for the chromosomes
3. Sort and collapse
4. Run assign-reads-to-genes.py to deal with exon v intron
**/

process remove_dups_assign_genes {
    cache 'lenient'
    memory '3 GB'
    module 'java/latest:modules:modules-init:modules-gs:bedtools/2.26.0:python/3.6.4:coreutils/8.24'

    input:
        set key, val(gtf_path) from split_bams

    output:
        set key, file("*.bed"), file("*.txt") into remove_dup_part_out

    """
    rmdup.py --bam $in_bam --output_bam out.bam

    bedtools bamtobed -i out.bam -split \
            | sort -k1,1 -k2,2n -k3,3n -S 5G \
            > "${in_bam}.bed"

    bedtools map \
        -a "${in_bam}.bed" \
        -b "${gtf_path}/latest.exons.bed" \
        -nonamecheck -s -f 0.95 -c 7 -o distinct -delim '|' \
    | bedtools map \
        -a - -b "${gtf_path}/latest.genes.bed" \
        -nonamecheck -s -f 0.95 -c 4 -o distinct -delim '|' \
    | sort -k4,4 -k2,2n -k3,3n -S 5G\
    | datamash \
        -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8 \
    | assign-reads-to-genes.py "${gtf_path}/latest.genes.bed" \
    | awk '\$3 == "exonic" || \$3 == "intronic" {{
            split(\$1, arr, "|")
            printf "%s|%s_%s_%s\t%s\t%s\\n", arr[2], arr[3], arr[4], arr[5], \$2, \$3
    }}' \
    | sort -k2,2 -k1,1 -S 5G > "${in_bam}.txt"

    """

}

remove_dup_part_out
    .groupTuple()
    .join(bam_and_log)
    .set { for_cat_dups }

/**
make cell, gene table for each read
sort
count instances of the gene for each read
**/

process umi_rollup {
    cache 'lenient'
    memory '1 GB'
    module 'java/latest:modules:modules-init:modules-gs:coreutils/8.24'

    input:
        set key, file(in_bed), file(in_assign), file(logfile) from for_cat_dups

    output:
        set key, file("*.gz"), file("*_ga.txt"), file("*.bed"), file("umi_rollup.log") into umi_rollup_out


    """
    cat ${logfile} > umi_rollup.log

    cat $in_bed > "${key}.bed"
    sort -m -k2,2 -k1,1 $in_assign > "${key}_ga.txt"
    datamash -g 1,2 count 2 < "${key}_ga.txt" \
    | gzip > "${key}.gz"

    printf "
        remove_dups ending reads  : \$(wc -l ${key}.bed | awk '{print \$1;}')\n\n
        Read assignments:\n\$(awk '{count[\$3]++} END {for (word in count) { printf "            %-20s %10i\\n", word, count[word]}}' $in_bed)\n\n" >> umi_rollup.log

    printf "** End processes 'remove duplicates, assign_genes, umi_rollup' at: \$(date)\n\n" >> umi_rollup.log
    
    """
}

save_umi_per_cell = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.txt/ + "/umis_per_cell_barcode.txt"}
save_umi_per_int = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.intronic.txt/ + "/intronic_umis_per_cell_barcode.txt"}

/**
Count intronic and total umis per cell
**/

process umi_by_sample_summary {
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'
    cache 'lenient'
    memory '8 GB'

    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_int, pattern: "*intronic.txt", mode: 'copy' 
    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_cell, pattern: "*barcode.txt", mode: 'copy'  
  
    input:
        set val(key), file(umi_rollup), file(gene_assignments_file), file(input_bed), file(logfile) from umi_rollup_out  

    output:
        set key, file(umi_rollup), file(input_bed), file("umi_by_sample_summary.log") into ubss_out
        set key, file("*UMIs.per.cell.barcode.txt") into umis_per_cell_barcode
        file "*UMIs.per.cell.barcode.intronic.txt" into umi_per_cell_intronic

    """
    cat ${logfile} > umi_by_sample_summary.log
    printf "** Start process 'umi_by_sample_summary' at: \$(date)\n\n" >> umi_by_sample_summary.log
    printf "    Process versions: 
        \$(python --version)\n\n" >> umi_by_sample_summary.log
    printf "    Process command:  
        tabulate_per_cell_counts.py 
            --gene_assignment_files "$gene_assignments_file" 
            --all_counts_file "${key}.UMIs.per.cell.barcode.txt" 
            --intron_counts_file "${key}.UMIs.per.cell.barcode.intronic.txt"\n\n"      >> umi_by_sample_summary.log

    tabulate_per_cell_counts.py \
        --gene_assignment_files "$gene_assignments_file" \
        --all_counts_file "${key}.UMIs.per.cell.barcode.txt" \
        --intron_counts_file "${key}.UMIs.per.cell.barcode.intronic.txt"


    printf "    Process stats:
        Total cells                            : \$(wc -l ${key}.UMIs.per.cell.barcode.txt | awk '{print \$1;}')
        Total cells > 100 reads                : \$(awk '\$3>100{c++} END{print c+0}' ${key}.UMIs.per.cell.barcode.txt)
        Total cells > 1000 reads               : \$(awk '\$3>1000{c++} END{print c+0}' ${key}.UMIs.per.cell.barcode.txt)
        Total reads in cells with > 100 reads  : \$(awk '\$3>100{c=c+\$3} END{print c+0}' ${key}.UMIs.per.cell.barcode.txt)\n\n" >> umi_by_sample_summary.log

    printf "** End process 'umi_by_sample_summary' at: \$(date)\n\n" >> umi_by_sample_summary.log
    """
}

save_umi = {params.output_dir + "/" + it - ~/.umi_counts.mtx/ + "/umi_counts.mtx"}
save_cell_anno = {params.output_dir + "/" + it - ~/.cell_annotations.txt/ + "/cell_annotations.txt"}
save_gene_anno = {params.output_dir + "/" + it - ~/.gene_annotations.txt/ + "/gene_annotations.txt"}

ubss_out.join(gtf_info2).set{make_matrix_prepped}
/**
sum up total assigned reads per cell and keep only those above the cutoff
make the number matrix
**/

process make_matrix {
    cache 'lenient'
    memory '15 GB'
    publishDir path: "${params.output_dir}/", saveAs: save_umi, pattern: "*umi_counts.mtx", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_anno, pattern: "*cell_annotations.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_gene_anno, pattern: "*gene_annotations.txt", mode: 'copy'

    input:
        set key, file(umi_rollup_file), file(input_bed), file(logfile), val(gtf_path) from make_matrix_prepped

    output:
        set key, file("*cell_annotations.txt"), file("*umi_counts.mtx"), file("*gene_annotations.txt"), val(gtf_path), file(input_bed), file("make_matrix.log") into mat_output

    """
    cat ${logfile} > make_matrix.log
    printf "** Start process 'make_matrix' at: \$(date)\n\n" >> make_matrix.log

    echo '    Process command:
        make_matrix.py <(zcat $umi_rollup_file) --gene_annotation "${gtf_path}/latest.gene.annotations" --key "$key"   
        cat ${gtf_path}/latest.gene.annotations > "${key}.gene_annotations.txt"  ' >> make_matrix.log
   
    make_matrix.py <(zcat $umi_rollup_file) --gene_annotation "${gtf_path}/latest.gene.annotations" --key "$key" 
    cat "${gtf_path}/latest.gene.annotations" > "${key}.gene_annotations.txt"

    printf "\n** End process 'make_matrix' at: \$(date)\n\n" >> make_matrix.log
    """

}

process make_cds {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:gcc/8.1.0:R/3.6.1'
    memory '15 GB'

    input:
        set key, file(cell_data), file(umi_matrix), file(gene_data), val(gtf_path), file(input_bed), file(logfile) from mat_output

    output:
        set key, file("*for_scrub.mtx"), file("*.RDS"), file("*cell_qc.csv"), file(input_bed), file("make_cds.log") into for_scrub
        file("*cell_qc.csv") into cell_qcs

"""
    cat ${logfile} > make_cds.log
    printf "** Start process 'make_cds' at: \$(date)\n\n" >> make_cds.log
    printf "    Process versions: 
        \$(R --version | grep 'R version')
            monocle3 version \$(Rscript -e 'packageVersion("monocle3")')\n\n" >> make_cds.log
    echo '    Process command:  
        make_cds.R \
            "$umi_matrix"\
            "$cell_data"\
            "$gene_data"\
            "${gtf_path}/latest.genes.bed"\
            "$key" ' >> make_cds.log

    make_cds.R \
        "$umi_matrix"\
        "$cell_data"\
        "$gene_data"\
        "${gtf_path}/latest.genes.bed"\
        "$key"

    printf "** End process 'make_cds' at: \$(date)\n\n" >> make_cds.log
"""

}

save_hist = {params.output_dir + "/" + it - ~/_scrublet_hist.png/ + "/" + it}

process run_scrublet {
    publishDir path: "${params.output_dir}/", saveAs: save_hist, pattern: "*png", mode: 'copy'
    cache 'lenient'
    module 'modules:java/latest:modules-init:modules-gs:python/3.6.4'
    memory '10 GB'
    input:
        set key, file(scrub_mat), file(cds), file(cell_qc), file(input_bed), file(logfile) from for_scrub
    output:
        set key, file("*scrublet_out.csv"), file(cds), file(cell_qc), file(input_bed), file("run_scrublet.log") into scrublet_out
        file ("*.png") into scrub_pngs


"""
    cat ${logfile} > run_scrublet.log
    printf "** Start process 'run_scrublet' at: \$(date)\n\n" >> run_scrublet.log
    printf "    Process versions: 
        \$(python --version)
            \$(pip freeze | grep scrublet | tr '==' ' ')\n\n" >> run_scrublet.log
    echo '    Process command:  
        run_scrublet.py --key $key --mat $scrub_mat\n'  >> run_scrublet.log

    run_scrublet.py --key $key --mat $scrub_mat

    printf "\n** End process 'run_scrublet' at: \$(date)\n\n" >> run_scrublet.log
"""

}


scrublet_out.join(read_count).set{umi_by_sample_in}

process umi_by_sample {
    cache 'lenient'
    memory '20 GB'    
    module 'java/latest:modules:modules-init:modules-gs:samtools/1.4:coreutils/8.24'

    input:
        set key, file(scrublet_outs), file(cds), file(cell_qc), file(input_bed), file(logfile), file(read_count) from umi_by_sample_in

    output:
        set file("*.UMI_count.txt"), file("*.read_count.txt") into summarize_dup_out
        set key, file(scrublet_outs), file(cds), file(cell_qc), file("*duplication_rate_stats.txt"), file("umi_by_sample.log") into duplication_rate_out

    """

    cat ${logfile} > umi_by_sample.log
    printf "** Start process 'umi_by_sample' at: \$(date)\n\n" >> umi_by_sample.log
    printf "    Process versions: 
        \$(samtools --version | tr '\n' ' ')\n\n" >> umi_by_sample.log
    echo '    Process command:      
    awk "{{ split(\$4, arr, "|")
            if (!seen[arr[1]]) {{
                seen[arr[1]] = 1; count[arr[2]]++;
            }}
            }} END {{
                for (sample in count) {{
                print sample "\\t" count[sample]
                }}
            }}" "$input_bed" \
    | sort -k1,1 -S 5G\
    >"${key}.UMI_count.txt"

    cat ${key}.UMI_count.txt \
    | join - "${key}.read_count.txt" \
    | awk "{{ 
            printf "%-18s   %10d    %10d    %7.1f%\\n",
                \$1, \$3, \$2, 100 * (1 - \$2/\$3);
    }}" \
    >"${key}.duplication_rate_stats.txt" \n'  >> umi_by_sample.log


       awk '{{ split(\$4, arr, "|")
            if (!seen[arr[1]]) {{
                seen[arr[1]] = 1; count[arr[2]]++;
            }}
            }} END {{
                for (sample in count) {{
                print sample "\\t" count[sample]
                }}
            }}' "$input_bed" \
    | sort -k1,1 -S 5G\
    >"${key}.UMI_count.txt"


    cat ${key}.UMI_count.txt \
    | join - "$read_count" \
    | awk '{{ 
            printf "%-18s   %10d    %10d    %7.1f%\\n",
                \$1, \$3, \$2, 100 * (1 - \$2/\$3);
    }}' \
    >"${key}.duplication_rate_stats.txt"

    printf "\n** End process 'umi_by_sample' at: \$(date)\n\n" >> umi_by_sample.log

    printf "** Start processes to generate cds, qc metrics and dashboard at: \$(date)\n\n" >> umi_by_sample.log
    """
}

save_cds = {params.output_dir + "/" + it - ~/_cds.RDS/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_cell_qc = {params.output_dir + "/" + it - ~/_cell_qc.csv/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_samp_stats = {params.output_dir + "/" + it - ~/_sample_stats.csv/ + "/" + it}

process reformat_scrub {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4:gcc/8.1.0:R/3.6.1'
    memory '10 GB'
    publishDir path: "${params.output_dir}/", saveAs: save_cds, pattern: "temp_fold/*cds.RDS", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_qc, pattern: "temp_fold/*cell_qc.csv", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_samp_stats, pattern: "*sample_stats.csv", mode: 'copy'

    input:
        set key, file(scrublet_outs), file(cds), file(cell_qc), file(dup_stats), file(logfile) from duplication_rate_out

    output: 
        set key, file("temp_fold/*.RDS"),  file("temp_fold/*.csv") into rscrub_out
        file("*sample_stats.csv") into sample_stats
        file("*collision.txt") into barn_collision
        set key, file(logfile) into pipe_log

"""
#!/usr/bin/env Rscript
library(monocle3)
dir.create("temp_fold")
cds <- readRDS("$cds")
cell_qc <- read.csv("$cell_qc")
if(nrow(pData(cds)) > 0) {
scrublet_out <- read.csv("$scrublet_outs", header=F)
pData(cds)\$scrublet_score <- scrublet_out\$V1
pData(cds)\$scrublet_call <- ifelse(scrublet_out\$V2 == 1, "Doublet", "Singlet")
cell_qc\$scrublet_score <- scrublet_out\$V1
cell_qc\$scrublet_call <- ifelse(scrublet_out\$V2 == 1, "Doublet", "Singlet")
}
write.csv(cell_qc, quote=FALSE, file="temp_fold/$cell_qc")

dup_stats <- read.table(paste0("$key", ".duplication_rate_stats.txt"))

df <- data.frame(sample="$key", n.reads = dup_stats\$V2, n.umi = dup_stats\$V3, duplication_rate = dup_stats\$V4,
                 doublet_count = sum(cell_qc\$scrublet_call == "Doublet", na.rm=TRUE),
                 doublet_perc = paste0(round(sum(cell_qc\$scrublet_call == "Doublet", na.rm=TRUE)/nrow(cell_qc) * 100, 1), "%"),
                 doublet_NAs=sum(is.na(cell_qc\$scrublet_call)))

write.csv(df, file=paste0("$key", "_sample_stats.csv"), quote=FALSE, row.names=FALSE)
saveRDS(cds, file="temp_fold/$cds")

if ("$key" == "Barnyard") {
  fData(cds)\$mouse <- grepl("ENSMUSG", fData(cds)\$id)
  fData(cds)\$human <- grepl("ENSG", fData(cds)\$id)

  pData(cds)\$mouse_reads <- Matrix::colSums(exprs(cds)[fData(cds)\$mouse,])
  pData(cds)\$human_reads <- Matrix::colSums(exprs(cds)[fData(cds)\$human,])
  pData(cds)\$total_reads <- pData(cds)\$mouse_reads + pData(cds)\$human_reads
  pData(cds)\$human_perc <- pData(cds)\$human_reads/pData(cds)\$total_reads
  pData(cds)\$mouse_perc <- pData(cds)\$mouse_reads/pData(cds)\$total_reads
  pData(cds)\$collision <- ifelse(pData(cds)\$human_perc >= .9 | pData(cds)\$mouse_perc >= .9, FALSE, TRUE)

  collision_rate <- round(sum(pData(cds)\$collision/nrow(pData(cds))) * 200, 1)
  fileConn<-file("Barn_collision.txt")
  writeLines(paste0("$key", "\t", collision_rate, "%"), fileConn)
  close(fileConn)

} else {
  fileConn<-file("no_collision.txt")
  writeLines(paste0("$key", "\t", "NA", fileConn)
  close(fileConn)
}
"""

}

for_gen_qc = rscrub_out.join(umis_per_cell_barcode)
save_knee = {params.output_dir + "/" + it - ~/_knee_plot.png/ + "/" + it}
save_umap = {params.output_dir + "/" + it - ~/_UMAP.png/ + "/" + it}
save_cellqc = {params.output_dir + "/" + it - ~/_cell_qc.png/ + "/" + it}

process generate_qc_metrics {
    module 'java/latest:modules:modules-init:modules-gs:gcc/8.1.0:R/3.6.1'
    memory '10 GB'
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_umap, pattern: "*UMAP.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_knee, pattern: "*knee_plot.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cellqc, pattern: "*cell_qc.png", mode: 'copy'

    input:
        set key, file(cds), file(cell_qc), file(umis_per_cell) from for_gen_qc

    output:
        file("*.png") into qc_plots
        file("*.txt") into cutoff
"""
mkdir temp2
generate_qc.R\
    $cds $umis_per_cell $key \
    --specify_cutoff $params.umi_cutoff\
"""

}

process zip_up_duplication {
    cache 'lenient'
    publishDir = [path: "${params.output_dir}/", pattern: "all_sample_stats.csv", mode: 'copy']

    input:
        file files from sample_stats.collect()
    output:
        file "*ll_sample_stats.csv" into all_dups
    """
     sed -s 1d $files > all_sample_stats.csv
    """      
}

process calc_cell_totals {
    module 'java/latest:modules:modules-init:modules-gs'
    memory '1 GB'
    cache 'lenient'

    input:
        file cell_qcs from cell_qcs.collect()

    output:
        file "*.txt" into cell_counts

"""
    for f in *.csv
    do
      awk 'BEGIN {FS=","}; \$2>100{c++} END{print FILENAME, "100", c-1}' \$f >> cell_counts.txt
      awk 'BEGIN {FS=","}; \$2>500{c++} END{print FILENAME, "500", c-1}' \$f >> cell_counts.txt
      awk 'BEGIN {FS=","}; \$2>1000{c++} END{print FILENAME, "1000", c-1}' \$f >> cell_counts.txt
    done
"""

}

process collapse_collision {
    module 'java/latest:modules:modules-init:modules-gs'
    memory '1 GB'
    cache 'lenient'

    input:
        file col_file from barn_collision.collect()

    output:
        file "*.txt" into all_collision

"""
cat $col_file > all_collision.txt


"""

}

process generate_dash_info {
    module 'java/latest:modules:modules-init:modules-gs:gcc/8.1.0:R/3.6.1'
    cache 'lenient'
    memory '1 GB'

    input:
        file(dup_file) from all_dups
        file(cell_counts) from cell_counts
        file(barn_col) from all_collision       
    output:
        file("*.js") into run_data

"""
#!/usr/bin/env Rscript

library(jsonlite)

all_dups <- read.csv("$dup_file", header=FALSE, stringsAsFactors=FALSE)
output_folder <- "$params.output_dir"
count_info <- read.table("$cell_counts", stringsAsFactors=FALSE)

project_name <- unlist(stringr::str_split(output_folder, "/"))
project_name <- project_name[[length(project_name)]]

c100 <- sum(count_info[count_info\$V2 == 100,]\$V3)
c500 <- sum(count_info[count_info\$V2 == 500,]\$V3)
c1000 <- sum(count_info[count_info\$V2 == 1000,]\$V3)

count_info_tab <- count_info

count_info_tab\$V1 <- gsub("_cell_qc.csv", "", count_info_tab\$V1)
ct100 <- count_info_tab[count_info_tab\$V2 == 100,]
ct500 <- count_info_tab[count_info_tab\$V2 == 500,]
ct1000 <- count_info_tab[count_info_tab\$V2 == 1000,]
row.names(ct100) <- ct100\$V1
row.names(ct500) <- ct500\$V1
row.names(ct1000) <- ct1000\$V1

all_dups\$V1 <- as.character(all_dups\$V1)
all_dups\$c100 <- ct100[all_dups\$V1,"V3"]
all_dups\$c1000 <- ct1000[all_dups\$V1,"V3"]

all_dups\$V5[all_dups\$V7 > 0] <- "Fail"
all_dups\$V6[all_dups\$V7 > 0] <-  "Fail"
all_dups\$V6[all_dups\$V6 == "NaN%"] <-  "Fail"

row.names(all_dups) <- all_dups\$V1
names(all_dups) <- c("Sample", "Total_reads",
                     "Total_UMIs",
                     "Duplication_rate",
                     "Doublet_Number", 
                     "Doublet_Percent",
                     "Doublet_NAs",
                     "Cells_100_UMIs",
                     "Cells_1000_UMIs" 
                     )

all_dups\$Doublet_Number[is.na(all_dups\$Doublet_Number)] <- "Fail"
all_dups\$Doublet_Percent[is.na(all_dups\$Doublet_Percent)] <- "Fail"


all_dup_lst <- apply(all_dups, 1, as.list)
sample_list <- as.character(all_dups\$Sample)[order(as.character(all_dups\$Sample))]

if("Sentinel" %in% sample_list) {
  sample_list <- c("Sentinel", setdiff(sample_list, c("Sentinel")))
}
barn_collision <- NA
if("Barnyard" %in% sample_list) {
  sample_list <- c("Barnyard", setdiff(sample_list, c("Barnyard")))
  barn_collision <- read.table("$barn_col")
  barn_collision <- barn_collision[barn_collision$V1 == "Barnyard", barn_collision$V2]
}
sample_list <- as.list(sample_list)
json_info <- list("run_name" = project_name,
                  "cell_counts" = c(sum(all_dups\$Cells_100_UMIs), sum(all_dups\$Cells_1000_UMIs)),
                  "sample_list" = sample_list,
                  "barn_collision" = barn_collision,
                  "sample_stats" = all_dup_lst)
                  
fileConn<-file("data.js")
writeLines(c("const run_data =", toJSON(json_info, pretty=TRUE, auto_unbox=TRUE)), fileConn)
close(fileConn)

"""
}

process exp_dash {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:gcc/8.1.0:R/3.6.1'
    memory '1 GB'

    publishDir path: "${params.output_dir}/", pattern: "exp_dash", mode: 'copy'


    input:
        file plots from qc_plots.collect()
        file run_data
        file scrub_png from scrub_pngs.collect()

    output:
        file exp_dash into exp_dash_out

    """
    mkdir exp_dash
    cp -R $baseDir/bin/skeleton_dash/* exp_dash/
    mv *.png exp_dash/img/

    mv $run_data exp_dash/js/
    
    """
}

save_logs = {params.output_dir + "/" + it - ~/_read_metrics.log/ - ~/_full.log/ + "/" + it}
save_json = {params.output_dir + "/" + it - ~/_log_data.js/ + "/" + it}
process generate_summary_log {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_logs, pattern: "*.log", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_json, pattern: "*.js", mode: 'copy'
    
    input:
        set key, file(logfile) from pipe_log
        file exp_dash from exp_dash_out

    output:
        file("*_full.log") into full_log
        file("*_read_metrics.log") into sum_log
        file("*log_data.js") into log_json

    """
    head -n 2 ${logfile} > ${key}_full.log
    printf "Git Version, Commit ID, Session ID: $workflow.revision, $workflow.commitId, $workflow.sessionId\n" >> ${key}_full.log
    printf "Command:\n$workflow.commandLine\n\n" >> ${key}_full.log
    printf "***** PARAMETERS *****: \n\n" >> ${key}_full.log
    printf "    params.run_dir:               $params.run_dir\n" >> ${key}_full.log
    printf "    params.output_dir:            $params.output_dir\n" >> ${key}_full.log
    printf "    params.sample_sheet:          $params.sample_sheet\n" >> ${key}_full.log
    printf "    params.p7_rows:               $params.p7_rows\n" >> ${key}_full.log
    printf "    params.p5_cols:               $params.p5_cols\n" >> ${key}_full.log
    printf "    params.demux_out:             $params.demux_out\n" >> ${key}_full.log
    printf "    params.level:                 $params.level\n" >> ${key}_full.log
    printf "    params.max_cores:             $params.max_cores\n" >> ${key}_full.log
    printf "    params.samples:               $params.samples\n" >> ${key}_full.log
    printf "    params.star_file:             $params.star_file\n" >> ${key}_full.log
    printf "    params.gene_file:             $params.gene_file\n" >> ${key}_full.log
    printf "    params.umi_cutoff:            $params.umi_cutoff\n" >> ${key}_full.log
    printf "    params.rt_barcode_file:       $params.rt_barcode_file\n\n" >> ${key}_full.log
    printf "    params.hash_list:             $params.hash_list\n\n" >> ${key}_full.log
    printf "    params.max_wells_per_sample:  $params.max_wells_per_sample\n\n" >> ${key}_full.log

    tail -n +2 ${logfile} >> ${key}_full.log
    printf "\n** End processes generate cds, qc metrics and dashboard at: \$(date)\n\n" >> ${key}_full.log
    printf "***** END PIPELINE *****: \n\n" >> ${key}_full.log
    filename=${key}_full.log

    # Trimming:
    trim_start=`cat \$filename | grep 'sequences processed in total' | awk -F ' ' '{sum += \$1} END {print sum}'`
    trim_lost=`cat \$filename | grep 'Sequences removed because they became shorter' | awk -F ' ' '{sum += \$14} END {print sum}'`
    trim_end=\$((\$trim_start - \$trim_lost))

    # Alignment:
    align_start=`cat \$filename | grep 'Number of input reads' | awk -F '|' '{sum += \$2} END {print sum}'`
    align_mapped=`cat \$filename | grep 'Uniquely mapped reads number' | awk -F '|' '{sum += \$2} END {print sum}'`
    align_totals=(\$(cat \$filename | grep 'Number of input reads' | cut -d "|" -f 2 | awk '{print \$1}'))
    align_multimapped=`cat \$filename | grep 'Number of reads mapped to multiple loci' |  awk -F '|' '{sum += \$2} END {print sum}'`
    align_too_short_arr=(\$(cat \$filename | grep 'unmapped: too short' | cut -d "|" -f 2 | tr '%' ' ' | awk '{\$1=\$1/100;print}'))
    align_too_short=`a=0
    for i in \${align_too_short_arr[@]}
    do
        echo "\${align_too_short_arr[\$a]} * \${align_totals[\$a]}" | bc 
        a=\$((a+1))
    done | awk '{sum += \$1} END {print sum}'`

    # Sort and Filter:
    sf_start=`cat \$filename | grep 'sort_and_filter starting reads' | awk -F ':' '{sum += \$2} END {print sum}'`
    sf_end=`cat \$filename | grep 'sort_and_filter ending reads' | awk -F ':' '{sum += \$2} END {print sum}'`

    # Dups:
    dup_start=`cat \$filename | grep 'remove_dups starting reads' | awk -F ':' '{sum += \$2} END {print sum}'`
    dup_end=`cat \$filename | grep 'remove_dups ending reads' | awk -F ':' '{sum += \$2} END {print sum}'`

    # Assignment:
    assigned_exonic=`cat \$filename | grep '    exonic     ' | awk -F ' ' '{sum += \$2} END {print sum}'`
    assigned_intronic=`cat \$filename | grep '    intronic     ' | awk -F ' ' '{sum += \$2} END {print sum}'`
    assigned_end=\$((\$assigned_exonic + \$assigned_intronic))

    # In real cells:
    reads_in_cells=`cat \$filename | grep 'Total reads in cells with > 100 reads' | awk -F ':' '{sum += \$2} END {print sum}'`
    
    printf "
        const log_data = { 
            "alignment_start" : \$align_start,
            "alignment_mapped" : \$align_mapped,
            "align_multimapped" : \$align_multimapped,
            "align_too_short" : \$align_too_short,
            "sf_start" : \$sf_start,
            "sf_end" : \$sf_end,
            "dup_start" : \$dup_start,
            "dup_end" : \$dup_end,
            "assigned_exonic" : \$assigned_exonic,
            "assigned_intronic" : \$assigned_intronic,
            "reads_in_cells" : \$reads_in_cells }
    " > ${key}_log_data.js


    printf "***** PIPELINE READ STATS *****: \n\n" >> ${key}_read_metrics.log

    printf "%20s %20s %20s %20s %20s\n" "Process" "Starting reads" "Ending reads" "% lost" "% of total lost" >> ${key}_read_metrics.log
    printf "========================================================================================================\n" >> ${key}_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f\n" "Trimming" \$trim_start \$trim_end \$(echo "(\$trim_start - \$trim_end)/\$trim_start * 100" | bc -l ) \$(echo "(\$trim_start - \$trim_end)/\$trim_start * 100" | bc -l ) >> ${key}_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f\n" "Alignment" \$align_start \$sf_start \$(echo "(\$align_start - \$sf_start)/\$align_start * 100" | bc -l ) \$(echo "(\$align_start - \$sf_start)/\$trim_start * 100" | bc -l ) >> ${key}_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f\n" "Filtering" \$sf_start \$sf_end \$(echo "(\$sf_start - \$sf_end)/\$sf_start * 100" | bc -l ) \$(echo "(\$sf_start - \$sf_end)/\$trim_start * 100" | bc -l ) >> ${key}_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f\n" "Deduplication" \$dup_start \$dup_end \$(echo "(\$dup_start - \$dup_end)/\$dup_start * 100" | bc -l ) \$(echo "(\$dup_start - \$dup_end)/\$trim_start * 100" | bc -l ) >> ${key}_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f\n" "Gene assignment" \$dup_end \$assigned_end \$(echo "(\$dup_end - \$assigned_end)/\$dup_end * 100" | bc -l ) \$(echo "(\$dup_end - \$assigned_end)/\$trim_start * 100" | bc -l ) >> ${key}_read_metrics.log
 
    printf "\nAlignment details: \n" >> ${key}_read_metrics.log
    printf "%25s %20s %20s\n" "" "Count" "Percent"  >> ${key}_read_metrics.log
    printf "========================================================================================================\n" >> ${key}_read_metrics.log
    printf "%25s %20s %20s\n" "Total reads processed:" \$align_start "" >> ${key}_read_metrics.log
    printf "%25s %20s %20.2f\n" "Reads uniquely mapped:" \$align_mapped \$(echo "(\$align_mapped)/\$align_start * 100" | bc -l ) >> ${key}_read_metrics.log
    printf "%25s %20s %20.2f\n" "Reads multi-mapped:" \$align_multimapped \$(echo "(\$align_multimapped)/\$align_start * 100" | bc -l )  >> ${key}_read_metrics.log
    printf "%25s %20s %20.2f\n" "Reads too short:" \$align_too_short \$(echo "(\$align_too_short)/\$align_start * 100" | bc -l ) >> ${key}_read_metrics.log

 
    cat ${key}_read_metrics.log >> ${key}_full.log

    """

}

workflow.onComplete {
	println ( workflow.success ? "Done! Saving output" : "Oops .. something went wrong" )
}

