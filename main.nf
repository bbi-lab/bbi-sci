/*
** Check that Nextflow version meets minimum version requirements.
*/
def minMajorVersion = 20
def minMinorVersion = 0
checkNextflowVersion( minMajorVersion, minMinorVersion )


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
params.garnett_file = false 

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
    log.info '    params.garnett_file = false                Path to a csv with two columns, first is the sample name, and second is a path to the Garnett classifier to be applied to that sample. Default is false - no classification.'
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

 Pass through:

 Summary:
    Check and process sample sheet - check_sample_sheet.py
    Start log

 Downstream:
    gather_info
    trim_fastqs
    combine_logs

 Published:

 Notes:

*************/

process check_sample_sheet {
    cache 'lenient'
    module 'modules:modules-init:modules-gs:python/3.6.4'

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
        check_sample_sheet.py 
            --sample_sheet $params.sample_sheet 
            --star_file $params.star_file
            --level $params.level --rt_barcode_file $params.rt_barcode_file
            --max_wells_per_samp $params.max_wells_per_sample\n\n" >> start.log


    check_sample_sheet.py --sample_sheet $params.sample_sheet --star_file $params.star_file \
        --level $params.level --rt_barcode_file $params.rt_barcode_file \
        --max_wells_per_samp $params.max_wells_per_sample


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
    samp_list = samp_list.intersect(params.samples.collect{"$it".replaceAll(/\s/, ".").replaceAll(/_/, ".").replaceAll(/-/, ".").replaceAll(/\\//, ".")})
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

 Pass through:
    input_fastq - all fastq files from params.demux_out folder

 Summary:
    Trim fastqs - trim_galore
    Continue log

 Downstream:
    gather_info
    process_hashes

 Published:

 Notes:
    Only moves forward if sample is in samp_list - where params.samples comes in

*************/

process trim_fastqs {
    cache 'lenient'
    memory '12G'
    module 'modules:modules-init:modules-gs:python/2.7.3:cutadapt/1.8.3:trim_galore/0.4.1'

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


    cat trim_out/*trimming_report.txt | sed '/Overview of/,/RUN/{//!d}' | sed 's/Overview of removed sequences//' >> piece.log
    printf "** End process 'trim_fastqs' at: \$(date)\n\n" >> piece.log
    cp piece.log ${name}_trim.txt
    cat piece.log >> trim.log
    """
}


/*************

Process: gather_info

 Inputs:
    good_sample_sheet - corrected csv sample sheet
    key - sample id
    params.star_file
    params.gene_file

 Outputs:
    key - sample id
    star_path - path to star index folder
    gtf_path - path to gtf info folder
    star_mem - GB needed for star alignment

 Pass through:
    name - file id (including lane and split info)
    trimmed_fastq - trimmed, gzipped fastq
    log_piece2 - piece of log to be concatenated for full log
    logfile - running log

 Summary:
    Gather the star and gtf paths and info for downstream

 Downstream:
    split_bam
    make_matrix
    align_reads

 Published:

 Notes:

*************/

process gather_info {
    cache 'lenient'
    memory '1G'
    module 'modules:modules-init:modules-gs:python/3.6.4'

    input:
        file good_sample_sheet
        set val(key), val(name), file(trimmed_fastq), file(logfile), file(log_piece2) from trimmed_fastqs

    output:
        set val(key), val(name), env(star_path), env(star_mem), file(trimmed_fastq), file(logfile), file(log_piece2) into align_prepped
        set val(key), env(gtf_path) into gtf_info
        set val(key), env(gtf_path) into gtf_info2

    """

    spec=`awk 'BEGIN {FS=",";OFS=","};{gsub(" ", ".", \$2);gsub("/", ".", \$2);gsub("-", ".", \$2);split(\$2,a,"_fq_part");gsub("_", ".", \$2);print(\$1, a[1], \$3)}' $good_sample_sheet | awk 'BEGIN {FS=","}; \$2=="$key" {print \$3}' | uniq`
    star_mem=`awk -v var="\$spec" '\$1==var {print \$3}' $params.star_file | uniq`
    star_path=`awk -v var="\$spec" '\$1==var {print \$2}' $params.star_file | uniq`
    gtf_path=`awk -v var="\$spec" '\$1==var {print \$2}' $params.gene_file | uniq`

    """
}


/*************

Process: process_hashes

 Inputs:
    key - sample id
    input_fastq - all fastq files from params.demux_out folder
    params.hash_list

 Outputs:
    hash_log - log of the hash information
    hash_mtx - MatrixMarket matrix of hash information
    hash_cell - Cell info for matrix of hash information
    hash_hash - Hash info for matrix of hash information

 Pass through:

 Summary:
    Collect and process hash barcodes - process_hashes.py

 Downstream:

 Published:
    hash_mtx - MatrixMarket matrix of hash information
    hash_cell - Cell info for matrix of hash information
    hash_hash - Hash info for matrix of hash information

 Notes:
    Only when params.hash = true

*************/

// Group fastqs for finding hash barcodes
fastqs_out
    .groupTuple()
    .set { for_hash }

save_hash_cell = {params.output_dir + "/" + it - ~/.hashumis_cells.txt/ + "/" + it}
save_hash_hash = {params.output_dir + "/" + it - ~/.hashumis_hashes.txt/ + "/" + it}
save_hash_mtx = {params.output_dir + "/" + it - ~/.hashumis.mtx/ + "/" + it}

process process_hashes {
    cache 'lenient'
    memory '12G'
    module 'modules:modules-init:modules-gs:python/3.6.4'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_cell, pattern: "*hashumis_cells.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_hash, pattern: "*hashumis_hashes.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_mtx, pattern: "*.mtx", mode: 'copy'

    input:
        set key, file(input_fastq) from for_hash

    output:
        file("*hash.log") into hash_logs
        set file("*mtx"), file("*hashumis_cells.txt"), file("*hashumis_hashes.txt") into hash_mats

    when:
        params.hash_list != false

    """
    process_hashes.py --hash_sheet $params.hash_list \
        --fastq <(zcat $input_fastq) --key $key

    """
}


/*************

Process: align_reads

 Inputs:
    name - file id (including lane and split info)
    star_path - path to star index folder
    star_mem - GB needed for star alignment
    trimmed_fastq - trimmed, gzipped fastq
    logfile - running log
    cores_align - number of cores to use

 Outputs:
    align_out - folder of all alignment output - stops here
    name - file id (including lane and split info)
    aligned_bam - bam output from star alignment
    logfile - running log
    log_piece3 - piece of log to be concatenated for full log

 Pass through:
    key - sample id
    log_piece2 - piece of log to be concatenated for full log

 Summary:
    Align reads to genome using STAR

 Downstream:
    sort_and_filter

 Published:

 Notes:

*************/

// Cores for alignment set at 8 unless limit is lower
if (params.max_cores < 8) {
    cores_align = params.max_cores
} else {
    cores_align = 8
}

process align_reads {
    cache 'lenient'
    module 'modules:modules-init:modules-gs:STAR/2.5.2b'
    memory { star_mem.toInteger()/cores_align + " GB" }
    penv 'serial'
    cpus cores_align

    input:
        set val(key), val(name), val(star_path), val(star_mem), file(trimmed_fastq), file(logfile), file(log_piece2) from align_prepped

    output:
        file "align_out" into align_output
        set val(key), val(name), file("align_out/*Aligned.out.bam"), file('align.log'), file(log_piece2), file("*align.txt") into aligned_bams

    """
    cat ${logfile} > align.log
    printf "** Start process 'align_reads' for $trimmed_fastq at: \$(date)\n\n" > piece.log
    printf "    Process versions:
        \$(STAR --version)\n\n" >> piece.log

    printf "    Process command:
        STAR --runThreadN $cores_align --genomeDir $star_path
            --readFilesIn $trimmed_fastq --readFilesCommand zcat 
            --outFileNamePrefix ./align_out/${name} --outSAMtype BAM Unsorted 
            --outSAMmultNmax 1 --outSAMstrandField intronMotif\n

    Reference genome information:
      \$(grep fastq_url $star_path/../*gsrc/record.out | awk '{\$1=\$2=""; print \$0}')
        FASTA download date: \$(grep fastq_download_date $star_path/../*gsrc/record.out | awk '{\$1=\$2=""; print \$0}')
        Non REF sequences removed.

      \$(grep gtf_url $star_path/../*gsrc/record.out | awk '{\$1=\$2=""; print \$0}')
        GTF download date: \$(grep gtf_download_date $star_path/../*gsrc/record.out | awk '{\$1=\$2=""; print \$0}')
        \$(grep gtf_include_biotypes $star_path/record.out | awk '{\$1=\$2=""; print \$0}')

    Process output:\n" >> piece.log


    mkdir align_out
    STAR \
        --runThreadN $cores_align \
        --genomeDir $star_path \
        --readFilesIn $trimmed_fastq \
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


/*************

Process: sort_and_filter

 Inputs:
    name - file id (including lane and split info)
    aligned_bam - bam output from star alignment
    logfile - running log
    cores_sf - number of cores to use

 Outputs:
    name - file id (including lane and split info)
    sorted_bam - sorted and quality filtered bam
    log_piece4 - piece of log to be concatenated for full log

 Pass through:
    key - sample id
    log_piece2 - piece of log to be concatenated for full log
    log_piece3 - piece of log to be concatenated for full log

 Summary:
    Use samtools to filter for read quality 30 and sort bam

 Downstream:
    merge_bams
    combine_logs

 Published:

 Notes:

*************/

// Cores for sort and filter set at 10 unless limit is lower
if (params.max_cores < 10) {
    cores_sf = params.max_cores
} else {
    cores_sf = 10
}

process sort_and_filter {
    cache 'lenient'
    module 'modules:modules-init:modules-gs:samtools/1.4'
    memory '1 GB'
    penv 'serial'
    cpus cores_sf

    input:
        set val(key), val(name), file(aligned_bam), file(logfile), file(log_piece2), file(log_piece3) from aligned_bams

    output:
        set val(key), file("*.bam") into sorted_bams
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

    """
}


/*************

Process: combine_logs

 Inputs:
    log_piece1 - piece of log to be concatenated for full log
    log_piece2 - piece of log to be concatenated for full log
    log_piece3 - piece of log to be concatenated for full log
    log_piece4 - piece of log to be concatenated for full log

 Outputs:
    logfile - concatenated running log

 Pass through:
    key - sample id

 Summary:
    Combine the log pieces from the first steps in the correct order

 Downstream:
    merge_bams

 Published:

 Notes:

*************/

log_pieces
    .groupTuple()
    .set { logs_to_combine }

process combine_logs {
    cache 'lenient'
    module 'modules:modules-init:modules-gs'
    memory '1 GB'

    input:
        file log_piece1
        set val(key), file(log_piece2), file(log_piece3), file(log_piece4) from logs_to_combine

    output:
        set val(key), file("*_pre.log") into log_premerge

    """

    cat $log_piece1 $log_piece2 $log_piece3 $log_piece4 > ${key}_pre.log

    """
}


/*************

Process: merge_bams

 Inputs:
    key - sample id
    logfile - running log
    sorted_bam - sorted and quality filtered bam

 Outputs:
    key - sample id
    merged_bam - sorted and quality filtered bam merged by sample
    read_count - file with the total reads listed
    logfile - running log

 Pass through:

 Summary:
    Use samtools to merge bams from the same sample
    Count the number of reads in the sample

 Downstream:
    split_bam
    calc_duplication_rate

 Published:
    merged_bam - sorted and quality filtered bam merged by sample

 Notes:

*************/

if (params.max_cores < 8) {
    cores_merge = params.max_cores
} else {
    cores_merge = 8
}

sorted_bams
    .groupTuple()
    .set { bams_to_merge }

log_premerge.join(bams_to_merge).set{for_merge_bams}

save_bam = {params.output_dir + "/" + it - ~/.bam/ + "/" + it}

process merge_bams {
    cache 'lenient'
    module 'modules:modules-init:modules-gs:samtools/1.4'
    publishDir path: "${params.output_dir}/", saveAs: save_bam, pattern: "*.bam", mode: 'copy'
    penv 'serial'
    cpus cores_merge

    input:
        set key, file(logfile), file(sorted_bam) from for_merge_bams

    output:
        set key, file("*.bam"), file("merge_bams.log") into sample_bams
        set key, file("*.read_count.txt") into read_count

    """
    cat ${logfile} > merge_bams.log
    printf "** Start process 'merge_bams' at: \$(date)\n\n" >> merge_bams.log
    printf "    Process versions:
        \$(samtools --version | tr '\n' ' ')\n\n" >> merge_bams.log
    printf "    Process command:
        samtools merge ${key}.bam $sorted_bam\n\n" >> merge_bams.log


    samtools merge -@ $cores_merge ${key}.bam $sorted_bam


    printf "${key}\t\$(samtools view -c ${key}.bam)" > ${key}.read_count.txt

    printf "** End process 'merge_bams' at: \$(date)\n\n" >> merge_bams.log
    """
}


/*************

Process: split_bam

 Inputs:
    merged_bam - sorted and quality filtered bam merged by sample
    logfile - running log

 Outputs:
    merged_bam - sorted and quality filtered bam merged by sample - stops here
    split_bam - bams split by reference (chromosome)

 Pass through:
    key - sample id
    gtf_path - path to gtf info folder
    logfile - running log

 Summary:
    Use bamtools split to split bam by reference (chromosome) for speed
    Combine the small non-chromosomal references to keep file number reasonable

 Downstream:
    merge_assignment
    remove_dups_assign_genes

 Published:

 Notes:
    Potential improvement: find a way to split bams into more evenly sized chunks

*************/

sample_bams.join(gtf_info).set{assign_prepped}

process split_bam {
    cache 'lenient'
    memory '5 GB'
    module 'modules:modules-init:modules-gs:samtools/1.4:bamtools/2.2.3:bedtools/2.26.0:python/3.6.4'

    input:
        set key, file(merged_bam), file(logfile), val(gtf_path) from assign_prepped

    output:
        set key, file("split_bams/*.bam"), val(gtf_path) into split_bams mode flatten
        set key, file("remove_dups.log") into split_bam_log
        file merged_bam into output

    """
    cat ${logfile} > remove_dups.log
    printf "** Start processes 'remove duplicates, assign_genes, merge_assignment' at: \$(date)\n\n" >> remove_dups.log
    printf "    Process versions:
        \$(bedtools --version)
        \$(samtools --version | tr '\n' ' ')
        \$(bamtools --version | grep bamtools)
        \$(python --version)\n\n" >> remove_dups.log

    echo '    Process command:
        mkdir split_bams
        bamtools split -in $merged_bam -reference -stub split_bams/split

        rmdup.py --bam in_bam --output_bam out.bam
    
        samtools view -c out.bam > split_bam_umi_count.txt

        bedtools bamtobed -i out.bam -split
                | sort -k1,1 -k2,2n -k3,3n -S 5G
                > "in_bam.bed"

        bedtools map
            -a in_bam.bed
            -b exon_index
            -nonamecheck -s -f 0.95 -c 7 -o distinct -delim "|"
        | bedtools map
            -a - -b gene_index
            -nonamecheck -s -f 0.95 -c 4 -o distinct -delim "|"
        | sort -k4,4 -k2,2n -k3,3n -S 5G
        | datamash
            -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8
        | assign-reads-to-genes.py gene_index
        | awk \$3 == "exonic" || \$3 == "intronic" {{
                split(\$1, arr, "|")
                printf "%s_%s_%s\t%s\t%s\\n", arr[3], arr[4], arr[5], \$2, \$3
        }}
        | sort -k2,2 -k1,1 -S 5G > in_bam.txt

        cat logfile > merge_assignment.log
        cat split_bed > key.bed
        sort -m -k1,1 -k2,2 split_gene_assign > key_ga.txt

        datamash -g 1,2 count 2 < key_ga.txt 
        | gzip > key.gz
        ' >> remove_dups.log

    printf "    Process stats:
        remove_dups starting reads: \$(samtools view -c $merged_bam)" >> remove_dups.log


    mkdir split_bams
    bamtools split -in $merged_bam -reference -stub split_bams/split
    cd split_bams
    if [[ \$(ls | grep "_[0-9A-Za-z\\.]\\{3,\\}.bam\$") ]]; then 
        ls | grep "_[0-9A-Za-z\\.]\\{3,\\}.bam\$" | samtools merge split.REFnonstand.bam -b -
        ls | grep "_[0-9A-Za-z\\.]\\{3,\\}.bam\$" | xargs -d"\\n" rm
        mv split.REFnonstand.bam split.REF_nonstand.bam
    fi
    """

}


/*************

Process: remove_dups_assign_genes

 Inputs:
    split_bam - bams split by reference (chromosome)
    gtf_path - path to gtf info folder
    split_umi_count - file with count of umis in split bam

 Outputs:
    split_gene_assign - text file with 3 columns: sample|cell, gene, gene type (exonic, intronic)
    split_bed - deduplicated sorted bed file

 Pass through:
    key - sample id

 Summary:
    1. Remove duplicate reads using rmdup.py
    2. Use bedtools map to map the dedupped bed file to all exons with options:
        -s forced strandedness
        -f 0.95 95% of read must overlap exon
        -c 7 map the name of the gene
        -o distinct concatenate list of gene names
        -delim "|" custom delimiter
        -nonamecheck Don't error if there are different naming conventions for the chromosomes
    3. Use bedtools map to map output to gene index
        -s forced strandedness
        -f 0.95 95% of read must overlap exon
        -c 4 map the name of the cell name
        -o distinct concatenate list of gene names
        -delim "|" custom delimiter
        -nonamecheck Don't error if there are different naming conventions for the chromosomes
    4. Sort and collapse
    5. Run assign-reads-to-genes.py to deal with exon v intron

 Downstream:
    merge_assignment

 Published:

 Notes:
    Potential speed up - remove non-genic reads before sort?

*************/

process remove_dups_assign_genes {
    cache 'lenient'
    memory '4 GB'
    module 'modules:modules-init:modules-gs:bedtools/2.26.0:python/3.6.4:coreutils/8.24:samtools/1.4'

    input:
        set key, file(split_bam), val(gtf_path) from split_bams

    output:
        set key, file("*.bed"), file("*_ga.txt"), file("*_umi_count.txt") into remove_dup_part_out

    """
    rmdup.py --bam $split_bam --output_bam out.bam

    samtools view -c out.bam > ${split_bam}_umi_count.txt

    bedtools bamtobed -i out.bam -split \
            | sort -k1,1 -k2,2n -k3,3n -S 5G \
            > "${split_bam}.bed"

    bedtools map \
        -a "${split_bam}.bed" \
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
            printf "%s_%s_%s\t%s\t%s\\n", arr[3], arr[4], arr[5], \$2, \$3
    }}' \
    | sort -k1,1 -k2,2 -S 5G > "${split_bam}_ga.txt"

    """

}

/*************

Process: merge_assignment

 Inputs:
    key - sample id
    split_gene_assign - text file with 3 columns: sample|cell, gene, gene type (exonic, intronic)
    split_bed - deduplicated sorted bed file
    logfile - running log
    split_umi_count - file with count of umis in split bam
    read_count - file with the total reads listed

 Outputs:
    key - sample id
    gene_assign - text file with 3 columns: sample|cell, gene, gene type (exonic, intronic)
    cell_gene_count - gzipped text file with a count of cell, gene pairs
    logfile - running log
    dup_stats - file with duplication rate information for the sample

 Pass through:

 Summary:
    merge bed files by sample
    merge gene assignment files by sample
    make cell gene count file
    calculate duplication rate

 Downstream:
    count_umis_by_sample
    reformat_scrub

 Published:

 Notes:

*************/

remove_dup_part_out
    .groupTuple()
    .join(split_bam_log)
    .join(read_count)
    .set { for_cat_dups }

process merge_assignment {
    cache 'lenient'
    memory '1 GB'
    module 'modules:modules-init:modules-gs:coreutils/8.24'

    input:
        set key, file(split_bed), file(split_gene_assign), file(split_umi_count), file(logfile), file(read_count) from for_cat_dups

    output:
        set key, file("*.gz"), file("*_ga.txt"), file("merge_assignment.log") into merge_assignment_out
        set val(key), file("*duplication_rate_stats.txt") into duplication_rate_out
        file "*.bed"  into temp_bed

    """
    cat ${logfile} > merge_assignment.log
    cat $split_bed > "${key}.bed"
    sort -m -k1,1 -k2,2 $split_gene_assign > "${key}_ga.txt"

    datamash -g 1,2 count 2 < "${key}_ga.txt" \
    | gzip > "${key}.gz"


    umi=`cat $split_umi_count | awk '{ sum += \$1 } END { print sum }'`
    read=`cut -f2 $read_count`
    perc=\$(echo "100.0 * (1 - \$umi/\$read)" | bc -l)
    printf "%-18s   %10d    %10d    %7.1f\\n" $key \$read \$umi \$perc \
    >"${key}.duplication_rate_stats.txt"

    printf "
        remove_dups ending reads  : \$(wc -l ${key}.bed | awk '{print \$1;}')\n\n
        Read assignments:\n\$(awk '{count[\$3]++} END {for (word in count) { printf "            %-20s %10i\\n", word, count[word]}}' ${key}_ga.txt)\n\n" >> merge_assignment.log

    printf "** End processes 'remove duplicates, assign_genes, merge_assignment' at: \$(date)\n\n" >> merge_assignment.log

    """
}


/*************

Process: count_umis_by_sample

 Inputs:
    key - sample id
    gene_assign - text file with 3 columns: sample|cell, gene, gene type (exonic, intronic)
    logfile - running log

 Outputs:
    key - sample id
    logfile - running log
    umis_per_cell - count of umis per cell
    umis_per_cell_intronic - count of umis per cell only from intronic reads - stops here

 Pass through:
    cell_gene_count - gzipped text file with a count of cell, gene pairs

 Summary:
    calculate umis per sample - tabulate_per_cell_counts.py

 Downstream:
    make_matrix
    generate_qc_metrics

 Published:
    umis_per_cell - count of umis per cell
    umis_per_cell_intronic - count of umis per cell only from intronic reads

 Notes:

*************/

save_umi_per_cell = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.txt/ + "/umis_per_cell_barcode.txt"}
save_umi_per_int = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.intronic.txt/ + "/intronic_umis_per_cell_barcode.txt"}

process count_umis_by_sample {
    module 'modules:modules-init:modules-gs:python/3.6.4'
    cache 'lenient'
    memory '8 GB'
    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_int, pattern: "*intronic.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_cell, pattern: "*barcode.txt", mode: 'copy'

    input:
        set val(key), file(cell_gene_count), file(gene_assign), file(logfile) from merge_assignment_out

    output:
        set key, file(cell_gene_count), file("count_umis_by_sample.log") into ubss_out
        set key, file("*UMIs.per.cell.barcode.txt") into umis_per_cell
        file "*UMIs.per.cell.barcode.intronic.txt" into umi_per_cell_intronic

    """
    cat ${logfile} > count_umis_by_sample.log
    printf "** Start process 'count_umis_by_sample' at: \$(date)\n\n" >> count_umis_by_sample.log
    printf "    Process versions:
        \$(python --version)\n\n" >> count_umis_by_sample.log
    printf "    Process command:
        tabulate_per_cell_counts.py
            --gene_assignment_files "$gene_assign"
            --all_counts_file "${key}.UMIs.per.cell.barcode.txt"
            --intron_counts_file "${key}.UMIs.per.cell.barcode.intronic.txt"\n\n"      >> count_umis_by_sample.log


    tabulate_per_cell_counts.py \
        --gene_assignment_files "$gene_assign" \
        --all_counts_file "${key}.UMIs.per.cell.barcode.txt" \
        --intron_counts_file "${key}.UMIs.per.cell.barcode.intronic.txt"


    printf "    Process stats:
        Total cells                            : \$(wc -l ${key}.UMIs.per.cell.barcode.txt | awk '{print \$1;}')
        Total cells > 100 reads                : \$(awk '\$2>100{c++} END{print c+0}' ${key}.UMIs.per.cell.barcode.txt)
        Total cells > 1000 reads               : \$(awk '\$2>1000{c++} END{print c+0}' ${key}.UMIs.per.cell.barcode.txt)
        Total reads in cells with > 100 reads  : \$(awk '\$2>100{c=c+\$2} END{print c+0}' ${key}.UMIs.per.cell.barcode.txt)\n\n" >> count_umis_by_sample.log

    printf "** End process 'count_umis_by_sample' at: \$(date)\n\n" >> count_umis_by_sample.log
    """
}


/*************

Process: make_matrix

 Inputs:
    key - sample id
    gtf_path - path to gtf info folder
    cell_gene_count - gzipped text file with a count of cell, gene pairs
    logfile - running log

 Outputs:
    key - sample id
    logfile - running log
    umi_matrix - MatrixMarket format matrix of cells by genes
    cell_anno - Cell annotations for umi_matrix
    gene_anno - Gene annotations for umi_matrix
    gtf_path - path to gtf info folder

 Pass through:

 Summary:
    Generate a matrix of cells by genes - make_matrix.py

 Downstream:
    make_cds

 Published:
    umi_matrix - MatrixMarket format matrix of cell by umi
    cell_anno - Cell annotations for umi_matrix
    gene_anno - Gene annotations for umi_matrix

 Notes:

*************/

ubss_out.join(gtf_info2).set{make_matrix_prepped}
save_umi = {params.output_dir + "/" + it - ~/.umi_counts.mtx/ + "/umi_counts.mtx"}
save_cell_anno = {params.output_dir + "/" + it - ~/.cell_annotations.txt/ + "/cell_annotations.txt"}
save_gene_anno = {params.output_dir + "/" + it - ~/.gene_annotations.txt/ + "/gene_annotations.txt"}

process make_matrix {
    cache 'lenient'
    memory '15 GB'
    publishDir path: "${params.output_dir}/", saveAs: save_umi, pattern: "*umi_counts.mtx", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_anno, pattern: "*cell_annotations.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_gene_anno, pattern: "*gene_annotations.txt", mode: 'copy'

    input:
        set key, file(cell_gene_count), file(logfile), val(gtf_path) from make_matrix_prepped

    output:
        set key, file("*cell_annotations.txt"), file("*umi_counts.mtx"), file("*gene_annotations.txt"), val(gtf_path), file("make_matrix.log") into mat_output

    """
    cat ${logfile} > make_matrix.log
    printf "** Start process 'make_matrix' at: \$(date)\n\n" >> make_matrix.log

    echo '    Process command:
        make_matrix.py <(zcat $cell_gene_count) 
            --gene_annotation "${gtf_path}/latest.gene.annotations" 
            --key "$key"
        cat ${gtf_path}/latest.gene.annotations > "${key}.gene_annotations.txt"  ' >> make_matrix.log


    make_matrix.py <(zcat $cell_gene_count) --gene_annotation "${gtf_path}/latest.gene.annotations" --key "$key"
    cat "${gtf_path}/latest.gene.annotations" > "${key}.gene_annotations.txt"


    printf "\n** End process 'make_matrix' at: \$(date)\n\n" >> make_matrix.log
    """

}


/*************

Process: make_cds

 Inputs:
    key - sample id
    umi_matrix - MatrixMarket format matrix of cells by genes
    cell_anno - Cell annotations for umi_matrix
    gene_anno - Gene annotations for umi_matrix
    gtf_path - path to gtf info folder
    logfile - running log
    params.umi_cutoff

 Outputs:
    key - sample id
    scrub_matrix - matrix of counts output in proper format for scrublet
    cds_object - cds object in RDS format
    cell_qc - csv of cell quality control information
    logfile - running log

 Pass through:

 Summary:
    Generate a monocle3 cds object - make_cds.R

 Downstream:
    run_scrublet
    calc_cell_totals

 Published:

 Notes:

*************/

process make_cds {
    cache 'lenient'
    module 'modules:modules-init:modules-gs:gcc/8.1.0:R/3.6.1'
    memory '15 GB'

    input:
        set key, file(cell_data), file(umi_matrix), file(gene_data), val(gtf_path), file(logfile) from mat_output

    output:
        set key, file("*for_scrub.mtx"), file("*.RDS"), file("*cell_qc.csv"), file("make_cds.log") into cds_out
        file("*cell_qc.csv") into cell_qcs

    """
    cat ${logfile} > make_cds.log
    printf "** Start process 'make_cds' at: \$(date)\n\n" >> make_cds.log
    printf "    Process versions:
        \$(R --version | grep 'R version')
            monocle3 version \$(Rscript -e 'packageVersion("monocle3")')\n\n" >> make_cds.log
    echo '    Process command:
        make_cds.R
            "$umi_matrix"
            "$cell_data"
            "$gene_data"
            "${gtf_path}/latest.genes.bed"
            "$key"
            "$params.umi_cutoff"\n' >> make_cds.log


    make_cds.R \
        "$umi_matrix"\
        "$cell_data"\
        "$gene_data"\
        "${gtf_path}/latest.genes.bed"\
        "$key"\
        "$params.umi_cutoff"


    printf "** End process 'make_cds' at: \$(date)\n\n" >> make_cds.log
    """
}


process apply_garnett {
    cache 'lenient'
    module 'modules:modules-init:modules-gs:gcc/8.1.0:R/3.6.1'
    memory '15 GB'

    input:
        set key, file(scrub_matrix), file(cds_object), file(cell_qc), file(logfile) from cds_out

    output:
        set key, file(scrub_matrix), file("new_cds/*.RDS"), file(cell_qc), file("apply_garnett.log") into for_scrub

"""
    cat ${logfile} > apply_garnett.log
    printf "** Start process 'apply_garnett' at: \$(date)\n\n" >> apply_garnett.log
    mkdir new_cds
    echo "No Garnett classifier provided for this sample" > garnett_error.txt
    if [ $params.garnett_file == 'false' ]
then
    cp $cds_object new_cds/
else
    apply_garnett.R $cds_object $params.garnett_file $key
fi

cat garnett_error.txt >> apply_garnett.log
printf "\n** End process 'apply_garnett' at: \$(date)\n\n" >> apply_garnett.log

"""


}


/*************

Process: run_scrublet

 Inputs:
    key - sample id
    scrub_matrix - matrix of counts output in proper format for scrublet
    logfile - running log

 Outputs:
    key - sample id
    scrub_csv - scrublet results csv
    scrublet_png - png histogram of scrublet scores
    logfile - running log

 Pass through:
    cds_object - cds object in RDS format
    cell_qc - csv of cell quality control information

 Summary:
    Run scrublet to generate doublet scores - run_scrublet.py

 Downstream:
    calc_duplication_rate
    generate_dashboard
    finish_log

 Published:
    scrublet_png - png histogram of scrublet scores

 Notes:

*************/

save_hist = {params.output_dir + "/" + it - ~/_scrublet_hist.png/ + "/" + it}

process run_scrublet {
    publishDir path: "${params.output_dir}/", saveAs: save_hist, pattern: "*png", mode: 'copy'
    cache 'lenient'
    module 'modules:java/latest:modules-init:modules-gs:python/3.6.4'
    memory '25 GB'

    input:
        set key, file(scrub_matrix), file(cds_object), file(cell_qc), file(logfile) from for_scrub

    output:
        set key, file("*scrublet_out.csv"), file(cds_object), file(cell_qc) into scrublet_out
        file ("*.png") into scrub_pngs
        set key, file("run_scrublet.log") into pipe_log

    """
    cat ${logfile} > run_scrublet.log
    printf "** Start process 'run_scrublet' at: \$(date)\n\n" >> run_scrublet.log
    printf "    Process versions:
        \$(python --version)
            \$(pip freeze | grep scrublet | tr '==' ' ')\n\n" >> run_scrublet.log
    echo '    Process command:
        run_scrublet.py --key $key --mat $scrub_matrix\n'  >> run_scrublet.log


    run_scrublet.py --key $key --mat $scrub_matrix


    printf "** End process 'run_scrublet' at: \$(date)\n\n" >> run_scrublet.log

    printf "** Start processes to generate qc metrics and dashboard at: \$(date)\n\n" >> run_scrublet.log
    """

}


/*************

Process: reformat_scrub

 Inputs:
    key - sample id
    cds_object - cds object in RDS format
    cell_qc - csv of cell quality control information
    scrub_csv - scrublet results csv
    dup_stats - file with duplication rate information for the sample

 Outputs:
    key - sample id
    cds_object - cds object in RDS format
    sample_stats - csv with sample-wise statistics
    cell_qc - csv of cell quality control information
    collision - file containing collision rate if barnyard sample

 Pass through:

 Summary:
    Add scrublet info to cell_qc and cds object
    Calculate collision rate for barnyard
    Calculate sample statistics

 Downstream:
    generate_qc_metrics
    zip_up_sample_stats
    collapse_collision

 Published:
    cds_object - cds object in RDS format
    sample_stats - csv with sample-wise statistics
    cell_qc - csv of cell quality control information

 Notes:

*************/

scrublet_out.join(duplication_rate_out).set{reformat_scrub_in}

save_cds = {params.output_dir + "/" + it - ~/_cds.RDS/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_cell_qc = {params.output_dir + "/" + it - ~/_cell_qc.csv/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_samp_stats = {params.output_dir + "/" + it - ~/_sample_stats.csv/ + "/" + it}

process reformat_scrub {
    cache 'lenient'
    module 'modules:modules-init:modules-gs:python/3.6.4:gcc/8.1.0:R/3.6.1'
    memory '10 GB'
    publishDir path: "${params.output_dir}/", saveAs: save_cds, pattern: "temp_fold/*cds.RDS", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_qc, pattern: "temp_fold/*cell_qc.csv", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_samp_stats, pattern: "*sample_stats.csv", mode: 'copy'

    input:
        set key, file(scrub_csv), file(cds_object), file(cell_qc), file(dup_stats) from reformat_scrub_in

    output:
        set key, file("temp_fold/*.RDS"), file("temp_fold/*.csv") into rscrub_out
        file("*sample_stats.csv") into sample_stats
        file("*collision.txt") into collision


    """
    #!/usr/bin/env Rscript

    library(monocle3)

    dir.create("temp_fold")
    cds <- readRDS("$cds_object")
    cell_qc <- read.csv("$cell_qc")

    if(nrow(pData(cds)) > 0) {
        scrublet_out <- read.csv("$scrub_csv", header=F)
        pData(cds)\$scrublet_score <- scrublet_out\$V1
        pData(cds)\$scrublet_call <- ifelse(scrublet_out\$V2 == 1, "Doublet", "Singlet")
        cell_qc\$scrublet_score <- scrublet_out\$V1
        cell_qc\$scrublet_call <- ifelse(scrublet_out\$V2 == 1, "Doublet", "Singlet")
    }

    write.csv(cell_qc, quote=FALSE, file="temp_fold/$cell_qc")

    dup_stats <- read.table(paste0("$key", ".duplication_rate_stats.txt"))

    df <- data.frame(sample="$key", n.reads = dup_stats\$V2, n.umi = dup_stats\$V3,
                     duplication_rate = dup_stats\$V4,
                     doublet_count = sum(cell_qc\$scrublet_call == "Doublet", na.rm=TRUE),
                     doublet_perc = paste0(round(sum(cell_qc\$scrublet_call == "Doublet",
                                                     na.rm=TRUE)/nrow(cell_qc) * 100, 1), "%"),
                     doublet_NAs=sum(is.na(cell_qc\$scrublet_call)))

    write.csv(df, file=paste0("$key", "_sample_stats.csv"), quote=FALSE, row.names=FALSE)
    saveRDS(cds, file="temp_fold/$cds_object")

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
        fileConn<-file("${key}_no_collision.txt")
        writeLines(paste0("$key", "\t", "NA"), fileConn)
        close(fileConn)
    }
    """

}


/*************

Process: generate_qc_metrics

 Inputs:
    key - sample id
    params.umi_cutoff
    cds_object - cds object in RDS format
    cell_qc - csv of cell quality control information
    umis_per_cell - count of umis per cell

 Outputs:
    cutoff - not currently used
    umap_png - png sample UMAP
    knee_png - png sample knee plot
    qc_png - png of cell qc stats

 Pass through:

 Summary:
    Generate a bunch of qc metrics and plots - generate_qc.R

 Downstream:
    generate_dashboard

 Published:
    umap_png - png sample UMAP
    knee_png - png sample knee plot
    qc_png - png of cell qc stats

 Notes:
    Need to test umi cutoff here and in cds function
    Need to either remove or use output cutoff

*************/

for_gen_qc = rscrub_out.join(umis_per_cell)
save_knee = {params.output_dir + "/" + it - ~/_knee_plot.png/ + "/" + it}
save_umap = {params.output_dir + "/" + it - ~/_UMAP.png/ + "/" + it}
save_cellqc = {params.output_dir + "/" + it - ~/_cell_qc.png/ + "/" + it}
save_garnett = {params.output_dir + "/" + it.split("_")[0] + "/" + it}

process generate_qc_metrics {
    module 'modules:modules-init:modules-gs:gcc/8.1.0:R/3.6.1'
    memory '10 GB'
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_umap, pattern: "*UMAP.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_knee, pattern: "*knee_plot.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cellqc, pattern: "*cell_qc.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_garnett, pattern: "*Garnett.png", mode: 'copy'

    input:
        set key, file(cds_object), file(cell_qc), file(umis_per_cell) from for_gen_qc

    output:
        file("*.png") into qc_plots
        file("*.txt") into cutoff

    """
    generate_qc.R\
        $cds_object $umis_per_cell $key \
        --specify_cutoff $params.umi_cutoff\

    """
}

/*************

Process: zip_up_sample_stats

 Inputs:
    sample_stats - csv with sample-wise statistics - collected

 Outputs:
    all_sample_stats - concatenated table of sample stats from all samples

 Pass through:

 Summary:
    Concatenate duplication information into one table

 Downstream:
    generate_dashboard

 Published:
    all_sample_stats - concatenated table of sample stats from all samples

 Notes:

*************/

process zip_up_sample_stats {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", pattern: "all_sample_stats.csv", mode: 'copy'

    input:
        file files from sample_stats.collect()

    output:
        file "*ll_sample_stats.csv" into all_sample_stats

    """

     sed -s 1d $files > all_sample_stats.csv

    """
}


/*************

Process: calc_cell_totals

 Inputs:
    cell_qc - csv of cell quality control information - collected

 Outputs:
    cell_counts - table cell totals above set UMI thresholds for all samples

 Pass through:

 Summary:
    Count cell totals above set UMI thresholds for all samples

 Downstream:
    generate_dashboard

 Published:

 Notes:

*************/

process calc_cell_totals {
    memory '1 GB'
    cache 'lenient'

    input:
        file cell_qc from cell_qcs.collect()

    output:
        file "*.txt" into cell_counts

    """

    for f in $cell_qc
    do
      awk 'BEGIN {FS=","}; \$2>100{c++} END{print FILENAME, "100", c-1}' \$f >> cell_counts.txt
      awk 'BEGIN {FS=","}; \$2>500{c++} END{print FILENAME, "500", c-1}' \$f >> cell_counts.txt
      awk 'BEGIN {FS=","}; \$2>1000{c++} END{print FILENAME, "1000", c-1}' \$f >> cell_counts.txt
    done

    """

}


/*************

Process: collapse_collision

 Inputs:
    collision - file containing collision rate if barnyard sample - collected

 Outputs:
    all_collision - concatenate collision values for all samples (all NA except Barnyard)

 Pass through:

 Summary:
    Concatenate collision values for all samples

 Downstream:
    generate_dashboard

 Published:

 Notes:

*************/

process collapse_collision {
    memory '1 GB'
    cache 'lenient'

    input:
        file col_file from collision.collect()

    output:
        file "*.txt" into all_collision

    """

    cat $col_file > all_collision.txt

    """
}


/*************

Process: generate_dashboard

 Inputs:
    cell_counts - table cell totals above set UMI thresholds for all samples
    all_collision - concatenate collision values for all samples (all NA except Barnyard)
    all_sample_stats - concatenated table of sample stats from all samples
    params.output_dir
    umap_png - png sample UMAP - combined as qc_plots
    knee_png - png sample knee plot - combined as qc_plots
    qc_png - png of cell qc stats - combined as qc_plots
    scrublet_png - png histogram of scrublet scores
    params.garnett_file

 Outputs:
    exp_dash - experimental dashboard

 Pass through:

 Summary:
    Collect plots and generate data file for experimental dashboard
    Assemble dashboard

 Downstream:
    generate_summary_log

 Published:
    exp_dash - experimental dashboard

 Notes:

*************/

process generate_dashboard {
    module 'modules:modules-init:modules-gs:gcc/8.1.0:R/3.6.1'
    cache 'lenient'
    memory '1 GB'
    publishDir path: "${params.output_dir}/", pattern: "exp_dash", mode: 'copy'

    input:
        file all_sample_stats
        file cell_counts
        file all_collision
        file plots from qc_plots.collect()
        file scrublet_png from scrub_pngs.collect()

    output:
        file exp_dash into exp_dash_out

    """
    generate_dash_data.R $all_sample_stats $params.output_dir $cell_counts $all_collision $params.garnett_file

    mkdir exp_dash
    cp -R $baseDir/bin/skeleton_dash/* exp_dash/
    mv *.png exp_dash/img/

    mv *.js exp_dash/js/

    """
}


/*************

Process: finish_log

 Inputs:
    key - sample id
    logfile - running log
    exp_dash - experimental dashboard - input so runs last

 Outputs:
    full_log - Final full pipeline log
    summary_log - Summary log
    log_data - Logging info for dashboards

 Pass through:

 Summary:
    Add parameter info to front of pipeline - allows restart when changing minor parameters
    Generate summary log
    Generate log info for dashboards

 Downstream:
    zip_up_log_data

 Published:
    full_log - Final full pipeline log
    summary_log - Summary log
    log_data - Logging info for dashboards

 Notes:

*************/

save_logs = {params.output_dir + "/" + it - ~/_read_metrics.log/ - ~/_full.log/ + "/" + it}
save_txt_for_wrap = {params.output_dir + "/" + it - ~/_log_data.txt/ + "/" + it}

process finish_log {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_logs, pattern: "*.log", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_txt_for_wrap, pattern: "*.txt", mode: 'copy'

    input:
        set key, file(logfile) from pipe_log
        file exp_dash from exp_dash_out

    output:
        file("*_full.log") into full_log
        file("*_read_metrics.log") into summary_log
        file("*log_data.txt") into log_txt_for_wrap

    """
    head -n 2 ${logfile} > ${key}_full.log
    printf "Nextflow version: $nextflow.version\n" >> ${key}_full.log
    printf "Pipeline version: $workflow.manifest.version\n" >> ${key}_full.log
    printf "Git Repository, Version, Commit ID, Session ID: $workflow.repository, $workflow.revision, $workflow.commitId, $workflow.sessionId\n\n" >> ${key}_full.log
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
    printf "    params.rt_barcode_file:       $params.rt_barcode_file\n" >> ${key}_full.log
    printf "    params.hash_list:             $params.hash_list\n" >> ${key}_full.log
    printf "    params.max_wells_per_sample:  $params.max_wells_per_sample\n\n" >> ${key}_full.log

    tail -n +2 ${logfile} >> ${key}_full.log
    printf "\n** End processes generate qc metrics and dashboard at: \$(date)\n\n" >> ${key}_full.log
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
    done | awk '{sum += \$1} END {printf "%1.0f", sum}'`

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
            \\"${key}\\": {
            \\"sample\\": \\"${key}\\",
            \\"alignment_start\\" : \\"\$align_start\\",
            \\"alignment_mapped\\" : \\"\$align_mapped\\",
            \\"align_multimapped\\" : \\"\$align_multimapped\\",
            \\"align_too_short\\" : \\"\$align_too_short\\",
            \\"sf_start\\" : \\"\$sf_start\\",
            \\"sf_end\\" : \\"\$sf_end\\",
            \\"dup_start\\" : \\"\$dup_start\\",
            \\"dup_end\\" : \\"\$dup_end\\",
            \\"assigned_exonic\\" : \\"\$assigned_exonic\\",
            \\"assigned_intronic\\" : \\"\$assigned_intronic\\",
            \\"reads_in_cells\\" : \\"\$reads_in_cells\\" }
      " > ${key}_log_data.txt


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

/*************

Process: zip_up_log_data

 Inputs:
    summary_log - collected summary log files
    full_log - collected full log files

 Outputs:
    log_js - Logging info in a js format for dashboards

 Pass through:

 Summary:
    Generate log data js file for dashboard

 Downstream:
    End

 Published:
    log_data.js - Logging info for dashboards

 Notes:

*************/

process zip_up_log_data {
    cache 'lenient'
    publishDir path: "${params.output_dir}/exp_dash/js/", pattern: "*.js", mode: 'copy'

    input:
        file summary_log from summary_log.collect()
        file full_log from full_log.collect()

    output:
        file "*.js" into log_js

    """

    echo 'const log_data = {' > log_data.js
    for file in $summary_log
    do
        samp_name=\$(basename \$file | sed 's/_read_metrics.log//')
        echo "\\"\$samp_name\\" :  \\`" >> log_data.js
        cat \$file >> log_data.js
        echo "\\`," >> log_data.js
    done
    sed -i '\$ s/,\$//' log_data.js
    echo '}' >> log_data.js

    echo 'const full_log_data = {' >> log_data.js
    for file in $full_log 
    do  
        samp_name=\$(basename \$file | sed 's/_full.log//')
        echo "\\"\$samp_name\\" :  \\`" >> log_data.js
        cat \$file >> log_data.js
        echo "\\`," >> log_data.js
    done
    sed -i '\$ s/,\$//' log_data.js
    echo '}' >> log_data.js

    """
}

workflow.onComplete {
	println ( workflow.success ? "Done! Saving output" : "Oops .. something went wrong" )
}


/*************
Groovy functions
*************/

def checkNextflowVersion( Integer minMajorVersion, Integer minMinorVersion )
{
  def sVersion = nextflow.version.toString()
  def aVersion = sVersion.split( /[.]/ )
  def majorVersion = aVersion[0].toInteger()
  def minorVersion = aVersion[1].toInteger()
  if( majorVersion < minMajorVersion || minorVersion < minMinorVersion )
  {
    def serr = "This pipeline requires Nextflow version at least %s.%s: you have version %s."
    println()
    println( '****  ' + String.format( serr, minMajorVersion, minMinorVersion, sVersion ) + '  ****' )
    println()
    System.exit( -1 )
    /*
    ** An exception produces an exceptionally verbose block of confusing text. I leave
    ** the command here in case the println() output is obscured by fancy Nextflow tables.
    **
    ** throw new Exception( String.format( serr, minMajorVersion, minMinorVersion, sVersion ) )
    */
  }
  return( 0 )
}

