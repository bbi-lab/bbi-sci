/*
** This pipeline is written for Nextflow DSL 1.
*/
nextflow.enable.dsl = 1

/*
** Check that Nextflow version meets minimum version requirements.
*/
def minMajorVersion = 20
def minMinorVersion = 07
checkNextflowVersion( minMajorVersion, minMinorVersion )


/*
** Check OS version.
** Notes:
**   o  works only for Linux systems
**   o  used to distinguish between CentOS 6 and CentOS 7
*/
( osName, osDistribution, osRelease ) = getOSInfo()


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
params.skip_doublet_detect = false
params.run_emptyDrops = true
params.hash_umi_cutoff = 5
params.hash_ratio = false
params.hash_dup = false // Default is false. Other options are "p5" or "pcr_plate". params.hash_list must also be true.
params.hash_rt_split = false // Default is false. If true, will process on hash samples demuxed by RT barcodes.


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
    log.info '    params.skip_doublet_detect = false         Whether to skip doublet detection, i.e. scrublet - useful for very large datasets.'
    log.info '    params.hash_umi_cutoff = 5                 The hash umi cutoff to determine hash vs background in top to second best hash oligo. Default is 5'
    log.info '    params.hash_ratio = false                  The min hash umi ratio for top to second best. Default is false and not filtered'
    log.info '    params.hash_dup = false                    Whether to run hash PCR duplication rate per plate and what indicates a plate. Options are "p5" and "pcr_plate". params.hash_list also needs to be set to true. Default is false.'
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

    output:
        file "*.csv" into good_sample_sheet
        file '*.log' into log_check_sample
        file 'start.txt' into log_piece1

    """
    # bash watch for errors
    set -ueo pipefail

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
    # bash watch for errors
    set -ueo pipefail

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
   o  the 'spec' variable uses the awk split() function to remove
      '_fq_part', which is added when very large fastq files are
      split. The gsub() function removes unacceptable characters
      from the sample names.

*************/

process gather_info {
    cache 'lenient'

    input:
        file good_sample_sheet
        set val(key), val(name), file(trimmed_fastq), file(logfile), file(log_piece2) from trimmed_fastqs

    output:
        set val(key), val(name), env(star_path), env(star_mem), file(trimmed_fastq), file(logfile), file(log_piece2) into align_prepped
        set val(key), env(gtf_path) into gtf_info
        set val(key), env(gtf_path) into gtf_info2
        file good_sample_sheet into sample_sheet_for_rt_split

    """
    # bash watch for errors
    set -ueo pipefail

    spec=`sed 's/ *\$//g' good_sample_sheet.csv | awk 'BEGIN {FS=",";OFS=","}{split(\$2,a,"_fq_part");gsub("[_ /-]", ".", a[1]);print(\$1, a[1], \$3)}' | awk 'BEGIN {FS=","}; \$2=="$key" {print \$3}' | uniq`
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
    hash_reads_per_cell.txt
    hash_umis_per_cell.txt
    hash_assigned_table.txt
    hash_dup_per_cell.txt

 Pass through:

 Summary:
    Collect and process hash barcodes - process_hashes.py

 Downstream:

 Published:
    hash_mtx - MatrixMarket matrix of hash information
    hash_cell - Cell info for matrix of hash information
    hash_hash - Hash info for matrix of hash information
    hash_reads_per_cell.txt
    hash_umis_per_cell.txt
    hash_assigned_table.txt
    hash_dup_per_cell.txt

 Notes:
    Only when params.hash = true

*************/

//  Group fastqs for finding hash barcodes
fastqs_out
    .groupTuple()
    .set { fastqs_for_hash }

save_hash_cell = {params.output_dir + "/" + it - ~/.hashumis_cells.txt/ + "/" + it}
save_hash_hash = {params.output_dir + "/" + it - ~/.hashumis_hashes.txt/ + "/" + it}
save_hash_mtx = {params.output_dir + "/" + it - ~/.hashumis.mtx/ + "/" + it}

save_hash_reads = {params.output_dir + "/" + it - ~/_hash_reads_per_cell.txt/ + "/" + it}
save_hash_umis = {params.output_dir + "/" + it - ~/_hash_umis_per_cell.txt/ + "/" + it}
save_hash_table = {params.output_dir + "/" + it - ~/_hash_assigned_table.txt/ + "/" + it}
save_hash_dup = {params.output_dir + "/" + it - ~/_hash_dup_per_cell.txt/ + "/" + it}

// save_sorted_hash = {params.output_dir + "/" + it - ~/_sorted_hash_combined.gz/ + "/" + it}

/*
** We no longer need or make the _sorted_hash_combined file.
*/
process process_hashes {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_cell, pattern: "*hashumis_cells.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_hash, pattern: "*hashumis_hashes.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_mtx, pattern: "*.mtx", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_reads, pattern: "*hash_reads_per_cell.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_umis, pattern: "*hash_umis_per_cell.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_table, pattern: "*hash_assigned_table.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_dup, pattern: "*hash_dup_per_cell.txt", mode: 'copy'
//    publishDir path: "${params.output_dir}/", saveAs: save_sorted_hash, pattern: "*sorted_hash_combined.gz", mode: 'copy'

    input:
        set key, file(input_fastq) from fastqs_for_hash

    output:
        file("*hash.log") into hash_logs
        set key, file("*mtx"), file("*hashumis_cells.txt"), file("*hashumis_hashes.txt") into hash_mats
        set key, file("*hash_reads_per_cell.txt"), file("*hash_umis_per_cell.txt"), file("*hash_assigned_table.txt") into hash_results
        set key, file("*hash_dup_per_cell.txt") into for_hash_calc
//        set key, file("*_sorted_hash_combined.gz") into sorted_hash_combined

    when:
        params.hash_list != false

    script:
        def sample_name = key.split(/\.P[0-9]\.[A-H][0-9]{2}/)[0]

    """
    # bash watch for errors
    set -ueo pipefail

    input_key="${key}"

    if [ ${params.hash_rt_split} == true ]; then
        input_key="${sample_name}"
    fi

   echo "Input key: \$input_key"

    awk 'NF < 4 {print "ERROR: Hash sample sheet contains rows with fewer than 4 columns"; exit 1}' $params.hash_list

    process_hashes --hash_sheet <(awk -v search="\$input_key" '\$1 ~ search {print \$2 "\t" \$3}' $params.hash_list) \
        --fastq $input_fastq --key $key --sample_name $key


#    process_hashes --hash_sheet $params.hash_list \
#        --fastq $input_fastq --key $key --sample_name $key

#    LL_ALL=C sort ${key}_hash_combined -S 50G -T /tmp/ -k2,2 -k4,4 -k3,3 --parallel=8 > ${key}_sorted_hash_combined
#    pigz -p 8 "${key}_sorted_hash_combined"

    """
}



hash_mats.into{hash_mats_for_cat; hash_mats_for_assign_hash}

concat_hash_in = hash_mats_for_cat
    .map { sample, mtx, cells, hashes ->
        def sample_name = sample.split(/\.P[0-9]\.[A-H][0-9]{2}/)[0] // "GENE3"
        tuple(sample_name, mtx, cells, hashes)
    }
    .groupTuple()

// Concatenate all rt split hash files (matrices, )


process concat_hash_matrices {
    cache 'lenient'

    input:
        tuple val(sample_name), path(mtx), path(cells), path(hashes) from concat_hash_in

    output:
        set val(sample_name), file("temp_fold/*hashumis.mtx"), file("temp_fold/*hashumis_hashes.txt"), file("temp_fold/*hashumis_cells.txt") into concat_hash_out

    publishDir (
        path: "${params.output_dir}",
        mode: 'copy',
        saveAs: { fileName -> 
            // Remove the temporary folder prefix and publish into a folder named after the sample
            def baseName = fileName.replaceFirst('^temp_fold/', '') 
            def sampleDir = baseName.replaceFirst('_hashumis.*', '') 
            return "${sample_name}/${baseName}"
        }
    )

    when:
        params.hash_rt_split != false 

    script:

    """
    # bash watch for errors
    set -ueo pipefail
    
    mkdir temp_fold

    cat_sparse_matrix.py \
        -i $mtx \
        -o temp_fold/$sample_name \
        -c .hashumis_cells.txt \
        -m .hashumis.mtx \
        -f .hashumis_hashes.txt

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
cores_align = params.max_cores < 8 ? params.max_cores : 8

process align_reads {
    cache 'lenient'
    memory { star_mem.toInteger()/cores_align + " GB" }
    cpus cores_align

    input:
        set val(key), val(name), val(star_path), val(star_mem), file(trimmed_fastq), file(logfile), file(log_piece2) from align_prepped

    output:
        file "align_out" into align_output
        set val(key), val(name), file("align_out/*Aligned.out.bam"), file('align.log'), file(log_piece2), file("*align.txt") into aligned_bams

    """
    # bash watch for errors
    set -ueo pipefail

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
cores_sf = params.max_cores < 10 ? params.max_cores : 10

process sort_and_filter {
    cache 'lenient'
    cpus cores_sf

    input:
        set val(key), val(name), file(aligned_bam), file(logfile), file(log_piece2), file(log_piece3) from aligned_bams

    output:
        set val(key), file("*.bam") into sorted_bams
        set val(key), file(log_piece2), file(log_piece3), file("*_sf.txt") into log_pieces

    """
    # bash watch for errors
    set -ueo pipefail

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

    input:
        file log_piece1
        set val(key), file(log_piece2), file(log_piece3), file(log_piece4) from logs_to_combine

    output:
        set val(key), file("*_pre.log") into log_premerge

    """
    # bash watch for errors
    set -ueo pipefail

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

cores_merge = params.max_cores < 8 ? params.max_cores : 8

sorted_bams
    .groupTuple()
    .set { bams_to_merge }

log_premerge.join(bams_to_merge).set{for_merge_bams}

save_bam = {params.output_dir + "/" + it - ~/.bam/ + "/" + it}

process merge_bams {
    cache 'lenient'
    cpus cores_merge
    publishDir path: "${params.output_dir}/", saveAs: save_bam, pattern: "*.bam", mode: 'copy'

    input:
        set key, file(logfile), file(sorted_bam) from for_merge_bams

    output:
        set key, file("*.bam"), file("merge_bams.log") into sample_bams
        set key, file("*.bam") into rt_bams
        set key, file("*.read_count.txt") into read_count
        set key, file("*.reads_per_cell.txt") into reads_per_cell

    """
    # bash watch for errors
    set -ueo pipefail

    cat ${logfile} > merge_bams.log
    printf "** Start process 'merge_bams' at: \$(date)\n\n" >> merge_bams.log
    printf "    Process versions:
        \$(samtools --version | tr '\n' ' ')\n\n" >> merge_bams.log
    printf "    Process command:
        samtools merge ${key}.bam $sorted_bam\n\n" >> merge_bams.log


    samtools merge -@ $cores_merge ${key}.bam $sorted_bam


    printf "${key}\t\$(samtools view -c ${key}.bam)" > ${key}.read_count.txt

    # Get reads per cell
    samtools view ${key}.bam | awk '{
    n = split(\$1, a, "|");
    cell = a[3] "_" a[4] "_" a[5];
    print cell
    }' | sort | uniq -c | awk '{print \$2"\t"\$1}' > ${key}.reads_per_cell.txt


    printf "** End process 'merge_bams' at: \$(date)\n\n" >> merge_bams.log
    """
}

rt_bams_in = rt_bams
    .map { sample, bam ->
        def sample_name = sample.split(/\.P[0-9]\.[A-H][0-9]{2}/)[0] 
        tuple(sample_name, bam)
    }
    .groupTuple()

// save_rt_bam = {params.output_dir + "/" + it - ~/.bam/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
// save_rt_bai = {params.output_dir + "/" + it - ~/.bai/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
// save_bam_count = {params.output_dir + "/" + it - ~/count.txt/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}

process merge_rt_bams {
    cache 'lenient'
    // publishDir path: "${params.output_dir}/", saveAs: save_rt_bam, pattern: "temp_fold/*.bam", mode: 'copy' 
    // publishDir path: "${params.output_dir}/", saveAs: save_rt_bam, pattern: "temp_fold/*.bai", mode: 'copy' 
    // publishDir path: "${params.output_dir}/", saveAs: save_bam_count, pattern: "temp_fold/*count.txt", mode: 'copy'   
    publishDir (
        path: "${params.output_dir}",
        mode: 'copy',
        saveAs: { fileName -> 
            // Remove the temporary folder prefix and publish into a folder named after the sample
            def baseName = fileName.replaceFirst('^temp_fold/', '') 
            return "${sample_name}/${baseName}"
        }
    )

    input:
        tuple val (sample_name), path(bam) from rt_bams_in

    output: 
        set val(sample_name), file("temp_fold/*") into rt_bams_out

    when: 
        params.hash_rt_split != false && params.hash_list != false

    """        
    # bash watch for errors
    set -ueo pipefail

    module load sambamba/0.6.8
    
    mkdir -p temp_fold

    bam_list=(${bam})
    if [ \${#bam_list[@]} -gt 1 ]
    then
        sambamba merge -t 8 "temp_fold/${sample_name}.bam" ${bam}
        sum=0
        for f in ${bam}
        do  
            # echo "bam: \${bam}"
            count=\$(sambamba view \${f} | wc -l)
            sum=\$((\$sum+\$count))
            # echo "count: \${count}"
            # echo "sum: \${sum}"
        done
        
        echo "${sample_name}_total: \${sum}" > "temp_fold/merged_bam_count.txt"

        # echo "${sample_name}_total: \${sum}"

        og=\$(sambamba view "temp_fold/${sample_name}.bam" | wc -l)

        echo "${sample_name}_merged: \${og}" >> "temp_fold/merged_bam_count.txt" 
        # echo "${sample_name}_merged: \${og}" 

    else
        cp ${bam} temp_fold
        echo "No additional bams to merge" > "temp_fold/merged_bam_count.txt"
        
    fi

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

    input:
        set key, file(merged_bam), file(logfile), val(gtf_path) from assign_prepped

    output:
        set key, file("split_bams/*.bam"), val(gtf_path) into split_bams mode flatten
        set key, file("remove_dups.log") into split_bam_log
        file merged_bam into output

    """
    # bash watch for errors
    set -ueo pipefail

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

        samtools index -c in_bam
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
    if [[ \$(ls *.bam | egrep -v '_(chr)?([0-9]?[0-9]?[0-9]|M|Mt|MT|MtDNA|mitochondrion_genome|X|Y)\\.bam\$') ]]; then
        ls *.bam | egrep -v '_(chr)?([0-9]?[0-9]?[0-9]|M|Mt|MT|MtDNA|mitochondrion_genome|X|Y)\\.bam\$' | samtools merge split.REFnonstand.bam.xxx -b -
        ls *.bam | egrep -v '_(chr)?([0-9]?[0-9]?[0-9]|M|Mt|MT|MtDNA|mitochondrion_genome|X|Y)\\.bam\$' | xargs -d"\\n" rm
        mv split.REFnonstand.bam.xxx split.REF_nonstand.bam
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

    input:
        set key, file(split_bam), val(gtf_path) from split_bams

    output:
        set key, file("*.bed"), file("*_ga.txt"), file("*_umi_count.txt") into remove_dup_part_out

    """
    # bash watch for errors
    set -ueo pipefail

    samtools index -c $split_bam
    rmdup.py --bam $split_bam --output_bam out.bam

    samtools view -c out.bam > ${split_bam}_umi_count.txt

# bedtools complains now about the sorting. I imagine that this
# may be a result of something like a LOCALE variable. I use
# bedtools to sort instead.
#    bedtools bamtobed -i out.bam -split \
#            | sort -k1,1 -k2,2n -k3,3n -S 5G \
#            > "${split_bam}.bed"

    bedtools bamtobed -i out.bam -split > "${split_bam}.bed.unsorted"
    bedtools sort -i "${split_bam}.bed.unsorted" > "${split_bam}.bed"
    rm ${split_bam}.bed.unsorted

    bedtools sort -i "${gtf_path}/latest.exons.bed" > latest.exons.bed.sort
    bedtools sort -i "${gtf_path}/latest.genes.bed" > latest.genes.bed.sort

    bedtools map \
        -a "${split_bam}.bed" \
        -b "latest.exons.bed.sort" \
        -nonamecheck -s -f 0.95 -c 7 -o distinct -delim '|' \
    | bedtools map \
        -a - -b "latest.genes.bed.sort" \
        -nonamecheck -s -f 0.95 -c 4 -o distinct -delim '|' \
    | sort -k4,4 -k2,2n -k3,3n -S 5G\
    | datamash \
        -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8 \
    | assign-reads-to-genes.py "latest.genes.bed.sort" \
    | awk '\$3 == "exonic" || \$3 == "intronic" {{
            split(\$1, arr, "|")
            printf "%s_%s_%s\t%s\t%s\\n", arr[3], arr[4], arr[5], \$2, \$3
    }}' \
    | sort -k1,1 -k2,2 -S 5G > "${split_bam}_ga.txt"

    rm latest.exons.bed.sort latest.genes.bed.sort
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
    reformat_qc

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

    input:
        set key, file(split_bed), file(split_gene_assign), file(split_umi_count), file(logfile), file(read_count) from for_cat_dups

    output:
        set key, file("*.gz"), file("*_ga.txt"), file("merge_assignment.log") into merge_assignment_out
        set val(key), file("*duplication_rate_stats.txt") into duplication_rate_out
        file "*.bed"  into temp_bed

    """
    # bash watch for errors
    set -ueo pipefail

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
    fraction_per_cell_intronic - fraction of barcode UMIs that are intronic

 Pass through:
    cell_gene_count - gzipped text file with a count of cell, gene pairs

 Summary:
    calculate umis per sample - tabulate_per_cell_counts.py

 Downstream:
    make_matrix
    generate_qc_metrics
    assign_hash (if true)

 Published:
    umis_per_cell - count of umis per cell
    umis_per_cell_intronic - count of umis per cell only from intronic reads

 Notes:

*************/

save_umi_per_cell = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.txt/ + "/umis_per_cell_barcode.txt"}
save_umi_per_int = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.intronic.txt/ + "/intronic_umis_per_cell_barcode.txt"}

process count_umis_by_sample {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_int, pattern: "*intronic.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_cell, pattern: "*UMIs.per.cell.barcode.txt", mode: 'copy'

    input:
        set val(key), file(cell_gene_count), file(gene_assign), file(logfile) from merge_assignment_out

    output:
        set key, file(cell_gene_count), file("count_umis_by_sample.log") into ubss_out
        set key, file("*UMIs.per.cell.barcode.txt") into umis_per_cell
        set key, file("*fraction_intron_barcode.txt") into fraction_per_cell_intronic
        file "*UMIs.per.cell.barcode.intronic.txt" into umi_per_cell_intronic
        set key, file("*UMIs.per.cell.barcode.txt") into for_assign_hash_umis

    """
    # bash watch for errors
    set -ueo pipefail

    cat ${logfile} > count_umis_by_sample.log
    printf "** Start process 'count_umis_by_sample' at: \$(date)\n\n" >> count_umis_by_sample.log
    printf "    Process versions:
        \$(python --version)\n\n" >> count_umis_by_sample.log
    printf "    Process command:
        tabulate_per_cell_counts.py
            --gene_assignment_files "$gene_assign"
            --all_counts_file "${key}.UMIs.per.cell.barcode.txt"
            --intron_counts_file "${key}.UMIs.per.cell.barcode.intronic.txt"
            --intron_fraction_file "${key}.fraction_intron_barcode.txt"\n\n"            >> count_umis_by_sample.log


    tabulate_per_cell_counts.py \
        --gene_assignment_files "$gene_assign" \
        --all_counts_file "${key}.UMIs.per.cell.barcode.txt" \
        --intron_counts_file "${key}.UMIs.per.cell.barcode.intronic.txt" \
        --intron_fraction_file "${key}.fraction_intron_barcode.txt"


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
    run_emptyDrops

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
    publishDir path: "${params.output_dir}/", saveAs: save_umi, pattern: "*umi_counts.mtx", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_anno, pattern: "*cell_annotations.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_gene_anno, pattern: "*gene_annotations.txt", mode: 'copy'

    input:
        set key, file(cell_gene_count), file(logfile), val(gtf_path) from make_matrix_prepped
 
    output:
        set key, file("*cell_annotations.txt"), file("*umi_counts.mtx"), file("*gene_annotations.txt"), val(gtf_path), file("make_matrix.log") into mat_output
        set key, file("*cell_annotations.txt"), file("*umi_counts.mtx"), file("*gene_annotations.txt") into for_rt_cell_by_gene

    """
    # bash watch for errors
    set -ueo pipefail

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

Process: run_emptyDrops

 Inputs:
    key - sample id
    logfile - running log
    umi_matrix - MatrixMarket format matrix of cells by genes
    cell_anno - Cell annotations for umi_matrix
    gene_anno - Gene annotations for umi_matrix
    gtf_path - path to gtf info folder

 Outputs:
    <sample_name>_emptyDrops.RDS

 Pass through:
   cell_data
   umi_matrix
   gene_data
   gtf_path
   logfile (modified)

 Summary:
    Run the emptyDrops utility on the umi_counts.mtx matrix.

 Downstream:
    make_cds

 Published:
    *_emptyDrops.RDS

 Notes:

*************/


save_empty_drops = {params.output_dir + "/" + it - ~/_emptyDrops.RDS/ + "/" + it}

process run_emptyDrops {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_empty_drops, pattern: "*_emptyDrops.RDS", mode: 'copy'

    input:
        set key, file(cell_data), file(umi_matrix), file(gene_data), val(gtf_path), file(logfile) from mat_output

    output:
        set key, file(cell_data), file(umi_matrix), file(gene_data), file("*_emptyDrops.RDS"), val(gtf_path), file("run_emptyDrops.log") into emptyDrops_output
        set key, file("*_emptyDrops.RDS") into for_gen_qc_emptyDrops
        file("*_emptyDrops.RDS") into for_combine_eds

"""
    # bash watch for errors
    set -ueo pipefail

    echo "test 2342324234"

    output_file="${key}_emptyDrops.RDS"

    cat ${logfile} > run_emptyDrops.log
    printf "** Start process 'run_emptyDrops' at: \$(date)\n\n" >> run_emptyDrops.log
    
    if [ "$params.run_emptyDrops" == 'true' ]
    then
      printf "    Process versions:
          \$(R --version | grep 'R version')
              emptyDrops version \$(Rscript -e 'packageVersion("DropletUtils")')\n\n" >> run_emptyDrops.log
      echo '    Process command:
          run_emptyDrops.R
              "$umi_matrix"
              "$cell_data"
              "$gene_data"
              "$key"
              "${key}_emptyDrops.RDS"\n' >> run_emptyDrops.log

      run_emptyDrops.R \
          "$umi_matrix" \
          "$cell_data" \
          "$gene_data" \
          "$key" \
          "${key}_emptyDrops.RDS"
    else
      # make an empty emptyDrops.RDS file
      Rscript -e 'note <- "emptyDrops was skipped"; saveRDS(note, file="${key}_emptyDrops.RDS")'
      printf "    emptyDrops skipped by request\n\n" >> run_emptyDrops.log
    fi

    printf "** End process 'run_emptyDrops' at: \$(date)\n\n" >> run_emptyDrops.log
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
    fraction_per_cell_intronic - fraction of barcode UMIs that are intronic

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
    generate_qc

 Published:

 Notes:

*************/

emptyDrops_output.combine(fraction_per_cell_intronic, by:0).set{emptyDrops_fraction_intronic}
// emptyDrops_fraction_intronic.into{emptyDrops_fraction_intronic_copy01; emptyDrops_fraction_intronic_copy02}

/*
** Diagnostic.
emptyDrops_output.combine(fraction_per_cell_intronic, by:0).into{emptyDrops_fraction_intronic; emptyDrops_fraction_intronic_tmp}
emptyDrops_fraction_intronic_tmp.view()
*/


// save_cds = {params.output_dir + "/" + it - ~/_cds.RDS/ + "/" + it}
// save_cell_qc = {params.output_dir + "/" + it - ~/_cell_qc.csv/ + "/" + it}


make_cds_in = emptyDrops_fraction_intronic.join(reads_per_cell)
process make_cds {
    cache 'lenient'
    // publishDir path: "${params.output_dir}/", saveAs: save_cds, pattern: "*cds.RDS", mode: 'copy'
    // publishDir path: "${params.output_dir}/", saveAs: save_cell_qc, pattern: "*cell_qc.csv", mode: 'copy'


    input:
      set key, file(cell_data), file(umi_matrix), file(gene_data), file(emptyDrops), val(gtf_path), file(logfile), file(fraction_intron_barcode), file(reads_per_cell) from make_cds_in

    output:
        set key, file("*for_scrub.mtx"), file("*_cds.mobs"), file("*cell_qc.csv"), file("make_cds.log") into cds_out
        file("*cell_emptyDrops.csv") into cell_eds
        file("*cell_emptyDrops.csv") into for_combine_cell_counts

    """
    # bash watch for errors
    set -ueo pipefail

    echo "testing"

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
            "$emptyDrops"
            "$fraction_intron_barcode"
            "$key"
            "$params.umi_cutoff"
            "$reads_per_cell"\n' >> make_cds.log

    make_cds.R \
        "$umi_matrix"\
        "$cell_data"\
        "$gene_data"\
        "${gtf_path}/latest.genes.bed"\
        "$emptyDrops"\
        "$fraction_intron_barcode"\
        "$key"\
        "$params.umi_cutoff"\
        "$reads_per_cell"

    printf "** End process 'make_cds' at: \$(date)\n\n" >> make_cds.log
    """
}


process apply_garnett {
    cache 'lenient'

    input:
        set key, file(scrub_matrix), file(cds_object), file(cell_qc), file(logfile) from cds_out

    output:
        set key, file(scrub_matrix), file("new_cds/*.mobs"), file(cell_qc), file("apply_garnett.log") into for_scrub

"""
    # bash watch for errors
    set -ueo pipefail

    cat ${logfile} > apply_garnett.log
    printf "** Start process 'apply_garnett' at: \$(date)\n\n" >> apply_garnett.log
    mkdir new_cds
    echo "No Garnett classifier provided for this sample" > garnett_error.txt
    if [ $params.garnett_file == 'false' ]
    then
        cp -r $cds_object new_cds/
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
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_hist, pattern: "*png", mode: 'copy'

    input:
        set key, file(scrub_matrix), file(cds_object), file(cell_qc), file(logfile) from for_scrub

    output:
        set key, file("*scrublet_out.csv"), file(cds_object), file(cell_qc) into scrublet_out
        file ("*.png") into scrub_pngs
        set key, file("run_scrublet.log") into pipe_log

    """
    # bash watch for errors
    set -ueo pipefail

    cat ${logfile} > run_scrublet.log
    printf "** Start process 'run_scrublet' at: \$(date)\n\n" >> run_scrublet.log
    printf "    Process versions:
        \$(python --version)
            \$(pip freeze | grep scrublet | tr '==' ' ')\n\n" >> run_scrublet.log

    if [ $params.skip_doublet_detect == 'false' ]
    then
        run_scrublet.py --key $key --mat $scrub_matrix
        echo '    Process command:
        run_scrublet.py --key $key --mat $scrub_matrix\n'  >> run_scrublet.log
    else
        run_scrublet.py --key $key --mat $scrub_matrix --skip
        echo '    Process command:
        run_scrublet.py --key $key --mat $scrub_matrix --skip\n'  >> run_scrublet.log
        printf "    Scrublet skipped by request\n\n" >> run_scrublet.log
    fi

    printf "** End process 'run_scrublet' at: \$(date)\n\n" >> run_scrublet.log

    printf "** Start processes to generate qc metrics and dashboard at: \$(date)\n\n" >> run_scrublet.log
    """

}


/*************

Process: reformat_qc

 Inputs:
    key - sample id
    cds_object - cds object in RDS format
    cell_qc - csv of cell quality control information
    scrub_csv - scrublet results csv
    dup_stats - file with duplication rate information for the sample

 Outputs:
    key - sample id
    cds_object - cds object in RDS format
    for_hash_cds_dir - directory to cds object for assigning hash if true
    sample_stats - csv with sample-wise statistics
    cell_qc - csv of cell quality control information
    collision - file containing collision rate if barnyard sample
    cds_dir - temp directory containing cds object and cell qc csv for assigning hash

 Pass through:

 Summary:
    Add scrublet info to cell_qc and cds object
    Calculate sample statistics

 Downstream:
    generate_qc_metrics
    zip_up_sample_stats
    assign_hash (if true)
    collapse_collision

 Published:
    sample_stats - csv with sample-wise statistics

 Notes:
    Nextflow dsl1 doesn't allow for conditional channel output or publishing, "temp_dir" is used 
    as a work around to publish one final cds object when hash assignment is true 
    without modifying the input cds object. 

*************/

scrublet_out.join(duplication_rate_out).set{reformat_qc_in}

save_cds = {params.output_dir + "/" + it - ~/_cds.mobs/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_cell_qc = {params.output_dir + "/" + it - ~/_cell_qc.csv/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_samp_stats = {params.output_dir + "/" + it - ~/_sample_stats.csv/ + "/" + it}

process reformat_qc {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_cds, pattern: "temp_fold/*cds.mobs", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_qc, pattern: "temp_fold/*cell_qc.csv", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_samp_stats, pattern: "*sample_stats.csv", mode: 'copy'

    input:
        set key, file(scrub_csv), file(cds_object), file(cell_qc), file(dup_stats) from reformat_qc_in

    output:
        // set key, file("temp_fold/*.mobs"), file("temp_fold/*.csv") into rscrub_out
        set key, file("temp_fold/*.mobs") into rscrub_out
        file("temp_fold/*.mobs") into for_combine_cds
        file("*sample_stats.csv") into sample_stats
        // file("*collision.txt") into collision
        set key, file("temp_fold") into temp_dir

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(monocle3)
        library(BPCells)
    })

    dir.create("temp_fold")
    cell_qc <- read.csv("$cell_qc")
    dup_stats <- read.table(paste0("$key", ".duplication_rate_stats.txt"))

    if (file.size("$cds_object/bpcells_matrix_dir/col_names") == 0) {
        file.copy("$cds_object", "temp_fold/", recursive=TRUE)
        file.copy("$cell_qc", "temp_fold/")

        df <- data.frame(sample="$key", n.reads = dup_stats\$V2, n.umi = dup_stats\$V3,
                    median_umis = "NA",
                    median_perc_mito_umis = "NA",
                    duplication_rate = dup_stats\$V4)

        write.csv(df, file=paste0("$key", "_sample_stats.csv"), quote=FALSE, row.names=FALSE)
        quit(save="no", status=0)
    } 

    cds <- load_monocle_objects("$cds_object")
    cell_qc <- read.csv("$cell_qc")
    dup_stats <- read.table(paste0("$key", ".duplication_rate_stats.txt"))


    df <- data.frame()

    if(nrow(pData(cds)) > 0) {
        if("$params.skip_doublet_detect" == 'false') {
            scrublet_out <- read.csv("$scrub_csv", header=F)
            pData(cds)\$scrublet_score <- scrublet_out\$V1
            cell_qc\$scrublet_score <- scrublet_out\$V1
        }
        df <- data.frame(sample="$key", n.reads = dup_stats\$V2, n.umi = dup_stats\$V3,
                median_umis = median(pData(cds)\$n.umi),
                median_perc_mito_umis = round(median(pData(cds)\$perc_mitochondrial_umis), 1),
                duplication_rate = dup_stats\$V4)
    } else {
        df <- data.frame(sample="$key", n.reads = dup_stats\$V2, n.umi = dup_stats\$V3,
                    median_umis = "NA",
                    median_perc_mito_umis = "NA",
                    duplication_rate = dup_stats\$V4)
    }

    write.csv(cell_qc, quote=FALSE, file="temp_fold/$cell_qc")
    write.csv(df, file=paste0("$key", "_sample_stats.csv"), quote=FALSE, row.names=FALSE)

    suppressMessages(
    save_monocle_objects(cds, directory_path="temp_fold/$cds_object", archive_control=list(archive_type='none', archive_compression='none'))
    )

    """
}

// Temp dir is used in "assign_hash",  "publish_cds_and_cell_qc", and "combine_cds" process block.
// It's copied here because these blocks are conditional and temp dir is required 
// to publish cds object and cell qc. 

temp_dir.into{temp_dir_copy01; temp_dir_copy02}


// collect all cds objects
// find sample name by splitting the cds object file name 
// set name as key to group the cds objects ; split by sample_rt 

/*************

Process: combine_cds

 Inputs:
    for_combine_cds - temp directory containing rt split cds objects

 Outputs:
    combined_cds_input - combined rt split cds objects by sample

 Pass through:

 Summary:
    Collect all rt split cds objects. Find sample name by splitting the cds object 
    file name. Set name as key to group cds objects by sample name. Then combine 
    all cds objects for assign_hash process. 

 Downstream:
    
 Published:

 Notes:
    runs only when params.hash_rt_split = true
    
*************/

combine_cds_input = for_combine_cds
    .map { file ->
        def fname = file.getName() // e.g. "GENE3.P1.A02_cds.RDS"
        def base = fname.replaceFirst(/_cds\.mobs$/, '') // "GENE3.P1.A02"
        def sample_name = base.split(/\.P[0-9]\.[A-H][0-9]{2}/)[0] // "GENE3"
        tuple(sample_name, file)
    }
    .groupTuple()
    // .view()
                        
process combine_cds {
    cache 'lenient'

    input:
        // set key, file(input_cds) from rscrub_out_for_combine_cds.collect()
        tuple val(sample_name), path(cds_list) from combine_cds_input


    output:
        set val(sample_name), file("*combined_cds.mobs") into combined_cds_out 
        set val(sample_name), file("*combined_cds.mobs") into for_gen_qc_cds 

    when: 
        params.hash_rt_split != false
    
    script: 

    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages({
        library(monocle3)
        library(data.table)
    })

    new_cds_list = strsplit("$cds_list", " ")[[1]]

    if (length(new_cds_list) < 2) {
        dir.create(paste0("$sample_name", "_combined_cds.mobs"))
        file.copy(from = list.files(new_cds_list[1], full.names = TRUE),
          to = paste0("$sample_name", "_combined_cds.mobs"),
          recursive = TRUE)
        quit(save="no", status=0)
    }

    temp_cds_list <- list()

    for(cds in new_cds_list) {
        temp_cds <- load_monocle_objects(cds) 
        temp_cds_list[[cds]] <- temp_cds
    }

    cds <- combine_cds(temp_cds_list) 
    cds\$sample <- "$sample_name"
    rownames(colData(cds)) <- cds\$cell
    
   # saveRDS(cds, paste0("$sample_name", "_combined_cds.RDS"))
    save_monocle_objects(cds,directory_path=paste0("$sample_name", '_combined_cds.mobs'), archive_control=list(archive_type='none', archive_compression='none'))
    

    """
}


/*************

Process: collect_stats

 Inputs:
    sample_stat - csv with rt-level statistics - flattened

 Outputs:
    combined_sample_stats - sample-level statistics

 Pass through:

 Summary:
    Sum up rt-level stats into sample-level stats 

 Downstream:
    zip_up_sample_stats
    generate_dashboard

 Published:

 Notes:

*************/

sample_stats_in = sample_stats.flatten()
    .map { stats ->
        def fname = stats.getName()
        // def sample_name = sample.split(/\.P[0-9]{1,2}\.[A-H][0-9]{1,2}/)[0] 
        def sample_name = fname.split(/\.P[0-9]{1,2}\.[A-H][0-9]{1,2}/)[0] 
        tuple(sample_name, stats)
    }
    .groupTuple()

process collect_stats {
    cache 'lenient' 

    input: 
        tuple val (sample_name), path(stats) from sample_stats_in 

    output: 
        file ("*_sample_stats.csv") into combined_sample_stats

    when:
        params.hash_rt_split != false
    
    """
    # bash watch for errors
    set -ueo pipefail
    
    echo "test"

    for file in ${stats}; do
        echo "sample,n.reads,n.umi,median_umis,median_perc_mito_umis,duplication_rate" > "${sample_name}_sample_stats.csv"
         awk -F"," 'BEGIN {OFS=","}  NR>1 {sum_reads+=\$2; sum_umis +=\$3} END {print ${sample_name}, sum_reads, sum_umis, \$4, \$5, \$6}' "\$file" 
    done >> "${sample_name}_sample_stats.csv"
    
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


zip_up_sample_stats_in = null

if (params.hash_rt_split != false) {
   zip_up_sample_stats_in = combined_sample_stats.collect()
} else {
   zip_up_sample_stats_in =  sample_stats.collect()
}

process zip_up_sample_stats {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", pattern: "all_sample_stats.csv", mode: 'copy'

    input:
        // file files from sample_stats.collect()
        // file files from combined_sample_stats.collect()
        file files from zip_up_sample_stats_in

    output:
        file "*ll_sample_stats.csv" into all_sample_stats

    """
    # bash watch for errors
    set -ueo pipefail

    echo "testingasdkjf"

    sed -s 1d $files > all_sample_stats.csv

    """
}



/*************

Process: combine_cell_counts

 Inputs:
    for_combine_cell_counts - text file containing emptyDrops cell counts at rt-level

 Outputs:
    combined_cell_counts - combined emptyDrops cell counts at sample-level

 Pass through:
    

 Summary:
    Collect all rt split emptyDrops cell counts. Find sample name by splitting cell counts
    file name. Set name as key to group cell counts by sample name. Then combine 
    all cell counts at the sample-level for generate_dashboard

 Downstream:
    generate_dashboard

 Published:

 Notes:
    runs only when params.hash_rt_split = true
    
*************/

combine_cell_counts_in = for_combine_cell_counts.flatten()
    .map { counts ->
        def fname = counts.getName()
        // def sample_name = sample.split(/\.P[0-9]{1,2}\.[A-H][0-9]{1,2}/)[0] 
        def sample_name = fname.split(/\.P[0-9]{1,2}\.[A-H][0-9]{1,2}/)[0] 
        tuple(sample_name, counts)
    }
    .groupTuple()


process combine_cell_counts {
    cache 'lenient'
    
    input: 
        tuple val (sample_name), path(counts) from combine_cell_counts_in

    output:
        file("*.csv") into combined_cell_eds
        
    when: 
        params.hash_rt_split != false

    """
    # bash watch for errors
    set -ueo pipefail

    cat ${counts} >> "${sample_name}_cell_emptyDrops.csv"

    """
    
}

/*************

Process: calc_cell_totals

 Inputs:
    cell_ed - csv of cell quality control information with emptyDrops FDR values - collected

 Outputs:
    cell_counts - table cell totals above set UMI thresholds for all samples

 Pass through:

 Summary:
    Count cell totals above set UMI thresholds for all samples

 Downstream:
    generate_dashboard

 Published:
    cell_counts - table cell totals above set UMI thresholds for all samples
 Notes:

*************/

calc_cell_total_in = cell_eds.collect()

if (params.hash_rt_split) {
    calc_cell_total_in = combined_cell_eds.collect()
}


process calc_cell_totals {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", pattern: "cell_counts.txt", mode: 'copy'

    input:
        // file(cell_ed) from cell_eds.collect()
        file(cell_ed) from calc_cell_total_in

    output:
        file "*.txt" into cell_counts

    """
    # bash watch for errors
    set -ueo pipefail
    
    rm -f cell_counts.txt
    for f in $cell_ed
    do
      awk 'BEGIN {FS=","; counter=0} {if(\$2 ~ /^[0-9]+\$/ && \$2 > 100)  {counter++}} END{print FILENAME, "100", counter}' \$f >> cell_counts.txt
      awk 'BEGIN {FS=","; counter=0} {if(\$2 ~ /^[0-9]+\$/ && \$2 > 500)  {counter++}} END{print FILENAME, "500", counter}' \$f >> cell_counts.txt
      awk 'BEGIN {FS=","; counter=0} {if(\$2 ~ /^[0-9]+\$/ && \$2 > 1000) {counter++}} END{print FILENAME, "1000", counter}' \$f >> cell_counts.txt

      awk 'BEGIN {FS=","; counter=0} {if(NF >= 3) {if(\$2 ~ /^[0-9]+\$/ && \$2 > 100 && \$3 <= 0.01) {counter++}}else{counter="-"}} END{print FILENAME, "FDR_p01", counter}' \$f >> cell_counts.txt
      awk 'BEGIN {FS=","; counter=0} {if(NF >= 3) {if(\$2 ~ /^[0-9]+\$/ && \$2 > 100 && \$3 <= 0.001) {counter++}}else{counter="-"}} END{print FILENAME, "FDR_p001", counter}' \$f >> cell_counts.txt
    done

    """

}

/*************

Process: combine_eds



*************/


combine_eds_input = for_combine_eds
    .map { file ->
        def fname = file.getName() // e.g. "GENE3.P1.A02_emptyDrops.RDS"
        def base = fname.replaceFirst(/_emptyDrops\.RDS$/, '') // "GENE3.P1.A02"
        def sample_name = base.split(/\.P[0-9]\.[A-H][0-9]{2}/)[0] // "GENE3"
        tuple(sample_name, file)
    }
    .groupTuple()


save_combined_eds = {params.output_dir + "/" + it - ~/_emptyDrops.RDS/ + "/" + it}
    

process combine_eds {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_combined_eds, pattern: "*_emptyDrops.RDS", mode: 'copy'

    input:
        tuple val(sample_name), path(eds_list) from combine_eds_input


    output:
        set val(sample_name), file("*.RDS") into combined_eds_out 

    when: 
        params.hash_rt_split != false
    
    script: 

     """
    #!/usr/bin/env Rscript

    new_eds_list = strsplit("$eds_list", " ")[[1]]
    new_eds_list <- lapply(new_eds_list, readRDS)
    combined_df <- do.call(rbind, new_eds_list)

    saveRDS(combined_df, "${sample_name}_emptyDrops.RDS")
    """

}


/*************

Process: assign_hash

 Inputs:
    cds_dir - temp directory containing cds object and cell qc csv
    hash_cell - text file with list of cell names with hash umis 
    hash_list - text file with list of hash names 
    hash_mtx - sparse matrix with hash umi counts for each cell 
    umis_per_cell - txt file with number of umis per cell barcode

 Outputs:
    corrected_hash_table - csv file with data frame of cells and hash stats 
    hash_cds - cds object with hash info
    cell_qc - csv of cell quality control information

 Pass through:

 Summary:
    Assign hash to cells and find top hash oligo for each cell 

 Downstream:
    
 Published:
    hash_table - csv file with data frame of cells and hash stats 
    cds - cds object with hash info in RDS format
    cell_qc - csv of cell quality control information

 Notes:
    runs only when params.hash_list = true
    
*************/

for_assign_hash_umis.into{for_assign_hash_umis_copy01; for_assign_hash_umis_copy02}
make_hash_cds = temp_dir_copy01.join(hash_mats_for_assign_hash).join(for_assign_hash_umis_copy01)

save_hash_cds = {params.output_dir + "/" + it - ~/_cds.RDS/ + "/" + it}
// save_cell_qc = {params.output_dir + "/" + it - ~/_cell_qc.csv/ + "/" + it}
save_hash_table = {params.output_dir + "/" + it - ~/_hash_table.csv/ + "/" + it}

process assign_hash {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_cds, pattern: "*_cds.RDS", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_table, pattern: "*_hash_table.csv", mode: 'copy'

    input:
        set key, file(cds_dir), file(hash_mtx), file(hash_cell), file(hash_hash), file(umis_per_cell) from make_hash_cds

    output:
        file("*cds.RDS") into hash_cds
        file("*hash_table.csv") into hash_table

    when: 
        params.hash_list != false && params.hash_rt_split == false

    """
    # bash watch for errors
    set -ueo pipefail

    cp ${cds_dir}/*.csv .

    assign_hash.R \
        $key \
        $hash_mtx \
        $hash_cell \
        $hash_hash \
        ${cds_dir}/*cds.RDS \
        $umis_per_cell \
        $params.hash_umi_cutoff \
        $params.hash_ratio

    """
}


// concat all rt split cell by gene files 


concat_cell_by_gene = for_rt_cell_by_gene.join(for_assign_hash_umis_copy02)
    .map {sample, cell_anno, umi_mtx, gene_anno, umi_per_cell ->
        def sample_name = sample.split(/\.P[0-9]\.[A-H][0-9]{2}/)[0] // "GENE3"
        tuple(sample_name, cell_anno, umi_mtx, gene_anno, umi_per_cell)
    }
    .groupTuple()

save_merged_umi = {params.output_dir + "/" + it - ~/UMIs.per.cell.barcode.txt/  - ~/temp_fold/  + "/umis_per_cell_barcode.txt"}
save_merged_features = {params.output_dir + "/" + it - ~/gene_annotations.txt/  - ~/temp_fold/  + "/gene_annotations.txt"}
save_merged_counts = {params.output_dir + "/" + it - ~/umi_counts.mtx/  - ~/temp_fold/  + "/umi_counts.mtx"}
save_merged_cells = {params.output_dir + "/" + it - ~/cell_annotations.txt/  - ~/temp_fold/  + "/cell_annotations.txt"}


process merge_rt_cell_by_gene {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_merged_umi, pattern: "temp_fold/*/*UMIs.per.cell.barcode.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_merged_features, pattern: "temp_fold/*/*gene_annotations.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_merged_counts, pattern: "temp_fold/*/*umi_counts.mtx", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_merged_cells, pattern: "temp_fold/*/*cell_annotations.txt", mode: 'copy'

    input:
        tuple val(sample_name), path(cell_anno), path(umi_mtx), path(gene_anno), path(umi_per_cell) from concat_cell_by_gene

    output:
        set val(sample_name), file("temp_fold/*/*UMIs.per.cell.barcode.txt") into concat_umi_out
        set val(sample_name), file("temp_fold/*/*UMIs.per.cell.barcode.txt") into for_gen_qc_umi
        set val(sample_name), file("temp_fold/*/*.mtx"), file("temp_fold/*/*.txt") into merge_rt_umis_out

    when:
        params.hash_rt_split != false && params.hash_list != false

    script:

    """
    # bash watch for errors
    set -ueo pipefail
    
    mkdir -p temp_fold/${sample_name}

    cat_sparse_matrix.py \
        -i ${umi_mtx} \
        -o temp_fold/${sample_name}/ \
        -c UMIs.per.cell.barcode.txt \
        -m umi_counts.mtx \
        -f gene_annotations.txt \
        -a cell_annotations.txt

    """
}

/*************

Process: assign_hash_rt_split

 Inputs:
    cds_dir - temp directory containing cds object and cell qc csv
    hash_cell - text file with list of cell names with hash umis 
    hash_list - text file with list of hash names 
    hash_mtx - sparse matrix with hash umi counts for each cell 
    umis_per_cell - txt file with number of umis per cell barcode

 Outputs:
    corrected_hash_table - csv file with data frame of cells and hash stats 
    hash_cds - cds object with hash info
    cell_qc - csv of cell quality control information

 Pass through:

 Summary:
    Assign hash to cells and find top hash oligo for each cell 

 Downstream:
    
 Published:
    hash_table - csv file with data frame of cells and hash stats 
    cds - cds object with hash info in RDS format
    cell_qc - csv of cell quality control information

 Notes:
    runs only when params.hash_list = true
    
*************/

// need to combine hash matrices, hash cell barcodes, hash names and umis per-cell 


assign_hash_rt_in = combined_cds_out.join(concat_hash_out).join(concat_umi_out)
save_combined_hash_cds = {params.output_dir + "/" + it - ~/_cds.mobs/ + "/" + it}
save_combined_hash_table = {params.output_dir + "/" + it - ~/_hash_table.csv/ + "/" + it}

process assign_hash_rt_split {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_combined_hash_cds, pattern: "*_cds.mobs", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_combined_hash_table, pattern: "*_hash_table.csv", mode: 'copy'

    input:
        set key, file(cds), file(hash_mtx), file(hash_hash), file(hash_cell), file(umis_per_cell) from assign_hash_rt_in
    output:
        file("*_cds.mobs") into combined_hash_cds
        file("*hash_table.csv") into combined_hash_table
        set key, file("*_cds.mobs") into combined_hash_cds_for_qc

    when: 
        params.hash_rt_split != false 

    """
    # bash watch for errors
    set -ueo pipefail

    assign_hash.R \
        $key \
        $hash_mtx \
        $hash_cell \
        $hash_hash \
        $cds \
        $umis_per_cell \
        $params.hash_umi_cutoff \
        $params.hash_ratio

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
    wellcheck_png - png of RT barcode qc 
    rt_stats - csv of RT barcode stats

 Pass through:

 Summary:
    Generate a bunch of qc metrics and plots - generate_qc.R

 Downstream:
    generate_dashboard

 Published:
    umap_png - png sample UMAP
    knee_png - png sample knee plot
    qc_png - png of cell qc stats
    wellcheck_png - png of RT barcode qc 

 Notes:

*************/


if (params.hash_rt_split != false) {
    for_gen_qc = combined_hash_cds_for_qc.join(for_gen_qc_umi).join(combined_eds_out)
} else {
    for_gen_qc = rscrub_out.join(umis_per_cell).join(for_gen_qc_emptyDrops)
}

save_knee = {params.output_dir + "/" + it - ~/_knee_plot.png/ + "/" + it}
save_umap = {params.output_dir + "/" + it - ~/_UMAP.png/ + "/" + it}
save_cellqc = {params.output_dir + "/" + it - ~/_cell_qc.png/ + "/" + it}
save_garnett = {params.output_dir + "/" + it.split("_")[0] + "/" + it}
save_umi_rt_stats = {params.output_dir + "/" + it - ~/_umi_rt_stats.csv/ + "/" + it}
save_mito_rt_stats = {params.output_dir + "/" + it - ~/_mito_rt_stats.csv/ + "/" + it}
save_wellcheck_combo = {params.output_dir + "/" + it - ~/_wellcheck.png/ + "/" + it}
// save_empty_hash_plot = {params.output_dir + "/" + it - ~/_hash_knee_plot.png/ + "/" + it}

process generate_qc_metrics {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_umap, pattern: "*UMAP.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_knee, pattern: "*knee_plot.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cellqc, pattern: "*cell_qc.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_garnett, pattern: "*Garnett.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_umi_rt_stats, pattern: "*_umi_rt_stats.csv", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_mito_rt_stats, pattern: "*_mito_rt_stats.csv", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_wellcheck_combo, pattern: "*_wellcheck.png", mode: 'copy'
    // publishDir path: "${params.output_dir}/", saveAs: save_empty_hash_plot pattern: "*_hash_knee_plot.png", mode: 'copy'

    input:
        set key, file(cds_object), file(umis_per_cell), file(emptydrops) from for_gen_qc
         
    output:
        file("*.png") into qc_plots
        file("*.txt") into cutoff
        file("*collision.txt") into collision


    """
    # bash watch for errors
    set -ueo pipefail

    echo "test123"
    generate_qc.R\
        $cds_object $umis_per_cell $key $emptydrops $params.hash_list \
        --specify_cutoff $params.umi_cutoff

    """
}



/*************

Process: calc_hash_dup

 Inputs:
    sorted_hash_combined - all sorted hash umis combined

 Outputs:
    hash_knee - unique hash umi by RT barcode knee plot in .png format


 Pass through:

 Summary:
    Takes in a combined, sorted hash umi file and calculates the hash duplication rate per cell.
    First, sums up the total hash reads per cell (non-unique). Second, sums up unique hash umis per cell. 
    Third, output a hash table with unique hash umis found for each cell. 
    Fourth, produces a knee plot with unique hash umis by RT barcode.
    Fifth, calculates the hash duplication rate per by (1-(unique hash umis/ total hash read for thatcell)).

 Downstream:
    
 Published:
    hash_knee - unique hash umi by RT barcode knee plot in .png format
    hash_results - .txt files with hash reads per cell, unique hash umis per cell, and a table of hash umis for each cell


 Notes:
    runs only when params.hash_list != false and params.hash_dup!= false
    
*************/

save_hash_knee = {params.output_dir + "/" + it - ~/_hash_knee_plot.png/ + "/" + it}

process calc_hash_dup_cell {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_hash_knee, pattern: "*hash_knee_plot.png", mode: 'copy'

    input:
        set key, file(hash_reads_per_cell), file(hash_umis_per_cell), file(hash_assigned_table)from hash_results

    output:
        file("*.png") into hash_knee

    when: 
        params.hash_list != false && params.hash_dup != false
    
    script:

    """
    # bash watch for errors
    set -ueo pipefail

    knee-plot.R \
    "$hash_umis_per_cell" \
    $key 

    """
}


/*************

Process: calc_tot_hash_dup

 Inputs:
    hash_dup_per_cell - hash duplication rate per cell 

 Outputs:
    total_hash_dup - total hash duplication rate per sample by pcr plate

 Pass through:

 Summary:
    Calculates total hash duplication rate for each PCR plate. 

 Downstream:
    
 Published:
    total_hash_dup - total hash duplication rate per sample by pcr plate

 Notes:
    runs only when params.hash_list != false and params.hash_dup!= false

*************/


for_hash_calc.into{ for_hash_calc_copy01;for_hash_calc_copy02}

save_total_hash_dup = {params.output_dir + "/" + it - ~/_total_hash_dup_rate.csv/ + "/" + it}

process calc_tot_hash_dup {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_total_hash_dup, pattern: "*_total_hash_dup_rate.csv", mode: 'copy'

    input:
        set key, file(hash_dup) from for_hash_calc_copy01

    output:
        file("*.csv") into total_hash_dup

    when: 
        params.hash_list != false && params.hash_dup != false
    
    
    script:

    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(tidyverse)

    dup = fread("${key}_hash_dup_per_cell.txt", header = FALSE,
                data.table = F,
                col.names = c("Expt", "Cell", "V4", "V5", "V6"))
     
    dup_rate = NULL
    
    if (dim(dup)[1] != 0) {
        dup = dup %>%
            separate(Cell, into = c("p5", "p7", "rt_plate_well", "lig_well"), sep = "_") %>%
            mutate(pcr_plate = paste(str_sub(p7, start = 1, end = 1), str_sub(p5, start = 2, end = 3), sep = ""))

        if ("$params.hash_dup" == 'pcr_plate') {
            dup = dup %>% group_by(pcr_plate) 
        } else if("$params.hash_dup" == 'p5') {
            dup = dup %>% group_by(p5) 
        } else {
            stop("params.hash_dup must be either 'pcr_plate' or 'p5'.")     
        }

        dup_rate = dup %>% summarize(dup_rate = 1-(sum(V4)/sum(V5))) %>% data.frame()
    }

    out <- file(paste0("$key", "_total_hash_dup_rate.csv"))
    write.csv(dup_rate, file = out, row.names = FALSE, quote = FALSE)

    """
}

// calculate total hash duplication rate for combined RT split sample

calc_hash_in = for_hash_calc_copy02
    .map { sample, hash_dup ->
        def sample_name = sample.split(/\.P[0-9]\.[A-H][0-9]{2}/)[0] // "GENE3"
        tuple(sample_name, hash_dup)
    }
    .groupTuple()

save_combined_hash_dup = {params.output_dir + "/" + it - ~/_total_hash_dup_rate.csv/ + "/" + it}

process calc_tot_hash_dup_combined {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_combined_hash_dup, pattern: "*_total_hash_dup_rate.csv", mode: 'copy'

    input:
        tuple(sample_name), path(hash_dup) from calc_hash_in

    output:
        file("*.csv") into combined_total_hash_dup

    when: 
        params.hash_rt_split != false && params.hash_dup != false
    
    
    script:

    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(tidyverse)

    hash_files = strsplit("$hash_dup", " ")[[1]]

    dup = data.frame()

    for (f in hash_files) {
        temp_dup = fread(f, header = FALSE,
                data.table = F,
                col.names = c("Expt", "Cell", "V4", "V5", "V6"))
        dup = rbind(dup, temp_dup)
    }

    dup_rate = NULL
    
    if (dim(dup)[1] != 0) {
        dup = dup %>%
            separate(Cell, into = c("p5", "p7", "rt_plate_well", "lig_well"), sep = "_") %>%
            mutate(pcr_plate = paste(str_sub(p7, start = 1, end = 1), str_sub(p5, start = 2, end = 3), sep = ""))

        if ("$params.hash_dup" == 'pcr_plate') {
            dup = dup %>% group_by(pcr_plate) 
        } else if("$params.hash_dup" == 'p5') {
            dup = dup %>% group_by(p5) 
        } else {
            stop("params.hash_dup must be either 'pcr_plate' or 'p5'.")     
        }

        dup_rate = dup %>% summarize(dup_rate = 1-(sum(V4)/sum(V5))) %>% data.frame()
    }

    out <- file(paste0("$sample_name", "_total_hash_dup_rate.csv"))
    write.csv(dup_rate, file = out, row.names = FALSE, quote = FALSE)

    """
}

/*************

Process: publish_cds_and_cell_qc

 Inputs: 
    temp_dir - directory containing monocle cds object and cell qc summary 

 Outputs: 
    cds_obj - monocle cds object in RDS format 
    cell_qc - csv of cell quality control information

 Summary: 
    Publish monocle cds object and cell qc when not a hash experiment

 Pass through: 

 Downstream: 

 Published:
    cds - cds object in RDS format
    cell_qc - csv of cell quality control information

 Notes: 
    Runs only when params.hash_list == false

    Nextflow dsl1 doesn't allow for conditional channel output or publishing, "temp_dir" is used 
    as a work around to publish one final cds object when hash assignment is true 
    without modifying the input cds object. This conditional process is used to publish the initial
    cds object (without hash info) and cell qc when the hash processes are skipped.

*************/

save_cds = {params.output_dir + "/" + it - ~/_cds.RDS/ + "/" + it}
save_cell_qc = {params.output_dir + "/" + it - ~/_cell_qc.csv/ + "/" + it}

process publish_cds_and_cell_qc {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_cds, pattern: "*cds.RDS", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_qc, pattern: "*cell_qc.csv", mode: 'copy'

    input:
        set key,file(cds_dir) from temp_dir_copy02

    output:
        file("*cds.RDS") into pub_cds
        file("*cell_qc.csv") into pub_cell_qc
    
    when:
        params.hash_list == false

    """
    # bash watch for errors
    set -ueo pipefail
 
    cp $cds_dir/*.RDS . 
    cp $cds_dir/*.csv . 

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
    cache 'lenient'

    input:
        file col_file from collision.collect()

    output:
        file "*.txt" into all_collision

    """
    # bash watch for errors
    set -ueo pipefail

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
    wellcheck_png - png of RT barcode qc stats - combined as qc_plots
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


// include a condition to use emptyDrops original cell counts file or rt-level mod cell counts

// cell_counts_in = cell_counts
// if (params.hash_rt_split) {
//     cell_counts_in = combined_cell_counts
// }


process generate_dashboard {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", pattern: "exp_dash", mode: 'copy'

    input:
        file all_sample_stats
        file cell_counts
        // file cell_counts_in
        file all_collision
        file plots from qc_plots.collect()
        file scrublet_png from scrub_pngs.collect()

    output:
        file exp_dash into exp_dash_out

    """
    # bash watch for errors
    set -ueo pipefail

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
    # bash watch for errors
    set -ueo pipefail

    head -n 2 ${logfile} > ${key}_full.log
    printf "Nextflow version: $nextflow.version\n" >> ${key}_full.log
    printf "Pipeline version: $workflow.manifest.version\n" >> ${key}_full.log
    printf "Git Repository, Version, Commit ID, Session ID: $workflow.repository, $workflow.revision, $workflow.commitId, $workflow.sessionId\n\n" >> ${key}_full.log
    printf "Command:\n$workflow.commandLine\n\n" >> ${key}_full.log
    printf "***** PARAMETERS *****: \n\n" >> ${key}_full.log
    printf "    params.run_dir:               $params.run_dir\n" >> ${key}_full.log
    printf "    params.output_dir:            $params.output_dir\n" >> ${key}_full.log
    printf "    params.sample_sheet:          $params.sample_sheet\n" >> ${key}_full.log
    printf "    params.demux_out:             $params.demux_out\n" >> ${key}_full.log
    printf "    params.level:                 $params.level\n" >> ${key}_full.log
    printf "    params.max_cores:             $params.max_cores\n" >> ${key}_full.log
    printf "    params.samples:               $params.samples\n" >> ${key}_full.log
    printf "    params.star_file:             $params.star_file\n" >> ${key}_full.log
    printf "    params.gene_file:             $params.gene_file\n" >> ${key}_full.log
    printf "    params.umi_cutoff:            $params.umi_cutoff\n" >> ${key}_full.log
    printf "    params.rt_barcode_file:       $params.rt_barcode_file\n" >> ${key}_full.log
    printf "    params.hash_list:             $params.hash_list\n" >> ${key}_full.log
    printf "    params.hash_umi_cutoff:       $params.hash_umi_cutoff\n" >> ${key}_full.log
    printf "    params.hash_ratio:            $params.hash_ratio\n" >> ${key}_full.log
    printf "    params.max_wells_per_sample:  $params.max_wells_per_sample\n\n" >> ${key}_full.log
    printf "    params.garnett_file:          $params.garnett_file\n\n" >> ${key}_full.log
    printf "    params.skip_doublet_detect:   $params.skip_doublet_detect\n\n" >> ${key}_full.log

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
    # bash watch for errors
    set -ueo pipefail

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
  if( majorVersion < minMajorVersion || ( majorVersion == minMajorVersion && minorVersion < minMinorVersion ) )
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


/*
** getOSInfo()
**
** Purpose: get information about the operating system.
**
** Returns:
**    list of strings with OS name, OS distribution, OS distribution release
**
** Notes:
**   o  limited to Linux operating systems at this time
*/
def getOSInfo()
{
  def osName = System.properties['os.name']
  def osDistribution
  def osRelease
  if( osName == 'Linux' )
  {
    def proc
    proc = "lsb_release -a".execute() | ['awk', 'BEGIN{FS=":"}{if($1=="Distributor ID"){print($2)}}'].execute()
    proc.waitFor()
    osDistribution = proc.text.trim()
    proc = "lsb_release -a".execute() | ['awk', 'BEGIN{FS=":"}{if($1=="Release"){print($2)}}'].execute()
    proc.waitFor()
    osRelease = proc.text.trim()
  }
  return( [ osName, osDistribution, osRelease ] )
}
