
// Parse input parameters
params.help = false
params.samples = false
params.star_file = "$baseDir/bin/star_file.txt"
params.gene_file = "$baseDir/bin/gene_file.txt"
params.umi_cutoff = 100
params.level = 3
params.align_mem = 80

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
    log.info '    --help                              Show this message and exit.'
    log.info ''
    log.info 'Required parameters (specify in your config file):'
    log.info '    params.run_dir = RUN_DIRECTORY             Path to the sequencer output.'
    log.info '    params.output_dir OUTPUT DIRECTORY         Output directory.'
    log.info '    params.sample_sheet = SAMPLE_SHEET_PATH    Sample sheet of the format described in the README.'
    log.info '    params.p7_rows = "A B C"                   The PCR rows used - must match order of params.p5_cols.'
    log.info '    params.p5_cols = "1 2 3"                   The PCR columns used - must match order of params.p7_rows.'
    log.info '    params.demux_out = DEMUX OUTPUT DIR        Path to the demux_out folder from the bbi-dmux run.'
    log.info '    params.level = 3                           2 or 3 level sci?'
    log.info ''
    log.info 'Optional parameters (specify in your config file):'
    log.info '    params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.'
    log.info '    process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.'
    log.info '    process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted. '
    log.info '    params.samples = [sample1, sample2]        Add to only run certain samples from trimming on. Default is to run all.'
    log.info '    params.star_file = PATH/TO/FILE            File with the genome to star maps, similar to the one included with the package.'
    log.info '    params.gene_file = PATH/TO/FILE            File with the genome to gene model maps, similar to the one included with the package.'
    log.info '    params.umi_cutoff = 100                    The umi cutoff to be called a cell in matrix output.'
    log.info '    params.align_mem = 80                      Gigs of memory to use for alignment. Default is 80.'
    log.info ''
    log.info 'Issues? Contact hpliner@uw.edu'
    exit 1
}

// check required options
if (!params.run_dir || !params.output_dir || !params.sample_sheet || !params.p7_rows || !params.p5_cols) {
    exit 1, "Must include config file using -c CONFIG_FILE.config that includes output_dir, sample_sheet, run_dir, p7_rows and p5_cols"
}

star_file = file(params.star_file)
gene_file = file(params.gene_file)

process check_sample_sheet {
    cache 'lenient'
    module 'modules:java/latest:modules-init:modules-gs:python/3.6.4'


    output:
        file "*.csv" into good_sample_sheet
        file '*.log' into log_check_sample

    """

    printf "BBI bbi-sci Pipeline Log\n\n" > start.log
    printf "Nextflow version: $nextflow.version\n" >> start.log
    printf "Pipeline ID: $workflow.scriptId\n" >> start.log
    printf "Git Version: $workflow.revision, $workflow.commitId\n" >> start.log
    printf "Run started at: \$(date)\n\n" >> start.log
    printf "Command:\n$workflow.commandLine\n\n" >> start.log
    printf "***** PARAMETERS *****: \n\n" >> start.log
    printf "    params.run_dir:      $params.run_dir\n" >> start.log
    printf "    params.output_dir:   $params.output_dir\n" >> start.log
    printf "    params.sample_sheet: $params.sample_sheet\n" >> start.log
    printf "    params.p7_rows:      $params.p7_rows\n" >> start.log
    printf "    params.p5_cols:      $params.p5_cols\n" >> start.log
    printf "    params.demux_out:    $params.demux_out\n" >> start.log
    printf "    params.level:        $params.level\n" >> start.log
    printf "    params.max_cores:    $params.max_cores\n" >> start.log
    printf "    params.samples:      $params.samples\n" >> start.log
    printf "    params.star_file:    $params.star_file\n" >> start.log
    printf "    params.gene_file:    $params.gene_file\n" >> start.log
    printf "    params.umi_cutoff:   $params.umi_cutoff\n" >> start.log
    printf "    params.align_mem:    $params.align_mem\n\n" >> start.log

    printf "***** BEGIN PIPELINE *****: \n\n" >> start.log
    printf "** Start process 'check_sample_sheet' at: \$(date)\n" >> start.log
    printf "    Process versions: \$(python --version)\n" >> start.log
    printf "    Process command: check_sample_sheet.py --sample_sheet $params.sample_sheet --star_file $star_file --level $params.level\n" >> start.log

    check_sample_sheet.py --sample_sheet $params.sample_sheet --star_file $star_file --level $params.level

    printf "** End process 'check_sample_sheet' at: \$(date)\n\n" >> start.log
    """
}

good_sample_sheet.into { sample_sheet_file1; sample_sheet_file2; sample_sheet_file3 }

process trim_fastqs {
    cache 'lenient'
    memory '12G'
    module 'java/latest:modules:modules-init:modules-gs:python/2.7.3:cutadapt/1.8.3:trim_galore/0.4.1'

    input:
        file input_fastq from Channel.fromPath("${params.demux_out}/*.fastq.gz")
        file logfile from log_check_sample

    output:
        file "trim_out" into trim_output
        set file("trim_out/*.fq.gz"), val("${input_fastq.baseName - ~/.fastq/}"), file('*.log') into trimmed_fastqs mode flatten
        file input_fastq into sample_fastqs
        file 'trim_out/*trimming_report.txt' into trimgalore_results
	    stdout trim_stdout

    when:
	!((input_fastq.name - ~/-L00\d.fastq.gz/) in "Undetermined") && (!params.samples || ((input_fastq.name - ~/-L00\d.fastq.gz/) in params.samples) || ((input_fastq.name  - ~/-L00\d.fastq.gz/) in params.samples.collect{"$it".replaceAll(/\s/, ".").replaceAll(/_/, ".").replaceAll(/-/, ".").replaceAll(/\\//, ".")}))
    """
    cat ${logfile} > trim.log
    printf "** Start process 'trim_fastqs' at: \$(date)\n" >> trim.log
    printf "    Process command: trim_galore $input_fastq -a AAAAAAAA --three_prime_clip_R1 1 \
        --gzip -o ./trim_out/" >> trim.log
    mkdir trim_out
    trim_galore $input_fastq \
        -a AAAAAAAA \
        --three_prime_clip_R1 1 \
        --gzip \
        -o ./trim_out/
    cat trim_out/*trimming_report.txt >> trim.log
    printf "** End process 'trim_fastqs' at: \$(date)\n\n" >> trim.log
    cutadapt --version &> tool_2.txt
    trim_galore --version &> tool_3.txt
    """
}

if (params.max_cores < 8) {
    cores_align = params.max_cores
} else {
    cores_align = 8
}

process prep_align {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'

    input:
        file sample_sheet_file1
        set file(trimmed_fastq), val(name), file(logfile) from trimmed_fastqs

    output:
        set file(trimmed_fastq), file('info.txt'), val(name), file(logfile), stdout  into align_prepped

    """
#!/usr/bin/env python
def quick_parse(file_path):
    # a copy of only the relevant lines from easygrid.read_delim
    fh = open(file_path)
    columns = next(fh).strip().split(",")
    for line_number, line in enumerate(fh):
        entries = line.strip().split(",")
        if len(entries) != len(columns):
            raise ValueError('Length of entries: %s, not equal to length of columns %s, in file: %s, line number %s' % (len(entries), len(columns), file_path, line_number))
        entries_dict = dict(zip(columns, entries))
        yield entries_dict
lookup = {}
for rt_well in quick_parse("$sample_sheet_file1"):
    lookup[rt_well['Sample ID'].replace('(', '.').replace(')', '.').replace(' ', '.').replace('-', '.').replace('_', '.').replace('/', '.')] = rt_well['Reference Genome']

STAR_INDICES = {}
MEM = {}
with open("$star_file", 'r') as f:
    for line in f:
        items = line.strip().split()
        key, values, mem = items[0], items[1],  items[2]
        STAR_INDICES[key] = values
        MEM[key] = mem
samp = "${trimmed_fastq}".split('-')[0]
samp_name = "${trimmed_fastq}".replace('_trimmed.fq.gz', '.')
star_index = STAR_INDICES[lookup[samp]]
print(MEM[lookup[samp]])
prefix = "./align_out/" + samp_name
f = open("info.txt", 'w')
f.write(star_index + '\\n' + prefix)
f.close()
    """

}

process align_reads {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:STAR/2.5.2b'
    memory { mem.toInteger()/cores_align + " GB" }
    penv 'serial'
    cpus cores_align

    input:
        set file(input_file), file(info), val(orig_name), file(logfile), val(mem) from align_prepped

    output:
        file "align_out" into align_output
        set file("align_out/*Aligned.out.bam"), val(orig_name), file('*.log') into aligned_bams mode flatten
        file 'v_star.txt' into tool_ch4
        file 'align_out/*Log.final.out' into star_results

    """
    cat ${logfile} > align.log
    printf "** Start process 'align_reads' at: \$(date)\n" >> align.log
    printf "    Process versions: \$(STAR --version)\n" >> align.log

    mkdir align_out
    info1=`head -n 1 $info`
    info2=`head -2 $info | tail -1`

    printf "    Process command: STAR --runThreadN $cores_align --genomeDir \$info1 \
        --readFilesIn $input_file --readFilesCommand zcat --outFileNamePrefix \$info2 \
        --outSAMtype BAM Unsorted --outSAMmultNmax 2 --outSAMstrandField intronMotif" >> align.log

    STAR \
        --runThreadN $cores_align \
        --genomeDir \$info1 \
        --readFilesIn $input_file \
        --readFilesCommand zcat \
        --outFileNamePrefix \$info2 \
        --outSAMtype BAM Unsorted \
        --outSAMmultNmax 2 \
        --outSAMstrandField intronMotif

    cat align_out/*Log.final.out >> align.log

    printf "** End process 'align_reads' at: \$(date)\n\n" >> align.log

    STAR --version &> v_star.txt

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
        set file(aligned_bam), val(orig_name), file(logfile) from aligned_bams

    output:
        file "*.bam" into sorted_bams
        file 'v_samtools.txt' into tool_ch5
        file '*.log' into log_sort_and_filter

    """
    cat ${logfile} > ${orig_name}.log
    printf "** Start process 'sort_and_filter' at: \$(date)\n" >> ${orig_name}.log
    printf "    Process versions: \$(samtools --version)\n" >> ${orig_name}.log
    printf "    Process command: samtools view -bh -q 30 -F 4 "$aligned_bam" \
        | samtools sort -@ $cores_sf - > "${orig_name}.bam""

    samtools view -bh -q 30 -F 4 "$aligned_bam" \
        | samtools sort -@ $cores_sf - \
        > "${orig_name}.bam" >> ${orig_name}.log

    printf "** End process 'sort_and_filter' at: \$(date)\n\n" >> ${orig_name}.log
    samtools --version &> v_samtools.txt
    """
}

sorted_bams
    .map { file ->
        def key = file.name.toString().split(/-L[0-9]{3}/)[0]
        return tuple(key, file)
    }
    .groupTuple()
    .set { bams_to_merge }

log_sort_and_filter
    .map { file ->
        def key = file.name.toString().split(/-L[0-9]{3}/)[0]
        return tuple(key, file)
    }
    .groupTuple()
    .set { logs_to_merge }

logs_to_merge.join(bams_to_merge).set{for_merge_bams}

save_bam = {params.output_dir + "/" + it - ~/.bam/ + "/" + it}

process merge_bams {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:samtools/1.4'
    publishDir = [path: "${params.output_dir}/", saveAs: save_bam, pattern: "*.bam", mode: 'copy' ]

    input:
        set key, file(logfile), file(bam_set) from for_merge_bams

    output:
        set key, file("*.bam"), file("*.log") into sample_bams

    """
    cat ${logfile} > merge_bams.log
    printf "** Start process 'merge_bams' at: \$(date)\n" >> merge_bams.log
    printf "    Process versions: \$(samtools --version)\n" >> merge_bams.log
    printf "    Process command: samtools merge ${key}.bam $bam_set" >> merge_bams.log

    samtools merge ${key}.bam $bam_set

    printf "** End process 'merge_bams' at: \$(date)\n\n" >> merge_bams.log
    """

}

process remove_dups {
    cache 'lenient'
    memory '20 GB'
    module 'java/latest:modules:modules-init:modules-gs:samtools/1.4:bedtools/2.26.0:python/3.6.4:coreutils/8.24'

    input:
        set key, file(merged_bam), file(logfile) from sample_bams


    output:
        set key, file("*.bed"), file(merged_bam), file("*.log") into remove_dup_out
        file 'v_bedtools.txt' into tool_ch6

    """
    cat ${logfile} > remove_dups.log
    printf "** Start process 'remove_dups' at: \$(date)\n" >> remove_dups.log
    printf "    Process versions: \$(bedtools --version), \$(samtools --version), \$(python --version)\n" >> remove_dups.log
    printf "    Process command:     samtools view -h "$merged_bam" \
            | rmdup.py --bam - \
            | samtools view -bh \
            | bedtools bamtobed -i - -split \
            | sort -k1,1 -k2,2n -k3,3n -S 5G \
            > "${key}.bed"" >> remove_dups.log

    export LC_ALL=C

    samtools view -h "$merged_bam" \
            | rmdup.py --bam - \
            | samtools view -bh \
            | bedtools bamtobed -i - -split \
            | sort -k1,1 -k2,2n -k3,3n -S 5G \
            > "${key}.bed"

    printf "** End process 'remove_dups' at: \$(date)\n\n" >> remove_dups.log
    bedtools --version &> v_bedtools.txt
    """
}

remove_dup_out.into { for_prep_assign; for_umi_by_sample }

process prep_assign {
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'
    cache 'lenient'

    input:
        file sample_sheet_file2
        set key, file(sample_bed), file(merged_bam), file(logfile) from for_prep_assign

    output:
        set key, file(sample_bed), file('info.txt'), file(logfile) into assign_prepped

    """
#!/usr/bin/env python
def quick_parse(file_path):
    # a copy of only the relevant lines from easygrid.read_delim
    fh = open(file_path)
    columns = next(fh).strip().split(",")
    for line_number, line in enumerate(fh):
        entries = line.strip().split(",")
        if len(entries) != len(columns):
            raise ValueError('Length of entries: %s, not equal to length of columns %s, in file: %s, line number %s' % (len(entries), len(columns), file_path, line_number))
        entries_dict = dict(zip(columns, entries))
        yield entries_dict
lookup = {}
for rt_well in quick_parse("$sample_sheet_file2"):
    lookup[rt_well['Sample ID'].replace('(', '.').replace(')', '.').replace(' ', '.').replace('-', '.').replace('_', '.').replace('/', '.')] = rt_well['Reference Genome']
GENE_MODELS = {}
with open("$gene_file", 'r') as f:
    for line in f:
        items = line.strip().split()
        key, values = items[0], items[1]
        GENE_MODELS[key] = values
print(lookup)
samp = "${key}"
print(samp)
exon_index = GENE_MODELS[lookup[samp]] + "latest.exons.bed"
gene_index = GENE_MODELS[lookup[samp]] + "latest.genes.bed"
f = open("info.txt", 'w')
f.write(exon_index + '\\n' + gene_index + '\\n' + samp + ".txt")
f.close()
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

process assign_genes {
    cache 'lenient'
    memory '15 GB'
    module 'java/latest:modules:modules-init:modules-gs:bedtools/2.26.0:coreutils/8.24'

    input:
        set key, file(input_bed), file(info), file(logfile) from assign_prepped

    output:
        set key, file("*.txt"), file("*.log") into assign_genes_out

    """
    cat ${logfile} > assign_genes.log
    printf "** Start process 'assign_genes' at: \$(date)\n" >> assign_genes.log
    printf "    Process versions: \$(bedtools --version), \$(python --version)\n" >> assign_genes.log
    printf "    Process command:         bedtools map \
        -a "$input_bed" \
        -b \$exon_index \
        -nonamecheck -s -f 0.95 -c 7 -o distinct -delim '|' \
    | bedtools map \
        -a - -b \$gene_index \
        -nonamecheck -s -f 0.95 -c 4 -o distinct -delim '|' \
    | sort -k4,4 -k2,2n -k3,3n -S 5G\
    | datamash \
        -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8 \
    | assign-reads-to-genes.py \$gene_index \
    > "\$prefix"
    if [[ ! -s \$prefix ]]; then echo "File is empty"; exit 125; fi" >> assign_genes.log


    exon_index=`head -n 1 $info`
    gene_index=`head -2 $info | tail -1`
    prefix=`head -3 $info | tail -1`

    bedtools map \
        -a "$input_bed" \
        -b \$exon_index \
        -nonamecheck -s -f 0.95 -c 7 -o distinct -delim '|' \
    | bedtools map \
        -a - -b \$gene_index \
        -nonamecheck -s -f 0.95 -c 4 -o distinct -delim '|' \
    | sort -k4,4 -k2,2n -k3,3n -S 5G\
    | datamash \
        -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8 \
    | assign-reads-to-genes.py \$gene_index \
    > "\$prefix"
    if [[ ! -s \$prefix ]]; then echo "File is empty"; exit 125; fi

     printf "** End process 'assign_genes' at: \$(date)\n\n" >> assign_genes.log
    """

}

process umi_by_sample {
    cache 'lenient'
    memory '20 GB'
    module 'java/latest:modules:modules-init:modules-gs:samtools/1.4:coreutils/8.24'

    input:
        set key, file(input_bed), file(filtered_bam), file(logfile) from for_umi_by_sample

    output:
        set key, file("*.UMI_count.txt"), file("*.read_count.txt") into for_summarize_dup

    """
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
    samtools view "$filtered_bam" \
    | cut -d '|' -f 2 \
    | datamash -g 1 count 1 \
    | sort -k1,1 -S 5G \
    | datamash -g 1 sum 2 \
    > "${key}.read_count.txt"
    """
}

save_dup = {params.output_dir + "/" + it - ~/.duplication_rate_stats.txt/ + "/duplication_stats.txt"}

process summarize_duplication {
    cache 'lenient'
    memory '8 GB'
    publishDir = [path: "${params.output_dir}/", saveAs: save_dup, pattern: "*duplication_rate_stats.txt", mode: 'copy']

    input:
        set key, file(umi_count_file), file(read_count_file) from for_summarize_dup

    output:
        file("*duplication_rate_stats.txt") into duplication_rate_out

    """
    cat $umi_count_file \
        | join - "$read_count_file" \
        | awk 'BEGIN {{
            printf "%-18s    %10s    %10s    %8s\\n",
                "sample", "n.reads", "n.UMI", "dup.rate"
        }} {{
                printf "%-18s   %10d    %10d    %7.1f%\\n",
                    \$1, \$3, \$2, 100 * (1 - \$2/\$3);
        }}' \
        >"${key}.duplication_rate_stats.txt"


    """

}

process zip_up_duplication {
    cache 'lenient'
    publishDir = [path: "${params.output_dir}/", pattern: "all_duplication_rate.txt", mode: 'copy']

    input:
        file files from duplication_rate_out.collect()
    output:
        file "*ll_duplication_rate.txt" into all_dups

    """
    tail -1 files | grep '' *.duplication_rate_stats.txt > all_duplication_rate.txt
    """

}

assign_genes_out.into { for_umi_rollup; for_umi_by_sample_summary }

/**
make cell, gene table for each read
sort
count instances of the gene for each read
**/

process umi_rollup {
    cache 'lenient'
    memory '8 GB'
    module 'java/latest:modules:modules-init:modules-gs:coreutils/8.24'

    input:
        set key, file(gene_assignments_file), file(logfile) from for_umi_rollup

    output:
        set key, file("*.gz"), file(gene_assignments_file), file("*.log") into umi_rollup_out


    """
    cat ${logfile} > umi_rollup.log
    printf "** Start process 'umi_rollup' at: \$(date)\n" >> umi_rollup.log
    printf "    Process versions: \$(python --version)\n" >> umi_rollup.log
    printf "    Process command:  awk '\$3 == "exonic" || \$3 == "intronic" {{
            split(\$1, arr, "|")
            printf "%s|%s_%s_%s\t%s\\n", arr[2], arr[3], arr[4], arr[5], \$2
    }}' "$gene_assignments_file" \
    | sort -k1,1 -k2,2 -S 5G \
    | datamash -g 1,2 count 2 \
    | gzip > "${key}.gz"      >> umi_rollup.log

    awk '\$3 == "exonic" || \$3 == "intronic" {{
            split(\$1, arr, "|")
            printf "%s|%s_%s_%s\t%s\\n", arr[2], arr[3], arr[4], arr[5], \$2
    }}' "$gene_assignments_file" \
    | sort -k1,1 -k2,2 -S 5G \
    | datamash -g 1,2 count 2 \
    | gzip > "${key}.gz"


     printf "** End process 'umi_rollup' at: \$(date)\n\n" >> umi_rollup.log
    """
}

save_umi_per_cell = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.txt/ + "/umis_per_cell_barcode.txt"}
save_umi_per_int = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.intronic.txt/ + "/intronic_umis_per_cell_barcode.txt"}
save_plot = {params.output_dir + "/" + it - ~/.knee_plot.png/ + "/knee_plot.png"}


/**
Count intronic and total umis per cell and plot knee plot
**/

process umi_by_sample_summary {
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4:gcc/8.1.0:R/3.6.1'
    cache 'lenient'
    memory '8 GB'

    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_int, pattern: "*intronic.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_plot, pattern: "*.knee_plot.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_cell, pattern: "*barcode.txt", mode: 'copy'


    input:
        set key, file(umi_rollup), file(gene_assignments_file), file(logfile) from umi_rollup_out

    output:
        set key, file(umi_rollup), file(gene_assignments_file), file("*.log") into ubss_out
        file "*UMIs.per.cell.barcode.txt" into umis_per_cell_barcode
        file "*UMIs.per.cell.barcode.intronic.txt" into umi_per_cell_intronic
        file "*.knee_plot.png" into knee_plots
        file 'v_R.txt' into tool_ch7

    """
    cat ${logfile} > umi_by_sample_summary.log
    printf "** Start process 'umi_by_sample_summary' at: \$(date)\n" >> umi_by_sample_summary.log
    printf "    Process versions: \$(python --version), \$(R --version)\n" >> umi_by_sample_summary.log
    printf "    Process command:  tabulate_per_cell_counts.py \
        --gene_assignment_files "$gene_assignments_file" \
        --all_counts_file "${key}.UMIs.per.cell.barcode.txt" \
        --intron_counts_file "${key}.UMIs.per.cell.barcode.intronic.txt"

    knee-plot.R \
        "${key}.UMIs.per.cell.barcode.txt" \
        --knee_plot "${key}.knee_plot.png" \
        --specify_cutoff 100\
        --umi_count_threshold_file "${key}.umi_cutoff.txt""      >> umi_by_sample_summary.log


    tabulate_per_cell_counts.py \
        --gene_assignment_files "$gene_assignments_file" \
        --all_counts_file "${key}.UMIs.per.cell.barcode.txt" \
        --intron_counts_file "${key}.UMIs.per.cell.barcode.intronic.txt"

    knee-plot.R \
        "${key}.UMIs.per.cell.barcode.txt" \
        --knee_plot "${key}.knee_plot.png" \
        --specify_cutoff 100\
        --umi_count_threshold_file "${key}.umi_cutoff.txt"

    printf "** End process 'umi_rollup' at: \$(date)\n\n" >> umi_rollup.log
    R --version &> v_R.txt
    """
}

process prep_make_matrix {
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'
    cache 'lenient'

    input:
        file sample_sheet_file3
        set key, file(umi_rollup), file(gene_assignments_file), file(logfile) from ubss_out

    output:
        set key, file(umi_rollup), file(gene_assignments_file), stdout, file("*_info.txt"), file(logfile) into make_matrix_prepped

    """
#!/usr/bin/env python
def quick_parse(file_path):
    # a copy of only the relevant lines from easygrid.read_delim
    fh = open(file_path)
    columns = next(fh).strip().split(",")
    for line_number, line in enumerate(fh):
        entries = line.strip().split(",")
        if len(entries) != len(columns):
            raise ValueError('Length of entries: %s, not equal to length of columns %s, in file: %s, line number %s' % (len(entries), len(columns), file_path, line_number))
        entries_dict = dict(zip(columns, entries))
        yield entries_dict
lookup = {}
for rt_well in quick_parse("$sample_sheet_file3"):
    lookup[rt_well['Sample ID'].replace('(', '.').replace(')', '.').replace(' ', '.').replace('-', '.').replace('_', '.').replace('/', '.')] = rt_well['Reference Genome']
GENE_MODELS = {}
with open("$gene_file", 'r') as f:
    for line in f:
        items = line.strip().split()
        key, values = items[0], items[1]
        GENE_MODELS[key] = values
samp = "${key}"
exon_index = GENE_MODELS[lookup[samp]] + "latest.gene.annotations"
print(exon_index, end="")
with open("bed_info.txt", 'w') as f:
    f.write(GENE_MODELS[lookup[samp]] + "latest.genes.bed")
    """
}

save_umi = {params.output_dir + "/" + it - ~/.umi_counts.matrix/ + "/umi_counts.matrix"}
save_cell_anno = {params.output_dir + "/" + it - ~/.cell_annotations.txt/ + "/cell_annotations.txt"}
save_gene_anno = {params.output_dir + "/" + it - ~/.gene_annotations.txt/ + "/gene_annotations.txt"}


/**
sum up total assigned reads per cell and keep only those above the cutoff
make the number matrix
**/

process make_matrix {
    cache 'lenient'
    memory '15 GB'
    publishDir path: "${params.output_dir}/", saveAs: save_umi, pattern: "*umi_counts.matrix", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_anno, pattern: "*cell_annotations.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_gene_anno, pattern: "*gene_annotations.txt", mode: 'copy'

    input:
        set key, file(umi_rollup_file), file(gene_assignments_file), val(annotations_path), file(gene_bed), file(logfile) from make_matrix_prepped

    output:
        set key, file("*cell_annotations.txt"), file("*umi_counts.matrix"), file("*gene_annotations.txt"), file(gene_bed), file("*.log") into mat_output

    """
    output="${key}.cell_annotations.txt"
    UMI_PER_CELL_CUTOFF=$params.umi_cutoff
    gunzip < "$umi_rollup_file" \
    | datamash -g 1 sum 3 \
    | tr '|' '\t' \
    | awk '\$3 >= int( \$UMI_PER_CELL_CUTOFF ) {
        print \$2
    }'  - \
    | sort -k1,1 -S 5G \
    > "\$output"
    gunzip < "$umi_rollup_file" \
    | tr '|' '\t' \
    | awk '{ if (ARGIND == 1) {
                gene_idx[\$1] = FNR
            } else if (ARGIND == 2) {
                cell_idx[\$1] = FNR
            } else if (\$2 in cell_idx) {
                printf "%d\t%d\t%d\\n", gene_idx[\$3], cell_idx[\$2], \$4
            }
    }' $annotations_path "\$output" - \
    > "${key}.umi_counts.matrix"
    cat $annotations_path > "${key}.gene_annotations.txt"
    """

}

save_cds = {params.output_dir + "/" + it - ~/_cds.RDS/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_cell_qc = {params.output_dir + "/" + it - ~/_cell_qc.csv/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
process make_cds {
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4:gcc/8.1.0:R/3.6.1'
//#    publishDir path: "${params.output_dir}/", saveAs: save_cds, pattern: "*cds.RDS", mode: 'copy'
//#    publishDir path: "${params.output_dir}/", saveAs: save_cell_qc, pattern: "*cell_qc.csv", mode: 'copy'
    memory '15 GB'

    input:
        set key, file(cell_data), file(umi_matrix), file(gene_data), file(gene_bed) from mat_output

    output:
        set key, file("*for_scrub.mtx"), file("*.RDS"), file("*cell_qc.csv") into for_scrub

"""
    make_cds.R \
        "$umi_matrix"\
        "$cell_data"\
        "$gene_data"\
        "$gene_bed"\
        "$key"
"""

}

save_hist = {params.output_dir + "/" + it - ~/_scrublet_hist.png/ + "/" + it}
process run_scrublet {
    publishDir path: "${params.output_dir}/", saveAs: save_hist, pattern: "*png", mode: 'copy'

    module 'modules:java/latest:modules-init:modules-gs:python/3.6.4'
    memory '10 GB'
    input:
        set key, file(scrub_mat), file(cds), file(cell_qc) from for_scrub
    output:
        file "*.png" into scrublet_png
        set file("*scrublet_out.csv"), file(cds), file(cell_qc) into scrublet_out



"""
#!/usr/bin/env python
import scrublet as scr
import scipy.io
import numpy
import numpy.ma
from PIL import Image, ImageDraw, ImageFont
import os
counts_matrix = scipy.io.mmread("$scrub_mat").T.tocsc()
scrub = scr.Scrublet(counts_matrix)
try:
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    scrub.plot_histogram()[0].savefig("$key" + "_scrublet_hist.png")
    all_scores = numpy.vstack((doublet_scores, predicted_doublets))
    all_scores = numpy.transpose(all_scores)
    numpy.savetxt("$key" + "_scrublet_out.csv", all_scores, delimiter=",")
except (ZeroDivisionError, ValueError):
    temp = numpy.array(["NA"] * numpy.size(counts_matrix, 0))
    all_scores = numpy.vstack((temp, temp))
    all_scores = numpy.transpose(all_scores)
    filename = "$key" + "_scrublet_hist.png"
    image = Image.new(mode = "RGB", size = (250,50), color = "white")
    draw = ImageDraw.Draw(image)
    draw.text((10,10), "Scrublet failed. This is generally \\nbecause there aren't enough cells.", fill = "black")
    image.save(filename)
    numpy.savetxt("$key" + "_scrublet_out.csv", all_scores, fmt="%s", delimiter=",")
except (AttributeError):
    predicted_doublets = scrub.call_doublets(threshold=0.15)
    scrub.plot_histogram()[0].savefig("$key" + "_scrublet_hist.png")
    all_scores = numpy.vstack((doublet_scores, predicted_doublets))
    all_scores = numpy.transpose(all_scores)
    numpy.savetxt("$key" + "_scrublet_out.csv", all_scores, delimiter=",")
"""

}

process reformat_scrub {
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4:gcc/8.1.0:R/3.6.1'
    memory '10 GB'
    publishDir path: "${params.output_dir}/", saveAs: save_cds, pattern: "temp_fold/*cds.RDS", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_qc, pattern: "temp_fold/*cell_qc.csv", mode: 'copy'

    input:
        set file(scrublet_outs), file(cds), file(cell_qc) from scrublet_out
    output:
        file "temp_fold/*.RDS" into cdss
        file "temp_fold/*.csv" into cell_qcs
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
saveRDS(cds, file="temp_fold/$cds")
"""

}

process calc_cell_totals {
    module 'java/latest:modules:modules-init:modules-gs'
    memory '1 GB'

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

process exp_dash {
    module 'java/latest:modules:modules-init:modules-gs:gcc/8.1.0:R/3.6.1'
    memory '8 GB'

    publishDir path: "${params.output_dir}/", pattern: "exp_dash", mode: 'copy'


    input:
        file cds_file from cdss.collect()
        file dups from all_dups
        file counts from cell_counts
    output:
        file exp_dash

    """
    mkdir exp_dash
    mkdir exp_dash/img
    cp $baseDir/bin/bbi_icon.png exp_dash/img/
    generate_exp_dash.R \
        "$params.output_dir" --all_dups "$dups" --counts "$counts"
    """
}

/*
process get_software_versions {
    cache 'lenient'
    publishDir = [path: "${params.output_dir}/", pattern: "software_versions.txt", mode: 'copy']

    input:
        file "v_python.txt" from tool_ch1
        file "v_cutadapt.txt" from tool_ch2
        file "v_trim_galore.txt" from tool_ch3
        file "v_star.txt" from tool_ch4
        file "v_samtools.txt" from tool_ch5
        file "v_bedtools.txt" from tool_ch6
        file "v_R.txt" from tool_ch7

    output:
        file 'software_versions.txt' into software_versions

    """
    cat v_python.txt v_cutadapt.txt v_trim_galore.txt v_star.txt v_samtools.txt v_bedtools.txt v_R.txt > software_versions.txt
    """
}

process trim_log {
    cache 'lenient'
    publishDir = [path: "${params.output_dir}/", pattern: "all_trimming.txt", mode: 'copy']

    input:
        file files from trimgalore_results.collect()

    output:
        file 'all_trimming.txt' into all_trim

    """
    cat files *trimming_report.txt > all_trimming.txt
    """

}

process align_log {
    cache 'lenient'
    publishDir = [path: "${params.output_dir}/", pattern: "all_alignment.txt", mode: 'copy']

    input:
        file files from star_results.collect()

    output:
        file 'all_alignment.txt' into all_align

    """
    cat files *Log.final.out > all_alignment.txt
    """

}

process zip_up_logs {
    cache 'lenient'
    memory '8 GB'
    publishDir = [path: "${params.output_dir}/", pattern: "pipeline_logs.txt", mode: 'copy']

    input:
        file "all_trimming.txt" from all_trim
        file "all_alignment.txt" from all_align

    output:
        file 'pipeline_logs.txt' into pipeline_logs

    """
    cat all_trimming.txt all_alignment.txt > pipeline_logs.txt
    """
}

*/
workflow.onComplete {
	println ( workflow.success ? "Done! Saving output" : "Oops .. something went wrong" )
}

/** send mail
\
workflow.onComplete {
    def subject = 'indropSeq execution'
    def recipient = "${params.email}"
    def attachment = "${outputMultiQC}/multiqc_report.html"
    ['mail', '-s', subject, '-a', attachment, recipient].execute() << """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}
*/
