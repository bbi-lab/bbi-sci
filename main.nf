
// Parse input parameters
params.help = false
params.rerun = false

//print usage
if (params.help) {
    log.info ''
    log.info 'BBI 2-level sci-RNA-seq Pipeline'
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
    log.info ''
    log.info 'Optional parameters (specify in your config file):'
    log.info '    params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.'
    log.info '    process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.'
    log.info '    process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted. '
    log.info '    params.rerun = [sample1, sample2]          Add to only rerun certain samples from trimming on.'
    log.info ''
    log.info 'Issues? Contact hpliner@uw.edu'
    exit 1
}

// check required options
if (!params.run_dir || !params.output_dir || !params.sample_sheet || !params.p7_rows || !params.p5_cols) {
    exit 1, "Must include config file using -c CONFIG_FILE.config that includes output_dir, sample_sheet, run_dir, p7_rows and p5_cols"
}

process check_sample_sheet {
    module 'modules:java/latest:modules-init:modules-gs:python/3.6.4'

    input:
        val params.sample_sheet

    output:
        file "*.csv" into good_sample_sheet

    """
    check_sample_sheet.py --sample_sheet $params.sample_sheet
    """
}

sample_sheet_file = good_sample_sheet


process trim_fastqs {
    cache 'lenient'
    clusterOptions "-l mfree=8G"
    module 'java/latest:modules:modules-init:modules-gs:python/2.7.3:cutadapt/1.8.3:trim_galore/0.4.1'

    input:
        file input_fastq from Channel.fromPath("${params.demux_out}/*.fastq")
    
    output:
        file "trim_out" into trim_output
        file "trim_out/*.fq.gz" into trimmed_fastqs mode flatten
        file input_fastq into sample_fastqs
    """
    mkdir trim_out
    trim_galore $input_fastq \
        -a AAAAAAAA \
        --three_prime_clip_R1 1 \
        --no_report_file \
        --gzip \
        -o ./trim_out
        
    """
}

lane = ~/-L[0-9]{3}.fastq/

sample_fastqs
    .collectFile() { item ->
     [ "${(item.name - lane)}.fastq", item.toString() + '\n' ]
    }
    .set { fastqs_to_merge }

save_fastq = {params.output_dir + "/" + it - ~/.fastq/ + "/" + it - ~/.fastq/}


process save_sample_fastqs {
    cache = 'lenient'
    publishDir = [path: "${params.output_dir}/", saveAs: save_fq, pattern: "*.fq.gz", mode: 'copy' ]

    input:
       file fastq_set from fastqs_to_merge

    output:
       file "*.fq.gz" into fqs

    """
    cat $fastq_set | gzip > ${fastq_set}.fq.gz

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
        file sample_sheet_file
        file trimmed_fastq from trimmed_fastqs

    output:
        set file(trimmed_fastq), file('info.txt') into align_prepped

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
for rt_well in quick_parse("$sample_sheet_file"):
    lookup[rt_well['Sample ID'].replace('-', '.').replace('_', '.').replace(' ', '.')] = rt_well['Reference Genome']
    

STAR_INDICES = {
    'Human': '/net/trapnell/vol1/jspacker/STAR/human/GRCh38-primary-assembly',
    'Mouse': '/net/trapnell/vol1/jspacker/STAR/mouse/GRCm38-primary-assembly',
    'Barnyard': '/net/trapnell/vol1/jspacker/STAR/human-and-mouse/GRCh38-GRCm38-primary-assembly',
    'Barn': '/net/trapnell/vol1/jspacker/STAR/human-and-mouse/GRCh38-GRCm38-primary-assembly',
    'Celegans': '/net/trapnell/vol1/jspacker/STAR/c-elegans/WS260',
    'Camelina': None,
    'Arabidopsis': None,
    'Maize': None,
    'Rat': '/net/bbi/vol1/data/genomes/rat_star',
    'Macaque': '/net/bbi/vol1/data/genomes/macaque_star',
    'Zebrafish': '/net/bbi/vol1/data/genomes/zebrafish_star',
    'Drosophila': '/net/bbi/vol1/data/genomes/drosophila_star'
}
samp = "${trimmed_fastq}".split('-')[0]
samp_name = "${trimmed_fastq}".replace('_trimmed.fq.gz', '.')
star_index = STAR_INDICES[lookup[samp]]
prefix = "./align_out/" + samp_name
f = open("info.txt", 'w')
f.write(star_index + '\\n' + prefix)
f.close()
    """

}


process align_reads {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:STAR/2.5.2b'
    clusterOptions "-l mfree=25G -pe serial $cores_align"
    

    input:
        set file(input_file), file(info) from align_prepped

    output:
        file "align_out" into align_output
        file "align_out/*Aligned.out.bam" into aligned_bams mode flatten

    """
    mkdir align_out
    info1=`head -n 1 $info`
    info2=`head -2 $info | tail -1`
    STAR \
        --runThreadN $cores_align \
        --genomeDir \$info1 \
        --readFilesIn $input_file \
        --readFilesCommand zcat \
        --outFileNamePrefix \$info2 \
        --outSAMtype BAM Unsorted \
        --outSAMmultNmax 2 \
        --outSAMstrandField intronMotif

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
    clusterOptions "-l mfree=1G -pe serial $cores_sf"

    input:
        file aligned_bam from aligned_bams

    output:
        file "*.bam" into sorted_bams



    """
    firstString="$aligned_bam"
    secondString="sorted"

    output_bam=`echo "\${firstString/Aligned/\$secondString}"`

    samtools view -bh -q 30 -F 4 $aligned_bam \
        | samtools sort -@ $cores_sf - \
        > \$output_bam
    """


}

lane = ~/-L[0-9]{3}.sorted.out.bam/

sorted_bams
    .collectFile() { item ->
     [ "${(item.name - lane)}.txt", item.toString() + '\n' ]
    }
    .set { Bams_to_merge }



save_bam = {params.output_dir + "/" + it - ~/.txt.bam/ + "/" + it - ~/.txt/}

process merge_bams {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:samtools/1.4'
    publishDir = [path: "${params.output_dir}/", saveAs: save_bam, pattern: "*.bam", mode: 'copy' ]

    input:
        file bam_set from Bams_to_merge

    output:
        file "*.bam" into sample_bams

    """
    samtools merge -b $bam_set ${bam_set}.bam
    
    """

}

process remove_dups {
    cache 'lenient'
    clusterOptions "-l mfree=4G"
    module 'java/latest:modules:modules-init:modules-gs:samtools/1.4:bedtools/2.26.0:python/3.6.4'

    input:
        file merged_bam from sample_bams

    output:
        set file("*.bed"), file(merged_bam) into remove_dup_out

    """
    export LC_ALL=C
    samtools view -h $merged_bam \
            | rmdup.py --bam - \
            | samtools view -bh \
            | bedtools bamtobed -i - -split \
            | sort -k1,1 -k2,2n -k3,3n -S 3G \
            > ${merged_bam}.bed
    """
}

remove_dup_out.into { for_prep_assign; for_umi_by_sample }

process prep_assign {
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'
    cache 'lenient'
    
    input:
        file sample_sheet_file
        set file(sample_bed), file(merged_bam) from for_prep_assign

    output:
        set file(sample_bed), file('info.txt') into assign_prepped

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
for rt_well in quick_parse("$sample_sheet_file"):
    lookup[rt_well['Sample ID'].replace('-', '.').replace('_', '.').replace(' ', '.')] = rt_well['Reference Genome']

GENE_MODELS = {
    'Human': '/net/trapnell/vol1/jspacker/gene-models/human/',
    'Mouse': '/net/trapnell/vol1/jspacker/gene-models/mouse/',
    'Barnyard': '/net/trapnell/vol1/jspacker/gene-models/human-and-mouse/',
    'Barn': '/net/trapnell/vol1/jspacker/gene-models/human-and-mouse/',
    'Celegans': '/net/trapnell/vol1/jspacker/gene-models/c-elegans/',
    'Camelina': '/net/gs/vol1/home/ajh24/common_data/bbi/gene_models/camelina/',
    'Arabidopsis': '/net/gs/vol1/home/ajh24/common_data/bbi/gene_models/arabidopsis/',
    'Maize': '/net/gs/vol1/home/ajh24/common_data/bbi/gene_models/maize/',
    'Rat': '/net/bbi/vol1/data/genomes/rat/',
    'Macaque': '/net/bbi/vol1/data/genomes/macaque/',
    'Zebrafish': '/net/bbi/vol1/data/genomes/zebrafish/',
    'Drosophila': '/net/bbi/vol1/data/genomes/drosophila/',
    'Temp': 'NA'
}

samp = "${sample_bed}".replace(".txt.bam.bed", "")
samp_name = "${sample_bed}".replace('.txt.bam.bed', '.')
exon_index = GENE_MODELS[lookup[samp]] + "latest.exons.bed"
gene_index = GENE_MODELS[lookup[samp]] + "latest.genes.bed"
f = open("info.txt", 'w')
f.write(exon_index + '\\n' + gene_index + '\\n' + samp_name + ".txt")
f.close()
    """
}

process assign_genes {
    cache 'lenient'
    clusterOptions "-l mfree=6G"
    module 'java/latest:modules:modules-init:modules-gs:bedtools/2.26.0'

    input:
        set file(input_bed), file(info) from assign_prepped

    output:
        file "*.txt" into assign_genes_out

    """
    exon_index=`head -n 1 $info`
    gene_index=`head -2 $info | tail -1`
    prefix=`head -3 $info | tail -1`
    bedtools map \
        -a $input_bed \
        -b \$exon_index \
        -nonamecheck -s -f 0.95 -c 7 -o distinct -delim '|' \
    | bedtools map \
        -a - -b \$gene_index \
        -nonamecheck -s -f 0.95 -c 4 -o distinct -delim '|' \
    | sort -k4,4 -k2,2n -k3,3n -S 3G \
    | datamash \
        -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8 \
    | assign-reads-to-genes.py \$gene_index \
    > \$prefix

    """

}



process umi_by_sample {
    cache 'lenient'
    clusterOptions "-l mfree=8G"
    module 'java/latest:modules:modules-init:modules-gs:samtools/1.4'

    input:
        set file(input_bed), file(filtered_bam) from for_umi_by_sample

    output:
        set file("*.UMI_count.txt"), file("*.read_count.txt") into for_summarize_dup

    """
    awk '{{ split(\$4, arr, "|")
            if (!seen[arr[1]]) {{
                seen[arr[1]] = 1; count[arr[2]]++;
            }}
            }} END {{
                for (sample in count) {{
                print sample "\\t" count[sample]
                }}
            }}' $input_bed \
    | sort -k1,1 \
    >${input_bed}.UMI_count.txt

    samtools view $filtered_bam \
    | cut -d '|' -f 2 \
    | datamash -g 1 count 1 \
    | sort -k1,1 -S 2G \
    | datamash -g 1 sum 2 \
    > ${input_bed}.read_count.txt
    """

}

save_dup = {params.output_dir + "/" + it - ~/.txt.bam.bed.UMI_count.txt.duplication_rate_stats.txt/ + "/duplication_stats.txt"}

process summarize_duplication {
    cache 'lenient'
    clusterOptions "-l mfree=8G"
    publishDir = [path: "${params.output_dir}/", saveAs: save_dup, pattern: "*duplication_rate_stats.txt", mode: 'copy']

    input:
        set file(umi_count_file), file(read_count_file) from for_summarize_dup

    output:
        file "*duplication_rate_stats.txt" into duplication_rate_out

    """
    cat $umi_count_file \
        | join - $read_count_file \
        | awk 'BEGIN {{ 
            printf "%-18s    %10s    %10s    %8s\\n",
                "sample", "n.reads", "n.UMI", "dup.rate"
        }} {{ 
                printf "%-18s   %10d    %10d    %7.1f%\\n",
                    \$1, \$3, \$2, 100 * (1 - \$2/\$3);
        }}' \
        >${umi_count_file}.duplication_rate_stats.txt
    
    
    """

}


assign_genes_out.into { for_umi_rollup; for_umi_by_sample_summary }


process umi_rollup {
    cache 'lenient'
    clusterOptions "-l mfree=4G"

    input:
        file gene_assignments_file from for_umi_rollup

    output:
        set file("*.gz"), file(gene_assignments_file) into umi_rollup_out


    """
    awk '\$3 == "exonic" || \$3 == "intronic" {{
            split(\$1, arr, "|")
            printf "%s|%s_%s_%s\t%s\\n", arr[2], arr[3], arr[4], arr[5], \$2
    }}' $gene_assignments_file \
    | sort -k1,1 -k2,2 -S 2G \
    | datamash -g 1,2 count 2 \
    | gzip > ${gene_assignments_file}.gz
    """

}

/*


save_bam = {params.output_dir + "/" + it - ~/.txt.bam/ + "/" + it - ~/.txt/}

    publishDir = [path: "${params.output_dir}/", saveAs: save_bam, pattern: "*.bam", mode: 'copy']


*/

save_umi_per_cell = {params.output_dir + "/" + it - ~/..txt.UMIs.per.cell.barcode.txt/ + "/umis_per_cell_barcode.txt"}
save_umi_per_int = {params.output_dir + "/" + it - ~/..txt.UMIs.per.cell.barcode.intronic.txt/ + "/intronic_umis_per_cell_barcode.txt"}
save_plot = {params.output_dir + "/" + it - ~/..txt.knee_plot.png/ + "/knee_plot.png"}

process umi_by_sample_summary {
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4:gcc/8.1.0:R/3.5.2'
    cache 'lenient'
    clusterOptions "-l mfree=8G"

    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_int, pattern: "*intronic.txt", mode: 'copy' 
    publishDir path: "${params.output_dir}/", saveAs: save_plot, pattern: "*.knee_plot.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_cell, pattern: "*barcode.txt", mode: 'copy'  
  

    input:
        set file(umi_rollup), file(gene_assignments_file) from umi_rollup_out        

    output:
        set file(umi_rollup), file(gene_assignments_file), file("*umi_cutoff.txt") into ubss_out
        file "*UMIs.per.cell.barcode.txt" into umis_per_cell_barcode
        file "*UMIs.per.cell.barcode.intronic.txt" into umi_per_cell_intronic
        file "*.knee_plot.png" into knee_plots

    """
    tabulate_per_cell_counts.py \
        --gene_assignment_files $gene_assignments_file \
        --all_counts_file ${gene_assignments_file}.UMIs.per.cell.barcode.txt \
        --intron_counts_file ${gene_assignments_file}.UMIs.per.cell.barcode.intronic.txt
    
    knee-plot.R \
        ${gene_assignments_file}.UMIs.per.cell.barcode.txt \
        --knee_plot ${gene_assignments_file}.knee_plot.png \
        --umi_count_threshold_file ${gene_assignments_file}.umi_cutoff.txt

    """


}

process prep_make_matrix {
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'
    cache 'lenient'
    
    input:
        file sample_sheet_file
        set file(umi_rollup), file(gene_assignments_file), file(umi_cutoff) from ubss_out

    output:
        set file(umi_rollup), file(gene_assignments_file), file(umi_cutoff), stdout into make_matrix_prepped

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
for rt_well in quick_parse("$sample_sheet_file"):
    lookup[rt_well['Sample ID'].replace('-', '.').replace('_', '.').replace(' ', '.')] = rt_well['Reference Genome']

GENE_MODELS = {
    'Human': '/net/trapnell/vol1/jspacker/gene-models/human/',
    'Mouse': '/net/trapnell/vol1/jspacker/gene-models/mouse/',
    'Barnyard': '/net/trapnell/vol1/jspacker/gene-models/human-and-mouse/',
    'Barn': '/net/trapnell/vol1/jspacker/gene-models/human-and-mouse/',
    'Celegans': '/net/trapnell/vol1/jspacker/gene-models/c-elegans/',
    'Camelina': '/net/gs/vol1/home/ajh24/common_data/bbi/gene_models/camelina/',
    'Arabidopsis': '/net/gs/vol1/home/ajh24/common_data/bbi/gene_models/arabidopsis/',
    'Maize': '/net/gs/vol1/home/ajh24/common_data/bbi/gene_models/maize/',
    'Rat': '/net/bbi/vol1/data/genomes/rat/',
    'Macaque': '/net/bbi/vol1/data/genomes/macaque/',
    'Zebrafish': '/net/bbi/vol1/data/genomes/zebrafish/',
    'Drosophila': '/net/bbi/vol1/data/genomes/drosophila/'
}

samp = "${gene_assignments_file}".replace("..txt", "")
exon_index = GENE_MODELS[lookup[samp]] + "latest.gene.annotations"
print(exon_index, end="")
    """
}

save_umi = {params.output_dir + "/" + it - ~/..txt.umi_counts.matrix/ + "/umi_counts.matrix"}
save_cell_anno = {params.output_dir + "/" + it - ~/..txt.cell_annotations.txt/ + "/cell_annotations.txt"}
save_gene_anno = {params.output_dir + "/" + it - ~/..txt.gene_annotations.txt/ + "/gene_annotations.txt"}

process make_matrix {
    cache 'lenient'
    clusterOptions "-l mfree=4G"
    publishDir path: "${params.output_dir}/", saveAs: save_umi, pattern: "*umi_counts.matrix", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_anno, pattern: "*cell_annotations.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_gene_anno, pattern: "*gene_annotations.txt", mode: 'copy'

    input:
        set file(umi_rollup_file), file(gene_assignments_file), file(umi_cutoff_file), val(annotations_path) from make_matrix_prepped

    output:
        set file("*cell_annotations.txt"), file("*umi_counts.matrix"), file("*gene_annotations.txt") into mat_output

    """
    output=${gene_assignments_file}.cell_annotations.txt
    touch samples_to_exclude_file
    UMI_PER_CELL_CUTOFF=\$(cat $umi_cutoff_file)
    gunzip < $umi_rollup_file \
    | datamash -g 1 sum 3 \
    | tr '|' '\t' \
    | awk -v CUTOFF=\$UMI_PER_CELL_CUTOFF 'ARGIND == 1 {{
        exclude[\$1] = 1
    }} \$3 >= int( CUTOFF ) {{
        print \$2
    }}' samples_to_exclude_file - \
    | sort -k1,1 -S 4G \
    > \$output
    gunzip < $umi_rollup_file \
    | tr '|' '\t' \
    | awk '{{ if (ARGIND == 1) {{
                gene_idx[\$1] = FNR
            }} else if (ARGIND == 2) {{ 
                cell_idx[\$1] = FNR
            }} else if (\$2 in cell_idx) {{
                printf "%d\t%d\t%d\\n", gene_idx[\$3], cell_idx[\$2], \$4
            }} 
    }}' $annotations_path \$output - \
    > ${gene_assignments_file}.umi_counts.matrix

    cat $annotations_path > ${gene_assignments_file}.gene_annotations.txt
    
    rm samples_to_exclude_file
    """

}

save_cds = {params.output_dir + "/" + it - ~/_cds.RDS/ + "/" + it}

process make_cds {
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4:gcc/8.1.0:R/3.5.2'
    publishDir path: "${params.output_dir}/", saveAs: save_cds, pattern: "*cds.RDS", mode: 'copy'

    input:
        set file(cell_data), file(umi_matrix), file(gene_data) from mat_output

    output:
        file "*.RDS"


"""
    make_cds.R \
        $umi_matrix\
        $cell_data\
        $gene_data

"""

}


workflow.onComplete { 
	println ( workflow.success ? "Done! Saving output" : "Oops .. something went wrong" )
}

/*
process summarize_alignments {


}



* send mail
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

        force_symlink(duplication_stats_file, joindir(FINAL_OUTPUT, 'duplicaton_stats.txt'))
        force_symlink(umis_per_cell_barcode_file, joindir(FINAL_OUTPUT, 'umis_per_cell_barcode.txt'))
        force_symlink(intronic_umis_per_cell_barcode_file, joindir(FINAL_OUTPUT, 'intronic_umis_per_cell_barcode.txt'))
        force_symlink(knee_plot_file, joindir(FINAL_OUTPUT, 'knee_plot.png'))
        force_symlink(region_stats_output, joindir(FINAL_OUTPUT, 'region_stats.txt'))
        force_symlink(alignment_stats_output, joindir(FINAL_OUTPUT, 'alignment_stats.txt'))
*/
