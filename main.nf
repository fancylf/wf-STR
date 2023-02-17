#!/usr/bin/env nextflow

// Developer notes
// 
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { start_ping; end_ping } from './lib/ping'

log.info """\
C A L L I N G S  -  N F    v 2.1 
================================
out_dir : $params.out_dir
fastq : $params.fastq
reference_genome : $params.reference_genome 
target_STR_region : $params.target_STR_region
minimap2_addition_para : $params.minimap2_addition_para
samplename : $params.samplename
report_name : $params.report_name
disable_ping : $params.disable_ping
"""

process minimap2 {
    // concatenate fastq and fastq.gz in a dir

    label "alignment"
    
    cpus 48
    input:
        path(fastq)
        path(reference_genome)
        val(minimap2_addition_para)
        val(samplename)
    output:
        path "${samplename}.ont.minimap2.sort.bam"
    shell:
    """
        minimap2 -t 48 ${minimap2_addition_para} ${reference_genome} ${fastq} | samtools sort -@ 8 -m 8G -o ${samplename}.ont.minimap2.sort.bam > ${samplename}.minimap2.log_o_file 2> ${samplename}.minimap2.log_e_file
    """
}

process grandSTR {
    // concatenate fastq and fastq.gz in a dir

    label "STRdetection"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    cpus 32
    input:
        path(bam_file)
        path(reference_genome)
        val(minimap2_addition_para)
        val(samplename)
        path(target_STR_region)
        
    output:
        path "${bam_file}.stat.txt"
        path "${samplename}.str.bam"
        path "${samplename}.str.bam.bai"
        path "${samplename}.STR_infos"
    shell:
    """
        samtools index -@ 32 ${bam_file}
        samtools stats -@ 32 ${bam_file} > ${bam_file}.stat.txt
        /home/software/GrandSTR_v1.2.9/GrandSTR ${target_STR_region} ${samplename} -rf ${reference_genome} -bf ${bam_file} -rt ont -em 0 -mi 0.6
        awk '{if(\$6 != 0){print \$0}else{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t.\\t.\\t.\\t0,0\\t0\\tNA" }}' ${samplename}.STR_infos > ${samplename}.STR_infos.modify
        #get str bam
        samtools view -H ${samplename}.ont.minimap2.sort.bam >  ${samplename}.str.sam
        awk -F ',' '{print "samtools view ${samplename}.ont.minimap2.sort.bam "\$2":" \$3-1000 "-" \$4+1000}' |sh |samtools view -bS - >  ${samplename}.str.bam
        samtools index -@ 8 ${samplename}.str.bam
	
    """
        //awk -F ',' '{print \$2 "\\t" \$3-1000 "\\t" \$4+1000 "\\t" \$1}' ${target_STR_region} > Target_STR_region_refine.V3.pa_updown.bed
        //samtools view -bh -L Target_STR_region_refine.V3.pa_updown.bed ${samplename}.ont.minimap2.sort.bam > ${samplename}.str.bam
}

process getParams {
    label "wfSTR"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process getVersions_alignment {
    label "alignment"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    minimap2 --version |sed 's/^/minimap2,/' >> versions.txt
    samtools --version |grep samtools | sed 's/ /,/' >> versions.txt
    """
}

process getVersions_grandSTR {
    label "STRdetection"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    cpus 1
    input:
        path(version_file)    
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> ${version_file}
    /home/software/GrandSTR_v1.2.9/GrandSTR --version |grep GrandSTR |sed 's/ version//' |sed 's/ //' >> ${version_file}
    """
}

process annoSTR {
    label "wfSTR"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    cpus 1
    input:
        path(target_str_config)
        val(str_result)

    output:
        path "STR_infos.annot.xls"
    script:
    """
    python $baseDir/bin/add_str_annotation.py --info ${str_result} --pa ${target_str_config} --outfile STR_infos.annot.xls
    """
}

process makeReport {
    label "wfSTR"
    //label "STRdetection"
    //label "alignment"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    cpus 1
    input:
        val(samplename)
        val(params_report)
        val(all_version_report)
        val(str_result_report)
        val(bam_stat_report)

    output:
        path "${samplename}.map_stat.xls"
        path "STR_report.html"
    script:
    """
    python $baseDir/bin/map_stat_from_bamstat.py --sample ${samplename} --bamstat ${bam_stat_report}
    python $baseDir/bin/GrandSTR_report.py --mapstat ${samplename}.map_stat.xls --param  ${params_report} --software  ${all_version_report} --info ${str_result_report} --template $baseDir/bin/str_report_template.html --outdir "./"
    """
}



process configure_jbrowse {
    label "wfSTR"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        val(reference)
        val(alignments)
    output:
        path("jbrowse.json")
    script:

   """
    python $baseDir/bin/make_jbrowse_json.py --genome_fa ${reference} --location_chr 1 --location_start 57801669 --location_end 57809669 --bam ${alignments} --template $baseDir/bin/jbowse_template.json
    """
}
    // publish inputs to output directory

//process output {
//    label "wfSTR"
//    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    
//    input:
//        val(params_file)
//        val(all_version)
//        val(str_bam_file)
//        val(str_bam_index)
//        val(str_result)
//        val(str_report)
//        val(jbrowse_result)
//        val(out_dir)
        
//    output:
//        path 'result.list'
//        path(params_file)
//        path(all_version)
//        path(str_bam_file)
//        path(str_bam_index)
//        path(str_result)
//        path(str_report)
//        path(jbrowse_result)       
//    """
//        echo ${params_file} ${all_version} ${str_bam_file} ${str_bam_index} ${str_result} ${str_report} ${jbrowse_result} ${out_dir} > result.list
 //   """
//}

workflow {
    //start_ping()
    params_file = getParams()
    bam_file = minimap2(params.fastq, params.reference_genome, params.minimap2_addition_para, params.samplename)
    (bam_stat_file, str_bam_file, str_bam_index, str_result) = grandSTR(bam_file, params.reference_genome, params.minimap2_addition_para, params.samplename, params.target_STR_region)
    str_annot_result = annoSTR(params.target_STR_region, str_result)
    aligner_version = getVersions_alignment()
    all_version =  getVersions_grandSTR(aligner_version)
    str_report = makeReport(params.samplename, params_file, all_version, str_annot_result, bam_stat_file)
    jbrowse_result = configure_jbrowse(params.reference_genome, str_bam_file)
    //output(params_file, all_version, str_bam_file, str_bam_index, str_result, str_report, jbrowse_result, params.out_dir)
    
    //output(params_file, all_version, str_bam_file, str_bam_index, str_result, str_report)
    //end_ping(pipeline.out.telemetry)
}

