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
    cpus 4
    input:
        path(fastq)
        path(reference_genome)
        val(minimap2_addition_para)
        val(samplename)
    output:
        path "${samplename}.ont.minimap2.sort.bam"
        path "${samplename}.map_stat.xls"
    shell:
    """
        minimap2 -t 8 ${minimap2_addition_para} ${reference_genome} ${fastq} | samtools sort -@ 8 -o ${samplename}.ont.minimap2.sort.bam > ${samplename}.minimap2.log_o_file 2> ${samplename}.minimap2.log_e_file
        samtools index -@ 8 ${samplename}.ont.minimap2.sort.bam
        samtools stats -@ 8 ${samplename}.ont.minimap2.sort.bam > ${samplename}.ont.minimap2.sort.bam.stat.txt
        echo -e "Sample_ID\tPass_Reads\tMapped_Reads\tMapped_Reads_Rate(%)\tPass_Bases\tMapped_Bases\tMapped_Bases_Rate(%)\tDepth(X)" > ${samplename}.map_stat.xls
    """
}

process grandSTR {
    // concatenate fastq and fastq.gz in a dir

    label "STRdetection"
    cpus 4
    input:
        path(bam_file)
        path(reference_genome)
        val(minimap2_addition_para)
        val(samplename)
        path(target_STR_region)

    output:
        path "${samplename}.str.bam"
        path "${samplename}.STR_infos.modify"
    shell:
    """
        /home/software/GrandSTR_v1.2.9/GrandSTR ${target_STR_region} ${samplename} -rf ${reference_genome} -bf ${bam_file} -rt ont -em 0 -mi 0.6
        awk '{if(\$6 != 0){print \$0}else{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t.\\t.\\t.\\t0,0\\t0\\tNA" }}' ${samplename}.STR_infos > ${samplename}.STR_infos.modify
        #get str bam
        awk -F ',' '{print \$2 "\\t" \$3-1000 "\\t" \$4+1000 "\\t" \$1}' ${target_STR_region} > Target_STR_region_refine.V3.pa_updown.bed
        samtools view -bh -L Target_STR_region_refine.V3.pa_updown.bed ${samplename}.ont.minimap2.sort.bam > ${samplename}.str.bam
        samtools index ${samplename}.str.bam
	
    """
}




process getParams {
    label "wfSTR"
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

process output {
    // publish inputs to output directory
    label "wfSTR"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files."
    """
}



workflow {
    getParams()

    (bam_file, bam_stat_file) = minimap2(params.fastq, params.reference_genome, params.minimap2_addition_para, params.samplename)

    grandSTR(bam_file, params.reference_genome, params.minimap2_addition_para, params.samplename, params.target_STR_region)

    //pipeline(samples)
    //output(pipeline.out.results)

}
