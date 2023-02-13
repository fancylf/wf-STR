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



process minimap2 {
    // concatenate fastq and fastq.gz in a dir

    label "aligment"
    cpus 8
    input:
        tuple path(fastq), path(reference_genome), val(minimap2_addition_para), val(samplename)
    output:
        tuple path "${samplename}.ont.minimap2.sort.bam", path "${samplename}.map_stat.xls"
    shell:
    """
        minimap2 -t 8 ${minimap2_addition_para} ${reference_genome} ${fastq} | samtools sort -@ 8 -o ${samplename}.ont.minimap2.sort.bam > ${samplename}.minimap2.log_o_file 2> ${samplename}.minimap2.log_e_file
        samtools index -@ 8 ${samplename}.ont.minimap2.sort.bam
        samtools stats -@ 8 ${samplename}.ont.minimap2.sort.bam > ${samplename}.ont.minimap2.sort.bam.stat.txt
		echo -e "Sample_ID\tPass_Reads\tMapped_Reads\tMapped_Reads_Rate(%)\tPass_Bases\tMapped_Bases\tMapped_Bases_Rate(%)\tDepth(X)" > ${samplename}.map_stat.xls
		Pass_Reads=`grep ^SN ${samplename}.ont.minimap2.sort.bam.stat.txt |cut -f 2-|sed -n '1,1p'|awk '{print \$4}'`
		Mapped_Reads=`grep ^SN ${samplename}.ont.minimap2.sort.bam.stat.txt |cut -f 2-|sed -n '7,7p'|awk '{print \$3}'`
		Mapped_Reads_Rate=`awk 'BEGIN{printf "%.2f\n",('\$Mapped_Reads'/'\$Pass_Reads')*100}'`
		Pass_Bases=`grep ^SN ${samplename}.ont.minimap2.sort.bam.stat.txt |cut -f 2-|sed -n '16,16p'|awk '{print \$3}'`
		Mapped_Bases=`grep ^SN ${samplename}.ont.minimap2.sort.bam.stat.txt |cut -f 2-|sed -n '19,19p'|awk '{print \$3}'`
		Mapped_Bases_Rate=`awk 'BEGIN{printf "%.2f\n",('\$Mapped_Bases'/'\$Pass_Bases')*100}'`
		Depth=`awk 'BEGIN{printf "%.2f\n",'\$Mapped_Bases'/2858658097}'`
		echo -e "${samplename}\t\$Pass_Reads\t\$Mapped_Reads\t\$Mapped_Reads_Rate\t\$Pass_Bases\t\$Mapped_Bases\t\$Mapped_Bases_Rate\t\$Depth" >> ${sample}.map_stat.xls
    """
}

process grandSTR {
    // concatenate fastq and fastq.gz in a dir

    label "STRdetection"
    cpus 8
    input:
        tuple path(fastq), path(reference_genome), val(minimap2_addition_para), val(samplename),path(target_STR_region)

    output:
        path "${samplename}.str.bam", path "${samplename}.STR_infos.modify"
    shell:
    """
    /home/software/GrandSTR_v1.2.9/GrandSTR ${target_STR_region} ${samplename} -rf $reference_genome -bf ${samplename}.ont.minimap2.sort.bam -rt ont -em 0 -mi 0.6
	awk '{if($6 != 0){print $0}else{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t.\t.\t.\t0,0\t0\tNA" }}' ${samplename}.STR_infos > ${samplename}.STR_infos.modify
	#get str bam
	awk -F ',' '{print $2 "\t" $3-1000 "\t" $4+1000 "\t" $1}' ${target_STR_region} > Target_STR_region_refine.V3.pa_updown.bed
	samtools view -bh -L Target_STR_region_refine.V3.pa_updown.bed ${bam} > ${sample}.str.bam
	samtools index ${samplename}.str.bam
	
    """
}




process getParams {
    label "wf-STR"
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
    label "wf-STR"
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
    minimap2()
    grandSTR()

    //pipeline(samples)
    //output(pipeline.out.results)

}
