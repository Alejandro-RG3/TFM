nextflow.enable.dsl=2

params.fastqs     = "../fastq/*_R{1,2}_*.fastq.gz" 
params.vcfs       = "../vcf_original/*.vcf.gz"
params.ref        = "../referencia/Homo_sapiens_assembly38.fasta"
params.clinvar    = "../bdd/clinvar.vcf.gz"
params.dbnsfp     = "../bdd/dbNSFP5.3.1a_grch38.gz"
params.snpeff_db  = "hg38"
params.outdir     = "resultados"


// Procesa FASTQ a VCF
process FASTQ_TO_VCF {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads)
    val ref

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

    script:
    // reads[0] es R1, y reads[1] es R2. BWA los lee a la vez.
    """
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" ${ref} ${reads[0]} ${reads[1]} | \\
    samtools sort -@ ${task.cpus} -o ${sample_id}.bam -
    samtools index -@ ${task.cpus} ${sample_id}.bam

    bcftools mpileup -a FORMAT/DP -Ou -f ${ref} ${sample_id}.bam | \\
    bcftools call -a GQ -mv -Oz -o ${sample_id}.vcf.gz
    """
}

// Filtros de calidad y normalizacion
process QC_VCF {
    tag "${sample_id}" 
    
    input:
    tuple val(sample_id), path(vcf)
    path ref

    output:
    path "${sample_id}_limpio.vcf.gz", emit: vcf_clean
    path "${sample_id}_limpio.vcf.gz.tbi", emit: vcf_idx

    script:
    """
    bcftools annotate -x INFO/ANN ${vcf} -Ou | \
    bcftools norm -m -any -f ${ref} -Ou | \
    bcftools view -f "PASS,." -i 'FORMAT/DP>=20 & FORMAT/GQ>=30' -Oz -o ${sample_id}_limpio.vcf.gz

    bcftools index -t ${sample_id}_limpio.vcf.gz
    """
}

// Junta todo en un solo VCF
process MERGE_COHORTE {
    publishDir "${params.outdir}/vcf_cohorte", mode: 'copy'

    input:
    path vcfs
    path tbis

    output:
    path "cohorte.vcf.gz", emit: vcf_merged
    path "cohorte.vcf.gz.tbi", emit: vcf_idx
    path "chroms.txt", emit: chr_list

    script:
    """
    bcftools merge -m all -O z -o cohorte.vcf.gz ${vcfs}
    bcftools index -t cohorte.vcf.gz

    # Extrae la lista de cromosomas ignorando algunos
    bcftools index -s cohorte.vcf.gz | cut -f1 | grep -v -i -E "decoy|hla|un|random|alt|chrm|mt" > chroms.txt
    """
}

// Anotacion dividida por cromosoma
process ANNOTATE_CHR {
    tag "${chr}"

    input:
    val chr
    path vcf
    path idx
    path clinvar
    path clinvar_tbi
    path dbnsfp
    path dbnsfp_tbi

    output:
    path "${chr}_anotado.vcf.gz"

    script:
    """
    set -e -o pipefail
    
    # Anotar variantes con clinvar, dbnsfp, gnomad ypredictores in silico
    bcftools view -r ${chr} ${vcf} | \
    snpEff -Xmx${task.memory.toGiga()}g -v ${params.snpeff_db} - | \
    SnpSift filter "( (ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'MODIFIER') )" | \
    SnpSift -Xmx4g annotate -name CLINVAR_ ${clinvar} - | \
    SnpSift -Xmx${task.memory.toGiga()}g dbnsfp -v -db ${dbnsfp} \
      -f CADD_phred,REVEL_score,gnomAD4.1_joint_AF,gnomAD4.1_joint_NFE_AF - | \
    bgzip -c > ${chr}_anotado.vcf.gz
    """
}

// Concatenar cromosomas y exportar tablas finales
process EXPORT_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path vcfs_anotados

    output:
    path "cohorte_final.vcf.gz"
    path "matriz_anotaciones.tsv"
    path "matriz_genotipos.tsv"

    script:
    """
    # Filtrar para ignorar los decoys, HLAs y cromosomas raros y M
    ls *_anotado.vcf.gz | grep -v -E "decoy|HLA|Un|random|alt|chrM|MT" | sort -V > lista_vcfs.txt
    bcftools concat -f lista_vcfs.txt -O z -o cohorte_final.vcf.gz
    tabix -p vcf cohorte_final.vcf.gz

    SnpSift -Xmx${task.memory.toGiga()}g extractFields -s "," cohorte_final.vcf.gz \
        CHROM POS REF ALT \
        "ANN[*].GENE" "ANN[*].IMPACT" "ANN[*].EFFECT" \
        "ANN[*].HGVS_C" "ANN[*].HGVS_P" \
        CLINVAR_CLNSIG \
        dbNSFP_gnomAD4.1_joint_AF dbNSFP_gnomAD4.1_joint_NFE_AF \
        dbNSFP_CADD_phred dbNSFP_REVEL_score \
        > matriz_anotaciones.tsv

    echo -e "CHROM\tPOS\tREF\tALT\t\$(bcftools query -l cohorte_final.vcf.gz | tr '\n' '\t' | sed 's/\t\$//')" > matriz_genotipos.tsv
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' cohorte_final.vcf.gz >> matriz_genotipos.tsv
    """
}


// Workflow Principal
workflow {
    ch_ref     = Channel.fromPath(params.ref).first()
    ch_clinvar = Channel.fromPath(params.clinvar).first()
    ch_dbnsfp  = Channel.fromPath(params.dbnsfp).first()
    ch_clinvar_tbi = Channel.fromPath("${params.clinvar}.tbi").first()
    ch_dbnsfp_tbi  = Channel.fromPath("${params.dbnsfp}.tbi").first()

    // 1. Cargar VCFs
    ch_vcf = Channel.fromPath(params.vcfs)
        | map { file -> tuple(file.baseName.replaceAll(/.vcf|.vcf.gz/, ''), file) }

    // 2. Cargar FASTQs
    ch_fastqs = Channel.fromFilePairs(params.fastqs, size: -1)
    ch_vcf_nuevos = FASTQ_TO_VCF(ch_fastqs, ch_ref)

    // 3. Juntar todo en el mismo pipeline
    ch_todos_vcfs = ch_vcf.mix(ch_vcf_nuevos)

    // 4. Limpieza
    ch_vcfs_limpios = QC_VCF(ch_todos_vcfs, ch_ref)

    // 5. Fusion
    ch_fusion = MERGE_COHORTE(ch_vcfs_limpios.vcf_clean.collect(),
    ch_vcfs_limpios.vcf_idx.collect())

    // 6. Paralelizar anotacion por cromosoma
    ch_chroms = ch_fusion.chr_list.splitText().map { it.trim() }
    ch_anotados = ANNOTATE_CHR(ch_chroms, ch_fusion.vcf_merged, ch_fusion.vcf_idx, ch_clinvar, ch_clinvar_tbi, ch_dbnsfp, ch_dbnsfp_tbi)

    // 7. Generar las matrices finales
    EXPORT_RESULTS(ch_anotados.collect())
}