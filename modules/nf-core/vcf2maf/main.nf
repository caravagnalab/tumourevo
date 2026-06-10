process VCF2MAF {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7cbf9421f0bee23a93a35c5d0c7166ac1e89a40008d8e474cecfddb93226bf65/data':
        'community.wave.seqera.io/library/ensembl-vep_vcf2maf:2d40b60b4834af73' }"

    input:
    tuple val(meta), path(vcf) // Use an uncompressed VCF file!
    tuple val(meta2), path(fasta)
    // path fasta                 // Required

    output:
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def VERSION       = '1.6.22' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """

    vcf2maf.pl \\
        $args \\
        --ref-fasta $fasta \\
        --input-vcf $vcf \\
        --tumor-id "${meta.tumour_sample}" \\
        --normal-id "${meta.normal_sample}" \\
        --output-maf ${prefix}.maf \\
        --inhibit-vep

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.6.22' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """

    touch ${prefix}.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: $VERSION
    END_VERSIONS
    """
}
