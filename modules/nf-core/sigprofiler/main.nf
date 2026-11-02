process SIGPROFILER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5f/5fb22903af2da5fd1747c02f514afcd9878809760f9f32543b11a0f233d8aa83/data' :
        'community.wave.seqera.io/library/python_pip_sigprofilerassignment_sigprofilerextractor_pruned:ba7b67176e1074c6'}"

    input:
    tuple val(meta), path(tsv_list, stageAs: '*.tsv')
    val(genome)                  // genome version
    path(genome_installed_path)  //optional

    output:
    tuple val(meta), path("results/*")    , emit: results_sigprofiler
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template "main_script.py"

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def context_types = task.ext.context_type?.split(',') ?: ['96','DINUC','ID']

    // Map context to signature type
    def context_map = [
        '96'    : 'SBS96',
        'DINUC' : 'DBS78',
        'ID'    : 'ID83'
    ]
    def signatures = context_types.collect { context_map[it] }

    """
    mkdir -p results/input
    touch results/input/input_data.txt

    # Create per-context outputs

    for sig in ${signatures.join(" ")}; do

        mkdir -p results/\$sig/\$sig/Suggested_Solution/COSMIC_\${sig}_Decomposed_Solution/Signatures/
        touch    results/\$sig/\$sig/Samples.txt
        touch    results/\$sig/\$sig/Suggested_Solution/COSMIC_\${sig}_Decomposed_Solution/Signatures/COSMIC_\${sig}_Signatures.txt

        type=\$(echo "\$sig" | sed 's/[0-9]*//')   # e.g. SBS96 → SBS, ID83 → ID, DBS78 → DBS
        mkdir -p results/output/\$type
        touch    results/output/\$type/${prefix}.\${sig}.all
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
        SigProfilerAssignment: \$(python3 -c "import SigProfilerAssignment as sig; print(sig.__version__)")
        SigProfilerExtractor: \$(python3 -c "import SigProfilerExtractor as sig; print(sig.__version__)")
        SigProfilerMatrixGenerator: \$(python3 -c "import SigProfilerMatrixGenerator as matGen; print(matGen.__version__)")
        pandas: \$(python3 -c "import pandas as pd; print(pd.__version__)")
        sigProfilerPlotting: \$(python3 -c "import sigProfilerPlotting as sig; print(sig.__version__)")
    END_VERSIONS
    """
}
