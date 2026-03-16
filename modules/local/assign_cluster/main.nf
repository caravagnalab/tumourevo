process ASSIGN_CLUSTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3e160064566f2529f87874cfc606b160d9f58c606fbb1842ca46023da2afe8d3/data':
        'community.wave.seqera.io/library/pip_sigprofilerassignment_sigprofilerextractor_sigprofilermatrixgenerator_pruned:02a3f95da35d8c9a' }"

    input:
    tuple val(meta), path(table), val(meta2), path(results), path(genome_installed_path)
    val(tool)

    output:
    tuple val(meta), path("${meta.dataset}*")  , emit: results_sigprofiler
    path "versions.yml"                         , emit: versions

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python
    import os
    import shlex
    import shutil
    import subprocess
    import sys
    from importlib.metadata import version

    import pandas as pd
    from SigProfilerAssignment import Analyzer as Analyze

    opt = dict()
    prefix = "${task.ext.prefix ?: meta.id}"
    tool = "${tool}"
    dataset = "${meta.dataset}"
    patient = "${meta.patient}"

    if tool == 'mobster':
        sample = "${meta.tumour_sample}"

    opt.update(
    {
        "genome": "${params.genome}",
        "input_type": "matrix",
        "context_type": "96,DINUC,ID",
        "seeds": "random",
        "volume": "./",
        "make_decomposition_plots": True,
        "download_genome_sigprofiler": True,
        "genome_installed_path": "${genome_installed_path}",
    }
    )

    genome = opt["genome"]

    name = f'{prefix}_{tool}'
    out = f'{prefix}_{tool}'
    if not os.path.exists(out):
        os.makedirs(out, exist_ok=True)

    table_path = os.path.join(out, "${table}")
    os.rename("${table}", table_path)

    if opt.get("download_genome_sigprofiler", True):
        print(f"Installing genome {genome} via SigProfilerMatrixGenerator...")
        install_genome = f"SigProfilerMatrixGenerator install {genome} -v {opt['volume']}"
        subprocess.run(install_genome, shell=True)
    else:
        if not opt.get("genome_installed_path"):
            raise ValueError("download_genome_sigprofiler is False but no genome_installed_path was provided.")
        print(f"Using pre-installed genome at: {opt['genome_installed_path']}")

    # Mutation counts matrix generation
    generate_count_matrix = (
        f"SigProfilerMatrixGenerator matrix_generator {name} {genome} {out} --volume {opt['volume']}"
    )
    subprocess.run(generate_count_matrix, shell=True)

    catalog_sbs = 'SBS96/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Signatures/COSMIC_SBS96_Signatures.txt'
    catalog_id = 'ID83/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/Signatures/COSMIC_ID83_Signatures.txt'

    data_sbs = os.path.join(out, f'output/SBS/{name}.SBS96.all')
    data_id = os.path.join(out, f'output/ID/{name}.ID83.all')

    os.mkdir(os.path.join(out, 'SBS/'))
    os.mkdir(os.path.join(out, 'ID/'))

    Analyze.cosmic_fit(data_sbs,
                    os.path.join(out, 'SBS/'),
                    input_type="matrix",
                    context_type="96",
                    cosmic_version=3.3,
                    collapse_to_SBS96 = False,
                    exome=False,
                    genome_build=genome,
                    signature_database=catalog_sbs,
                    export_probabilities=True,
                    make_plots=True,
                    verbose=True)

    Analyze.cosmic_fit(data_id,
                    os.path.join(out, 'ID/'),
                    input_type="matrix",
                    context_type="83",
                    cosmic_version=3.3,
                    collapse_to_SBS96 = False,
                    exome=False,
                    genome_build=genome,
                    signature_database=catalog_id,
                    export_probabilities=True,
                    make_plots=True,
                    verbose=True)

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    SigProfilerAssignment: {version("SigProfilerAssignment")}\\n')
        f.write(f'    SigProfilerMatrixGenerator: {version("SigProfilerMatrixGenerator")}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

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
