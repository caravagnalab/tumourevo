process SPARSESIGNATURE_ASSIGN {
    tag "$meta.id"
    label "process_single"
    label "error_retry"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/25/252387028f8b9db88bfe20235617fd2ff8712fb23b848dadae4dc4cff35ae89b/data':
        'community.wave.seqera.io/library/r-mutationalpatterns:0.2b--3b9858739b22ee16' }"

    input:
    tuple val(meta), path(signatures_nmfOut_rds)
    tuple val(meta), path(signatures_mutCounts_rds)
    val(genome)
    
    output:
    tuple val(meta), path("*.rds"), emit: sparsesignature_assigned
    path "versions.yml",            emit: versions

    script:
    def args   = task.ext.args   ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(dplyr)
    
    url_reference <- paste0("https://raw.githubusercontent.com/SigProfilerSuite/SigProfilerAssignment/main/SigProfilerAssignment/data/Reference_Signatures/","${genome}","/COSMIC_v3.5_SBS_","${genome}",".txt")
    destfile_path <- paste0("COSMIC_v3.5_SBS_","${genome}",".txt")
    download.file(
     url = url_reference,
     destfile = destfile_path,
     mode = "wb"
    )
    
    cosmic_path <-destfile_path
    
    cos_sim <- function(x, y) {
      res <- x %*% y / (sqrt(x %*% x) * sqrt(y %*% y))
      # coerce matrix to numeric
      res <- as.numeric(res)
      return(res)
    }
    
    cos_sim_matrix <- function (mut_matrix1, mut_matrix2) {
            if (!all(apply(mut_matrix1, 2, is.numeric))) {
                stop("The first input contains non-numeric columns, while all columns should be numeric.")
            }
            if (!all(apply(mut_matrix2, 2, is.numeric))) {
                stop("The second input contains non-numeric columns, while all columns should be numeric.")
            }
            n_samples1 <- ncol(mut_matrix1)
            n_samples2 <- ncol(mut_matrix2)
            res_matrix <- matrix(nrow = n_samples1, ncol = n_samples2)
            for (s in seq_len(n_samples1)) {
                signal1 <- mut_matrix1[, s, drop = TRUE]
                cos_sim_vector <- c()
                for (i in seq_len(n_samples2)) {
                    signal2 <- mut_matrix2[, i, drop = TRUE]
                    cos_sim_vector[i] <- cos_sim(signal1, signal2)
                }
                res_matrix[s, ] <- cos_sim_vector
            }
            rownames(res_matrix) <- colnames(mut_matrix1)
            colnames(res_matrix) <- colnames(mut_matrix2)
            return(res_matrix)
        }
    map_sparsesig_to_cosmic <- function(
        sparsesig_out,
        sparsesig_count,
        cosmic_path,
        threshold = 0.8,
        top_n = 2,
        min_weight = 0.1
    ) {
    
    # Check inputs
    samples <- rownames(sparsesig_count)
    
    if (is.null(sparsesig_out[["beta"]]) || is.null(sparsesig_out[["alpha"]])) {
        stop("sparsesig_out must contain 'beta' and 'alpha'.")
    }
    
    # Load COSMIC reference
    cosmic_signatures <- read.delim(cosmic_path, check.names = FALSE) |>
        tibble::column_to_rownames("Type") |>
        as.matrix()
    
    # Extract de novo signatures and exposures
    de_novo_signatures <- t(sparsesig_out[["beta"]]) |> as.matrix()
    de_novo_exposures  <- as.matrix(sparsesig_out[["alpha"]])
    rownames(de_novo_exposures) <- samples
    
    # Check that mutation types overlap
    common_types <- intersect(rownames(cosmic_signatures), rownames(de_novo_signatures))
    
    if (length(common_types) == 0) {
        stop("No shared mutation types between COSMIC reference and de novo signatures.")
    }
    
    # Restrict both matrices to shared mutation types in identical order
    cosmic_signatures  <- cosmic_signatures[common_types, , drop = FALSE]
    de_novo_signatures <- de_novo_signatures[common_types, , drop = FALSE]
    
    # Compute cosine similarity matrix
    similarity_matrix <- cos_sim_matrix(
        de_novo_signatures,
        cosmic_signatures
    )
    
    # Initialize output objects
    # Store decomposition weights for each de novo signature
    # Example: S1 -> SBS5 = 0.7, SBS44 = 0.3
    mapping_weights <- list()
    reconstruction_quality <- data.frame()
    
    # Decompose each DeNovo signature
    
    for (de_novo_sig in rownames(similarity_matrix)) {
        similarities <- similarity_matrix[de_novo_sig, ]
        
        # SparseSignatures background is treated as SBS5
        if (de_novo_sig == "Background") {
        if (!"SBS5" %in% colnames(cosmic_signatures)) {
            stop("Background is mapped to SBS5, but SBS5 is not present in COSMIC reference.")
        }
        
        mapping_weights[[de_novo_sig]] <- c(SBS5 = 1)
        
        reconstruction_quality <- rbind(
            reconstruction_quality,
            data.frame(
            de_novo_signature = de_novo_sig,
            reconstruction_cosine = similarities["SBS5"],
            n_cosmic_components = 1
            )
        )
        next
        }
        
        # Keep signatures above cosine similarity threshold
        candidates <- names(similarities[similarities >= threshold])
        
        if (length(candidates) == 0) {
        candidates <- names(sort(similarities, decreasing = TRUE))[seq_len(min(top_n, length(similarities)))]
        } else {
        candidates <- names(sort(similarities[candidates], decreasing = TRUE))[seq_len(min(top_n, length(candidates)))]
        }
        
        target <- de_novo_signatures[, de_novo_sig]
        reference <- cosmic_signatures[, candidates, drop = FALSE]
        
        # Objective function: minimise squared reconstruction error
        objective <- function(w) {
        reconstructed <- as.vector(reference %*% w)
        sum((target - reconstructed)^2)
        }
        
        # Fit non-negative weights using constrained optimization
        fit <- stats::optim(
        par = rep(1 / length(candidates), length(candidates)),
        fn = objective,
        method = "L-BFGS-B",
        lower = rep(0, length(candidates))
        )
        
        weights <- fit[["par"]]
        
        if (sum(weights) == 0) {
        mapping_weights[[de_novo_sig]] <- NA_real_
        next
        }
        
        # Normalize weights to sum to 1
        weights <- weights / sum(weights)
        names(weights) <- candidates
        # Remove weak contributors
        weights <- weights[weights >= min_weight]
        weights <- weights / sum(weights)
        
        mapping_weights[[de_novo_sig]] <- weights
        
        reconstructed <- as.vector(reference[, names(weights), drop = FALSE] %*% weights)
        
        reconstruction_quality <- rbind(
        reconstruction_quality,
        data.frame(
            de_novo_signature = de_novo_sig,
            # cosine similarity between reconstructed and original de novo signature
            reconstruction_cosine = sum(target * reconstructed) /
            sqrt(sum(target^2) * sum(reconstructed^2)),
            n_cosmic_components = length(weights)
        )
        )
    }
    
    # Build remapped exposure matrix
    cosmic_sigs <- unique(unlist(lapply(mapping_weights, names)))
    cosmic_sigs <- cosmic_sigs[!is.na(cosmic_sigs)]
    
    remapped_exposures <- matrix(
        0,
        nrow = nrow(de_novo_exposures),
        ncol = length(cosmic_sigs),
        dimnames = list(rownames(de_novo_exposures), cosmic_sigs)
    )
    
    # Redistribute de novo exposures to COSMIC signatures
    for (de_novo_sig in colnames(de_novo_exposures)) {
        mapping <- mapping_weights[[de_novo_sig]]
        
        if (length(mapping) == 1 && is.na(mapping)) {
        warning(sprintf("No COSMIC decomposition found for %s", de_novo_sig))
        next
        }
        
        for (cosmic_sig in names(mapping)) {
        remapped_exposures[, cosmic_sig] <-
            remapped_exposures[, cosmic_sig] +
            de_novo_exposures[, de_novo_sig] * mapping[cosmic_sig]
        }
    }
    
    # Convert to proportions
    row_totals <- rowSums(remapped_exposures)
    
    remapped_exposures_prop <- remapped_exposures
    nonzero_rows <- row_totals > 0
    
    remapped_exposures_prop[nonzero_rows, ] <-
        remapped_exposures[nonzero_rows, , drop = FALSE] / row_totals[nonzero_rows]
    
    remapped_exposures_prop[!nonzero_rows, ] <- 0
    
    return(list(
        remapped_exposures = remapped_exposures,
        remapped_exposures_prop = remapped_exposures_prop,
        mapping_weights = mapping_weights,
        similarity_matrix = similarity_matrix,
        reconstruction_quality = reconstruction_quality
    ))
    }
    nmf_out <- readRDS("${signatures_nmfOut_rds}")
    c_matrix <- readRDS("${signatures_mutCounts_rds}")
    
    assign_cosimic <- map_sparsesig_to_cosmic(sparsesig_out=nmf_out, sparsesig_count=c_matrix,
      cosmic_path = cosmic_path, threshold = 0.8)
    saveRDS(object=assign_cosimic, file=paste0("$prefix", "_cosmic_assigned.rds"))  
    # Version export #####
    f = file("versions.yml","w")
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    close(f)
    
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}