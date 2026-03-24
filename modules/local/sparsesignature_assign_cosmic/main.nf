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
    
    download.file(
     url = "https://raw.githubusercontent.com/SigProfilerSuite/SigProfilerAssignment/main/SigProfilerAssignment/data/Reference_Signatures/GRCh38/COSMIC_v3.5_SBS_GRCh38.txt",
     destfile = "COSMIC_v3.5_SBS_GRCh38.txt",
     mode = "wb"
    )
    
    cosmic_path <- "COSMIC_v3.5_SBS_GRCh38.txt"
    
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
    map_sparsesig_to_cosmic <- function(sparsesig_out, sparsesig_count,cosmic_path, threshold = 0.8) {
      # Check inputs
      
      samples <- rownames(sparsesig_count) 
      if (is.null(sparsesig_out[["beta"]]) || is.null(sparsesig_out[["alpha"]])) {
        stop("sparsesig_out must contain 'beta' and 'alpha'.")
      }
      
      # Load COSMIC reference
      cosmic_signatures <- read.delim(cosmic_path, check.names = FALSE) %>%
        tibble::column_to_rownames("Type") %>%
        as.matrix()
    
      # Extract de novo signatures and exposures
      de_novo_signatures <- t(sparsesig_out[["beta"]]) %>%  as.matrix()
      de_novo_exposures  <- as.matrix(sparsesig_out[["alpha"]])
    
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
      
      # Build mapping list: one entry per de novo signature
      # Each entry is a named numeric vector of cosine similarities
      sim_matrix_all <- vector("list", length = nrow(similarity_matrix))
      names(sim_matrix_all) <- rownames(similarity_matrix)
      
      for (de_novo_sig in rownames(similarity_matrix)) {
        similarities <- similarity_matrix[de_novo_sig, ]
        
        if (de_novo_sig == "Background") {
          # Force Background to SBS5
          if (!"SBS5" %in% colnames(similarity_matrix)) {
            stop("Background is mapped to SBS5, but SBS5 is not present in COSMIC reference.")
          }
          sim_matrix_all[[de_novo_sig]] <- similarities["SBS5"]
        } else {
          above_threshold <- similarities[similarities >= threshold]
          above_threshold <- sort(above_threshold, decreasing = TRUE)[1]
          
          if (length(above_threshold) == 0) {
            sim_matrix_all[[de_novo_sig]] <- NA_real_
          } else {
            sim_matrix_all[[de_novo_sig]] <- above_threshold
          }
        }
      }
      
      # Determine all COSMIC signatures that appear in mappings
      cosmic_sigs <- unique(unlist(
        lapply(sim_matrix_all, function(x) {
          if (length(x) == 1 && is.na(x)) {
            return(NULL)
          }
          names(x)
        })
      ))
      
      # Ensure SBS5 exists if Background is present
      if ("Background" %in% names(sim_matrix_all) && !"SBS5" %in% cosmic_sigs) {
        cosmic_sigs <- c(cosmic_sigs, "SBS5")
      }
      
      # Create remapped exposure matrix
      samples <- rownames(de_novo_exposures)
      remapped_exposures <- matrix(
        0,
        nrow = length(samples),
        ncol = length(cosmic_sigs),
        dimnames = list(samples, cosmic_sigs)
      )
      
      # Remap exposures
      for (de_novo_sig in colnames(de_novo_exposures)) {
        if (!de_novo_sig %in% names(sim_matrix_all)) {
          warning(sprintf("No similarity entry found for de novo signature: %s", de_novo_sig))
          next
        }
        
        mapping <- sim_matrix_all[[de_novo_sig]]
        
        # Skip signatures with no match above threshold
        if (length(mapping) == 1 && is.na(mapping)) {
          warning(sprintf(
            "No COSMIC match above threshold for de novo signature: %s",
            de_novo_sig
          ))
          next
        }
        
        # Normalize similarity weights
        sim_weights <- mapping / sum(mapping)
        
        # Add weighted exposure to mapped COSMIC signatures
        for (sig in names(sim_weights)) {
          tmp <- remapped_exposures[, sig] +
            de_novo_exposures[, de_novo_sig] * sim_weights[sig]
          remapped_exposures[, sig] <- round(tmp,0)
        }
      }
      
      # Convert to proportions
      row_totals <- rowSums(remapped_exposures)
      remapped_exposures_prop <- remapped_exposures
    
      nonzero_rows <- row_totals > 0
      remapped_exposures_prop[nonzero_rows, ] <-
        remapped_exposures[nonzero_rows, , drop = FALSE] / row_totals[nonzero_rows]
      
      remapped_exposures_prop[!nonzero_rows, ] <- 0
      rownames(remapped_exposures_prop) <- samples
      return(list(
        remapped_exposures_prop = remapped_exposures_prop,
        sim_matrix_all = sim_matrix_all,
        similarity_matrix = similarity_matrix
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
