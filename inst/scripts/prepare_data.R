library(curatedMetagenomicData)
library(ape)
library(mia)

# Load data from curatedMetagenomicData
pattern <- "^2021-03-31\\.GuptaA_2019\\.(pathway_abundance|relative_abundance)$"
res <- curatedMetagenomicData(pattern = pattern, dryrun = FALSE, counts = TRUE, rownames = "short")

# If vignette/data directory do not exist, create one
dir <- file.path("vignettes", "data")
if( !dir.exists(dir) ){
    dir.create(dir, recursive = TRUE)
}

# Save taxonomy data
phylo <- res[[2]] |> rowTree()
matches <- match(rowLinks(res[[2]])[[1L]], phylo$tip.label)
phylo$tip.label[ matches ] <- rownames(res[[2]])
write.tree(phylo, file.path(dir, "phylogeny.tree"))
phylo <- read.tree(file.path(dir, "phylogeny.tree"))
rownames(res[[2]]) <- phylo$tip.label[ matches ]
write.tree(phylo, file.path(dir, "phylogeny.tree"))
tab <- res[[2]] |> rowData()
write.csv(tab, file.path(dir, "taxonomy_table.csv"))
tab <- res[[2]] |> colData()
write.csv(tab, file.path(dir, "sample_metadata.csv"))
tab <- res[[2]] |> assay()
write.csv(tab, file.path(dir, "taxonomy_abundance.csv"))

addRowSds <- function(x, assay.type = "counts", name = "sd"){
    sds <- getRowSds(x, assay.type)
    rowData(x)[[name]] <- sds
    return(x)
}

getRowSds <- function(x, assay.type = "counts"){
    sds <- assay(x, assay.type) |> rowSds()
    return(sds)
}

subsetByVariance <- function(x, assay.type = "counts", threshold = NULL, top = NULL) {

    sds <- getRowSds(x, assay.type)
    keep <- rep(TRUE, nrow(x))

    if (!is.null(threshold)) {
        keep <- keep & sds >= threshold
    }

    if (!is.null(top)) {
        top <- pmin(top, nrow(x))
        eligible <- which(keep)
        keep <- eligible[order(sds[eligible], decreasing = TRUE)[seq_len(top)]]
    }

    x <- x[keep, ]
    return(x)
}


# Save pathway data
tab <- res[[1]] |> rowData()
tab[["taxonomy"]] <- rownames(res[[1]])
tab <- mia:::.parse_taxonomy(
    tab, column_name = "taxonomy", sep = "\\.")
tab[["pathway"]] <- sub("\\|.*$", "", rownames(res[[1]]))
tab[["pathway"]] <- tab[["pathway"]] |> sub("^[^:]+:\\s*", "", x = _)
rownames(tab) <- rownames(res[[1]])
rowData(res[[1]]) <- tab
res[[1]] <- agglomerateByVariable(res[[1]], by = 1L, group = "pathway")
res[[1]] <- res[[1]][!grepl("UNINTEGRATED|UNMAPPED", rownames(res[[1]]), ignore.case = TRUE), ]
res[[1]] <- subsetByVariance(res[[1]], assay.type = "pathway_abundance", top = 50)
tab <- res[[1]] |> assay()
write.csv(tab, file.path(dir, "pathway_abundance.csv"))
