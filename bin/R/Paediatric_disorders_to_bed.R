# Run separately

library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Read the TSV file
gene_data <- read.delim(
    "data/references/Paediatric disorders.tsv", 
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
)


ensembl_ids <- unique(gene_data[, "EnsemblId.GRch37."])

ensembl_ids <- ensembl_ids[ensembl_ids != ""]

# Retrieve gene location info
gene_locations <- getBM(
    attributes = c(
        "ensembl_gene_id", "external_gene_name",
        "chromosome_name", "start_position",
        "end_position", "strand"
    ),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
)


bed <- gene_locations[, c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id", "external_gene_name")]

# Writte the gene locations to a BED file
write.table(
    bed, 
    file = "data/referencesPaediatric_disorders_new.bed",
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
