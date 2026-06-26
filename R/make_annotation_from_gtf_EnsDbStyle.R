#' Create Annotation from GTF File in EnsDb Style
#'
#' This function creates a gene annotation data frame from a GTF file, similar to
#' the format used in EnsDb packages. It extracts gene information including
#' gene_id, gene_name, chromosome, strand, and coordinates.
#'
#' @param gtf_file Character string specifying the path to the GTF file.
#' @param gene_biotype Character string specifying which gene biotype to filter for.
#'                     Default is "protein_coding". Use NULL to include all biotypes.
#' @param include_non_coding Logical value indicating whether to include non-coding genes.
#'                           Default is FALSE.
#'
#' @return A data frame with gene annotation information in EnsDb-style format.
#' @export
#'
#' @examples
#' # Create annotation from GTF file
#' annotation <- make_annotation_from_gtf_EnsDbStyle("hg38.gtf")
#'
#' # Get only protein coding genes
#' protein_coding <- make_annotation_from_gtf_EnsDbStyle("hg38.gtf",
#'                                                     gene_biotype = "protein_coding")
#'
#' # Include non-coding genes
#' all_genes <- make_annotation_from_gtf_EnsDbStyle("hg38.gtf",
#'                                                include_non_coding = TRUE)
make_annotation_from_gtf_EnsDbStyle <- function(gtf_file,
                                               gene_biotype = "protein_coding",
                                               include_non_coding = FALSE) {
  # Input validation
  if (!file.exists(gtf_file)) {
    stop(paste("GTF file not found:", gtf_file))
  }

  # Read GTF file
  gtf_data <- read.delim(gtf_file, header = FALSE, comment.char = "#",
                        stringsAsFactors = FALSE)

  # Extract relevant columns (GTF format: seqname, source, feature, start, end, score, strand, frame, attribute)
  colnames(gtf_data) <- c("seqname", "source", "feature", "start", "end",
                         "score", "strand", "frame", "attribute")

  # Filter for gene lines
  gene_lines <- gtf_data[gtf_data$feature == "gene", ]

  # Parse gene_id and gene_name from attribute column
  gene_info <- data.frame()
  for (i in 1:nrow(gene_lines)) {
    attrs <- strsplit(gene_lines$attribute[i], "; ")[[1]]
    gene_id <- NA
    gene_name <- NA
    gene_biotype <- NA

    for (attr in attrs) {
      if (grepl("gene_id", attr)) {
        gene_id <- sub("gene_id \"([^\"]+)\".*", "\\1", attr)
      } else if (grepl("gene_name", attr)) {
        gene_name <- sub("gene_name \"([^\"]+)\".*", "\\1", attr)
      } else if (grepl("gene_biotype", attr)) {
        gene_biotype <- sub("gene_biotype \"([^\"]+)\".*", "\\1", attr)
      }
    }

    # Only keep protein coding genes if specified
    if (!is.na(gene_biotype)) {
      if (!include_non_coding && gene_biotype != "protein_coding") {
        next
      }
      if (!is.null(gene_biotype) && gene_biotype != "protein_coding") {
        next
      }
    }

    new_row <- data.frame(
      gene_id = gene_id,
      gene_name = gene_name,
      seqname = gene_lines$seqname[i],
      start = as.integer(gene_lines$start[i]),
      end = as.integer(gene_lines$end[i]),
      strand = gene_lines$strand[i],
      stringsAsFactors = FALSE
    )

    gene_info <- rbind(gene_info, new_row)
  }

  # Remove rows with missing gene_id
  gene_info <- gene_info[!is.na(gene_info$gene_id), ]

  # Return the annotation data frame
  return(gene_info)
}
