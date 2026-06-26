## Function -- MakeAnnotationFromGTF_EnsDbStyle
#' Generate EnsDb-style GRanges from a GTF file, optionally with seqinfo from .fai
#'
#' @param gtf_file Path to GTF file (supports .gtf or .gtf.gz)
#' @param genome_name Genome name string (e.g., "Mmul10", "T2T-MMU8v2.0")
#' @param fai_file Optional path to genome .fai file for seqinfo (chromosome lengths)
#' @return GRanges object compatible with Signac/Seurat
#' @examples
#' # Example usage:
#' # gtf <- MakeAnnotationFromGTF_EnsDbStyle("path/to/genome.gtf", "hg38")
#' @export
MakeAnnotationFromGTF_EnsDbStyle <- function(gtf_file, genome_name, fai_file = NULL) {
  stopifnot(file.exists(gtf_file))
  suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(S4Vectors)
    library(GenomeInfoDb)
  })

  # Import GTF
  gtf <- rtracklayer::import(gtf_file)

  # Fix annot
  # 修复逻辑：
  # 1. 优先取 gene_biotype
  # 2. 如果没有，取 transcript_biotype
  # 3. 如果还没有，根据 gbkey 判断 (mRNA -> protein_coding, lnc_RNA -> lncRNA)
  mcols(gtf)$gene_biotype <- ifelse(
    !is.na(mcols(gtf)$gene_biotype),
    mcols(gtf)$gene_biotype,
    ifelse(!is.na(mcols(gtf)$transcript_biotype),
           mcols(gtf)$transcript_biotype,
           mcols(gtf)$gbkey) # gbkey 往往是 NCBI 的保底标签
  )
  # 规范化名称 (可选)
  mcols(gtf)$gene_biotype[mcols(gtf)$gene_biotype == "mRNA"] <- "protein_coding"

  # Keep only features of interest: exon, CDS, five_prime_UTR, three_prime_UTR
  #keep_types <- c("exon", "CDS", "five_prime_UTR", "three_prime_UTR")
  #gtf <- gtf[tolower(gtf$type) %in% tolower(keep_types)]

  # Standardize type column: unify 5' and 3' UTR as "utr"
  gtf$type <- tolower(gtf$type)
  #gtf$type[gtf$type %in% c("five_prime_utr", "three_prime_utr")] <- "utr"

  # Safely generate gene_name vector
  gene_name <- rep(NA_character_, length(gtf))
  if (!is.null(mcols(gtf)$gene_name)) {
    gene_name <- mcols(gtf)$gene_name
  }
  gene_id <- mcols(gtf)$gene_id
  gene_name[is.na(gene_name)] <- gene_id[is.na(gene_name)]

  # Set rownames for GRanges: use exon_id if available, otherwise transcript_id
  range_names <- rep(NA_character_, length(gtf))
  if (!is.null(mcols(gtf)$exon_id)) range_names <- mcols(gtf)$exon_id
  missing_idx <- is.na(range_names) & !is.null(mcols(gtf)$transcript_id)
  range_names[missing_idx] <- mcols(gtf)$transcript_id[missing_idx]
  range_names[is.na(range_names)] <- paste0("range_", which(is.na(range_names)))
  names(gtf) <- range_names

  # Simplify elementMetadata to EnsDb-style
  mcols(gtf) <- DataFrame(
    tx_id        = mcols(gtf)$transcript_id,
    gene_name    = gene_name,
    gene_id      = gene_id,
    gene_biotype = mcols(gtf)$gene_biotype,
    type         = gtf$type
  )

  # Set genome name
  genome(gtf) <- genome_name
  # seqlevelsStyle(gtf) <- "Ensembl"
  seqlevelsStyle(gtf) <- "UCSC"

  # Optionally add seqinfo from .fai
  if (!is.null(fai_file)) {
    fai <- read.table(fai_file, stringsAsFactors = FALSE)
    fai_seqinfo <- Seqinfo(seqnames = fai[,1], seqlengths = fai[,2], genome = genome_name)

    # Use seqnames intersection to safely subset
    common_seq <- intersect(seqlevels(gtf), seqlevels(fai_seqinfo))
    if (length(common_seq) < length(seqlevels(gtf))) {
      warning("Some seqlevels in GTF not found in .fai; they will be dropped from seqinfo")
    }
    seqinfo(gtf) <- fai_seqinfo[common_seq]
  }

  return(gtf)
}
