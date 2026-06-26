\name{MakeAnnotationFromGTF_EnsDbStyle}
\alias{MakeAnnotationFromGTF_EnsDbStyle}
\title{Generate EnsDb-style GRanges from a GTF file, optionally with seqinfo from .fai}
\description{
  This function creates a GRanges object from a GTF file that is compatible with Signac/Seurat workflows.
  It handles gene annotation, biotype classification, and optional chromosome length information from a .fai file.
}
\usage{
MakeAnnotationFromGTF_EnsDbStyle(gtf_file, genome_name, fai_file = NULL)
}
\arguments{
  \item{gtf_file}{Path to GTF file (supports .gtf or .gtf.gz formats).}
  \item{genome_name}{Genome name string (e.g., "Mmul10", "T2T-MMU8v2.0", "hg38").}
  \item{fai_file}{Optional path to genome .fai file for seqinfo (chromosome lengths).}
}
\details{
  The function:
  \itemize{
    \item Imports GTF file using rtracklayer
    \item Standardizes gene biotype information (prioritizing gene_biotype, then transcript_biotype, then gbkey)
    \item Normalizes mRNA biotype to "protein_coding"
    \item Creates proper gene names and IDs
    \item Sets appropriate row names for GRanges
    \item Simplifies element metadata to EnsDb-style format
    \item Sets genome name and seqlevels style
    \item Optionally adds seqinfo from .fai file
  }
}
\value{
  A GRanges object compatible with Signac/Seurat containing:
  \itemize{
    \item tx_id: transcript ID
    \item gene_name: gene name
    \item gene_id: gene ID
    \item gene_biotype: gene biotype classification
    \item type: feature type (exon, CDS, UTR, etc.)
  }
}
\examples{
# Example usage:
# gtf <- MakeAnnotationFromGTF_EnsDbStyle("path/to/genome.gtf", "hg38")
# gtf <- MakeAnnotationFromGTF_EnsDbStyle("path/to/genome.gtf.gz", "mm10", "path/to/genome.fai")
}
\seealso{
  \code{\link{rtracklayer::import}}, \code{\link{GenomicRanges::GRanges}}
}
\keyword{annotation}
\keyword{genomics}
