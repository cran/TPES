% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TPES_purity.R
\name{TPES_purity}
\alias{TPES_purity}
\title{Tumor Purity Estimation using SNVs}
\usage{
TPES_purity(ID, SEGfile, SNVsReadCountsFile, ploidy, RMB = 0.47,
  maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)
}
\arguments{
\item{ID}{Sample ID. Must be the same ID as in SEGfile, SNVsReadCountsFile and ploidy.}

\item{SEGfile}{A standard SEG file (segmented data). It is a data frame object that lists
loci and associated numeric values. The header must be compatible with the standard format
defined by the Broad Institute. For more information please visit
\href{https://software.broadinstitute.org/software/igv/SEG}{SEG file format}.}

\item{SNVsReadCountsFile}{A standard MAF (Mutation Annotation Format) file.
It is a data frame object containing the read counts data of somatic
single nucleotide variants (SNVs) loci. The header must contains at least informations about the chromosme
that harbors the SNV ("chr" column), the position of the SNV (defined by the "start" and "end" columns),
the sample ID ("sample" column) and finally the informations about the reference and alternative base
counts ("ref.count" and    "alt.count" columns, respectively). For more information please visit
\href{https://software.broadinstitute.org/software/igv/MutationAnnotationFormat}{MAF file format}.}

\item{ploidy}{A data frame containing the ploidy status of a sample. It must contain at
least the sample ID ("sample" column) and the ploidy status ("ploidy" column).}

\item{RMB}{The Reference Mapping Bias value. The reference genome contains only one allele
at any given locus, so reads that carry a non-reference allele are less likely to be mapped
during alignment; this causes a shift from 0.5. It can be
estimated as: \eqn{1 - medAF}, where medAF is the median value of the allelic fraction of the
sample's germline heterozygous SNPs. Default is set to 0.47. For more informations see: PMID: 19808877.}

\item{maxAF}{The filter on the allelic fraction (AF) distribution of SNVs. This is necessary to be sure to keep
only heterozygous SNVs. Clonal and subclonal SNVs, which have an AF greater than maxAF, will be removed.}

\item{minCov}{The minimum coverage for a SNV to be retained.}

\item{minAltReads}{The minimum coverage for the alternative base of a SNV to be retained.}

\item{minSNVs}{The minimum number of SNVs required to make a purity call.}
}
\value{
TPES returns a data.frame object with one row per sample and the following columns:
\item{sample}{The sample ID;}
\item{purity}{The sample purity estimated by TPES;}
\item{purity.min}{The sample minimum purity estimated by TPES;}
\item{purity.max}{The sample maximum purity estimated by TPES;}
\item{n.segs}{The number of copy number neutral segments used by TPES;}
\item{n.SNVs}{The number of SNVs used by TPES;}
\item{RMB}{The Reference Mapping Bias value used to estimate the tumor purity;}
\item{BandWidth}{The smoothing bandwidth value of the \code{\link{density}} function chosen by TPES.}
\item{log}{Reports if the run was successful; otherwise provides debugging information.}
}
\description{
TPES_purity function estimates tumor purity.
}
\examples{
## Compute tumor purity for samples "TCGA-A8-A0A7" and "TCGA-HT-8564"
## https://cancergenome.nih.gov/
## Please copy and paste the following lines:
library(TPES)
TPES_purity(ID = "TCGA-A8-A0A7", SEGfile = TCGA_A8_A0A7_seg,
SNVsReadCountsFile = TCGA_A8_A0A7_maf, ploidy = TCGA_A8_A0A7_ploidy,
RMB = 0.47, maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)

TPES_purity(ID = "TCGA-HT-8564", SEGfile = TCGA_HT_8564_seg,
SNVsReadCountsFile = TCGA_HT_8564_maf, ploidy = TCGA_HT_8564_ploidy,
RMB = 0.47, maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)

}
