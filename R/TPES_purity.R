#' Tumor Purity Estimation using SNVs
#'
#' TPES_purity function estimates tumor purity.
#'
#' @param ID Sample ID. Must be the same ID as in SEGfile, SNVsReadCountsFile and ploidy.
#' @param SEGfile A standard SEG file (segmented data). It is a data frame object that lists
#' loci and associated numeric values. The header must be compatible with the standard format
#' defined by the Broad Institute. For more information please visit
#' \href{https://software.broadinstitute.org/software/igv/SEG}{SEG file format}.
#' @param SNVsReadCountsFile A standard MAF (Mutation Annotation Format) file.
#' It is a data frame object containing the read counts data of somatic
#' single nucleotide variants (SNVs) loci. The header must contains at least informations about the chromosme
#' that harbors the SNV ("chr" column), the position of the SNV (defined by the "start" and "end" columns),
#' the sample ID ("sample" column) and finally the informations about the reference and alternative base
#' counts ("ref.count" and	"alt.count" columns, respectively). For more information please visit
#' \href{https://software.broadinstitute.org/software/igv/MutationAnnotationFormat}{MAF file format}.
#' @param ploidy A data frame containing the ploidy status of a sample. It must contain at
#' least the sample ID ("sample" column) and the ploidy status ("ploidy" column).
#' @param RMB The Reference Mapping Bias value. The reference genome contains only one allele
#' at any given locus, so reads that carry a non-reference allele are less likely to be mapped
#' during alignment; this causes a shift from 0.5. It can be
#' estimated as: \eqn{1 - medAF}, where medAF is the median value of the allelic fraction of the
#' sample's germline heterozygous SNPs. Default is set to 0.47. For more informations see: PMID: 19808877.
#' @param maxAF The filter on the allelic fraction (AF) distribution of SNVs. This is necessary to be sure to keep
#' only heterozygous SNVs. Clonal and subclonal SNVs, which have an AF greater than maxAF, will be removed.
#' @param minCov The minimum coverage for a SNV to be retained.
#' @param minAltReads The minimum coverage for the alternative base of a SNV to be retained.
#' @param minSNVs The minimum number of SNVs required to make a purity call.
#'
#' @return TPES returns a data.frame object with one row per sample and the following columns:
#' \item{sample}{The sample ID;}
#' \item{purity}{The sample purity estimated by TPES;}
#' \item{purity.min}{The sample minimum purity estimated by TPES;}
#' \item{purity.max}{The sample maximum purity estimated by TPES;}
#' \item{n.segs}{The number of copy number neutral segments used by TPES;}
#' \item{n.SNVs}{The number of SNVs used by TPES;}
#' \item{RMB}{The Reference Mapping Bias value used to estimate the tumor purity;}
#' \item{BandWidth}{The smoothing bandwidth value of the \code{\link{density}} function chosen by TPES.}
#' \item{log}{Reports if the run was successful; otherwise provides debugging information.}
#'
#' @examples
#' ## Compute tumor purity for samples "TCGA-A8-A0A7" and "TCGA-HT-8564"
#' ## https://cancergenome.nih.gov/
#' ## Please copy and paste the following lines:
#' library(TPES)
#' TPES_purity(ID = "TCGA-A8-A0A7", SEGfile = TCGA_A8_A0A7_seg,
#' SNVsReadCountsFile = TCGA_A8_A0A7_maf, ploidy = TCGA_A8_A0A7_ploidy,
#' RMB = 0.47, maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)
#'
#' TPES_purity(ID = "TCGA-HT-8564", SEGfile = TCGA_HT_8564_seg,
#' SNVsReadCountsFile = TCGA_HT_8564_maf, ploidy = TCGA_HT_8564_ploidy,
#' RMB = 0.47, maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)
#'
#' @export
################################################################################
TPES_purity <- function(ID,
                        SEGfile,
                        SNVsReadCountsFile,
                        ploidy,
                        RMB         = 0.47,
                        maxAF       = 0.55,
                        minCov      = 10,
                        minAltReads = 5,
                        minSNVs     = 10)

  {


  ##############
  ## Get sample data
  colnames(SEGfile) <- c('ID','chrom','loc.start','loc.end','num.mark','seg.mean')
  if(ID %in% SEGfile$ID){
    SEGfile <- SEGfile[which(SEGfile$ID == ID),]
  } else {stop(paste0(ID, ' is not in the SEGfile\n'))}

  if(ID %in% SNVsReadCountsFile$sample){
    SNVsReadCountsFile <- SNVsReadCountsFile[which(SNVsReadCountsFile$sample == ID),]
  } else {stop(paste0(ID, ' is not in the SNVsReadCountsFile\n'))}

  if(ID %in% ploidy$sample){
    ploidy <- ploidy[which(ploidy$sample == ID),]$ploidy
  } else {stop(paste0(ID, ' is not in the ploidy file\n'))}


  ##############
  ## Apply ploidy correction
  log2shift <- round(-log2(ploidy/2),3)
  SEGfile$log2.plCorr <- SEGfile$seg.mean - log2shift


  ##############
  ## Init Purity Table
  purityTable            <- data.frame(row.names = ID, stringsAsFactors = F)
  purityTable$sample     <- ID
  purityTable$purity     <- NA
  purityTable$purity.min <- NA
  purityTable$purity.max <- NA
  purityTable$n.segs     <- NA
  purityTable$n.SNVs     <- NA
  purityTable$RMB        <- NA
  purityTable$BandWidth  <- NA
  purityTable$log        <- NA


  ##############
  ## Remove chromosomes X and Y from SEGfile
  SEGfile <- SEGfile[which(gsub('chr', '', SEGfile$chr) %in% seq(1,22,1)), ]

  ## Keep only segments with log2R corrected by ploidy in [-0.1, 0.1]
  SEGfile <- SEGfile[which(SEGfile$log2.plCorr >= -0.1 &
                             SEGfile$log2.plCorr <= 0.1), ]

  if (nrow(SEGfile) == 0){
    purityTable$log <- 'no available segments'
    warning(paste0(ID, ' no CN neutral segments are available'), call. = F)
    return(purityTable)}


  ##############
  ## Compute SNVs coverage
  SNVsReadCountsFile$cov <- SNVsReadCountsFile$ref.count + SNVsReadCountsFile$alt.count
  # Compute allelic fraction
  SNVsReadCountsFile$AF <- SNVsReadCountsFile$alt.count / SNVsReadCountsFile$cov
  # Filter SNVs
  SNVsReadCountsFile <- SNVsReadCountsFile[which(SNVsReadCountsFile$cov >= minCov &
                                                   SNVsReadCountsFile$alt.count >= minAltReads), ]
  ## Remove chromosomes X and Y from SNVsReadCountsFile
  SNVsReadCountsFile <- SNVsReadCountsFile[which(gsub('chr', '', SNVsReadCountsFile$chr) %in% seq(1,22,1)), ]

  if (nrow(SNVsReadCountsFile) == 0){
    purityTable$n.SNVs <- 0
    purityTable$log    <- 'no SNVs available'
    warning(paste0(ID, ': no SNVs available\n'), call. = F)
    return(purityTable)}

  ## Keep only SNVs within CN neutral segments
  SNVsReadCountsFile <- do.call(rbind,lapply(seq(1,nrow(SNVsReadCountsFile),1), function(j){
    segPos <- getSegmentsPos(SNVsReadCountsFile$chr[j],
                             SNVsReadCountsFile$start[j],
                             SNVsReadCountsFile$end[j],
                             SEGfile[,c('chrom', 'loc.start', 'loc.end')])
    if (length(segPos) == 1){
      return(SNVsReadCountsFile[j,,drop=F])}}))

  if (is.null(SNVsReadCountsFile)){
    purityTable$n.SNVs <- 0
    purityTable$log    <- 'no SNVs available in CN neutral segments'
    warning(paste0(ID, ': no SNVs available in CN neutral segments\n'), call. = F)
    return(purityTable)}

  ## Filter for maxAF
  SNVsReadCountsFile <- SNVsReadCountsFile[which(SNVsReadCountsFile$AF <= maxAF),]

  if (nrow(SNVsReadCountsFile) == 0){
    purityTable$n.SNVs <- 0
    purityTable$log    <- paste0('no SNVs available with maxAF = ', maxAF)
    warning(paste0(ID, ': no SNVs available with maxAF = ', maxAF, '\n'), call. = F)
    return(purityTable)}


  ##############
  ## Compute the number of modes of the AF distribution
  bwRange <- seq(0.010, 0.080, 0.001)
  nModes  <- sapply(bwRange, function(bw){
    return(nr.modes(stats::density(SNVsReadCountsFile$AF, bw=bw, na.rm=T)$y))})


  ##############
  ## Select optimal bandwidth
  if (length(which(nModes == 2)) >= 2 ){
    bwOptim <- bwRange[min(which(nModes == 2))]
  }else if (1 %in% nModes){
    bwOptim <- bwRange[min(which(nModes == 1))]
  }else{
    purityTable$log <- 'no suitable bandwidth'
    warning(paste0(ID, ': no suitable bandwidth was found\n'), call. = F)
    return(purityTable)}


  ##############
  ## Filter for minSNVs
  if (length(SNVsReadCountsFile$AF) > minSNVs){
    purityTable$n.SNVs <- length(SNVsReadCountsFile$AF)
  } else {
    purityTable$log <- 'no SNVs available'
    warning(paste0(ID, ': no SNVs available; required ', minSNVs, "\n"), call. = F)
    return(purityTable)}


  ##############
  ## Compute the number of peaks
  Peaks <- findPeaks(AF = SNVsReadCountsFile$AF, bw = bwOptim)

  ## If no peaks are found, return purityTable, else get the putative clonalPeak
  if (nrow(Peaks) == 0){
    purityTable$log <- 'no peaks found'
    warning(paste0(ID, ': no peaks found'), call. = F)
    return(purityTable)
  } else {
    clonalPeak <- unlist(Peaks[nrow(Peaks),])[1]
  }

  ## If the putative clonalPeak < RMB, find the minimum in the distribution
  if (clonalPeak < RMB){
    min <- findMin(AF = SNVsReadCountsFile$AF, bw = bwOptim)
  } else {
    # a. Remove SNVs whose AF > RMB;
    # b. Try to recompute clonalPeak
    SNVsReadCountsFile <- SNVsReadCountsFile[which(SNVsReadCountsFile$AF <= RMB),]

    if (nrow(SNVsReadCountsFile) == 0){
      purityTable$n.SNVs <- 0
      purityTable$log    <- 'no SNVs available below RMB value'
      warning(paste0(ID, ': no SNVs available below RMB value\n'), call. = F)
      return(purityTable)}

    nModes <- sapply(bwRange, function(bw){
      return(nr.modes(stats::density(SNVsReadCountsFile$AF, bw = bw, na.rm = T)$y))})

    if (length(which(nModes == 2)) >= 2){
      bwOptim <- bwRange[min(which(nModes == 2))]
    } else if (1 %in% nModes){
      bwOptim <- bwRange[min(which(nModes == 1))]
    } else{
      purityTable$log <- 'no suitable bandwidth'
      warning(paste0(ID, 'no suitable bandwidth\n'), call. = F)
      return(purityTable)}

    if (length(SNVsReadCountsFile$AF) > minSNVs){
      purityTable$n.SNVs <- length(SNVsReadCountsFile$AF)
    } else {
      purityTable$log <- 'no SNVs available'
      warning(paste0(ID, ': no SNVs available; required ', minSNVs, '\n'), call. = F)
      return(purityTable)}

    Peaks <- findPeaks(AF = SNVsReadCountsFile$AF, bw = bwOptim)

    if (nrow(Peaks) == 0){
      warning(paste0(ID, ' no peaks found'), call. = F)
      return(purityTable)
    } else {
      clonalPeak <- unlist(Peaks[nrow(Peaks),])[1]
    }

    if (clonalPeak < RMB){
      min <- findMin(AF = SNVsReadCountsFile$AF, bw = bwOptim)
    } else {
      # Set purity to 1
      purityTable$purity     <- 1
      purityTable$purity.min <- 1
      purityTable$purity.max <- 1
      purityTable$n.segs     <- nrow(SEGfile)
      purityTable$n.SNVs     <- nrow(SNVsReadCountsFile)
      purityTable$RMB        <- RMB
      purityTable$BandWidth  <- bwOptim
      purityTable$log        <- 'computation ok'
      return(purityTable)
    }

  }


  ##############
  ## Isolate the putative clonal AF
  if (length(colnames(min)) == 0){
    clonalAF <- SNVsReadCountsFile$AF
  } else{
    clonalAF <- SNVsReadCountsFile$AF[which(SNVsReadCountsFile$AF > min[1,])]
    if (length(clonalAF) < minSNVs){
      purityTable$n.SNVs <- length(clonalAF)
      purityTable$log    <- 'no putative clonal SNVs available'
      warning(paste0(ID, ': ', length(clonalAF), ' putative clonal SNVs available; required ',
                     minSNVs, '\n'), call. = F)
      return(purityTable)
    }
  }


  ##############
  ## Compute tumor purity
  purityTable$purity     <- round((clonalPeak/RMB),2)
  purityTable$purity.min <- round(((stats::quantile(stats::ecdf(clonalAF))[2])/RMB),2)
  purity.max             <- round(((stats::quantile(stats::ecdf(clonalAF))[4])/RMB),2)
  if(purity.max > 1){
    purityTable$purity.max <- 1
  } else {
    purityTable$purity.max <- purity.max
  }
  purityTable$n.segs     <- nrow(SEGfile)
  purityTable$n.SNVs     <- length(clonalAF)
  purityTable$RMB        <- RMB
  purityTable$BandWidth  <- bwOptim
  purityTable$log        <- 'computation ok'


  ##############
  # Return purityTable
  return(purityTable)



}
