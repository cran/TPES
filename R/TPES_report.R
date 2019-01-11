#' Tumor Purity Estimation using SNVs
#'
#' TPES_report function produces a graphical report regarding the allelic fraction
#' values of the putative clonal SNVs used by TPES_purity and the density function(s) computed
#' by TPES_purity.
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
#' @return A plot with:
#' \item{histogram}{Represents the allelic fraction distribution of putative clonal and subclonal
#' (if presents) SNVs within copy number neutral segments and the peak(s) detected by TPES;}
#' \item{density plot}{Represents how the density function varies according to
#' different bandwidth values (for more information see \code{\link{density}});
#' only the bandwidth values that result in at most 2 peaks are considered.}
#'
#' @examples
#' ## Generate TPES report for samples "TCGA-A8-A0A7" and "TCGA-HT-8564"
#' ## https://cancergenome.nih.gov/
#' ## Please copy and paste the following lines:
#' library(TPES)
#' TPES_report(ID = "TCGA-A8-A0A7", SEGfile = TCGA_A8_A0A7_seg,
#' SNVsReadCountsFile = TCGA_A8_A0A7_maf, ploidy = TCGA_A8_A0A7_ploidy,
#' RMB = 0.47, maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)
#'
#' TPES_report(ID = "TCGA-HT-8564", SEGfile = TCGA_HT_8564_seg,
#' SNVsReadCountsFile = TCGA_HT_8564_maf, ploidy = TCGA_HT_8564_ploidy,
#' RMB = 0.47, maxAF = 0.55, minCov = 10, minAltReads = 5, minSNVs = 10)
#'
#' @export
################################################################################
TPES_report <- function(ID,
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
  ## Define bandwidth range for the density function
  bwRange <- seq(0.010, 0.080, 0.001)


  ##############
  ## Apply ploidy correction
  log2shift <- round(-log2(ploidy/2),3)
  SEGfile$log2.plCorr <- SEGfile$seg.mean - log2shift


  ##############
  ## Init peaksTable
  peaksTable <- expand.grid(bwRange, stringsAsFactors = F)

  ##############
  ## Remove chromosomes X and Y from SEGfile
  SEGfile <- SEGfile[which(gsub('chr', '', SEGfile$chr) %in% seq(1,22,1)), ]

  ## Keep only segments with log2R corrected by ploidy in [-0.1, 0.1]
  SEGfile <- SEGfile[which(SEGfile$log2.plCorr >= -0.1 &
                             SEGfile$log2.plCorr <= 0.1), ]

  if(nrow(SEGfile) == 0){
    warning(paste0(ID, ' no CN neutral segments are available'), call. = F)
    return()}

  ##############
  ## Compute SNVs coverage
  SNVsReadCountsFile$cov <- SNVsReadCountsFile$ref.count + SNVsReadCountsFile$alt.count
  ## Compute allelic fraction
  SNVsReadCountsFile$AF <- SNVsReadCountsFile$alt.count / SNVsReadCountsFile$cov
  ## Filter SNVs
  SNVsReadCountsFile <- SNVsReadCountsFile[which(SNVsReadCountsFile$cov >= minCov &
                                                   SNVsReadCountsFile$alt.count >= minAltReads), ]
  ## Remove chromosomes X and Y from SNVsReadCountsFile
  SNVsReadCountsFile <- SNVsReadCountsFile[which(gsub('chr', '', SNVsReadCountsFile$chr) %in% seq(1,22,1)), ]

  if(nrow(SNVsReadCountsFile) == 0){
    warning(paste0(ID, ': no SNVs available\n'), call. = F)
    return()}

  ## Keep only SNVs within CN neutral segments
  SNVsReadCountsFile <- do.call(rbind,lapply(seq(1,nrow(SNVsReadCountsFile),1), function(j){
    segPos <- getSegmentsPos(SNVsReadCountsFile$chr[j],
                             SNVsReadCountsFile$start[j],
                             SNVsReadCountsFile$end[j],
                             SEGfile[,c("chrom","loc.start","loc.end")])
    if (length(segPos) == 1){
      return(SNVsReadCountsFile[j,,drop=F])}}))

  if(is.null(SNVsReadCountsFile)){
    warning(paste0(ID, ': no SNVs available in CN neutral segments\n'), call. = F)
    return()}

  ## Filter for maxAF
  SNVsReadCountsFile <- SNVsReadCountsFile[which(SNVsReadCountsFile$AF <= maxAF),]

  if (nrow(SNVsReadCountsFile) == 0){
    warning(paste0(ID, ': no SNVs available in CN neutral segments\n'), call. = F)
    return()}

  ##############
  ## Compute the number of modes of the AF distribution
  nModes  <- sapply(bwRange, function(bw){
    return(nr.modes(stats::density(SNVsReadCountsFile$AF, bw=bw, na.rm=T)$y))})

  peaksTable$nModes <- nModes
  peaksTable <- peaksTable[which(peaksTable$nModes < 3),]

  ##############
  ## Select optimal bandwidth value
  if (length(which(nModes == 2)) >= 2 ){
    bwOptim <- bwRange[min(which(nModes == 2))]
  }else if (1 %in% nModes){
    bwOptim <- bwRange[min(which(nModes == 1))]
  }else{
    warning(paste0(ID, ' no suitable bandwidth was found\n'), call. = F)
    return()}

  ##############
  ## Filter for minSNVs
  if (length(SNVsReadCountsFile$AF) < minSNVs){
    warning(paste0(ID, ': no SNVs available\n'), call. = F)
    return()}

  ##############
  ## Compute the number of peaks
  Peaks <- findPeaks(AF = SNVsReadCountsFile$AF, bw = bwOptim)


  ## If no peaks are found, return purityTable, else get the putative clonalPeak
  if (nrow(Peaks) == 0){
    warning(paste0(ID, ': no peaks found'), call. = F)
    return()
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
      warning(paste0(ID, ': no SNVs available below RMB value\n'), call. = F)
      return()}

    nModes <- sapply(bwRange, function(bw){
      return(nr.modes(stats::density(SNVsReadCountsFile$AF, bw = bw, na.rm = T)$y))})

    if (length(which(nModes == 2)) >= 2){
      bwOptim <- bwRange[min(which(nModes == 2))]
    } else if (1 %in% nModes){
      bwOptim <- bwRange[min(which(nModes == 1))]
    } else{
      warning(paste0(ID, 'no suitable bandwidth\n'), call. = F)
      return()}

    if (length(SNVsReadCountsFile$AF) < minSNVs){
      warning(paste0(ID, ': no SNVs available; required ', minSNVs, '\n'), call. = F)
      return()}

    Peaks <- findPeaks(AF = SNVsReadCountsFile$AF, bw = bwOptim)

    if (nrow(Peaks) == 0){
      warning(paste0(ID, ' no peaks found'), call. = F)
      return()
    } else {
      clonalPeak <- unlist(Peaks[nrow(Peaks),])[1]
    }

    if (clonalPeak < RMB){
      min <- findMin(AF = SNVsReadCountsFile$AF, bw = bwOptim)
    } else {
      purity <- 1
    }

  }


  ##############
  ## Isolate the putative clonal AF
  if (length(colnames(min)) == 0){
    clonalAF <- SNVsReadCountsFile$AF
  } else{
    clonalAF <- SNVsReadCountsFile$AF[which(SNVsReadCountsFile$AF > min[1,])]
    if (length(clonalAF) < minSNVs){
      warning(paste0(ID, ': ', length(clonalAF), ' putative clonal SNVs available; required ',
                     minSNVs, '\n'), call. = F)
      return()
    }
  }


  ##############
  ## Compute tumor purity
  purity <- round((clonalPeak/RMB),2)


  ##############
  # Generate report
  graphics::par(mfrow=c(1,2))

  ## FIRST PLOT
  dens <- stats::density(SNVsReadCountsFile$AF, bw = bwOptim, na.rm = TRUE)
  mylim <- max(graphics::hist(SNVsReadCountsFile$AF, plot = FALSE)$density)
  h <- graphics::hist(SNVsReadCountsFile$AF, col = 'grey60', border = 'grey60',
                      breaks = seq(0, 1, by = 0.02), main = "", xlab = "AF",
                      freq = FALSE, ylim = c(0, mylim+1))
  graphics::title("AF histogram")
  graphics::abline(v = 0.5, lty = 3, lwd = 3, col = "gray40")
  graphics::lines(dens, col = 'tomato3', lwd = 2)
  apply(Peaks, MARGIN = 1, function(x){
    graphics::abline(v = x$xID, col = "tomato4", lwd = 2)
    graphics::text(x = x$xID, y = x$yID, labels = round(x$xID, digits = 3), cex = 1.1, pos = 3)})

  ## SECOND PLOT
  e <- sapply(by(peaksTable, peaksTable$nModes, function(x){ min(x) }), identity)
  bw2 <- min(e)
  newDens <- stats::density(SNVsReadCountsFile$AF, bw = bw2, na.rm = T)
  graphics::plot(NA, type = "l", xlab = "AF", ylab = "", main = "Density function variation",
                 ylim = c(0, max(newDens$y))+1, xlim = c(0, 1), bty="n")
  graphics::abline(v = 0.5, lty = 3, lwd = 3, col = "gray40")
  colori <- grDevices::topo.colors(length(e))
  sapply(seq(1, length(e), 1), function(i, snvs, e, colori){
    graphics::lines(stats::density(SNVsReadCountsFile$AF, bw = e[[i]], na.rm = T), col = colori[i], lwd = 2)
  }, SNVsReadCountsFile, e, colori)
  graphics::legend("topright", legend = e, fill = colori, border = colori,
                   title = expression(bold("bw")), bty = "n")

  titleMain <- paste0(ID, '; Number of SNVs=', length(SNVsReadCountsFile$AF),
                      '; Bandwidth=', bwOptim, '; Purity=', purity)
  graphics::mtext(titleMain, side = 3, outer = TRUE, line = -1.2, cex = 1.3)


}
