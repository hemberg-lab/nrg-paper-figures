BrenneckeGetVariableGenes_ggplot <- function (expr_mat, spikes = NA, suppress.plot = FALSE, fdr = 0.1, 
                                              minBiolDisp = 0.5, font_size = 10) 
{
  rowVars <- function(x) {
    unlist(apply(x, 1, var, na.rm = TRUE))
  }
  colGenes <- "black"
  colSp <- "blue"
  fullCountTable <- expr_mat
  if (is.character(spikes)) {
    sp <- rownames(fullCountTable) %in% spikes
    countsSp <- fullCountTable[sp, ]
    countsGenes <- fullCountTable[!sp, ]
  }
  else if (is.numeric(spikes)) {
    countsSp <- fullCountTable[spikes, ]
    countsGenes <- fullCountTable[-spikes, ]
  }
  else {
    countsSp <- fullCountTable
    countsGenes <- fullCountTable
  }
  meansSp <- rowMeans(countsSp, na.rm = TRUE)
  varsSp <- rowVars(countsSp)
  cv2Sp <- varsSp/meansSp^2
  meansGenes <- rowMeans(countsGenes, na.rm = TRUE)
  varsGenes <- rowVars(countsGenes)
  cv2Genes <- varsGenes/meansGenes^2
  minMeanForFit <- unname(quantile(meansSp[which(cv2Sp > 0.3)], 
                                   0.8))
  useForFit <- meansSp >= minMeanForFit
  if (sum(useForFit, na.rm = TRUE) < 20) {
    warning("Too few spike-ins exceed minMeanForFit, recomputing using all genes.")
    meansAll <- c(meansGenes, meansSp)
    cv2All <- c(cv2Genes, cv2Sp)
    minMeanForFit <- unname(quantile(meansAll[which(cv2All > 
                                                      0.3)], 0.8))
    useForFit <- meansSp >= minMeanForFit
  }
  if (sum(useForFit, na.rm = TRUE) < 30) {
    warning(paste("Only", sum(useForFit), "spike-ins to be used in fitting, may result in poor fit."))
  }
  fit <- glmgam.fit(cbind(a0 = 1, a1tilde = 1/meansSp[useForFit]), 
                    cv2Sp[useForFit])
  a0 <- unname(fit$coefficients["a0"])
  a1 <- unname(fit$coefficients["a1tilde"])
  psia1theta <- a1
  minBiolDisp <- minBiolDisp^2
  m <- ncol(countsSp)
  cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
  testDenom <- (meansGenes * psia1theta + meansGenes^2 * cv2th)/(1 + 
                                                                   cv2th/m)
  p <- 1 - pchisq(varsGenes * (m - 1)/testDenom, m - 1)
  padj <- p.adjust(p, "BH")
  sig <- padj < fdr
  sig[is.na(sig)] <- FALSE
  
  dat <- data.frame(meansGenes = meansGenes, cv2Genes = cv2Genes,
                    padj = padj < 0.1)
  xg <- 10^seq(-1, 6, length.out = 1000)
  line1 <- data.frame(x = xg, y = (a1)/xg + a0)
  line2 <- data.frame(x = xg, y = psia1theta/xg + a0 + minBiolDisp)
  if (!suppress.plot) {
    p <- ggplot(dat, aes(x = meansGenes, y = cv2Genes, color = padj)) +
      geom_point(size = 0.1) +
      scale_x_log10() +
      scale_y_log10() +
      scale_color_manual(values = c("black", "magenta", "black")) +
      guides(color=FALSE) +
      geom_line(data = line1, aes(x = x, y = y), color = "red") +
      geom_line(data = line2, aes(x = x, y = y), color = "magenta") +
      xlab("Reads") +
      ylab(expression(CV^{2})) +
      theme_classic(base_size=font_size)
  }
  return(p)
}
