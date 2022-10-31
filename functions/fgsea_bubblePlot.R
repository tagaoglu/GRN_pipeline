#' @description Creates bubble plot for fgsea results
#' @param fgseaRes fgsea result
#' @param fdrcut defines FDR cut-off to use as output for significant signatures, 0.05 by default
#'  
#' @output The function returns a data.frame object with integrated centrality metric with corresponding p-values for each node


fgsea_bubblePlot <- function(fgseaRes, fdrcut=0.05)
{
  dat <- data.frame(fgseaRes)
  
  ### Settings
  #fdrcut <- 0.05 # FDR cut-off to use as output for significant signatures
  dencol_neg <- "blue" # bubble plot color for negative ES
  dencol_pos <- "red" # bubble plot color for positive ES
  signnamelength <- 4 # set to remove prefix from signature names (2 for "GO", 4 for "KEGG", 8 for "REACTOME")
  asp <- 3 # aspect ratio of bubble plot
  charcut <- 100 # cut signature name in heatmap to this nr of characters
  
  ### Make signature names more readable
  a <- as.character(dat$pathway) # 'a' is a great variable name to substitute row names with something more readable
  # for (j in 1:length(a)){
  #   a[j] <- substr(a[j], signnamelength+2, nchar(a[j]))
  # }
  a <- tolower(a) # convert to lower case (you may want to comment this out, it really depends on what signatures you are looking at, c6 signatures contain gene names, and converting those to lower case may be confusing)
  # for (j in 1:length(a)){
  #   if(nchar(a[j])>charcut) { a[j] <- paste(substr(a[j], 1, charcut), "...", sep=" ")}
  # } # cut signature names that have more characters than charcut, and add "..."
  a <- gsub("_", " ", a)
  
  dat$NAME <- a
  
  ### Determine what signatures to plot (based on FDR cut)
  dat2 <- subset(dat, dat$padj<fdrcut) ## Subset to pathways with FDR < 0.05
  dat2 <- dat2[order(dat2$padj),] #, decreasing = TRUE
  dat2$signature <- factor(dat2$NAME, rev(as.character(dat2$NAME)))
  
  
  ### Determine what labels to color
  sign_neg <- which(dat2[,"NES"]<0)
  sign_pos <- which(dat2[,"NES"]>0)
  # Color labels
  signcol <- rep(NA, length(dat2$signature))
  signcol[sign_neg] <- dencol_neg # text color of negative signatures
  signcol[sign_pos] <- dencol_pos # text color of positive signatures
  signcol <- rev(signcol) # need to revert vector of colors, because ggplot starts plotting these from below
  
  ### Plot bubble plot
  g<-ggplot(dat2, aes(x=padj,y=signature,size=size)) +
    geom_point(aes(fill=NES), shape=21, colour="white") +
    theme_bw() + # white background, needs to be placed before the "signcol" line
    xlim(0,fdrcut) +
    scale_size_area(max_size=10,guide="none") +
    scale_fill_gradient2(low=dencol_neg, high=dencol_pos) +
    theme(axis.text.y = element_text(colour=signcol)) +
    theme(aspect.ratio=asp, axis.title.y=element_blank()) # test aspect.ratio
  
  return(g)
}



