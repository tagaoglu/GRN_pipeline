#From Andrej's work

save_png_pdf <- function(p=NULL, basename="plot", height=5, width=5, res=200) {
  png(paste0(basename,".png"), units= "in", height = height, width = width, res=res)
  print(p)
  dev.off()
  pdf(paste0(basename,".pdf"), height = height, width = width)
  print(p)
  dev.off()
}