setwd("/Users/pmundra/Documents/Mucosal_segments/")

files1 <- list.files(pattern = "MELA*")
files2 <- list.files(pattern = "Sample*")
files3 <- list.files(pattern = "PT*")
files <- c(files1,files3,files2)




FractGenomeAlter <- matrix(data=NA,nrow=length(files),ncol=2)
for(i in 1:length(files)){
  X <- read.table(files[i],header=TRUE)
  barscn <- data.frame(size = X$end.pos - X$start.pos,
                       CNt = X$CNt)
  cn.sizes_chr0 <- sum(barscn$size[which(barscn$CNt<2)])
  cn.sizes_chr3 <- sum(barscn$size[which(barscn$CNt>2)])
  cn.sizes_chrm <- split(barscn$size, barscn$CNt)
  cn.sizes <- sapply(cn.sizes_chrm, "sum")
  print(files[i])
  FractGenomeAlter[i,1] <- cn.sizes_chr0/sum(cn.sizes)
  FractGenomeAlter[i,2] <- cn.sizes_chr3/sum(cn.sizes)
  
}
rownames(FractGenomeAlter) <- files
colnames(FractGenomeAlter) <- c("Loss","Gain")

write.csv(FractGenomeAlter,"FractionOfGenomeAltered_GainLoss.csv")