args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	inputFile1="plink.hwe"
}else{
	inputFile1=args[1]	
}

if (length(args)<2) {
	inputFile2="plinkzoomhwe.hwe"
}else{
	inputFile2=args[2]	
}

hwe<-read.table (file=inputFile1, header=TRUE)
pdf("histhwe.pdf")
hist(hwe[,9],main="Histogram HWE")
dev.off()

hwe_zoom<-read.table (file=inputFile2, header=TRUE)
pdf("histhwe_below_threshold.pdf")
hist(hwe_zoom[,9],main="Histogram HWE: strongly deviating SNPs only")
dev.off()
