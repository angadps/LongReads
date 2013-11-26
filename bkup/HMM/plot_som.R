
snvr=commandArgs(TRUE)[1]
errr=commandArgs(TRUE)[2]
cov=commandArgs(TRUE)[3]
readlen=commandArgs(TRUE)[4]
prev=commandArgs(TRUE)[5]
base=paste(snvr,errr,cov,readlen,21,sep="_")
inp=paste(base,".c",sep="")
op=paste(base,"_som.jpeg",sep="")

data = read.table(inp, header=TRUE)

# log.data = data
# log.data$count = log(log.data$count)

jpeg(
  file      = op,
  width     = 2.36,
  height    = 2.36,
  units     = "in",
  res       = 200,
  pointsize = 7
)

plot(data$NB_Specificity, data$NB_Sensitivity)
plot(data$Sam_Specificity, data$Sam_Sensitivity, type="l", xlab="1 - Specificity", ylab="Sensitivity", xlim=c(0,1), ylim=c(0,1), col="red")
lines(data$NB_Specificity, data$NB_Sensitivity,col="green")
text(0.4,0.6,paste("coverage =",cov), adj = c(0,NA))
text(0.4,0.5,paste("Avg Prevalence =",prev), adj = c(0,NA))
legend(0.4,0.4, c("Samtools","HapMut"), lty=c(1,1), lwd=c(1.0,1.0),col=c("red","green"))

