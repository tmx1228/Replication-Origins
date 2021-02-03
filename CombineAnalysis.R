###env
BiocManager::install("magrittr")
library(magrittr)
library(reshape2)
library(ggplot2)
library(reshape2)
library(devtools)
library(dplyr) #mutate
library(ggrepel) #Avoid overlapping labels
library(VennDiagram)
#library(easyGgplot2)
library(grid)
library(gridExtra)





setwd("/Users/mt9tz/home/220_DNA_Replication/combine_analysis/GC_count")
# totalPeak <- read.table("totalPeak.GCcount")[,5]
# cover.over24.peaks <- read.table("cover.over24.peaks.GCcount")[,5]
# shared24peaks.must.have.3 <- read.table("shared24peaks.must.have.3.GCcount")[,5]
# shared24peaks.must.have.4 <- read.table("shared24peaks.must.have.4.GCcount")[,5]
# cover.over48.peaks <- read.table("cover.over48.peaks.GCcount")[,5]
all.gc <- read.table("all.GCcount", header = 1)[,5]
b.gc <- read.table("Bubble.GCcount", header = 1)[,5]
ns.gc <- read.table("NS.GCcount", header = 1)[,5]
orc.gc <- read.table("ORC.GCcount", header = 1)[,5]
rep.gc <- read.table("Repli.GCcount", header = 1)[,5]
re.gc <- read.table("Rerep.GCcount", header = 1)[,5]
shared.gc <- read.table("shared_origins.GCcount", header = 1)[,5]

###### GC content historgam for shared peaks
hist(shared.gc, main = "GC content distribution of shared origins", xlab = "")

###### GC content box plot with diff selections
head(shared.gc)
df <- data.frame(c(totalPeak,cover.over24.peaks,shared24peaks.must.have.3,shared24peaks.must.have.4,cover.over48.peaks),c(rep("1", length(totalPeak)),rep("2", length(cover.over24.peaks)),rep("3", length(shared24peaks.must.have.3)),rep("4", length(shared24peaks.must.have.4)),rep("5", length(cover.over48.peaks))))
colnames(df) <- c("GC_content","Peak_selection_method")
mylabel <- c("totalPeak","cover over24 samples","cover over 24 samples and 3 data types", "cover over 24 samples and 4 data types", "cover over 48 samples")
mycol <- c("gray60","gray50","gray40","gray30","gray90")
mylegendtitle <- "Peak selection"
p <- ggplot(df,aes(x=Peak_selection_method,y=GC_content,color=Peak_selection_method)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x=element_blank())+
  scale_color_manual(name=mylegendtitle,
                     labels=mylabel,
                     values=mycol) 
p

###### GC content box plot with diff data types
df <- data.frame(c(all.gc, ns.gc, rep.gc, re.gc,b.gc),c(rep("1", length(all.gc)),rep("2", length(ns.gc)),rep("3", length(rep.gc)),rep("4", length(re.gc)),rep("5", length(b.gc))))
colnames(df) <- c("GC_content","Data type")
mylegendtitle <- "Data type"
mycol <- c("#BB9D00", "#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
head(df)
p <- ggplot(df,aes(x=`Data type`,y=GC_content,color=`Data type`)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x=element_blank()) +
  scale_color_manual(name=mylegendtitle,
                     labels=group.list,
                     values=mycol) 
p









setwd("~/home/220_DNA_Replication/combine_analysis/genomic_annotation")
mat.all <- read.table("genomic_coverage.all_totalPeak.matrix")  
mat.b <- read.table("genomic_coverage.Bubble_totalPeak.matrix")
mat.ns <- read.table("genomic_coverage.NS_old_totalPeak.matrix")
mat.repli <- read.table("genomic_coverage.Repli_totalPeak.matrix") 
mat.rerep <- read.table("genomic_coverage.Rerep_totalPeak.matrix") 
mat.shared <- read.table("genomic_coverage.shared_origins.matrix")
mat.orc <- read.table("genomic_coverage.ORC_totalPeak.matrix")   

mat.list <- c("mat.all","mat.ns","mat.repli","mat.rerep","mat.b","mat.shared","mat.orc")
#test <- read.table("genomic_coverage.test.matrix") #test <- test[,c(1,2,3,7,8,9,10,11,12)]
group.list <- c("All Technics","NS-seq","Repli-seq","Rerep-seq","Bubble-seq","Shared Origins","ORC ChIP-seq Peaks")
names <- c("chr","start","end","promoter","exon","intron","ctcf","hotspot","desert","pk","vally")

########## genome coverage bar distribution
## calculate sample coverage
mat=data.frame()
df=data.frame()
for (i in 1:5){
  i = 6  ## for shared origins
  i = 7  ## for ORC
  mymat <- get(mat.list[i])
  colnames(mymat) <- names
  tmp <- mymat %>% mutate(genome.group = ifelse(promoter > 0, "Promoter", ifelse(exon > 0, "Exon", ifelse(intron > 0, "Intron", "Intergenic"))))
  genome <- tmp$genome.group
  #count <- tapply(genome, genome, function(x) length(x))
  prep <- tapply(genome, genome, function(x) length(x)/length(genome))
  mat <- rbind(mat, prep)
  #df <- rbind(df, cbind(group=rep(group.list[i],length(genome)), genome))
}
## draw percentage barplot of data types
colnames(mat) <- c("Exon","Intergenic","Intron","Promoter")
rownames(mat) <- c("All Technics","NS-seq","Repli-seq","Rerep-seq","Bubble-seq")
head(mat)
df <- melt(t(mat)[c(4,1,3,2),])
colnames(df) <- c("Annotation","Data type", "Percentage")
head(df)
#df$`Data type`
theme <- theme(panel.background = element_blank(), panel.border=element_rect(fill=NA), 
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               strip.background=element_blank(), axis.text.x=element_text(angle = 90, hjust = 1, colour="black"), 
               axis.text.y=element_text(colour="black"), axis.ticks=element_line(colour="black"), 
               plot.margin=unit(c(1,1,1,1),"line"))

p <- ggplot(df, aes(x=`Data type`, y=Percentage, group=`Data type`)) +
  geom_bar(stat="identity", position="fill", aes(fill=Annotation)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label=round(Percentage,3)*100), position=position_fill(vjust = 0.5)) +
  theme +
  scale_fill_manual(values=c("#98FB98", "#6495ED", "#E0E0E0", "#FF7F50"))
p
## draw pie plot of shared peaks and ORC peaks
mat
rownames(mat) <- "Percentage"
df <- data.frame(t(mat))
df$type <- rownames(df)
df
theme <- theme(panel.background = element_blank(), panel.border=element_rect(fill=NA), 
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               strip.background=element_blank(), axis.text.x=element_text(angle = 90, hjust = 1, colour="black"), 
               axis.text.y=element_text(colour="black"), axis.ticks=element_line(colour="black"), 
               plot.margin=unit(c(1,1,1,1),"line"))
library(scales)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
Percentage <- df$Percentage*100
pie + scale_fill_grey() +  blank_theme +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(y = Percentage/4 + c(0, cumsum(Percentage)[-length(Percentage)]), 
                label = percent(Percentage/1)), size=5)


p <- ggplot(df, aes(x="", y=Percentage, fill=type)) +
  geom_bar(stat="identity", width = 1) +
  coord_polar("y", start=0) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  #geom_text(aes(label=round(Percentage,3)*100), position=position_fill(vjust = 0.5)) +
  theme +
  scale_fill_manual(values=c( "#6495ED","#FF7F50", "#E0E0E0", "#98FB98")) +
  xlab("") +
  theme(axis.text.x=element_blank())
  geom_text(aes(y = Percentage/4 + c(0, cumsum(Percentage)[-length(round(Percentage*100,4))]), 
                label = round(Percentage*100,4), size=5))
p

####### overlap with constitutive CTCF (with or without fisher test)
mat=data.frame()
for (i in 1:5){
  mymat <- get(mat.list[i])
  colnames(mymat) <- c("chr","start","end","tss","exon","intron","ctcf","hotspot","desert")
  tmp <- mymat %>% mutate(ctcf.group = ifelse(ctcf > 0, "Cover with CTCF", "Not cover with CTCF"))
  ctcf <- tmp$ctcf.group
  count <- tapply(ctcf, ctcf, function(x) length(x))
  prep <- tapply(ctcf, ctcf, function(x) length(x)/length(ctcf))
  mat <- rbind(mat, c(count, prep))
  #df <- rbind(df, cbind(group=rep(group.list[i],length(genome)), genome))
}
colnames(mat) <- c("Cover with CTCF","Not cover with CTCF","Cover with CTCF rate","Not cover with CTCF rate")
rownames(mat) <- c("All Technics","NS-seq","Repli-seq","Rerep-seq","Bubble-seq")
head(mat)
t=data.frame()

#fisher test
for (i in 2:5){
  mymat <- mat[i,]
  allmat <- mat[1,]
  TeaTasting <-
    matrix(c(mymat$`Cover with CTCF`, mymat$`Not cover with CTCF`, allmat$`Cover with CTCF`, allmat$`Not cover with CTCF`),
           nrow = 2,
           dimnames = list(Guess = c("cover with CTCF", "not cover"),
                           Truth = c("cover with CTCF", "not selected")))
  fisher <- fisher.test(TeaTasting, alternative = "greater")
  t <- rbind(t,fisher$p.value)
  print(TeaTasting)
}

##with p-value
#df <- as.data.frame(cbind(`Cover with CTCF rate` = mat[2:5,3], `p_value` = t[,1]))
#df$group <- c("NS-seq","Repli-seq","Rerep-seq","Bubble-seq")

##without p-value
df <- as.data.frame(mat[,3])
colnames(df) <- "Cover with CTCF rate"
df$group <- c("All technics","NS-seq","Repli-seq","Rerep-seq","Bubble-seq")

df
p <- ggplot(df, aes(x=group, y=`Cover with CTCF rate`, group = group, fill = group)) +
  geom_bar(stat="identity") +
  scale_x_discrete(limits=c("All technics", "NS-seq", "Repli-seq", "Rerep-seq","Bubble-seq")) +
  scale_fill_manual(values=c("#BB9D00", "#C77CFF", "#F8766D", "#7CAE00", "#00BFC4")) +
  scale_y_continuous(labels = scales::percent) +
  xlab(NULL)+
  ylab("Cover with constitutive CTCF percentage") +
  geom_text(aes(label=round(`Cover with CTCF rate`,3)*100), vjust = -0.5) +
  theme
p


######### Hotspot & desert coverage
mat=data.frame()
for (i in 1:5){
  mymat <- get(mat.list[i])
  colnames(mymat) <- c("chr","start","end","tss","exon","intron","ctcf","hotspot","desert")
  tmp <- mymat %>% mutate(ac.group = ifelse(hotspot > 0, "Cover with Hotspot", ifelse(desert > 0, "Cover with Desert", "Not cover")))
  ac <- tmp$ac.group
  count <- tapply(ac, ac, function(x) length(x))
  prep <- tapply(ac, ac, function(x) length(x)/length(ac))
  mat <- rbind(mat, c(count, prep))
  #df <- rbind(df, cbind(group=rep(group.list[i],length(genome)), genome))
}
colnames(mat) <- c("Cover with Hotspot","Cover with Desert","No overlap","Cover with Hotspot rate","Cover with Desert rate","No overlap rate")
rownames(mat) <- c("All Technics","NS-seq","Repli-seq","Rerep-seq","Bubble-seq")
head(mat)

df <- melt(t(mat[,4:5]))
df
colnames(df) <- c("Annotation","Data type", "Percentage")
head(df)
#df$`Data type`
theme <- theme(panel.background = element_blank(), panel.border=element_rect(fill=NA), 
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               strip.background=element_blank(), axis.text.x=element_text(angle = 90, hjust = 1, colour="black"), 
               axis.text.y=element_text(colour="black"), axis.ticks=element_line(colour="black"), 
               plot.margin=unit(c(1,1,1,1),"line"))

p <- ggplot(df, aes(x=`Data type`, y=Percentage, group=`Data type`)) +
  geom_bar(stat="identity",  aes(fill=Annotation)) +  #position="fill",
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label=round(Percentage,5)*100)) +   #, position=position_fill(vjust = 0.5)
  theme +
  scale_fill_manual(values=c("#FF8C00", "#E0E0E0", "#FFFFFF"))
p



######### Replication timing initiation sites(pk) & vally coverage

names
mat=data.frame()
for (i in 1:7){
  mymat <- get(mat.list[i])
  colnames(mymat) <- names
  tmp <- mymat %>% mutate(it.group = ifelse(pk > 0, "Replication initiation zones", ifelse(vally > 0, "Replication termination zones", "Not cover")))
  it <- tmp$it.group
  count <- tapply(it, it, function(x) length(x))
  prep <- tapply(it, it, function(x) length(x)/length(it))
  mat <- rbind(mat, c(count, prep))
  #df <- rbind(df, cbind(group=rep(group.list[i],length(genome)), genome))
}
colnames(mat) <- c("Other","Replication initiation zones","Replication termination zones","No overlap rate","Cover with initiation zones rate","Cover with termination zones rate")
rownames(mat) <- c("All Technics","NS-seq","Repli-seq","Rerep-seq","Bubble-seq","Shared origins","ORC ChIP-seq peaks")
mat

df <- melt(t(mat[,5:6]))
df
colnames(df) <- c("Annotation","Data type", "Percentage")
head(df)
#df$`Data type`
theme <- theme(panel.background = element_blank(), panel.border=element_rect(fill=NA), 
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               strip.background=element_blank(), axis.text.x=element_text(angle = 90, hjust = 1, colour="black"), 
               axis.text.y=element_text(colour="black"), axis.ticks=element_line(colour="black"), 
               plot.margin=unit(c(1,1,1,1),"line"))

p <- ggplot(df, aes(x=`Data type`, y=Percentage, group=`Data type`)) +
  geom_bar(stat="identity",  aes(fill=Annotation)) +  #position="fill",
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label=round(Percentage,5)*100), position=position_stack(vjust = 0.5)) + #position=position_fill(vjust = 0.5/100)
  theme +
  scale_fill_manual(values=c("olivedrab1", "khaki3", "#FFFFFF"))
p






#########shared peak bart2 result dotplot
bart.share <- read.table("/Users/mt9tz/home/220_DNA_Replication/combine_analysis/bart2/shared_origins_sampleCoverage/shared_origins_sampleCoverage_bart_results.txt", header = 1)
bart.shared.3 <- read.table("/Users/mt9tz/home/220_DNA_Replication/combine_analysis/bart2/totalPeak_sampleCoverage/totalPeak_sampleCoverage_bart_results.txt", header = 1)
bart.shared.4 <- read.table("/Users/mt9tz/home/220_DNA_Replication/combine_analysis/bart2/totalPeak_withoutORC_sampleCoverage/totalPeak_withoutORC_sampleCoverage_bart_results.txt", header = 1)

setwd("/Users/mt9tz/home/220_DNA_Replication/combine_analysis/bart2")

datalist.bart <- c("bart.shared.3","bart.shared.4", "bart.share")
namelist.bart <- c("All shared peaks bart2 result","Shared peaks(without ORC) bart2 result","Shared peaks")

# bart RNA-seq top20 dot plot
for (i in 1:length(datalist.bart)){
  #i=3  ##shared peaks
  usedata1 <- get(datalist.bart[i])
  #select top20
  df <-  mutate(usedata1, sig=ifelse(usedata1$re_rank <= usedata1[10,"re_rank"], "top 10",  "others")) #Will have different colors depending on significance
  df$index <- as.integer(rownames(df))
  myOutput <- paste("Rplot_bart2result_",i,".png",sep="")
  myName <- namelist.bart[i]
  p <- ggplot(df,aes(index,-log10(irwin_hall_pvalue),label=TF)) + 
    geom_point(aes(col=sig)) + #add points colored by significance
    scale_color_manual(values=c("black", "red")) + #blue for neg
    ggtitle(myName) +#e.g. 'Volcanoplot DESeq2'
    geom_text_repel(data = df[which(df$sig == "top 10"),], 
                    size = 10,
                    nudge_x = 700, 
                    direction = "y",
                    angle        = 0,
                    hjust        = 1,
                    segment.size = 0.2) +
    xlab("TFs") +
    theme_classic(base_size = 25) 
  png(myOutput, width = 520, height = 480)
  print(p)
  dev.off() 
}
















setwd("/Users/mt9tz/home/220_DNA_Replication/combine_analysis")
mpeaks <- read.table("totalPeak_intersectionCount.matrix",  header = 1, sep = "\t")
ns <- read.table("NS_totalPeak.bed")
rep <- read.table("Repli_totalPeak.bed")
re <- read.table("Rerep_totalPeak.bed")
b <- read.table("Bubble_totalPeak.bed")
#orc <- read.table("ORC_totalPeak.bed")
setwd("/Users/mt9tz/home/220_DNA_Replication/combine_analysis/R_output")


################################ peak length distribution of technics
pdf("Rplot_peak_len.pdf")
par(mfrow=c(3,2))
colors=c("#F8766D","#7CAE00","#00BFC4","#C77CFF","#BB9D00")
hist(ns[,3]-ns[,2],xlim=c(0,2000),breaks=300,main="NS-seq",col=colors[1], xlab = "origin lenth")
hist(rep[,3]-rep[,2],xlim=c(0,2000),breaks=450,main="Repli-seq",col=colors[2], xlab = "origin lenth")
hist(re[,3]-re[,2],xlim=c(0,2000),breaks=15,main="Rerep-seq",col=colors[3], xlab = "origin lenth")
hist(b[,3]-b[,2],xlim=c(0,2000),breaks=30,main="Bubble-seq",col=colors[4], xlab = "origin lenth")
hist(mpeaks[,3]-mpeaks[,2],xlim=c(0,2000),breaks=500,main="All data type",col=colors[5], xlab = "origin lenth")
dev.off()

############################# peak length violin plot
df.all.peaks <- as.data.frame(rbind(all[,1:3],ns,rep,orc,b))
df.all.peaks$group <- c(rep("All",nrow(all)),rep("NS-seq",nrow(ns)),rep("Repli-seq",nrow(rep)),rep("ORC ChIP-seq",nrow(orc)),rep("Bubble-seq",nrow(b)))
df.all.peaks$group <- as.factor(df.all.peaks$group)
p<-ggplot(df.all.peaks, aes(x = group, y = V3-V2)) 
pdf("Rplot_peak_len_violin.pdf")
p+geom_violin(aes(fill=group))+theme + ylab("peak length") + ylim(0,2000)
dev.off()













setwd("/Users/mt9tz/home/220_DNA_Replication/combine_analysis")
mpeaks <- read.table("totalPeak_intersectionCount.matrix",  header = 1, sep = "\t")
ns.mpeaks <- read.table("NS_totalPeak_intersectionCount.matrix")
rep.mpeaks <- read.table("Repli_totalPeak_intersectionCount.matrix")
re.mpeaks <- read.table("Rerep_totalPeak_intersectionCount.matrix")
b.mpeaks <- read.table("Bubble_totalPeak_intersectionCount.matrix")
orc.mpeaks <- read.table("ORC_totalPeak_intersectionCount.matrix")
setwd("/Users/mt9tz/home/220_DNA_Replication/combine_analysis/R_output")

mylist <- c("ns.mpeaks","rep.mpeaks","re.mpeaks","b.mpeaks","mpeaks")
mytitle <- c("NS-seq","Repli-seq","Rerep-seq","Bubble-seq","All data type")


############################################################### origins sample coverage
pdf("Rplot_sample_coverage_hist.pdf",width=5, height=5)
theme <- theme(panel.background = element_blank(), panel.border=element_rect(fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               strip.background=element_blank(), 
               plot.title=element_text(size=20,face="bold",colour="black"),
               axis.text.x=element_text(size=16,colour="black"), 
               axis.text.y=element_text(size=16,colour="black"), 
               axis.ticks=element_line(colour="black"), 
               plot.margin=unit(c(1,1,1,1),"line"))
#df.count <- data.frame()
for (i in 1:5){
	group <- mylist[i]
	m <- get(group)
	clean.peaks <- m[,4:ncol(m)]
	clean.peaks[clean.peaks>1] <- 1
	sample.count <- data.frame(sample.count=apply(clean.peaks, 1, sum),group=rep(group,nrow(clean.peaks)))
	#table(sample.count)
	p <- ggplot(sample.count, aes(x = sample.count)) +
	  theme +
	#theme_grey(base_size = 20) +
	  geom_histogram(binwidth = 1,alpha = 0.5, fill = colors[i]) +
 	  #scale_y_log10() + 
	  labs(title=mytitle[i], x="Sample count", y="Sites(in log10 scale)")
	print(p)
	#df.count <- rbind(df.count,sample.count)
}
dev.off()
p


######################################## ORC peak sample coverage for model fitting
clean.peaks <- orc.mpeaks[,4:ncol(orc.mpeaks)]
head(clean.peaks)
clean.peaks[clean.peaks>1] <- 1
clean.peaks.orc1 <- clean.peaks[,1:3]  #3
clean.peaks.orc2 <- clean.peaks[,4:5]   #2

sample.count.orc <- apply(clean.peaks, 1, sum)
sample.count.orc1 <- apply(clean.peaks.orc1, 1, sum)
sample.count.orc2 <- apply(clean.peaks.orc2, 1, sum)
output <- table(sample.count.orc)
write.table(output, "peaks_cover_sample_count_orc.csv", sep = "\t", quote = F, row.names = F)
output <- table(sample.count.orc2)
write.table(output, "peaks_cover_sample_count_orc2.csv", sep = "\t", quote = F, row.names = F)
output <- table(sample.count.orc1)
write.table(output, "peaks_cover_sample_count_orc1.csv", sep = "\t", quote = F, row.names = F)



######################################### Total peak sample coverage
setwd("/Users/mt9tz/home/220_DNA_Replication/combine_analysis")
mpeaks <- read.table("totalPeak_intersectionCount.matrix",  header = 1, sep = "\t")


mpeaks[1:5,1:5]
#rownames(mpeaks) <- paste0("peak",seq(nrow(mpeaks)))
clean.peaks <- mpeaks[,4:ncol(mpeaks)]
dim(clean.peaks)
clean.peaks[clean.peaks>1] <- 1
clean.peaks.ns <- clean.peaks[,1:59]  #59
clean.peaks.repli <- clean.peaks[,60:86]   #27
clean.peaks.re <- clean.peaks[,87:95] #9
clean.peaks.b <- clean.peaks[,96:98]  #3

#count.per.sample <- apply(mpeaks[,4:ncol(mpeaks)], 2, sum)
#mpeaks.over300 <- mpeaks[,-c(19,20,46,48,56,57,58,59,60,62,79,84,85,87,94,95)]
#clean.peaks <- mpeaks.over300[,4:ncol(mpeaks.over300)]
clean.peaks[1:6,1:5]
table(clean.peaks[,1])
table(clean.peaks[,6])

sample.count <- apply(clean.peaks, 1, sum)
sample.count.ns <- apply(clean.peaks.ns, 1, sum)
sample.count.repli <- apply(clean.peaks.repli, 1, sum)
sample.count.re <- apply(clean.peaks.re, 1, sum)
sample.count.b <- apply(clean.peaks.b, 1, sum)

score.avg.sample.count <- sample.count.ns/59 + sample.count.repli/27 + sample.count.re/9 + sample.count.b/3

#clean.peaks <- cbind(mpeaks, sample.count, sample.count.ns, sample.count.repli, sample.count.re, sample.count.b)
#dim(clean.peaks)

all.sample.count <- cbind(mpeaks[,1:3], sample.count, sample.count.ns, sample.count.repli, sample.count.re, sample.count.b)
head(all.sample.count)
write.table(all.sample.count, "all.sample.count.txt", sep = "\t", quote = F, row.names = T, col.names = T)
shared <- subset(all.sample.count, sample.count >= 25 & sample.count.ns >= 16 & sample.count.repli >= 4 & sample.count.re >= 1 & sample.count.b >= 1)
dim(shared)[1]
mymat = data.frame()
for(i in 21:30){
  mydata <- subset(all.sample.count, sample.count >= i & sample.count.ns >= 16 & sample.count.repli >= 4 & sample.count.re >= 1 & sample.count.b >= 1)
  mymat <- rbind(mymat,c(i,dim(mydata)[1]))
}

write.table(shared[,1:3], "shared_origins.bed", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(shared[,1:4], "shared_origins_sampleCoverage.bed", sep = "\t", quote = F, row.names = F, col.names = F)

#mydata <- subset(use.peaks, sample.count > 24 & sample.count.ns > 0 & sample.count.repli > 0 & sample.count.b > 0)
#mydata <- subset(use.peaks, sample.count > 24 & sample.count.ns > 0 & sample.count.repli > 0 & sample.count.b > 0 & sample.count.orc > 0)

dim(mydata)
len <- mydata[,3] - mydata[,2]
df <- cbind(mydata, len)
head(df)

p <- ggplot(df, aes(x = sample.count)) +
  theme_grey(base_size = 20) +
  geom_histogram(binwidth = 1,alpha = 0.1, color="salmon") +
  scale_y_log10() + 
  labs(title="Sample coverage level", x="Sample count", y="Sites(in log10 scale)")
print(p)


p <- ggplot(df, aes(x=sample.count)) + 
  geom_histogram(scale=log10)
p

clean.peaks[1:5,1:5]
dim(clean.peaks)


