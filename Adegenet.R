install.packages("adegenet", dep=TRUE)
install.packages("remotes")
remotes::install_github("nspope/r2vcftools")
install.packages("adegenet")
install.packages("vcfR")
install.packages("poppr")
install.packages("vcftools")
install.packages("dartR")
install.packages("diveRsity")
install.packages("ape", dependencies = T)
install.packages("gplots")
install.packages('rmetasim')
install.packages("VennDiagram")
install.packages("phangorn")
install.packages("LEA")
install.packages("gsubfn")
install.packages("gfortran")
install.packages("gsubfn")
install.packages("terra")
install.packages("rgdal")
install.packages("metasim")
install.packages("gdsfmt")
install.packages("PopGenReport")
install.packages("knitr")
install.packages("KRIS")
install.packages("ggtext")
install.packages("HardyWeinberg")
install.packages("mice")
devtools::install_github("pievos101/PopGenome")
BiocManager::install("SNPRelate")
devtools::install_github("softloud/metasim")
BiocManager::install("gdsfmt")
if (!require('devtools')) install.packages('devtools')
# install from GitHub


library(metasim)
library(knitr)
library("genefilter")
library("adegenet")
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library('vcfR')
library('poppr')
library('hierfstat')
library("dartR")
library("Biostrings")
library("ggtree")
library("ggtext")
library("gridtext")
library("ggplot2")
library("VennDiagram")
library("phangorn")
library(terra)
library(rgdal)
library(SNPRelate)
library(gdsfmt)
library(KRIS)
library(PopGenome)
library(HardyWeinberg)
library(PopGenReport)
library(mice)

source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.09.txt")
getpackages(myrepos = 'http://cran.us.r-project.org')
setwd('/Users/luisfuentes/Desktop/plink')
importdata(inputprefix = 'sambarFormat', samplefile = NULL, sumstatsfile = FALSE, 
           depthfile = FALSE, nchroms = NULL)


heatmap.2

FgenlightWithR
genlight2sambar(genlight_object = 'FgenlightWithR', do_confirm = TRUE)
filterdata(indmiss=0.25,snpmiss=0.1,min_mac=2,dohefilter=TRUE,min_spacing=500,nchroms=NULL,TsTvfilter=NULL)
mygenlight
snps
inds$popcol
pidf
mysambar$populations <- c("Ba1","Ba2","Ba3","Is4","Is5","Is6")
mysambar
calcdiversity()  
calcpi()
findstructure(Kmax=6,add_legend=TRUE,legend_pos="right",legend_cex=2,pop_order=NULL)
calcdistance()
calckinship()
inferdemography()


devtools::install_github("hadley/devtools")
install.packages("devtools")
options(download.file.method = "libcurl")
devtools::install_github("uqrmaie1/admixtools")
remotes::install_github("uqrmaie1/admixtools")
library("devtools")
install_github("uqrmaie1/admixtools")
utils::download.file("https://github.com/uqrmaie1/admixtoolss", destfile = tempfile(fileext=".zip"), method="libcurl")
file.exists(tempfile(fileext=".zip"))
library("admixtools")

##### FINDSTRUCTURE #####
### correspondence analyses:
mysambar$res.ca <<- CA(mydf,ncp=5,graph=FALSE)
mysambar$res.ca <<- dudi.coa(mydf,scannf=FALSE,nf=5)
### PCOA analyses:
neimatrix <- stamppNeisD(mygenlighttemp[indstemp$filter,snps$filter],pop=FALSE)
p <- pcoa(neimatrix, correction='none', rn=NULL)
hammingmatrix<- bitwise.dist(mygenlight[indselection,snpselection],mat=TRUE)
p <- pcoa(hammingmatrix, correction='none', rn=NULL)
calcpi(pi_per_pop=FALSE,myinput=mygenlight,indselection=inds$filter,snpselection=snps$filter,popnames=mysambar$populations,corrected=TRUE)
pidf<- mysambar$tajdlist_meta[[1]]
pimatrix<- pidf2matrix(pidf=mysambar$tajdlist_meta[[1]],myinds=inds$nr[inds$filter])
p <- pcoa(pimatrix, correction='none', rn=NULL)

devtools::install_github("pievos101/PopGenome")
library(PopGenome)

thetaw(fileWithR)


devtools::install_github(repo="zakrobinson/RLDNe")
library(RLDNe)

Fgenlightboth <- vcfR2genlight(fileWithR, n.cores = 1)
pop.databoth <- read.table("/Users/luisfuentes/Desktop/infoPopulation89.txt", sep="\t", header = T)
pop.databoth$Loc
#confirm that all samples are in the infoPopulation
file@gt
#to confirm if the names in vcf are the same with the population data file
all(colnames(file@gt)[-1] == pop.databoth$ID)
ploidy(Fgenlightboth) <- 2
pop(Fgenlightboth) <- pop.databoth$Loc
Fgenlightboth$pop

gl.LDNe(Fgenlightboth, outfile = "nePianguaboth.txt", outpath ="/Users/luisfuentes/Desktop",
        neest.path = "/Users/luisfuentes/Downloads/Zip Folder_64_bit_191125")
gl2genepop(Fgenlightboth, outfile ="both2.genepop", outpath = "/Users/luisfuentes/")

vcf2gpop()
tempdir()
getwd()
neest.path 

F_ST.stats()


require(strataG)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

install.packages("proto")
install.packages("KRIS")
install.packages("beast")
library("gsubfn")
library("proto")
library("KRIS")
library("ggtree")
library("beast")
library(P2C2M.SNAPP)
install.packages("/Users/luisfuentes/Downloads/P2C2M_SNAPP-master/P2C2M.SNAPP_1.0.0.tar.gz", repos = NULL, type = "source")
beast(fileWithR)
data("FungalGrowthDataset")

p2c2m.snapp()

gl2snapp(FgenlightWithR, outfile = "withReference.nex", outpath ="/Users/luisfuentes/Desktop" )
vcf2


run.p2c2m.snapp(path_to_xml ="/Users/luisfuentes/Desktop/" , 
                xml_file = "snapp.xml", 
                sum_stat = c("RF", "MLSD", "PFST"), num_sims = 100, dist_reps = 100, 
                sample_unif = TRUE, delete_sims = FALSE, save_graphs = TRUE, 
                save_output = TRUE, run_mode = 1)
FgenlightWithR
fileWithR

a <- read.csv("/Users/luisfuentes/Desktop/vcf/a.txt", header = FALSE)
b <- read.csv("/Users/luisfuentes/Desktop/vcf/b.txt", header = F)
GENOME.class <- readData("/Users/luisfuentes/Desktop/vcf2/", format = "VCF")
GENOME.class
GENOME.class <- set.populations(GENOME.class, list(
  c("genotype_10_1", "genotype_11", "genotype_12", "genotype_13", "genotype_14_1", "genotype_15", "genotype_16", "genotype_17_1", "genotype_1_1", "genotype_2", "genotype_3_1", "genotype_4","genotype_5", "genotype_6", "genotype_7", "genotype_8_1", "genotype_9"),
  c("genotype_18", "genotype_19_1", "genotype_20", "genotype_21", "genotype_22_1", "genotype_23", "genotype_24", "genotype_26", "genotype_27", "genotype_28", "genotype_30", "genotype_31", "genotype_32", "genotype_33", "genotype_34", "genotype_35"),
  c("genotype_36", "genotype_37", "genotype_38", "genotype_39", "genotype_40", "genotype_41", "genotype_42", "genotype_43", "genotype_44", "genotype_45", "genotype_46", "genotype_47", "genotype_48", "genotype_49",  "genotype_50", "genotype_51"), 
  c("genotype_52", "genotype_53", "genotype_54", "genotype_55", "genotype_56", "genotype_57", "genotype_58", "genotype_59", "genotype_60", "genotype_61", "genotype_62", "genotype_63", "genotype_64", "genotype_65"), 
  c("genotype_66", "genotype_68", "genotype_69", "genotype_71", "genotype_72", "genotype_73", "genotype_74", "genotype_75", "genotype_76", "genotype_77", "genotype_78", "genotype_79"), 
  c("genotype_80", "genotype_82", "genotype_83", "genotype_84", "genotype_85", "genotype_86", "genotype_87", "genotype_88", "genotype_89", "genotype_90", "genotype_91", "genotype_92", "genotype_93", "genotype_94")))
GENOME.class <- neutrality.stats(GENOME.class, FAST = F)
get.neutrality(GENOME.class)
GENOME.class <- F_ST.stats(GENOME.class)
GENOME.class <- F_ST.stats.2(GENOME.class)
GENOME.class <- diversity.stats(GENOME.class, pi='Nei')
GENOME.class@region.stats@nucleotide.diversity
GENOME.class@Hudson.K_ST
get.diversity(GENOME.class)[[2]]
GENOME.class@nuc.diversity.between
GENOME.class@genes

F_ST.stats(GENOME.class@populations[1:2])

slide <- sliding.window.transform(GENOME.class, width = 10000, jump = 10000)
#------------------------------------------------------#
#--------------  Aligned sequences --------------------#
#------------------------------------------------------#

alin <- read.dna("/Users/luisfuentes/Desktop/sequence.aln", format = "clustal")
alinBI <- read.dna("/Users/luisfuentes/Desktop/sequences_BMIs.aln", format = "clustal")

alinBI
alinGenid <- DNAbin2genind(alinBI)
alinGenid

amova(alinBI)
alleles2(alin, ploidy = 2, rownames = NULL, 
             population = NULL, phased = F)

mutations(netBI)
hap <- haplotype(alin)
hapBI <- haplotype(alinBI, trailingGapsAsN=FALSE, strict=FALSE)
hapBI
hap <- sort(hap, what = "labels")
hapBI <- sort(hapBI, what = "labels")
hapBI
net <- haploNet(hap)
netBI <- haploNet(hapBI)
hapFre<-haploFreq(alinBI, haplo = hapBI, what = 3 )
haploFreq(alin, haplo = hap, what = 3 )
hapFre
summary(hapFre)
length(hapFre)

range(1)

hapFre
sum = 0
for (i in 1:length(hapFre)){
  if(i > 49){
    if(hapFre[i] != 0){
      sum = sum + 1
    }
  }
}

for (i in 1:length(hapFre)){
  if(i){
    
  }
}
  
plot(netBI, fast=T)
summary(hap)
write.nexus.data(hapBI, file = "haplotypeNetwork.nex")
write.nexus.data(alin, file = "sequences.nex")
write.nexus.data(hapFre, file = "hapFreqBI.nex", gap = " ")


haploFreq(alinBI, haplo = h)

hap.div(alinBI, method="Nei", variance=T)
nuc.div(alinBI, variance = T)
tajima.test(alinBI)

as.network(net)
as.igraph(net)
plot(net, size=0.5)

install.packages("GGally")
install.packages("network")
install.packages("sna")
install.packages("scales")
library(GGally)
library(ggplot2)
library(network)
library(sna)
library(scales)
ggnet2(netBI, node.size=3)
netBI
print.default(netBI)
haploNet()

install.packages("/Users/luisfuentes/Downloads/Rfunctions.zip", repos = NULL, type = "source")
#------------------------------------------------------#
#--------  Population Genetic Analysis   --------------#
#------------------------------------------------------#

#to load VCF file
file <- read.vcfR("/Users/luisfuentes/Desktop/filtered_bi_89.vcf", verbose = FALSE)
fileTr <- read.vcfR("/Users/luisfuentes/Desktop/tranlationdeNovoGBS_89_0.01_bi.vcf", verbose = F)
fileWithR <- read.vcfR("/Users/luisfuentes/Desktop/withReference_89_0.01.vcf", verbose = F)

matrix <- as.matrix(extract.gt(fileWithR))
matrix
phy <- phyDat(matrix, type = "USER", levels = c("0", "1", "2"), as.character = TRUE)
phy
write.nexus.data(phy, file="WithReference.nex")


#to convert loaded file to genind object
Fgenind2 <- vcfR2genind(file, pop=NULL)
Fgenind2Tr<-vcfR2genind(fileTr, pop=NULL)
Fgenind2WithR<-vcfR2genind(fileWithR, pop=NULL)
Fgenind2WithR
nInd(Fgenind2)
indNames(Fgenind2)
nLoc(Fgenind2)
locNames(Fgenind2)

minorAllele(Fgenind2WithR)
#to convert loaded file to genlight object
Fgenlight <- vcfR2genlight(file, n.cores = 1)
FgenlightTr <- vcfR2genlight(fileTr, n.cores = 1)
FgenlightWithR <- vcfR2genlight(fileWithR, n.cores = 1)
FgenlightWithR
nInd(Fgenlight)
indNames(Fgenlight)
locNames(Fgenlight)

FgenlightWithR

allelic.richness(Fgenind2WithR)
pairwise.betas(Fgenind2WithR, diploid=TRUE)
fst.dosage(Fgenind2WithR, pop=Fgenind2WithR[,1])
fst.hudson(fileWithR)


#pop data
pop.data <- read.table("/Users/luisfuentes/Desktop/infoPopulations.txt", sep="\t", header = T)

pi<-read.table("/Users/luisfuentes/Desktop/UniAndes/pi_withReference.windowed_10.pi", header = T)
piBM<-read.table("/Users/luisfuentes/Desktop/UniAndes/pi_BM_withReference.windowed_10.pi", header = T)
piIs<-read.table("/Users/luisfuentes/Desktop/UniAndes/pi_Is_withReference.windowed_10.pi", header = T)
tajima<-read.table("/Users/luisfuentes/Desktop/Uniandes/tajima_withReference_10.Tajima.D", header = T)
tajimaBM<-read.table("/Users/luisfuentes/Desktop/UniAndes/tajima_BM_withReference_10.Tajima.D", header = T)
tajimaIs<-read.table("/Users/luisfuentes/Desktop/UniAndes/tajima_Is_withReference_10.Tajima.D", header = T)
tajima$BIN_START
pi$
hist(pi$tajima)
mean(piIs$PI)
mean(tajimaIs$TajimaD)
boxplot(pi$PI)
replace(stats, "NA", 0)
pi$CHROM
count=0
sum=0
boxplot(tajimaData)

mean(pi$PI)
mean(piBM$PI)
mean(piIs$PI)

mean(tajimaBM$TajimaD)

for (i in tajimaIs$TajimaD){
  if (i != "NA"){
    sum = sum + i
    count= count + 1
    
  }
}
sum/count

library(ggplot2)
fisher.test(pi$PI)

ggplot(pi, aes(PI)) + geom_histogram(bins = 80, color="white") + 
  geom_vline(xintercept=3.347463e-07, color="pink", linetype="dashed") +
  geom_vline(xintercept=3.53211e-07, color="#0972B5", linetype="dashed") +
  geom_vline(xintercept=3.543648e-07, color="#EE1b24", linetype="dashed") +
  labs(x="Diversidad nucleotídica", y="SNPs") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))

  
ggplot(tajima, aes(TajimaD)) + geom_histogram(bins = 80, color="white") + 
  geom_vline(xintercept=-0.52, color="red", linetype="dashed") +
  labs(x="Tajima D", y="SNPs") +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  geom_boxplot(outlier.colour = "red") + 
  stat_boxplot(geom="errorbar") + geom_hline(xintercept = 6.22e-5, color="red") +
  labs(x="contigs", y="") + 
  ggtitle("Diversitad Nucleotídica") +
  theme( legend.position="none",
         plot.title = element_text(size=16, hjust = 0.5, face="bold"))

ggplot(tajima, aes(color=CHROM, y=TajimaD)) + geom_boxplot(outlier.colour = "black") +
  geom_hline(yintercept = 0.14, color="red") + 
  theme( legend.position="none",
         plot.title = element_text(size=16, hjust = 0.5, face="bold"))+
  labs(x="contigs", y="") + ggtitle("Indice de Tajima") 

+
  labs(x="genoma", y="") +
  geom_jitter(size=0.4, alpha=0.9) +
  ggtitle("Indice de Tajima") +
  theme( legend.position="none",
    plot.title = element_text(size=16, hjust = 0.5, face="bold"))

ggplot(tajima, aes(fill=CHROM, y=TajimaD)) + geom_boxplot()

length(stats$TajimaD)
hist(stats$TajimaD)


#....He-Ho WithReference

HeHo <- read.table("/Users/luisfuentes/Desktop/He-Ho.txt", header = T, sep = "\t")
HeHoBa <- read.table("/Users/luisfuentes/Desktop/HeHoBahi.txt", header = T, sep = "\t")
HeHoIs <-read.table("/Users/luisfuentes/Desktop/HeHoIsc.txt", header = T, sep = "\t")

HeHoIs

sum(HeHoIs$HeHo < 0.1)
sum(HeHoBa$HeHo < 0.1)
ggplot(HeHoBa, aes(HeHo)) + geom_histogram( bins = 50, fill="#0972B5") + 
  scale_fill_manual(
                    labels=c("He-Ho")) + labs(subtitle = "6904 SNPs\n 4.18% missing data") + 
  theme(plot.subtitle = element_text(hjust = 0.9, vjust=-15, face="bold", size = 12)) +
  theme(legend.background = element_rect(fill=NA),  legend.position = c(0.9,1), text = element_text(face = "bold", size = 12)) +
  labs(x="Frequency of heterozygosity", y='SNPs', fill=" ") + ylim(0,2000) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggplot(HeHoIs, aes(HeHo)) + geom_histogram( bins = 50, fill="#EE1b24") + 
  scale_fill_manual(
    labels=c("He-Ho")) + labs(subtitle = "6904 SNPs\n 6.21% missing data") + 
  theme(plot.subtitle = element_text(hjust = 0.9, vjust=-15, face="bold", size = 12)) +
  theme(legend.background = element_rect(fill=NA),  legend.position = c(0.9,1), text = element_text(face = "bold", size = 12)) +
  labs(x="Frequency of heterozygosity", y='SNPs', fill=" ") + ylim(0,2000) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

require(reshape2)
hehoAll <- cbind(Ba=HeHoBa$HeHo, Isc=HeHoIs$HeHo)
hehoAll <- melt(hehoAll, id.vars=1:2)
hehoAll
sum(hehoAll$value == 0)

HeHoplot <- ggplot(hehoAll, aes(value, fill=Var2))  + geom_histogram(position="dodge", bins=50) +
  scale_fill_manual( values=c("#0972B5","#EE1b24"),
    labels=c("Bahía Málaga", "Iscuandé")) + 
  theme(plot.subtitle = element_text(hjust = 0.9, vjust=-15, face="bold", size = 12)) +
  theme(legend.background = element_rect(fill=NA),  legend.position = c(0.85,0.7), 
        text = element_text(face = "bold", size = 24, family = "Times New Roman")) +
  labs(x="Frequency of heterozygosity", y='SNPs', fill=" ") + ylim(0,2000) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("hehoplot.tiff", HeHoplot)

ggplot(HeHoBa, aes(HeHo)) + geom_histogram(bins=100)
ggplot(HeHoIs, aes(HeHo)) + geom_histogram(bins=100)

#confirm that all samples are in the infoPopulation
file@gt
#to confirm if the names in vcf are the same with the population data file
all(colnames(fileWithR@gt)[-1] == pop.data$ID)

poppr.amova(Fgenind2WithR, pop.data$Id)

#number os ploidy
ploidy(Fgenlight) <- 2
ploidy(FgenlightTr) <- 2
ploidy(FgenlightWithR) <- 2
#To agregate of location to samples
pop(Fgenlight) <- pop.data$Loc
pop(Fgenind2) <- pop.data$Loc
pop(FgenlightTr) <- pop.data$Loc
pop(Fgenind2Tr) <- pop.data$Loc
pop(FgenlightWithR) <- pop.data$Loc
pop(Fgenind2WithR) <- pop.data$Loc
Fgenind2WithR$ploidy
Fgenind2WithR$pop


library(HardyWeinberg)
#Hardy-Weinberg test
hw <- hw.test(Fgenind2WithR, B=1000)
fisher.test(fileWithR)

fisher.test(extract.gt(fileWithR))

inb <- inbreeding(Fgenind2WithR, N=89)
mean(inb$genotype_11)
inb$genotype_11

hwvcf <- read.table("/Users/luisfuentes/Desktop/UniAndes/hweAll.hwe", header = T, sep = "\t")
hwvcf$P_HWE

ggplot(hwvcf, aes(P_HWE)) + geom_histogram(bins = 100, fill="red") +
  theme(plot.subtitle = element_text(hjust = 0.9, vjust=-15, face="bold", size = 12)) +
  xlab("p-value")


a <- pairwise.t.test(hwvcf$P_HWE, hwvcf$CHR, p.adjust.method = "bonferroni")
aHW <- p.adjust(hwvcf$P_HWE, method = "bonferroni")

hist(aHW)

install.packages("genetics")
library(genetics)

hwe(FgenlightWithR, test="fisher")
HardyWeinberg()

sum(aHW <= 0.01)
boxplot(aHW)

sum=c()
count=0
aHW
for (i in aHW){
  count = count + 1
  if (i<0.01){
    print(count)
    print(aHW[i])
    sum=c(sum, (count)) 
  }
}
sum(aHW<0.01)
range <- 1:length(hwvcf$POS)
index <- range
hwvcf <- cbind(hwvcf, index)
length(hwvcf$POS)
range <- 1:length(hwvcf$POS)
length(range)
count = 0
sum[1]
sum

position <- c()
for(i in hwvcf$index){
  for(j in sum){
    if(i == j){
      print(hwvcf$POS[i])
      position <- c(position,hwvcf$POS[i])
    }
  }
}

fileWithR

position
write.csv(position, file="posicion.txt", row.names = TRUE)

for (x in sum){
}

sum[1527]
hwvcf$index

length(hw[,3])
mean(hw[,3])
#Genetic Diversity
gl.report.hwe(FgenlightWithR)
hw2 <- gl.hwe.pop(FgenlightWithR)
hw2
#Expected heterozigosity
Hs(Fgenind2WithR)
#Observed heterozigosity
Ho(Fgenind2WithR)
#MAF
Fgenind2WithR
gl.filter.maf(Fgenind2WithR, threshold = 0.01)

#Fst
geneDif <- genet.dist(Fgenind2WithR,diploid=TRUE,method="WC84")
basic.stats(Fgenind2WithR, diploid=TRUE,method="Nei87")
geneDif

genet.dist
pairwise.WCfst
wc
geneDifdata <- data.frame(
  "pop" = c("pop2/pop1","pop3/pop1","pop4/pop1","pop5/pop1","pop6/pop1","pop3/pop2","pop4/pop2","pop5/pop2","pop6/pop2","pop4/pop3","pop5/pop3","pop5/pop3",
             "pop5/pop4","pop6/pop4","pop6/pop5"),
  "Fst" = c(6.176496e-03,4.429830e-03,4.899819e-03,3.676016e-03,5.911310e-03,8.085121e-03,4.330717e-03,7.226372e-03,5.370682e-03,9.402285e-03,1.023518e-02,
            5.541726e-03,1.431965e-03,5.216164e-05,-2.274404e-04
    
  )
)

test.between.within(FgenlightWithR, pop.data)

geneDifdata 
names(geneDifdata)
levels(geneDifdata$Fst)


##... privates Alleles

privaBM <- "/Users/luisfuentes/Desktop/privateAlleles_BM_freq.txt"
privaBMTable <- read.table(text = gsub(":", "\t", readLines(privaBM)),skip=10)

privaIsc <- "/Users/luisfuentes/Desktop/privateAlleles_Isc_freq.txt"
privaIscTable <- read.table(text = gsub(":", "\t", readLines(privaIsc)),skip=10)


privaAll <- read.csv("/Users/luisfuentes/Desktop/privateAll.txt", header = T, sep="\t")


ggplot(privaAll, aes(MAF, fill=Loc)) + geom_histogram(position="dodge", bins = 50) + 
  scale_fill_manual(values = c("#0972B5","#EE1b24")) + 
  theme(plot.subtitle = element_text(hjust = 0.9, vjust=-15, face="bold", size = 12)) + ylim(0,500) + xlim(0,0.3) +
  theme(legend.background = element_rect(fill=NA),  legend.position = c(0.9,0.57), text = element_text(face = "bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
                    


privaBMTable
privaIscTable
hist(privaIscTable[,9])
##....permutation
install.packages("PermutationR")
library(PermutationR)
library(tidyr)
library(MASS)
library(readr)

read.csv("/Users/luisfuentes/Desktop/indices_F.csv", header = T, sep = "\t")

FindexData
permuANOVA(FindexData$Loci, FindexData$Fst, reps = 1000, perm.type="restricted")



edit(geneDif)
qt.test(geneDif)
geneDif <- genet.dist(Fgenind2WithR,diploid=TRUE,method="WC84")
Fst(as.loci(Fgenind2WithR))
x <- data.frame(as.loci(Fgenind2WithR))
order(a["population"])
test.between.within(x, pop.data )

inbreeding(Fgenind2WithR)

pop_2Loc <- read.csv("/Users/luisfuentes/Desktop/infoPopulation89.txt", header = F, sep = "\t")
pop_2Loc
length((pop_2Loc))
strata(FgenlightWithR) <- pop_2Loc
strata(FgenlightWithR)
addStrata(FgenlightWithR, pop_2Loc)
FgenlightWithR
fisher.test(x=geneDif, y=NULL)

FgenlightWithRFil <- mlg.filter(FgenlightWithR, missing = "zero")

a <- poppr.amova(FgenlightWithR, hier = ~V1/V2,  clonecorrect = TRUE)
a
randtest(a, nrepet=9999)
plot(randtest(a))
a$varcomp/sum(a$varcomp)
gl.fst.pop(FgenlightWithR, nboots = 10000)

wc(FgenlightWithR)

Findex <- Fst(as.loci(Fgenind2WithR))

Fst

mean(Findex["Fit"])

write.matrix(Findex, file="indices_F.csv", sep="\t")
write_delim(FindesData, file="indices_F.csv")
FindesData <- data.frame(Findex)
Loci <- c(1:6904)
FindexData <- cbind(FindesData, Loci)
FindexData
colnames(FindesData) = c("Loci", "Fit", "Fst", "Fis")
FindesData
separate(FindesData, col=Fit, into =c("Loci","Fit"))
shapiro.
shapiro.test(Findex[,1])
kruskal.test(Findex[,1])
chisq.test(absolute2)
hist(Findex[,3])

absolute2 <- 0

for (i in absolute){
  absolute2 <- c(absolute2, i)
}
absolute2

boot.ppfis(fileWithR,nboot=100,quant=c(0.025,0.975),diploid=TRUE,dig=4)
kruskal.test(genet.dist(Fgenind2WithR,diploid=TRUE,method="Fst"))

#Nucleotide diversity
nucDi <- read.table("/Users/luisfuentes/Desktop/deNovo_ND_allSam.windowed.pi", header = T)
hist(nucDi$PI, br=20)
plot(taj$N_SNPS,taj$TajimaD, xlab='position',ylab='TajimaD')

#Allele richness
Ar <- allel.rich(Fgenind2WithR)

pA<-private_alleles(Fgenind2WithR, level = "population", count.alleles=T,
                    report="table", form= alleles ~ pop)


fileWithR
vcfPopGenome <- readVCF("/Users/luisfuentes/Desktop/withReference_89_0.01.vcf", 6904,frompos = 1, topos = 6904)
neutrality.stats(fileWithR)
vcfPopGenome

Fgenind2WithR$pop
allele.count(Fgenind2WithR)
allelic.richness(Fgenind2WithR, diploid = T, min.n = NULL)

pi.dosage(dogasedeNovo,L=NULL)
pop.freq(Fgenind2WithR)

strata(FgenindDArT) <- data.frame(x=pop(FgenindDArT))
FgenindDArT$strata
other(Fgenind2WithR)$pop

gl <- 
gl2fasta(FgenlightWithR, method = 1, outfile = "withReference.fasta", 
         outpath = "/Users/luisfuentes/Desktop", verbose = NULL)


#TajimaD
TajimaD.dosage(dogasedeNovo)
theta.Watt.dosage(dogasedeNovo)
dogasedeNovo <- fstat2dos(Fgenind2, diploid = T)
dogaseWithR <- fstat2dos(Fgenind2WithR, diploid = T)
dogaseWithR


fstWR <- Fst(as.loci(Fgenind2WithR))
mean(fstWR[3])
fstWR <- data.frame(fstWR)
fstdN <- Fst(as.loci(Fgenind2))
sum(fstdN[,2]=="NaN")
average(fstWR[,2])

obs_fst = fstWR[2]
mean(obs_fst)  

mean(fstWR$Fis)

fstBa<-Fst(as.loci(FgenindBaMa))
fstIs<-Fst(as.loci(FgenindIsc))

fstBa

sum = 0
count =0
for(i in fstBa[,2]){
  if(i != "NaN"){
    sum = sum + i
    count = count + 1
    }
}

print(sum/count)

mean(fstBa[,1])
write.csv(fstWR[,3], file="Fis.csv")

sum = 0
conteo = 0
for (i in fstWR[,2]){
  if(i<1){
    sum = sum + i
    conteo = conteo + 1
  }
}

print(sum/conteo)

loci <- rownames(fstWR)
fstWR <- cbind(fstWR, loci)
fstWR <- subset(fstWR, select= c("Fit","Fst","Fis","loci"))

levels(fstWR)
hist(fstWR[,1])
names(fstWR)
mean(fstdN)
mean(fstdN[,2],na.rm=TRUE)
colMeans(fstdN,na.rm=TRUE)

Fst(dogaseWithR)
Fst(dogaseBaMa)
fst.dosage(dogaseWithR)
fstat2dos()

#------------------------------------------------------#
#--------------------   PCA   -------------------------#
#------------------------------------------------------#

pca <- glPca(FgenlightWithR, nf=NULL)
barplot(100*pca$eig/sum(pca$eig), col = heat.colors(50), main = 'PCA Eigenvalues')
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
pca$scores
plot(pca)
sum(pca$eig, na.rm = F)
inc <- function(x) {eval.parent(substitute(x <- x + i))}
sum = 0
for (i in pca$eig){
  total = inc(sum)
  if (total < 230){
    print(pca$eig[i])
  }
}

pc.pca.scores <- as.data.frame(pca$scores)
pop(FgenlightWithR)
pc.pca.scores$pop <- pop(Fgenlight)
varPC1<- round(var(pc.pca.scores$PC1, na.rm = FALSE), digits = 3)
varPC2<- round(var(pc.pca.scores$PC2, na.rm = FALSE), digits = 3)
varPC1
summary(pc.pca.scores)
# plot for PCA
p <- ggplot(pc.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw() + scale_color_manual(name = 'Populations', values = c("blue","red","green", "black", 'brown', 'pink'))
p<- p + ggtitle('PCA for populations') + theme(plot.title = element_text(hjust = 0.5))
p<- p + ylab(paste0('PC2'," ","(",varPC2, ')')) + xlab(paste0('PC1'," ","(",varPC1, ')'))
p


#------------------------------------------------------#
#--------------------  DAPC   -------------------------#
#------------------------------------------------------#


grp <- find.clusters(Fgenind2WithR ,max.n = 10, n.pca = 50)

grp
names(grp)
grp$size
grp$grp


testDAPC <- dapc(Fgenind2WithR, pop=grp$grp, n.pca=50, n.da=2)
scatter(testDAPC, posi.da="bottomleft",	bg="white",	pch=20,	
        scree.pca=TRUE,	posi.pca="bottomright",
        col = c("brown3", "darkviolet", "green4", "skyblue4", "maroon", "grey40"),
        cstar = 1, main="asjjadjsa")

mean(testDAPC$eig)

DAPCDefPop <- dapc(Fgenind2WithR, n.pca=50, n.da=2)
scatter.dapc(DAPCDefPop, posi.da="bottomleft",	bg="white",	pch=20,	
        scree.pca=TRUE,	posi.pca="bottomright", leg=T, 
        col=c("#0972B5", "#0972B5", "#0972B5","#EE1b24", "#EE1b24", "#EE1b24"),
        #label.inds = c("Ba", "Isc"))
        txt.leg = (c("[1] Caracas (BM)","[2] Dos peñas (BM)","[3] Corozal (BM)", 
                     "[4] Puya Puyal (Is)","[5] Punta de Caimanera (Is)","[6] Los Esterones (Is)")), 
        cstar = 1, cleg = 0.7, posi.leg= "topright")
DAPCDefPop$

mean(DAPCDefPop$var.contr[,1])
compoplot(DAPCDefPop, col=c("#0972B5", "#0972B5", "#0972B5","#EE1b24", "#EE1b24", "#EE1b24"))

a <- DAPCDefPop$eig/sum(DAPCDefPop$eig)*100
a
barplot(a, names.arg = round(a,2))

ind_coords <- as.data.frame(DAPCDefPop$ind.coord)
ind_coords$Ind <- indNames(FgenlightWithR)
ind_coords$pop <-FgenlightWithR$pop
ind_c

xvalDapc(FgenlightWithR,FgenlightWithR$pop)
FgenlightWithR$pop
var(DAPCDefPop$pca.eig)
var
mean(DAPCDefPop$eig)
plot(DAPCDefPop$posterior)

DAPCDefPop$var

posterior(DAPCDefPop)
count=1
sum = 0 
for (i in DAPCDefPop$pca.eig){
  if (i >= 9.543995){
    count = count + 1
    sum = sum + i
  }
}

sum/count

"#0972B5"
"#EE1b24"
help(scatter.dapc)

txt.leg = c("Puya Puyal [4]","Punta de Caimanera [5]","Los Esterones [6]")
txt.leg = c("Caracas [1]","Dos peñas [2]","Corozal [3]")

pop=Fgenind2WithR$pop


testDAPC$grp.coord
summary(testDAPC)
testDAPC$grp
scatter(testDAPC)

#------------------------------------------------------#
#------------   tree Neighbor joining  ----------------#
#------------------------------------------------------#

treeNJ <- aboot(Fgenind2, tree="nj", distance = "nei.dist", sample = 1000)
treeNJTr <- aboot(Fgenind2Tr, tree="nj", distance = "nei.dist", sample = 1000)
treeNJWithR <- aboot(Fgenind2WithR, tree="nj", distance = "nei.dist", sample = 1000)
write.tree(treeNJ, file = "dendrogram89deNovo.txt")
write.tree(treeNJTr, file = "dendrogram89Translated.txt")
write.tree(treeNJWithR, file = "dendrogram89WithRef.txt")


treeNJ$node.label
samplesBaMa
ggtree(treeNJ, layout="equal_angle", branch.length = "none") + 
  geom_tippoint(aes(color=(label %in% samplesBaMa))) + 
  scale_color_manual(values=c("firebrick","deepskyblue4"), 
                     labels=c("Iscuande", "Bahia_Malaga"),
                     name="Locations")  + 
  ggtitle("Dendrogram de Novo ") + 
  theme(plot.title =element_text(hjust = 8, size=16,face="bold", vjust = 1))



samplesBaMa <- treeNJ$tip.label[1:45]
samplesBaMa <- c(samplesBaMa, treeNJ$tip.label[54], treeNJ$tip.label[64], treeNJ$tip.label[83],
                 treeNJ$tip.label[84])
samplesIs <- c(treeNJ$tip.label[46:53], treeNJ$tip.label[55:63], treeNJ$tip.label[65:82],
               treeNJ$tip.label[85:89])


#------------------------------------------------------#
#---- Analysis with deNovo and ReferenceGenome --------#
#------------------------------------------------------#

hist(allData[,7])

#.....deNovo

#He and Ho
ggplot(hehoDeN, aes(value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) + 
  scale_fill_manual(values = c("gray40","forestgreen"),
                    labels=c("He", "Ho")) + labs(x="Frequency of Heterozygosity", y='SNPs', fill=" ") + 
  labs(subtitle = "3806 SNPs\n 6.3% missing data")+ 
  ##labs( title= "*de Novo* genotyping") + 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
                                      size=20, hjust =0.5)) + 
  theme(plot.subtitle = element_text(hjust = 0.9, vjust=-15, face="bold", size = 12)) + ylim(0,1000) +
  theme(legend.background = element_rect(fill=NA),  legend.position = c(0.85,0.57), 
        text = element_text(face = "bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#MAF
ggplot(allData, aes(V9)) + geom_histogram(bins = 100, fill="forestgreen", color="white") +
  #labs( title= "*de Novo* genotyping") + 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
                                      size=20, hjust =0.5))  +
  labs(x="Minor Allele Frequency", y='SNPs', fill=" ") + ylim(0,300) +
  theme(text = element_text(face = "bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#........withReference 

withRef <- "/Users/luisfuentes/Desktop/withReference_89_00.1_diversity"
withRefTable <- read.table(text = gsub(":", "\t", readLines(withRef)),skip=10)
head(withRefTable)
summary (withRefTable)
hist(withRefTable[,12])


sum((p.adjust(withRefTable[,12], method = "bonferroni")) <= 0.01)

require(reshape2)
length(transTable$V6)
hehoWith <- cbind(wiHo=withRefTable$V6, wiHe=withRefTable$V7)
hehoWith <- melt(hehoWith, id.vars=1:2)
hehoWith

#He and Ho
ggplot(hehoWith, aes(value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) + 
  scale_fill_manual(values = c("gray40","goldenrod3"),
                    labels=c("He", "Ho")) + labs(x="Frequency of Heterozygosity", y='SNPs', fill=" ") + 
  labs(subtitle = "6904 SNPs\n 5.9% missing data")+ 
  #labs( title= "Genotyping with reference genome") + 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
                                      size=20, hjust =0.5)) + 
  theme(plot.subtitle = element_text(hjust = 0.9, vjust=-15, face="bold", size = 12)) + ylim(0,1500) +
  theme(legend.background = element_rect(fill=NA),  legend.position = c(0.9,0.57), text = element_text(face = "bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
#MAF

withRefTable
length(withRefTable$V9)
sum(withRefTable$V9 < 0.1)

MAFWR <- ggplot(withRefTable, aes(V9)) + geom_histogram(bins = 60, fill="forestgreen", color="white")  +
  #labs( title= "Genotyping with reference genome") + 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
                                      size=24, hjust =0.5))  +
  labs(x="Minor Allele Frequency", y='SNPs', fill=" ") + ylim(0,800) + 
  theme(text = element_text(face = "bold", size = 26, family = "Times New Roman")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("MAFWR.tiff", MAFWR)

#.....translated 
trans <- "/Users/luisfuentes/Desktop/tranlationdeNovoGBS_89_0.01_bi_Diversity"
transTable <- read.table(text = gsub(":", "\t", readLines(trans)),skip=10)
head(transTable)
summary (transTable)
hist(transTable[,7])

require(reshape2)
length(transTable$V6)
hehoTr <- cbind(trHo=transTable$V7, trHe=transTable$V6)
hehoTr <- melt(hehoTr, id.vars=1:2)
hehoTr

#He and Ho
ggplot(hehoTr, aes(value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) + 
  scale_fill_manual(values = c("grey28", "deepskyblue4"),
                    labels=c("Ho", "He")) + labs(x="Frequency", y='', fill=" ") + ggtitle('Translated deNovo GBS') +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1000)

#MAF

ggplot(transTable, aes(V9)) + geom_histogram(bins = 150, fill="blue4") +
  ggtitle('Translated deNovo GBS - MAF') + 
  theme(plot.title = element_text(hjust = 0.5)) + labs(x="Frequency", y='', fill=" ") +
  ylim(0,400)


"#0972B5"
"#EE1b24"
#------------------------------------------------------#
#-----------------  Bahia Malaga    -------------------#
#------------------------------------------------------#

#to load VCF file
bahiaMal <- read.vcfR("/Users/luisfuentes/Desktop/withReference_89_0.01_Bahia.vcf", verbose = FALSE)

#to convert loaded file to genind object
FgenindBaMa <- vcfR2genind(bahiaMal, pop=NULL)
FgenindBaMa

#to convert loaded file to genlight object
FgenlightBaMa <- vcfR2genlight(bahiaMal, n.cores = 1)
FgenlightBaMa


allelic.richness(FgenindBaMa)
pairwise.betas(FgenindBaMa, diploid=TRUE)
fst.dosage(FgenindBaMa, pop=FgenindBaMa[,1])

#pop data
pop.dataBaMa <- read.table("/Users/luisfuentes/Desktop/infoPopulationsBaMa.txt", sep="\t", header = T)
pop.dataBaMa$Loc

#confirm that all samples are in the infoPopulation
bahiaMal@gt

#to confirm if the names in vcf are the same with the population data file
all(colnames(bahiaMal@gt)[-1] == pop.dataBaMa$ID)

  
#number os ploidy
ploidy(FgenlightBaMa) <- 2

#To agregate of location to samples
pop(FgenlightBaMa) <- pop.dataBaMa$Loc
pop(FgenindBaMa) <- pop.dataBaMa$Loc

FgenlightBaMa$pop
FgenindBaMa$pop


library(HardyWeinberg)
#Hardy-Weinberg test
hw.test(Fgenind2WithR, B=1000)

#Genetic Diversity

#Expected heterozigosity
Hs(FgenindBaMa)
#Observed heterozigosity
Ho(FgenindBaMa)
#MAF
Fgenind2WithR
gl.filter.maf(Fgenind2WithR, threshold = 0.01)

#Fst

genet.dist(FgenindBaMa,diploid=TRUE,method="Nei87")
basic.stats(FgenindBaMa, diploid=TRUE,digits=4)

genet.dist(Fgenind2WithR,diploid=TRUE,method="Nei87")
allele.count(Fgenind2WithR)

genet.dist(FgenindBaMa,diploid=TRUE,method="WC84")
Fst(as.loci(FgenindBaMa))

allelic.richness(FgenindBaMa, diploid = T, min.n = NULL)

# Dogase data
dogaseBaMa <- fstat2dos(FgenindBaMa, diploid = T)
dogaseBaMa

#TajimaD
TajimaD.dosage(dogaseBaMa)

#theta Watterson
theta.Watt.dosage(dogaseBaMa)

#Pi - nucleotipe diversity
pi.dosage(dogaseBaMa,L=NULL)

pop.freq(FgenindBaMa)

fstBM <- Fst(as.loci(FgenindBaMa))
fstBM[,2]
mean(fstWR[,2])
colMeans(fstBM, na.rm=TRUE)
##......plots......

BaMadiversity <- "/Users/luisfuentes/Desktop/withReference_Bahia_diversity"
BaMadiversity <- read.table(text = gsub(":", "\t", readLines(BaMadiversity)),skip=10)
head(BaMadiversity)
summary (BaMadiversity)
hist(BaMadiversity[,9])
histo()
sum(BaMadiversity[,9]==0)
require(reshape2)
length(transTable$V6)
hehoBaMa <- cbind(wiHo=BaMadiversity$V6, wiHe=BaMadiversity$V7)
hehoBaMa <- melt(hehoBaMa, id.vars=1:2)
hehoBaMa
scale_fill_manual(values = c("#1a5f77","goldenrod3"),
                  labels=c("He", "Ho")) + labs(x="Frequency of Heterozygosity", y='SNPs', fill=" ") + 
  labs(subtitle = "3806 SNPs\n 6.3% missing data")

#He and Ho
HeHoBm <- ggplot(hehoBaMa, aes(value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) + 
  scale_fill_manual(values = c("gray60","#0972B5"),
                    labels=c("He", "Ho")) + #labs(subtitle = "6904 SNPs\n 4.18% missing data") + 
  theme(plot.subtitle = element_text(hjust = 0.9, vjust=-15, face="bold", size = 18)) +
  #labs( title= "Bahía Málaga") + 
  #theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
   #                                   size=20, hjust =0.5)) + ylim(0,1500) +
  theme(legend.background = element_rect(fill=NA),  legend.position = "none", text = element_text(face = "bold", size = 18, family = "Times New Roman")) +
  labs(x="Frequency of heterozygosity", y='SNPs', fill=" ") + ylim(0,2000) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

HeHoBm

ggsave("HeHoBM.tiff", HeHoBm)

#MAF

ggplot(BaMadiversity, aes(V9)) + geom_histogram(bins = 150, fill="skyblue4") +
  labs( title= "Bahía Málaga") + 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
                                      size=20, hjust =0.5))  +
  labs(x="Frequency - MAF", y='', fill=" ") + ylim(0,500) + 
  theme(text = element_text(face = "bold", size = 12))
  


#----  DAPC   --

grpBM <- find.clusters(FgenindBaMa ,max.n = 10, n.pca = 100)

grpBM
names(grp)
grpBM$size
grpBM$grp


DAPCBM <- dapc(FgenindBaMa, pop=grpBM$grp, n.pca=20, n.da=2)
scatter.dapc(DAPCBM, posi.da="bottomleft",	bg="white",	pch=20,	
             scree.pca=TRUE,	posi.pca="bottomright",
             col = c("deepskyblue1", "deepskyblue3", "deepskyblue4"),
             cstar = 1, main="asjjadjsa")
deepskyblue4

col = c("#0000FF", "#4949FF", "#A3A3FF"),

DAPCDefBM <- dapc(FgenindBaMa, pop=FgenindBaMa$pop, n.pca=30, n.da=2)
scatter(DAPCDefBM, posi.da="bottomleft",	bg="white",	pch=20,	
        scree.pca=TRUE,	posi.pca="bottomright", leg=T, 
        col = c("deepskyblue1", "deepskyblue3", "deepskyblue4"),
        txt.leg = c("Caracas [1]","Dos peñas [2]","Corozal [3]"), cstar = 1)



#------------------------------------------------------#
#-------------------    Iscuande    -------------------#
#------------------------------------------------------#

#to load VCF file
iscua <- read.vcfR("/Users/luisfuentes/Desktop/withReference_89_0.01_Isc.vcf", verbose = FALSE)

#to convert loaded file to genind object
FgenindIsc <- vcfR2genind(iscua, pop=NULL)
FgenindIsc

#to convert loaded file to genlight object
FgenlightIsc <- vcfR2genlight(iscua, n.cores = 1)
FgenlightIsc


allelic.richness(FgenindIsc)
pairwise.betas(FgenindIsc, diploid=TRUE)
fst.dosage(FgenindIsc, pop=FgenindIsc[,1])

#pop data
pop.dataIsc <- read.table("/Users/luisfuentes/Desktop/infoPopulationsIsc.txt", sep="\t", header = T)
pop.dataIsc$Loc


#confirm that all samples are in the infoPopulation
iscua@gt

#to confirm if the names in vcf are the same with the population data file
all(colnames(iscua@gt)[-1] == pop.dataIsc$ID)


#number os ploidy
ploidy(FgenlightIsc) <- 2

#To agregate of location to samples
pop(FgenlightIsc) <- pop.dataIsc$Loc
pop(FgenindIsc) <- pop.dataIsc$Loc

FgenlightIsc$pop
FgenindIsc$pop


library(HardyWeinberg)
#Hardy-Weinberg test
hw.test(FgenindIsc, B=1000)

#Genetic Diversity

#Expected heterozigosity
Hs(FgenindIsc)
#Observed heterozigosity
Ho(FgenindIsc)

#MAF
gl.filter.maf(FgenlightIsc, threshold = 0.01)

#Fst
genet.dist(FgenindIsc,diploid=TRUE,method="Fst")
basic.stats(FgenindIsc, diploid=TRUE,digits=4)


allele.count(Fgenind2WithR)

allelic.richness(FgenindIsc, diploid = T, min.n = NULL)

# Dogase data
dogaseIsc <- fstat2dos(FgenindIsc, diploid = T)
dogaseBaMa
dogaseIsc
#TajimaD
TajimaD.dosage(dogaseIsc)

#theta Watterson
theta.Watt.dosage(dogaseIsc)

#Pi - nucleotipe diversity
pi.dosage(dogaseIsc,L=NULL)

pop.freq(FgenindIsc)

fstIs <- Fst(as.loci(FgenindIsc))
colMeans(fstIs, na.rm=TRUE)

##......plots......

Iscdiversity <- "/Users/luisfuentes/Desktop/withReference_89_0.01_Isc_diversity"
Iscdiversity <- read.table(text = gsub(":", "\t", readLines(Iscdiversity)),skip=10)
head(Iscdiversity)
sum(Iscdiversity[,9]==0)
summary (Iscdiversity)
hist(Iscdiversity[,7])

require(reshape2)
length(transTable$V6)
hehoIsc <- cbind(wiHo=Iscdiversity$V6, wiHe=Iscdiversity$V7)
hehoIsc <- melt(hehoIsc, id.vars=1:2)
hehoIsc



#He and Ho
ggplot(hehoIsc, aes(value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) + 
  scale_fill_manual(values = c("gray60","#EE1b24"),
                    labels=c("He", "Ho")) + labs(subtitle = "6904 SNPs\n 6.21% missing data") + 
  theme(plot.subtitle = element_text(hjust = 0.9, vjust=-15, face="bold", size = 12)) +
  #labs( title= "Iscuandé") + 
  #theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
   #                                   size=20, hjust =0.5)) + ylim(0,1500) +
  theme(legend.background = element_rect(fill=NA),  legend.position = c(0.9,0.57), text = element_text(face = "bold", size = 12)) +
  labs(x="Frequency of Heterozygosity", y='SNPs', fill=" ") + ylim(0,2000) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

HeHoI <- ggplot(hehoIsc, aes(value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) + 
  scale_fill_manual(values = c("gray60","#EE1b24"),
                    labels=c("He", "Ho")) + 
  theme(plot.subtitle = element_text(hjust = 0.9, vjust=-15, face="bold", size = 18)) +
  theme(legend.background = element_rect(fill=NA),  legend.position = "none", text = element_text(face = "bold", size = 18, family = "Times New Roman")) +
  labs(x="Frequency of Heterozygosity", y='SNPs', fill=" ") + ylim(0,2000) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

HeHoI
ggsave("HeHoIs.tiff", HeHoI)

#MAF

ggplot(Iscdiversity, aes(V9)) + geom_histogram(bins = 150, fill="skyblue4") +
  labs( title= "Iscuandé") + 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
                                      size=20, hjust =0.5))  +
  labs(x="Frequency - MAF", y='', fill=" ") + ylim(0,700) + 
  theme(text = element_text(face = "bold", size = 12))
Iscdiversity


#MAF

MAF<- ggplot(MAFLoc, aes(x=value, fill=Var2)) + geom_histogram(bins = 100, position="dodge") +
  #labs( title= "Bahía Málaga") + 
  scale_fill_manual(values = c("#0972B5", "#EE1b24"), labels=c("","")) +
  labs(x="Minor Allele Frequency", y='SNPs', fill=" ") + ylim(0,800) + 
  theme(text = element_text(face = "bold", size = 20, family = "Times New Roman")) + theme(legend.position = 'none') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold")) 


MAF

ggsave("MAF.tiff", MAF, device='tiff')


ggplot(df, aes(x=value, fill=Var2)) + geom_histogram(bins = 80, position="dodge") +
  scale_fill_manual(values = c("grey28", "deepskyblue4", "firebrick")) +
  labs(x="Samples", y='', fill=" ") + ggtitle('Genotyped Samples') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 1) ) + 
  labs(subtitle = "3806 SNPs")

#V6 column - He

MAFLoc<- cbind(BaMa=BaMadiversity$V9, Isc=Iscdiversity$V9)
MAFLoc <- melt(MAFLoc, id.vars=1:2)
length(MAFLoc$value)

for (i in MAFLoc$value){
  if(i == 0){
    suma = suma + 1
    print(sum(suma))
  }
}
if (MAFLoc$value == 0){
  "si"
}
  
#----  DAPC   --

grpIs <- find.clusters(FgenindIsc ,max.n = 10)
grpIs$size
grpIs$grp
grpIs

DAPCIs <- dapc(FgenindIsc, pop=grpIs$grp, n.pca=20, n.da=2)
scatter.dapc(DAPCIs, posi.da="bottomleft",	bg="white",	pch=20,	
             scree.pca=TRUE,	posi.pca="bottomright",
             col = c("#DC1C13", "#F07470", "#F6BDC0"),
             cstar = 1, inset.pca = 0.02)


DAPCDefIs <- dapc(FgenindIsc, pop=FgenindIsc$pop, n.pca=30, n.da=2)
scatter(DAPCDefIs, posi.da="bottomleft",	bg="white",	pch=20,	
        scree.pca=TRUE,	posi.pca="bottomright", leg=T, 
        col=c("#DC1C13", "#F07470", "#F6BDC0"),
        txt.leg = c("Puya Puyal [4]","Punta de Caimanera [5]","Los Esterones [6]"), cstar = 1)


#------------------------------------------------------#
#------------ Plots for He, Ho and MAF ----------------#
#------------------------------------------------------#

#all samples
allSamp <- "/Users/luisfuentes/Desktop/filtered_bi_89_diversity"
allData <- read.table(text = gsub(":", "\t", readLines(allSamp)),skip=10)
head(allData)
summary (allData)

require(reshape2)
hehoDeN<- cbind(heDN=allData$V6, hoDN=allData$V7 )
hehoDeN <- melt(hehoDeN, id.vars=1:2)
hehoDeN

#Population1 - Bahia Malaga
pop1 <- "/Users/luisfuentes/Desktop/filtered_bi_89_BaMa_diversity"
pop1Data <- read.table(text = gsub(":", "\t", readLines(pop1)),skip=10) 
head(pop1Data)
summary (pop1Data)
hist(pop1Data[,7])

#Population2-Iscuande
pop2 <- "/Users/luisfuentes/Desktop/filtered_bi_89_Is_diversity"
pop2Data <- read.table(text = gsub(":", "\t", readLines(pop2)),skip=10) 
head(pop2Data)
summary (pop2Data)
hist(pop2Data[,7])


corr <- cor.test(allData$V6, allData$V7)
corr
ggplot(allData, aes(V7, V6)) + geom_point()


#dataframe for:
#V5 column - genotyped samples
require(reshape2)
df<- cbind(allSamples=allData$V5, Bahía_Málaga=pop1Data$V5, Iscuandé=pop2Data$V5)
df <- melt(df, id.vars=1:2)

#V6 column - He
dfhe<- cbind(allSamples=allData$V6, Bahía_Málaga=pop1Data$V6, Iscuandé=pop2Data$V6)
dfhe <- melt(dfhe, id.vars=1:2)

#V7 column - Ho
dfho<- cbind(allSamples=allData$V7, Bahía_Málaga=pop1Data$V7, Iscuandé=pop2Data$V7)
dfho<- melt(dfho, id.vars=1:2)

#V9 column - MAF
maf<- cbind(allSamples=allData$V9, Bahía_Málaga=pop1Data$V9, Iscuandé=pop2Data$V9)
maf<- melt(maf, id.vars=1:2)

# plot for genotyped samples 

all <- ggplot(df, aes(x=value, fill=Var2)) + geom_histogram(bins = 80, colour='black', position="dodge") +
  scale_fill_manual(values = c("grey28", "deepskyblue4", "firebrick")) +
  labs(x="Samples", y='', fill=" ") + ggtitle('Genotyped Samples') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 1) ) + 
  labs(subtitle = "3806 SNPs")
all

#plot for Expected Heterozygosity

He<- ggplot(dfhe, aes(x=value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) +
  scale_fill_manual(values = c("grey28", "deepskyblue4", "firebrick")) +
  labs(x="Frequency", y='', fill=" ") + ggtitle('Expected Heterozygosity') +
  theme(plot.title = element_text(hjust = 0.5))
He

#plot for Observated Heterozygosity

Ho<-ggplot(dfho, aes(x=value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) +
  scale_fill_manual(values = c("grey28", "deepskyblue4", "firebrick")) +
  labs(x="Frequency", y='', fill=" ") + ggtitle('Observed Heterozygosity') +
  theme(plot.title = element_text(hjust = 0.5))
Ho

#plot for MAF

pmaf<-ggplot(maf, aes(x=value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) +
  scale_fill_manual(values = c("grey28", "deepskyblue4", "firebrick")) +
  labs(x="Frequency", y='', fill=" ") + ggtitle('Minor Allele Frequency') +
  theme(plot.title = element_text(hjust = 0.5))
pmaf



require(reshape2)
summary(allData)
summary(transTable)
summary(withRefTable)

max_length <- max(length(allData$V6), length(transTable$V6), length(withRefTable$V6))
max_length

#V6 column - He

cbind
rbind()
heSNPs <- cbind(deNovo=allData$V6,tr=transTable$V6, wR=withRefTable$V6)
heSNPs <- melt(heSNPs, id.vars=1:2)

#V7 column - Ho
hoSNPs<- cbind(deNovo=allData$V7, tr=transTable$V7, wR=withRefTable$V7)
dart7 <- melt(dart7, id.vars=1:2)

#V9 column - MAF
dart9<- cbind(allSamples=allData$V9, DArT=DArTData$V10)
dart9 <- melt(dart9, id.vars=1:2)


#------------------------------------------------------#
#------------ Analysis with DArT data ----------------#
#------------------------------------------------------#

# dendrogram 
DArTree <- read.tree('/Users/luisfuentes/Desktop/DArT_0.01_80_bi_NJ')
DArTree$tip.label
DArTplot <- ggtree(DArTree, layout="equal_angle", branch.length = 'none') + 
  geom_tippoint(aes(color=(label %in% names))) + 
  scale_color_manual(values=c("firebrick","deepskyblue4"), 
                     labels=c("Iscuande", "Bahia_Malaga"),
                     name="Locations") + 
  ggtitle("Dendrogram of DArT data") + theme(text=element_text(size=14))
DArTplot

#DArT samples
DArTSamp <- "/Users/luisfuentes/Desktop/DArT_0.01_80_bi_Diversity"
DArTData <- read.table(text = gsub(":", "\t", readLines(DArTSamp)),skip=10)
DArTData
hist(DArTData[,10])
#dataframe for:

#V5 column - genotyped samples
require(reshape2)
dart5<- cbind(He=DArTData$V7, Ho=DArTData$V8)
dart5 <- melt(dart5, id.vars=1:2)

#He and Ho
ggplot(dart5, aes(value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) + 
  scale_fill_manual(values = c("gray40","#FF66FF"),
                    labels=c("He", "Ho")) + 
  labs(subtitle = "1953 SNPs \n4.14% missing data")+
  theme(plot.subtitle = element_text(vjust= -15, face="bold", 
                                     size=12, hjust = 0.9)) +
  ##labs( title= "Frequency of Heterozygosity") + 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, -15, 0), face="bold", 
                                      size=20, hjust =0.5)) + ylim(0,1000) +
  theme(legend.background = element_rect(fill=NA),  legend.position = c(0.9,0.57), text = element_text(face = "bold", size = 12)) +
  labs(x="Frequency of Heterozygosity", y='SNPs', fill=" ") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#MAF

ggplot(DArTData, aes(V10)) + geom_histogram(bins = 100, fill="#FF66FF", color="white") +
 # labs( title= "DArT genotyping") + 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
                                      size=20, hjust =0.5))  +
  labs(x="Minor Allele Frequency", y='SNPs', fill=" ") + ylim(0,300) + 
  theme(text = element_text(face = "bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



require(reshape2)
mafNGDA <- read.table("/Users/luisfuentes/Desktop/MAF.csv", sep=";", header = T)
mafNGDA

ggplot(mafNGDA, aes(x=Value, fill=Var)) + geom_histogram(bins = 100, position = 'dodge') +
  scale_fill_manual(values = c("deeppink4","forestgreen"),
                    labels=c("DArT", "NGSEP")) +
  #labs( title= "*de Novo* genotyping") + 
  #theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
   #                                   size=20, hjust =0.5))  +
  labs(x="Minor Allele Frequency", y='SNPs', fill=" ") + ylim(0,300) +
  theme(text = element_text(face = "bold", size = 12)) +
  theme(legend.background = element_rect(fill=NA),  legend.position = c(0.9,0.57))

#......de Novo NGSEP .......
#He and Ho
ggplot(hehoDeN, aes(value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) + 
  scale_fill_manual(values = c("gray40","forestgreen"),
                    labels=c("He", "Ho")) + labs(x="Frequency", y='', fill=" ") + 
  labs(subtitle = "3806 SNPs \n6,3% missing data")+ 
  labs( title= "*de Novo* NGSEP genotyping") + 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, -15, 0), face="bold", 
                                      size=20, hjust =0.5)) + 
  theme(plot.subtitle = element_text(hjust = 1, vjust=-15, face="bold", size = 10)) + ylim(0,1000) +
  theme(legend.background = element_rect(fill=NA),  legend.position = c(0.85,0.57), text = element_text(face = "bold", size = 12))

#MAF
ggplot(allData, aes(V9)) + geom_histogram(bins = 150, fill="skyblue4") +
  labs( title= "*de Novo* NGSEP genotyping") + 
  theme(plot.title = element_markdown(margin = ggplot2::margin(5, 0, 5, 0), face="bold", 
                                      size=20, hjust =0.5))  +
  labs(x="Frequency - MAF", y='', fill=" ") + ylim(0,300) +
  theme(text = element_text(face = "bold", size = 12))

#......genetics Analysis........


#to load VCF file
DArTVCF <- read.vcfR("/Users/luisfuentes/Desktop/DArT_0.01_80_bi.vcf", verbose = FALSE)

#to convert loaded file to genind object
FgenindDArT <- vcfR2genind(DArTVCF, pop=NULL)
FgenindDArT

#to convert loaded file to genlight object
FgenlightDArT <- vcfR2genlight(DArTVCF, n.cores = 1)
FgenlightDArT$ind.names


allelic.richness(FgenindIsc)
pairwise.betas(FgenindIsc, diploid=TRUE)
fst.dosage(FgenindIsc, pop=FgenindIsc[,1])

#pop data
pop.dataDArT <- read.table("/Users/luisfuentes/Desktop/infoPopulationsDArT.txt", sep="\t", header = T)
pop.dataDArT$Loc


#confirm that all samples are in the infoPopulation
DArTVCF@gt

#to confirm if the names in vcf are the same with the population data file
all(colnames(DArTVCF@gt)[-1] == pop.dataIsc$ID)


#number os ploidy
ploidy(FgenlightIsc) <- 2

#To agregate of location to samples
pop(FgenlightDArT) <- pop.dataDArT$Loc
pop(FgenindDArT) <- pop.dataDArT$Loc

FgenlightDArT$pop
FgenindDArT$pop


library(HardyWeinberg)
#Hardy-Weinberg test
hw.test(FgenindIsc, B=1000)

#Genetic Diversity

#Expected heterozigosity
Hs(FgenindIsc)
#Observed heterozigosity
Ho(FgenindIsc)

#MAF
gl.filter.maf(FgenlightIsc, threshold = 0.01)

#Fst
genet.dist(FgenindIsc,diploid=TRUE,method="Fst")
basic.stats(FgenindDArT, diploid=TRUE,digits=4)


allele.count(Fgenind2WithR)

allelic.richness(FgenindIsc, diploid = T, min.n = NULL)

# Dogase data
dogaseDArT <- fstat2dos(FgenindDArT, diploid = T)
dogaseDArT

#TajimaD
TajimaD.dosage(dogaseDArT)

#theta Watterson
theta.Watt.dosage(dogaseDArT)

#Pi - nucleotipe diversity
pi.dosage(dogaseDArT,L=NULL)

pop.freq(FgenindIsc)

fstDArT <- Fst(as.loci(FgenindDArT))
colMeans(fstDArT, na.rm=TRUE)



#V6 column - He
dart6<- cbind(allSamples=allData$V6, DArT=DArTData$V7)
dart6 <- melt(dart6, id.vars=1:2)

#V7 column - Ho
dart7<- cbind(allSamples=allData$V7, DArT=DArTData$V8)
dart7 <- melt(dart7, id.vars=1:2)

#V9 column - MAF
dart9<- cbind(allSamples=allData$V9, DArT=DArTData$V10)
dart9 <- melt(dart9, id.vars=1:2)

# plot for genotyped samples 

DArTall <- ggplot(dart5, aes(x=value, fill=Var2)) + geom_histogram(bins = 80, colour='black') +
  scale_fill_manual(values = c("grey28", "deepskyblue4"),
                    labels=c("NESEP", "DArT")) +
  labs(x="Samples", y='', fill=" ") + ggtitle('Genotyped Samples') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle =element_text(hjust = 1) ) + 
  labs(subtitle = "1903 SNPs")
DArTall

#plot for Expected Heterozygosity

DArTHe<- ggplot(dart6, aes(x=value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) +
  scale_fill_manual(values = c("grey28", "deepskyblue4"),
                    labels=c("NESEP", "DArT")) +
  labs(x="Frequency", y='', fill=" ") + ggtitle('Expected Heterozygosity') +
  theme(plot.title = element_text(hjust = 0.5))
DArTHe

#plot for Observed Heterozygosity

DArTHo<-ggplot(dart7, aes(x=value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) +
  scale_fill_manual(values = c("grey28", "deepskyblue4"),
                    labels=c("NESEP", "DArT")) +
  labs(x="Frequency", y='', fill=" ") + ggtitle('Observed Heterozygosity') +
  theme(plot.title = element_text(hjust = 0.5))
DArTHo

#plot for MAF

DArTmaf<-ggplot(dart9, aes(x=value, fill=Var2)) + geom_histogram(position="dodge", bins = 50) +
  scale_fill_manual(values = c("grey28", "deepskyblue4"),
                    labels=c("NESEP", "DArT")) +
  labs(x="Frequency", y='', fill=" ") + ggtitle('Minor Allele Frequency') +
  theme(plot.title = element_text(hjust = 0.5))
DArTmaf




# 
# #DAPC
# dapc <- dapc(Fgenind2, pop=kfind$grp, center=T, scale=F, var.contrib=T, pca.info=T, truenames=T)
# dapc$grp
# Fgenind2$tab
# names(Fgenlight)
# x <- Fgenlight$gen
# x
# dapc1 <- dapc(Fgenlight)
# scatter(dapc$grp)
# scatter(dapc, scree.da=F, bg="white", pch=20, cell=0, cstar = 0,
#         col = c("darkblue","purple"), solid=.4, cex=3, clab=0, leg=T,
#         txt.leg = paste("Cluster", 1:6))
# 
# 
# 
# scatter(dapc, xax=1, yax=1, main="DAPC", bg="white", solid=0.5,
#         leg=T, txt.leg=c("Ba", "Is"), posi.leg="topright")
# assignplot(dapc)
# compoplot(dapc, xlab="individuals", leg=F)
# loadingplot(dapc$var.contr)
# barplot(t(dapc$posterior), las = 3, space = 0, border = NA, cex.names = 0.6)
# title(ylab = "Membership probability")
# 
# head(dapc.results, n=89)
# dapc$posterior
# 
# names(fileAde)
# dcpa1 <- dapc(fileAde$tab, dcpa$grp)
# scatter(dcpa)
# table(pop(dcpa), grp)
# 
# 
# co <- tab(Fgenind2WithR, freq=T, NA.method="mean")
# pco<- dudi.pco(dist(co), scannf = F, nf=3)
# pco
# scatter(pco, posi="bottomright", sub = NULL )
# cor(pco$li, pco$li)^2
# pco$li
# 
# install.packages("factoextra")# PCA visualization
# library("factoextra")
# 
# fviz_pca_ind(pco, col.ind = "cos2")
# ggplot(pco, aes(x=A1, y=A2)) + geom_point() + coord_fixed()
