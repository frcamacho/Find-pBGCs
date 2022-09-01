#####
#Summaries
#####
library(tidyverse)
library(scatterplot3d)
library(RColorBrewer)
library(stringr)
library(ggsci)

setwd('/Users/francine/Documents/Projects/find-pbgcs')
rm(list=ls())


writePDF=F

#read the log file and fix names
# all_log=read.table("all_log.tsv", sep="\t", header=T, comment.char = "")
all_log=read_tsv("all_log.tsv", col_names=T, comment = "")
names=gsub("#","",stringr::str_split_fixed(readLines("all_log.tsv",  n = 1), "\t",6))
colnames(all_log)=c("ID","genome_anti_hits","prophage_hits", "pBGC_hits", "pBGC_types","pBGC_coordinates","genome_length","tax_name")

#make species and genus names
all_log$species=apply(str_split_fixed(all_log$tax_name, " ", 8)[,2:3],MARGIN = 1, function(x) paste(x, collapse = " "))
all_log$genus=(gsub("\\[|\\]","",str_split_fixed(all_log$tax_name, " ", 8)[,2]))


#SUMMARIES
#%-age of genomes with prophage
100*length(which(all_log$prophage_hits>0))/NROW(all_log) #80.92214

#%-age of genomes without prophage
100-(100*length(which(all_log$prophage_hits>0))/NROW(all_log)) #19.07786

#%-age with pBGC
100*length(which(all_log$pBGC_hits>0))/NROW(all_log) #2.586306

#%-age of phage positives with pBGC
100*length(which(all_log$pBGC_hits>0))/length(which(all_log$prophage_hits>0)) #3.196043


##########
#SUBSETTING TO pBGC-positive
pBGC_indx=which(all_log$pBGC_hits>0)
only_pBGC=all_log[pBGC_indx,]

#WHAT ARE THE TYPES
sort(table(only_pBGC$pBGC_types))

#IN WHICH ARE pBGCs A MAJOR SOURCE OF BGC
sort(table(only_pBGC[only_pBGC$genome_anti_hits<2,]$species))

#CHECKING LISTERIA
# list_mono=subset(all_log, subset = species=="Listeria monocytogenes")
# plot(list_mono$pBGC_hits~jitter(list_mono$genome_anti_hits))

#CHECHING ENTEROCOCCUS
# ent_fae=subset(all_log, subset = species=="Enterococcus faecalis")
# plot(ent_fae$pBGC_hits~jitter(ent_fae$genome_anti_hits))

#CHECKING OTHERS
# subset(all_log, subset = pBGC_types=="other")
#
# mann=subset(all_log, subset = genus=="Mannheimia")
# sum(mann$prophage_hits)
#
# write.table(mann,file = "mannheimia_table.txt", row.names = F)

##########

#WHICH GENOMES ARE OVER REPRESENTED?
sort(table(all_log$species), decreasing = T)[1:50]

#WHICH pBGCs carriers ARE OVERREPRESENTED
sort(table(all_log$species[pBGC_indx]), decreasing = T)

all_log %>% filter(pBGC_hits > 0) %>% group_by(genus) %>% summarise(count = n()) %>%
  ggplot(.,aes(x = reorder(genus,(count)), y = count))+ geom_bar(stat="identity")+ theme_bw()+ coord_flip()+
  xlab("Genera")+ ylab('Species count with pBGC') +
ggtitle("Which genera are overrepresented pBGCs carriers?")
ggsave("pbgc_genera_count.png", width = 11, height = 8, units = "in")

# all_log %>% filter(species == 'Bacillus subtilis') %>% nrow()
# all_log %>% filter(genus == 'Bacillus') %>% nrow()
# all_log %>% filter(pBGC_hits > 0) %>% filter(species == 'Bacillus subtilis' ) %>% nrow()
# all_log %>% filter(pBGC_hits > 0) %>% filter(genus == "Bacillus") %>% nrow()

# all_log %>% filter(species == 'Enterococcus faecalis') %>% nrow()
# all_log %>% filter(genus == 'Enterococcus') %>% nrow()
# all_log %>% filter(pBGC_hits > 0) %>% filter(species == 'Enterococcus faecalis' ) %>% nrow()
# all_log %>% filter(pBGC_hits > 0) %>% filter(genus == "Enterococcus") %>% nrow()


#WHICH HAVE MULTIPLE pBGCs
table(all_log[which(all_log$pBGC_hits>1),]$species)
all_log %>% filter(pBGC_hits > 1) %>% select(ID, pBGC_hits, genus, species) %>%
  arrange(desc(pBGC_hits)) %>% write_tsv(., "genomes_morethanone_pBGC.tsv")



topNames=names(sort(table(all_log$genus[which(all_log$pBGC_hits>0)]), decreasing = T)[1:11])

all_log$mainNames=as.character(all_log$genus)
all_log$mainNames[which(is.na(match(all_log$mainNames, topNames)))]="others"
all_log$mainNames=factor(all_log$mainNames, levels=c(topNames,"others"))

cols=brewer.pal(12, "Paired")
all_log$cols=cols[as.numeric(all_log$mainNames)]

#Distribution of pBGC across genera
ggplot(all_log %>% filter(pBGC_hits  > 0) %>%group_by(genus,cols, mainNames) %>% tally() %>%
         mutate(prop = n / nrow(all_log %>% filter(pBGC_hits  > 0))), aes(x = "Prophages", y =prop, fill = mainNames)) +
  geom_col() + xlab('NCBI') + ylab('Fraction of pBGC hosts') +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual('Genera', values = cols) + theme_bw()
ggsave("distribution_bacteria_w_pBGC.png", width =4.06, height = 3.5, units = "in")

#gBGCs more abundant in pBGC hosts?
plot(jitter(c(rep(1, 672),rep(2,25311) ) ), (c(all_log$genome_anti_hits[pBGC_indx],all_log$genome_anti_hits[-pBGC_indx])))
t.test(log10(all_log$genome_anti_hits[pBGC_indx]+1),log10(all_log$genome_anti_hits[-pBGC_indx]+1))
wilcox.test(all_log$genome_anti_hits[pBGC_indx],all_log$genome_anti_hits[-pBGC_indx], paired = F)


GLM1=glm(c(all_log$genome_anti_hits[pBGC_indx],all_log$genome_anti_hits[-pBGC_indx]) ~ (c(rep("A", 672),rep("B",25311) ) ))
summary(GLM1)
plot(GLM1)

#prophages more abundant in pBGC hosts?
boxplot(list(all_log$prophage_hits[pBGC_indx],all_log$prophage_hits[-pBGC_indx]))
t.test(all_log$prophage_hits[pBGC_indx],all_log$prophage_hits[-pBGC_indx])
wilcox.test(all_log$prophage_hits[pBGC_indx],all_log$prophage_hits[-pBGC_indx], paired = F)


if(writePDF) pdf("pBGCvsPHAGEvsgBGC.pdf", onefile = T)

########
#Genomic BGCs as a function of phages
########

png("pBGCvsPHAGEvsgBGC.png", width=7.4, height = 4, units = "in", res=600)
par(mai=c(.5,0.5,0.1,0.1), mgp=c(1.4,0.4,0), font.lab=2, family="serif")
with(all_log[-pBGC_indx,],
     plot( jitter((genome_anti_hits),amount = 1) ~ jitter((prophage_hits), amount=1 ) ,cex= 0.5,
           xlab="#Prophages", ylab="#gBGCs" , xlim=c(-1,50),
           pch=21, bg=cols, col=ifelse(pBGC_hits>0,1,"grey50") )
)
with(all_log[pBGC_indx,],
     points( jitter((genome_anti_hits),amount = 1) ~ jitter((prophage_hits), amount=1 ) , cex=ifelse(pBGC_hits==1,1,1.5),
            pch=ifelse(pBGC_hits==1,24,23), bg=cols, col=1 )
)

legend("topright" , y.intersp = .8, legend = levels(all_log$mainNames), pt.bg=cols, pch=21, title="Genus")
legend("bottomright", y.intersp = 1,   legend = c("0","1",">1"), pt.bg="grey90", pch=c(21,24,23), title = "pBGC")
dev.off()
########
#pBGCs as a function of phages
########
png("pBGCvsPHAGE.png", width=7.4, height = 4, units = "in", res=600)
par(mai=c(.5,0.5,0.1,0.1), mgp=c(1.4,0.4,0), font.lab=2, family="serif")
plot( jitter((all_log$pBGC_hits),amount = .1) ~ jitter(((all_log$prophage_hits)),amount = .2), cex=ifelse(all_log$pBGC_hits>0,1,0.5),
      xlab="Prophages", ylab="#pBGCs" , pch=ifelse(all_log$pBGC_hits>0,24,21), bg=all_log$cols)
legend("topright" , legend = levels(all_log$mainNames), pt.bg=cols, pch=21)
dev.off()
########
#Genomic BGCs as a function of phages
########
png("pBGCvsgBGCs.png", width=7.4, height = 4, units = "in", res=600)
par(mai=c(.5,0.5,0.1,0.1), mgp=c(1.4,0.4,0), font.lab=2, family="serif")
plot( jitter((all_log$pBGC_hits),amount = .1) ~ jitter(((all_log$genome_anti_hits)),amount = .1), cex=ifelse(all_log$pBGC_hits>0,1,0.5),
      xlab="#gBGCs", ylab="#pBGCs" , pch=ifelse(all_log$pBGC_hits>0,24,21), bg=all_log$cols)
legend("topright" , legend = levels(all_log$mainNames), pt.bg=cols, pch=21)
dev.off()

# source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
# s3d=scatterplot3d(x = log1p(jitter(all_log$genome_anti_hits)), y=log1p(jitter(all_log$prophage_hits)),
#                   z=all_log$pBGC_hits,bg=all_log$cols, pch=" ",grid=TRUE, box=FALSE, xlab = "gBGCs", ylab = "Phages",zlab = "pBGCs")
#
# addgrids3d(x = log1p(jitter(all_log$genome_anti_hits)), y=log1p(jitter(all_log$prophage_hits)),
#            z=all_log$pBGC_hits, grid = c( "xz", "yz"))
# s3d$points3d(x = log1p(jitter(all_log$genome_anti_hits)), y=log1p(jitter(all_log$prophage_hits)),
#              z=all_log$pBGC_hits,bg=all_log$cols, pch=ifelse(all_log$pBGC_hits>0,24,21))
if(writePDF) dev.off()

#Distribution of pBGC types

only_pBGC %>% separate_rows(pBGC_types, sep=",\\s*") %>% group_by(pBGC_types) %>% summarise(count = n()) %>%
  ggplot(.,aes(x = reorder(pBGC_types,(count)), y = count))+
  geom_bar(stat="identity") +  geom_text(aes(label = count), vjust = 0, hjust= -0.5)+
  theme_bw()+ coord_flip()+
  xlab("Cluster Types")+ ylab('Count') +
  ggtitle("Which cluster type are overrepresented in pBGCs carriers?")
ggsave("distribution_bgc_types_bacteria_w_pBGC.png", width =11, height = 8, units = "in")

#Display percents for each cluster type
only_pBGC %>% separate_rows(pBGC_types, sep=",\\s*") %>%group_by(pBGC_types) %>% tally() %>%
  mutate(percent = (n / 702) *100) %>% arrange(desc(percent))
