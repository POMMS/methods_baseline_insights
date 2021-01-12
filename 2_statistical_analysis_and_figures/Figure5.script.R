#===================================================#
#========Figure 5 ==================================#
#========April 2020=================================#
#========McCann et al===============================#

#Library loading

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(ggthemes)
library(metagenomeSeq)
library(cowplot)
library(readr)

#===================================================#
#=======load data into phyloseq object==============#
#===================================================#

asv.taxon<- readr::read_tsv("Rawls-Rawls-20190118-R24_05_filterSample-table-taxon.tsv") %>% as.data.frame()
colnames(asv.taxon)[1]<- "ASV_ID" 
colnames(asv.taxon)[62:63]<- c("Taxonomy", "Confidence")
asv<- asv.taxon[2:61]
asv.matrix<- as.matrix(asv)
rownames(asv.matrix)<- asv.taxon$ASV_ID
otab<- otu_table(asv.matrix, taxa_are_rows = TRUE)

taxtab<- asv.taxon[62:63]
rownames(taxtab)<- asv.taxon$ASV_ID
taxtab<- taxtab %>% separate(Taxonomy, c("Kingdom", "Phylum", "Class", 
                                   "Order", "Family", "Genus", "Species"),
                       ";")
taxtab<-apply(taxtab[, 1:7], 1:2, function(x) strsplit(x, "__")[[1]][2]) %>% data.frame(stringsAsFactors = FALSE)

w<- which(!is.na(taxtab$Species))
if(length(w)>0){
  taxtab$Species[w]<- paste0(taxtab$Genus[w], " ", taxtab$Species[w], ";NA")
}

cols<- colnames(taxtab)[1:7]
taxtab$Combined <- do.call(paste, c(taxtab[cols], sep=";"))
taxtab$LastTax<- sapply(taxtab$Combined, function(x) strsplit(x, ";NA")[[1]][1])

find_last<- function(string, split_char){
  temp_str<- unlist(strsplit(string, split_char))
  l<- length(temp_str)
  return(temp_str[l])
}
taxtab$LastTax<-sapply(taxtab$LastTax, function(x) find_last(x, ";"))
taxtab<-as.matrix(taxtab)
rownames(taxtab)<-rownames(otab)
taxtab<- tax_table(taxtab)

tree<- read_tree("Rawls-Rawls-20190118-R24_09A_seq-aligned-mask-rooted-raxml.nwk")

library(Biostrings)
repSeq<- readDNAStringSet("sequences.fasta")
names(repSeq)<- rownames(otab)
m<- which(names(repSeq)%in%rownames(taxtab))
repSeq<- repSeq[m]

map<- read_tsv("Rawls-Rawls-20190118-R24_map.txt") 
map<- map[-1,]
colnames(map)[1]<- "SampleID"
map<- map %>% select(-BarcodeSequence, -LinkerPrimerSequence, 
                     -PrimerID, -Project, -Submission, 
                     -Description) %>% as.data.frame()
rownames(map)<- map$SampleID
map$Health<- dplyr::recode(map$Health, "healthy"="HWC",
                           "overweight"="OB")
map <- sample_data(map)

phy<- phyloseq(otab,map,taxtab, tree, repSeq)

##Remove samples that do not fit
phy<- subset_samples(phy, !SampleID %in% c("71T1", "78T1", "87T1", "105T1", "35T1"))

save(phy, file="phyloseq.RData")

#===================================================#
#=======Normalization of data=======================#
#===================================================#

load("phyloseq.RData")

library(Wrench)
met<- phyloseq_to_metagenomeSeq(phy)
met$Health<- factor(met$Health)
met<- wrenchNorm(met, condition=met$Health)

phy.norm<- phy
otu_table(phy.norm)<- otu_table(MRcounts(met, log=TRUE, norm=TRUE),
                                taxa_are_rows = TRUE)

phy.norm.genus<- tax_glom(phy.norm, taxrank = "Genus")

phy.genus<- phy
otu_table(phy.genus)<- otu_table(as(MRcounts(met, norm=TRUE, log=FALSE), "matrix"), taxa_are_rows = TRUE)
phy.genus<- tax_glom(phy.genus, taxrank = "Genus")

phy.genus.norm<- tax_glom(phy.norm, taxrank = "Genus")
samdat<- as(sample_data(phy.genus), "data.frame")

#===================================================#
#=========Clustering and Plotting===================#
#===================================================#

d<- phyloseq::distance(phy.norm.genus, method="wunifrac")
#d<- phyloseq::distance(phy.genus, method="unifrac")
h<- hclust(d)
plot(h)
k<- c(3:12)
ct.list<- sapply(k, function(x){
  ct<- ct<-cutree(h, k=x)
  cs<- chisq.test(samdat$Health, ct)
  list(k=x, ct=ct, chisq=cs, p=cs$p.value)
}, simplify=FALSE, USE.NAMES = TRUE)
ct.summary<-data.frame(do.call(rbind, ct.list))[c("k", "p")] %>%
  mutate(k=as.numeric(k),
         p=as.numeric(p))
ggplot(ct.summary, aes(k, p))+
  geom_line()+
  geom_point()
ggsave("chi sq different genus agglom cutoffs.png", dpi=300)

d<- phyloseq::distance(phy.norm.genus, method="wunifrac") #here we're not agglomerating at genus
h<- hclust(d)

k=5
ct<- cutree(h, k=k)

library(ggdendro)
df<- data.frame(SampleID=h$labels[h$order], Clusters=cutree(h, k=k))
df2<-as(sample_data(phy.norm.genus), "data.frame") %>% select(SampleID, Health)
df3<- merge(df, df2, by.x="SampleID", sort=FALSE) %>%
  mutate(SampleID=factor(.$SampleID, levels = .$SampleID))

p1<-ggdendrogram(h, rotate=FALSE, labels=FALSE, leaf_labels = FALSE, theme_dendro=FALSE) +
  theme_nothing()+
  theme(axis.text.x = element_blank())

p2<-ggplot(df3,aes(SampleID,y=1, fill=factor(Clusters)))+
  geom_tile()+
  scale_y_continuous(expand=c(0,0))+
  #geom_text(aes(label=SampleID))+
  scale_fill_tableau()+
  theme_nothing()+
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position="none")

p3<-ggplot(df3,aes(SampleID,y=1, fill=factor(Health)))+
  geom_tile()+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_tableau()+
  theme_nothing()+
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position="none")

library(grid)

p4<- plot_grid(p2, p3, ncol=1)

setEPS()
  postscript("Figure5A.eps")
  grid.newpage()
  print(p1, vp = viewport(x = 0.4, y = 0.53, width = 0.8, height = 0.8))
  print(p4, vp = viewport(x = 0.4, y = .1, width = 0.74, height = 0.1))
dev.off()

sample_data(phy.norm.genus)$ct<- ct[h$labels%in%sample_data(phy.norm.genus)$SampleID]
ps.total.prop<- psmelt(phy.norm.genus) %>%
  mutate(Patient=factor(Patient))%>%
  group_by(Patient) %>%
  mutate(Prop=100*Abundance/sum(Abundance))
ps.total.prop<- ps.total.prop[order(-ps.total.prop$Abundance),]
topg<-lapply(levels(ps.total.prop$Patient), function(x){
  head(unique(ps.total.prop$Genus[order(-ps.total.prop$Prop) & ps.total.prop$Patient==x]), 10)})
topg<-unique(unlist(topg))
topg.keep<- topg
ps.toplot<-ps.total.prop[ps.total.prop$Genus%in%topg,] %>%
  mutate(SampleID=factor(SampleID, levels=labels(dendro)),
         ct=factor(ct))
palette1<-c("#5d8aa8", "#e32636", "#ffbf00", "#9966cc", "#a4c639", 
            "#915c83", "#008000", "#9c3838", "#8db600", "#00ffff",
            "#7fffd4", "#4b5320",  "#87a96b", "#a52a2a", "#ff9966", 
            "#fdee00", "#007fff", "#a1caf1", "#ffe135", 
            "#21abcd", "#ff55a3", "#66ff00","#cc5500", "#006a4e",
            "#702963", "#536872", "#91a3b0", "#bf94e4", "#78866b", 
            "#b87333", "#ccff00","#e9d66b","#6e7f80","orchid","limegreen","royalblue4","blue4",
            "slategrey","forestgreen",
            "gray1", "blueviolet","grey89","orangered3","mediumturquoise")

new.palette<-c("blue","steelblue3","darkolivegreen4","red1",
               "springgreen4","olivedrab2","thistle","goldenrod2","midnightblue",
               "peru","yellow1","skyblue4","orangered1","chartreuse3","turquoise4",
               "seagreen3","wheat4","tomato2","cadetblue3","darkseagreen4",
               "orchid","limegreen","royalblue4","blue4","slategrey","forestgreen",
               "gray1", "blueviolet","grey89","orangered3","mediumturquoise")
ctp<-ggplot(ps.toplot, aes(x=SampleID, y=1, fill=factor(ct)))+
  geom_tile()+
  #geom_text(aes(label=ct))+
  scale_fill_tableau()+
  theme_nothing()+
  facet_wrap(~Health, scale="free_x")

taxap<-ggplot(ps.toplot, aes(x=factor(SampleID), y=Prop))+
  geom_bar(stat="identity", aes(fill=Genus))+
  scale_fill_manual(values=c(palette1, new.palette),
                    guide=guide_legend(ncol=3))+
  scale_y_continuous(limits=c(0,100), expand=c(0,0),)+
  xlab("Samples")+ylab("Proportion of Total")+
  facet_wrap(~Health, scale="free_x")+
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x.bottom = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "bottom")
tleg<- get_legend(taxap)
taxap<- taxap + theme(legend.position = "none")

