#######################################################################################
####Bacterial community Analysis for Drought Experiment using Bean and Switchgrass#####
#######################################################################################

#######################################################################################
#####Updated February 22, 2023#########################################################
#######################################################################################


#######################################################################################
####Navigate to HPCC remote address "/mnt/research/ShadeLab/WorkingSpace/Bandopadhyay_WorkingSpace/Bean_SW_drought_rhizobiome/FinalCode_Clean_ForGit"
#######################################################################################

#######################################################################################
#####Contaminant, mitochondria, chloroplast removal####################################
#######################################################################################


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(grid)
library(reshape2)
library(grid)
#install.packages("vctrs")
#update.packages("tidyverse")
#update.packages("scales")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library(decontam)

otu = read.csv("otu_table_dn_99.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("sample-metadata-decontam.csv", sep=",", row.names=1)
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged
sample_names(phyloseq_merged)

##check column names of taxonomy table
colnames(tax_table(phyloseq_merged))

##remove mito, chloroplast and archaea, eukarya
phyloseq_merged_clean <- phyloseq_merged %>%
  subset_taxa(
    Domain == "d__Bacteria" &
      Family   != " f__Chloroplast" &
      Family  != " f__Mitochondria"
  )
phyloseq_merged_clean

head(sample_data(phyloseq_merged_clean))

##inspect library sizes
df <- as.data.frame(sample_data(phyloseq_merged_clean)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phyloseq_merged_clean)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot<- ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point() ##Supplemental Figure S3

ggsave(filename = "librarySize.pdf", plot = plot,
       width = 16,
       height = 13, units = c("cm"),
       dpi = 300)


##identify contaminants by prevalence. in this method, the distribution 
#of the frequency of each sequence feature as a function of the 
#prevalence is used to identify contaminants. In this method, 
#the prevalence (presence/absence across samples) of each sequence feature 
#in true positive samples is compared to the prevalence in negative controls 
#to identify contaminants.

sample_data(phyloseq_merged_clean)$is.neg <- sample_data(phyloseq_merged_clean)$Sample_or_Control == "Control Sample"
contamdf.prev0.5 <- isContaminant(phyloseq_merged_clean, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev0.5$contaminant)
##FALSE  TRUE 
##83100    47 

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phyloseq_merged_clean, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev0.5$contaminant)
plot<- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") ##Supplemental Figure S3

ggsave(filename = "prevalence.pdf", plot = plot,
       width = 16,
       height = 13, units = c("cm"),
       dpi = 300)

write.csv(df.pa, "contaminant-table-dna-cdna.csv")

write.csv(contamdf.prev0.5, "contaminant-prev-0.5-dna-cdna.csv")

##removing contaminants from phyloseq object
df.pa <- read.csv("contaminant-prev-0.5-dna-cdna.csv")
View(df.pa)
subset.df <- subset(df.pa, contaminant== "FALSE")
View(subset.df)
keep.taxa <- as.vector(subset.df$X)
keep.taxa

phyloseq_merged_clean_decontam <- prune_taxa(keep.taxa, phyloseq_merged_clean)
phyloseq_merged_clean_decontam


##export otu table out of phyloseq object 

OTU1 = as(otu_table(phyloseq_merged_clean_decontam), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

write.csv(OTUdf, "OTU_clean.csv")

#remove neg controls from otu table-- might add more rows which are rowsums=0, read back in to rarefy reads to 15K

otu = read.csv("OTU_clean.csv", sep=",", row.names=1)#with no contaminants, mito, chloroplast, neg control 
tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("sample-metadata-decontam.csv", sep=",", row.names=1)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

#############################################
#####rarefaction curves-- all samples#######

otu1 = read.csv("OTU_clean.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
tax1 = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax1 = as.matrix(tax1)
metadata1 = read.csv("sample_metadata_dna_cdna.csv", sep=",", row.names=1)
OTU1 = otu_table(otu1, taxa_are_rows = TRUE)
TAX1 = tax_table(tax1)
meta1 = sample_data(metadata1)

phyloseq_merged1 = phyloseq(OTU1, TAX1, meta1)
phyloseq_merged1

col<- cyan4
rarecurve<-rarecurve(t(otu_table(phyloseq_merged1)), step=20, label = FALSE, col="cyan4") ## Supplemental Figure S3

###############################
sample_df <- data.frame(sample_data(phyloseq_merged))
View(sample_df)
write.csv(sample_df, "samples_pre_rarefied.csv")

library(vegan)
rarecurve(t(otu_table(phyloseq_merged)), step=50)
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged, sample.size = 15000, rngseed = TRUE, trimOTUs=FALSE)

##7 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: DNA_129DNA_146DNA_163DNA_164DNA_245

sample_df <- data.frame(sample_data(phyloseq_rarefy))

View(sample_df)
write.csv(sample_df, "samples_rarefied.csv")

##use anti-join to select non-present samples in rarefied object
#This join is like eg: for df1-df2, it selects all rows from df1 that are not present in df2.
df1<- read.csv("samples_pre_rarefied.csv")
df2 <- read.csv("samples_rarefied.csv")
df= df1 %>% anti_join(df2,by="sampleid")

View(df)

write.csv(df, "samples-lost-rarefaction.csv")

OTU2 = as(otu_table(phyloseq_rarefy), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_clean_noneg_rarefied15k.csv")

#######################################################################################
###### Divide into DNA and cDNA dataset for rna/dna ratio calculations ################
#######################################################################################
library(tidyverse)
require(dplyr) 
require(tibble)
dna_cdna<- read.csv("OTU_clean_noneg_rarefied15k_removedsamples.csv", row.names=1) ##new rarefied and cleaned OTU table, with 15k reads per sample, and exact match of DNA RNA samples
dna_cdna = dna_cdna[rowSums(dna_cdna[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(dna_cdna)
dna_cdna ##otus reduced to 48647
dna_cdna1 <-as_tibble(dna_cdna, rownames="OTUID") 
dna_cdna1
has_rownames(dna_cdna1)
rownames(dna_cdna1)
dna<- dna_cdna1 %>% select(OTUID, starts_with("DNA"))
View(dna)
cdna <- dna_cdna1 %>% select(OTUID, starts_with("cDNA"))
View(cdna)
write.csv(dna, "DNA.csv")
write.csv(cdna, "cDNA.csv")

##already removed all negative samples from DNA and CDNA dataset
##check that all rowsums are >1, remove any rowsums=0 

#subset DNA dataframe to check only those OTUs that appear in at least one sample
DNA<- read.csv("DNA.csv", row.names=1)
View(DNA)##48647 otus
dnanozero = DNA[rowSums(DNA[])>0, ,drop=FALSE]
dnanozero
View(dnanozero)##34168 otus
write.csv(dnanozero, "dna_noneg_nozero.csv")

dnanozero["Total"] <- rowSums(dnanozero)
dnanozero
dnanozero1 = dnanozero %>% filter(dnanozero$Total== 0)

View(dnanozero1)


##subset cDNA dataframe to check only those OTUs that appear in at least one sample
CDNA<- read.csv("cDNA.csv", row.names=1)
View(CDNA)##48647 otus
cdnanozero = CDNA[rowSums(CDNA[])>0, ,drop=FALSE]
cdnanozero
View(cdnanozero)##27012 otus
write.csv(cdnanozero, "cdna_noneg_nozero.csv")


##merge to see overlap in OTUs between dna and cdna datasets

A = read.csv("dna_noneg_nozero.csv")

B = read.csv("cdna_noneg_nozero.csv")

C<- inner_join(A, B) ##12533 taxa across 425 columns, join keeps both sets of data from A and B
y <- match_df(A, B) ##12533 taxa across 213 columns, match_df only keeps DNA columns (values from left table A) 
View(y)
View(C)

##if negative samples removed, then 12,533 taxa match between dna and cdna datasets.



#######################################################################################
###### Use DNA and cDNA dataset for rna/dna ratio calculations ########################
#######################################################################################

library(tidyverse)

dna<- read.csv("DNA_mat_noheader.csv", header=FALSE) ##remove OTU names (rownames) and header row with sample names to create matrix for calculation of ratios
cdna <- read.csv("cDNA_mat_noheader.csv", header=FALSE) ##remove OTU names (rownames) and header row with sample names to create matrix for calculation of ratios
View(dna)
View(cdna)



###Method 1 for ratio calculation
cdna.mat <- as.matrix(cdna)

View(cdna.mat)
dna.mat<- as.matrix(dna)

new.mat= ifelse (dna.mat == 0 & cdna.mat>0, 100, cdna.mat/dna.mat) ##this works
View(new.mat)

write.csv(new.mat, "cdna-dna-ratio-method1.csv")


##selecting ratios greater than equal to 1

new.mat[new.mat < 1] <- 0 
View(new.mat)

new.mat[!is.finite(new.mat)] <- 0
View(new.mat)

write.csv(new.mat, "cdna-dna-allfinite-greaterthanEqone_method1.csv")

#####Method 2 for calculating ratios

dna.mat<- as.matrix(dna)
dna.mat[dna.mat == 0] <- 1

View(dna.mat)

cdna.mat <- as.matrix(cdna)

View(cdna.mat)
merge.mat<- cdna.mat/dna.mat


View(merge.mat)
write.csv(dna.mat, "dna_mat_method2.csv")##need this for replacing with dna abun
write.csv(merge.mat, "cdna-dna-ratio-method2.csv")

##selecting ratios greater than equal to 1

merge.mat[merge.mat < 1] <- 0 
View(merge.mat)


merge.mat[!is.finite(merge.mat)] <- 0
View(merge.mat)

write.csv(merge.mat, "cdna-dna-allfinite-greaterthanEqone_method2.csv") 


##replacing ratios greater than equal to 1 with abund values in **dna table**, using method 2

cdna_ratio<- read.csv("cdna-dna-allfinite-greaterthanEqone_method2.csv")
View(cdna_ratio)
dna_abun <- read.csv("dna_mat_method2.csv") ##use this to filter to dna abun 
View(dna_abun)
cdna_ratio_mat<- as.matrix(cdna_ratio)
View(cdna_ratio_mat)
dna_abun_mat <- as.matrix(dna_abun)
View(dna_abun_mat)

mask <- cdna_ratio_mat > 0
mask

cdna_ratio_mat[mask] <- dna_abun_mat[mask]
cdna_ratio_mat
View(cdna_ratio_mat)


write.csv(cdna_ratio_mat, "finalactivecdna_method2_changing_to_dnaabun.csv")


##remove rowsums equals zero
cdna_active<- read.csv("finalactivecdna_method2_changing_to_dnaabun.csv", row.names=1)

activefinalnozero = cdna_active[rowSums(cdna_active[])>0, ,drop=FALSE]

View(activefinalnozero) 

## for ratio >= 1 26653 taxa across 212 samples method 1, 26653 taxa across 212 samples method 2 
#for ratio >1 17096 taxa across 212 samples for method 2, 26546 taxa across 212 samples for method 1  
##for ratio>1 this makes sense as method 1 will overestimate taxa due to ratio set as 100, which for a similar taxa in method 2 might be ratio 1 and hence not accounted for


write.csv(activefinalnozero, "activefinalnozero_method2_changing_to_dnaabun.csv")


########################################################################################################
##############Setting phantom taxa detection threshold for final active community count table ##########
########################################################################################################


#### Report the number and proportion of phantom taxa (DNA=0) that are detected in the unrarefied DNA dataset: detected as in has rowsums>0
DNA<- read.csv("DNA.csv", row.names=1)
View(DNA)##48647 otus


CDNA<- read.csv("cDNA.csv", row.names=1)
View(CDNA)


DNA["Total_DNA"] <- rowSums(DNA)
DNA

DNA1 = DNA %>% filter(DNA$Total_DNA== 0)

View(DNA1) ##14479 otus with rowsums == 0


DNA2<- DNA %>% select("Total_DNA")
View(DNA2)

CDNA["Total_CDNA"] <- rowSums(CDNA)
CDNA

CDNA1 = CDNA %>% filter(CDNA$Total_CDNA== 0)

View(CDNA1) ##21635 otus with rowsums == 0

CDNA2<- CDNA %>% select("Total_CDNA")
View(CDNA2)

dna_cdna <- merge(DNA2, CDNA2, by="row.names")

View(dna_cdna)

dna_cdna1<- dna_cdna %>% mutate(phantom = ifelse(Total_DNA == 0 & Total_CDNA > 0, "True", "False"))
dna_cdna1

dna_cdna2<- dna_cdna1  %>% filter(phantom == "True")

View(dna_cdna2) ##14479 phantoms
dna_cdna3<- dna_cdna2 %>% remove_rownames %>% column_to_rownames(var="Row.names")
View(dna_cdna3)


##plot distribution of ratios that equals to 1 across samples for the phantom taxa (14479 otus) in the rarefied dataset
library(tidyverse)

##using rarefied DNA and cDNA copies with 15k reads.
dna<- read.csv("DNAcopy.csv", header=TRUE, row.names = 1)
cdna <- read.csv("cDNAcopy.csv", header=TRUE, row.names = 1)

Active_DNA<- read.csv("activefinalnozero_method2_changing_to_dnaabun.csv", row.names=1)

View(Active_DNA)

View(dna)
View(cdna)


activeDNA_checkphantom <- Active_DNA[rownames(Active_DNA)%in%rownames(dna_cdna3),]
View(activeDNA_checkphantom) ##all 14479 phantoms are present in this dataset

dna_phantomtaxa <- dna[rownames(dna)%in%rownames(dna_cdna3),]
View(dna_phantomtaxa)

cdna_phantomtaxa <- cdna[rownames(cdna)%in%rownames(dna_cdna3),]
View(cdna_phantomtaxa)

dna1<- dna_phantomtaxa %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="dna", -ASV)#gather columns into key value pairs
View(dna1)
cdna1<- cdna_phantomtaxa %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="cdna", -ASV)
View(cdna1)
merge.df = full_join(dna1, cdna1, by=c("ASV", "SampleID"))
View(merge.df)##3069548 taxa which is 212 samples * 14479 taxa

merge.df[is.na(merge.df)] <-0
View(merge.df)
metadata<-read.csv("metadata_dna_cdna_new.csv")

newdf = full_join(merge.df, metadata, by=c("SampleID"))
View(newdf)

finaldf = newdf %>%  mutate(ratioMethod1= ifelse(dna == 0 & cdna > 0, 100, cdna/dna)) %>%  
  mutate(dna2= ifelse(dna == 0, 1,dna))  %>%  mutate(ratioMethod2= cdna/dna2) %>% mutate(ratio_nophantoms = cdna/dna)
View(finaldf)
finaldf[is.na(finaldf)] <-0
finaldf$ratioMethod1[!is.finite(finaldf$ratioMethod1)] <- 0 
finaldf$ratioMethod2[!is.finite(finaldf$ratioMethod2)] <- 0
finaldf$ratio_nophantoms[!is.finite(finaldf$ratio_nophantoms)] <- 0
View(finaldf) 

library(dplyr) ##use package version 1.0.7 , this is important to have the code below work


finaldf2 <- finaldf %>% group_by(ASV) %>% 
  mutate(ratio_eq_1_method2_abs = sum(ratioMethod2 == 1))  %>% 
  group_by(ASV) %>% 
  mutate(ratio_eq_1_method2_absoutoftotal_percent = (sum(ratioMethod2 == 1)/212)*100)  %>% 
  group_by(ASV) %>% 
  mutate(ratio_eq_1_method2_absoutofcdnapresent_percent = (sum(ratioMethod2 == 1)/sum(dna+cdna!=0)) *100) %>% 
  group_by(ASV) %>% 
  mutate(ratio_gr_eq_1_method2_absoutofcdnapresent_percent = (sum(ratioMethod2 >= 1)/sum(dna+cdna!=0)) *100)  %>% 
  group_by(ASV) %>% 
  mutate(ratio_gr_eq_1_method2_absoutoftotal_percent  = (sum(ratioMethod2 >= 1)/212) *100)  %>% 
  group_by(ASV) %>%
  mutate(ratio_gr_eq_1_method2_abs = sum(ratioMethod2 >= 1))  

#used this one instead of finaldf2 with mutate function
finaldf3 <- finaldf %>% group_by(ASV) %>% 
  summarise(ratio_eq_1_method2_abs = sum(ratioMethod2 == 1),
            ratio_eq_1_method2_absoutoftotal_percent = (sum(ratioMethod2 == 1)/212)*100,
            ratio_eq_1_method2_absoutofcdnapresent_percent = (sum(ratioMethod2 == 1)/sum(dna+cdna!=0)) *100,
            ratio_gr_eq_1_method2_absoutofcdnapresent_percent = (sum(ratioMethod2 >= 1)/sum(dna+cdna!=0)) *100,
            ratio_gr_eq_1_method2_absoutoftotal_percent  = (sum(ratioMethod2 >= 1)/212) *100,
            ratio_gr_eq_1_method2_abs = sum(ratioMethod2 >= 1))

View(finaldf2)
View(finaldf3)
write.csv(finaldf3, "phantoms_ratio_summary.csv")


##subset finaldf3 to include only those phantom otus which are represented in at least 5 or 10 % of samples
##filter based on  column "ratio_gr_eq_1_method2_absoutoftotal_percent"

finaldf4<- finaldf3 %>% filter(ratio_gr_eq_1_method2_absoutoftotal_percent>=10)
View(finaldf4) ##only 79 phantom otus remain 

finaldf4<- finaldf3 %>% filter(ratio_gr_eq_1_method2_absoutoftotal_percent>=5)
View(finaldf4) ##only 214 phantom otus remain 

##remove phantom taxa which do not belong to the 214 as selected above. This means deleting 14479-214 taxa == 14265
##use anti join : returns all rows from x without a match in y


dna_cdna_3_1 <- cbind(ASV = rownames(dna_cdna3), dna_cdna3)
View(dna_cdna_3_1)
phantoms_delete<- anti_join(dna_cdna_3_1, finaldf4, by='ASV')
View(phantoms_delete)##14265 remains, now antijoin this with the active dna dataset to get final active dna dataset with the threshold detection for phantoms set at 5%
Active_DNA_1 <- cbind(ASV = rownames(Active_DNA), Active_DNA)
View(Active_DNA_1)
active_final<- anti_join(Active_DNA_1, phantoms_delete, by='ASV')
View(active_final)
active_final<- active_final %>% select(!ASV)
View(active_final)
write.csv(active_final, "active_final_DNAabund_phantoms_threshold_detection5%.csv")



##plot distribution of singleton phantom otus
#finaldf6<- dna_cdna2 %>% filter(Total_DNA==0 & Total_CDNA==1)
#View(finaldf6)

##plot the results as distributions, x axis should be the percent ratios which meet a criteria(percent greater than 1, equals 1 etc)
#y axis should be the count of OTUs that meet that criteria
##black outline, black fill
plot<-ggplot(finaldf3, aes(x=ratio_eq_1_method2_abs)) + geom_histogram(binwidth=0.5)
plot

ggsave(filename = "ratio_eq1_method2_abs_countdata.tiff", plot = plot,
       width = 17,
       height = 15, units = c("cm"),
       dpi = 300)


# Draw with black outline, white fill
plot2<- ggplot(finaldf3, aes(x=ratio_eq_1_method2_absoutoftotal_percent)) + 
  geom_histogram(binwidth=0.1, colour="black", fill="white") + ylab("number of phantom OTUs") + xlab("percent detected in  total sample size(n=212)")
plot2
ggsave(filename = "ratio_eq_1_method2_absoutoftotal_percent_countdata_newaxes.tiff", plot = plot2,
       width = 17,
       height = 20, units = c("cm"),
       dpi = 300)



# Density curve
ggplot(finaldf3, aes(x=ratio_eq_1_method2_abs)) + geom_density()

# Histogram overlaid with kernel density curve
plot1<- ggplot(finaldf3, aes(x=ratio_eq_1_method2_abs)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=1,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
ggsave(filename = "ratio_eq1_method2_abs.tiff", plot = plot1,
       width = 17,
       height = 15, units = c("cm"),
       dpi = 300)





####analysis with non-rarefied dataset to see phantom taxa distrbution before rarefaction
otu_norarefied<- read.csv("OTU_clean.csv", row.names=1)
otu_norarefied_nozero = otu_norarefied[rowSums(otu_norarefied[])>0, ,drop=FALSE]
View(otu_norarefied) #83100
View(otu_norarefied_nozero) ##82911

#test<- read.csv("OTU_clean_noneg_rarefied15k_removedsamples.csv") ##checking for formatting
#View(test)

otu_norarefied_final <-as_tibble(otu_norarefied_nozero, rownames="OTUID") 
otu_norarefied_final
has_rownames(otu_norarefied_final)
rownames(otu_norarefied_final)
dna_otu_norarefied_final <- otu_norarefied_final %>% select(OTUID, starts_with("DNA"))
View(dna_otu_norarefied_final)
cdna_otu_norarefied_final <- otu_norarefied_final %>% select(OTUID, starts_with("cDNA"))
View(cdna_otu_norarefied_final)
write.csv(dna_otu_norarefied_final, "otu_norarefied_DNA.csv")
write.csv(cdna_otu_norarefied_final, "otu_norarefied_cDNA.csv")

##for non-rarefied dataset calculate phantom taxa 
dna_otu_norarefied_final <-as_tibble(dna_otu_norarefied_final, row.names=2) ##doesnt work
View(dna_otu_norarefied_final)

dna_otu_norarefied_final <- read.csv("otu_norarefied_DNA.csv", row.names=1)

View(dna_otu_norarefied_final)

cdna_otu_norarefied_final <- read.csv("otu_norarefied_cDNA.csv", row.names=1)
View(cdna_otu_norarefied_final)

dna_otu_norarefied_final["Total_DNA"] <- rowSums(dna_otu_norarefied_final)
dna_otu_norarefied_final

dna_otu_norarefied_final1 = dna_otu_norarefied_final %>% filter(dna_otu_norarefied_final$Total_DNA== 0)

View(dna_otu_norarefied_final1) ## 28459 otus with rowsums == 0


dna_otu_norarefied_final2<- dna_otu_norarefied_final %>% select("Total_DNA")
View(dna_otu_norarefied_final2)

cdna_otu_norarefied_final["Total_CDNA"] <- rowSums(cdna_otu_norarefied_final)
cdna_otu_norarefied_final

cdna_otu_norarefied_final1 = cdna_otu_norarefied_final %>% filter(cdna_otu_norarefied_final$Total_CDNA== 0)

View(cdna_otu_norarefied_final1) ## 39405 otus with rowsums == 0

cdna_otu_norarefied_final2<- cdna_otu_norarefied_final %>% select("Total_CDNA")
View(cdna_otu_norarefied_final2)

dna_cdna_norarefied <- merge(dna_otu_norarefied_final2, cdna_otu_norarefied_final2, by="row.names")

View(dna_cdna_norarefied)

dna_cdna_norarefied1<- dna_cdna_norarefied %>% mutate(phantom = ifelse(Total_DNA == 0 & Total_CDNA > 0, "True", "False"))
dna_cdna_norarefied1

dna_cdna_norarefied2<- dna_cdna_norarefied1  %>% filter(phantom == "True")
View(dna_cdna_norarefied2) ##28459 otus phantoms
dna_cdna_norarefied3<- dna_cdna_norarefied1  %>% filter(phantom == "False")
View(dna_cdna_norarefied3) ##54452 otus not phantoms

##we need to see if the phantoms detected in the rarefied data show up in the non rarefied dataset
##so lets merge the Total DNA dataset from unrarefied data with the phantom == true dataset from the rarefied data



write.csv(dna_cdna2, "dna_cdna2.csv")
View(dna_cdna2)
View(dna_otu_norarefied_final2)
write.csv(dna_otu_norarefied_final2, "dna_otu_norarefied_final2.csv")

dna_cdna2<- read.csv("dna_cdna2.csv", row.names=1)
dna_otu_norarefied_final2 <- read.csv("dna_otu_norarefied_final2.csv", row.names=1)
phantoms <- merge(dna_cdna2, dna_otu_norarefied_final2, by="row.names")
View(phantoms)

phantoms_select <- phantoms %>% filter (Total_DNA.x==0 & Total_DNA.y>0)
View(phantoms_select) ##922 taxa out of 14479 phantom taxa are detected in the unrarefied dataset. calculate proportions.

## 922/14479 *100= 6.4 %

View(phantoms_select)

write.csv(phantoms_select, "phantoms_rarefiedzero_unrarefiedgreaterthanzero.csv")


#### For those phantom taxa detected in unrarefied dataset, Plot distribution of DNA counts of no. sequence (meaning rowsums, should mostly be rare) from unrarefied dataset
##plot distribution of DNA counts (rowsums) of phantom taxa detected in unrarefied dataset

## These both result in the same output:

##black outline, black fill
plot<-ggplot(phantoms_select, aes(x=Total_DNA.y)) + geom_histogram(binwidth=0.5)
plot


ggsave(filename = "phantomDNAcount_in_unrarefied_countdata.tiff", plot = plot,
       width = 17,
       height = 15, units = c("cm"),
       dpi = 300)


# Draw with black outline, white fill
plot2<- ggplot(phantoms_select, aes(x=Total_DNA.y)) +
  geom_histogram(binwidth=1, colour="black", fill="white")
plot2
ggsave(filename = "phantomDNAcount_in_unrarefied_countdata.tiff", plot = plot2,
       width = 17,
       height = 15, units = c("cm"),
       dpi = 300)



# Density curve
ggplot(phatoms_select, aes(x=Total_DNA.y)) + geom_density()

# Histogram overlaid with kernel density curve
plot1<- ggplot(phantoms_select, aes(x=Total_DNA.y)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=1,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
ggsave(filename = "phantomDNAcount_in_unrarefied.tiff", plot = plot1,
       width = 17,
       height = 15, units = c("cm"),
       dpi = 300)


#### Report proportion of sequences attributed to phantoms in the unrarefied DNA dataset (out of total no. sequences): these are the phantoms from the rarefied data which may have sequences in unrarefied DNA dataset
View(dna_otu_norarefied_final2)
sum(dna_otu_norarefied_final2$Total_DNA) #19939133
View(phantoms_select)
sum(phantoms_select$Total_DNA.y) # 7576

## 0.04 % of reads attribute to phantom taxa out of the total DNA reads in unrarefied dataset 


#### Plot distribution of occupancy from unrarefied dataset: do the same as number 2 above but instead of rowsums, use occupancy from presence absence data
dna_otu_norarefied_final <- read.csv("otu_norarefied_DNA.csv")
View(dna_otu_norarefied_final)
dna.norarefied.long<-pivot_longer(dna_otu_norarefied_final, !OTUID, names_to = "sample_id", values_to = "Abundance")
View(dna.norarefied.long)

dna.norarefied.long.PA<- mutate(dna.norarefied.long, PA = if_else(Abundance>0, 1, 0))
View(dna.norarefied.long.PA)
write.csv(dna.norarefied.long.PA, 'dna.norarefied.long.PA.csv') 
##remove abundance column and read back in or do it using tidyverse directly
dna.norarefied.long.PA <- dna.norarefied.long.PA %>% select("OTUID","sample_id","PA")

#dna.norarefied.long.PA< read.csv("dna.norarefied.long.PA.csv")
dna.norarefied.long.PA.wide<-pivot_wider(dna.norarefied.long.PA, names_from = sample_id, values_from = PA, values_fill = 0)
View(dna.norarefied.long.PA.wide)
write.csv(dna.norarefied.long.PA.wide, 'dna.norarefied.long.PA.wide.csv')


##occupancy
##calculating occupancy of all taxa
otu_PA<-read.csv("dna.norarefied.long.PA.wide.csv") #removed first column with numbers and removed the first column name OTUID
otu_PA.df<-as.data.frame(otu_PA)
#checking if the factor is class or numeric
str(otu_PA)
#need to make numeric
otu_PA_transform<-transform(otu_PA.df, occupancy=rowSums(sapply(otu_PA.df[-1], as.numeric))/length(colnames(otu_PA.df[-1])))
View(otu_PA_transform)
#View(otu_PA_transform)          
write.csv(otu_PA_transform,"dna_norarefied_occupancy.csv")
otu_PA_transform<- read.csv("dna_norarefied_occupancy.csv") ##removed first columnn with number and added back column name OTUID
#dim(otu_PA_transform)

dna.norarefied.occ <- otu_PA_transform %>% select ("OTUID","occupancy")
View(dna.norarefied.occ)
dna.norarefied.occ <- dna.norarefied.occ %>% mutate(Row.names=OTUID)
View(dna.norarefied.occ)
View(phantoms_select)
phantoms_new <- merge(phantoms_select,dna.norarefied.occ, by="Row.names")
View(phantoms_new)

write.csv(phantoms_new, "phantom_taxa_occupancy_in_nonrarefied_dataset.csv")

phantoms_new<- phantoms_new %>% mutate(percent_occ=occupancy*100)

# Draw with black outline, white fill
plot1<- ggplot(phantoms_new, aes(x=percent_occ)) +
  geom_histogram(binwidth=0.5, colour="black", fill="white")
plot1
ggsave(filename = "phantom_occupancy_in_unrarefied_countdata.tiff", plot = plot1,
       width = 17,
       height = 15, units = c("cm"),
       dpi = 300)

# Histogram overlaid with kernel density curve
plot2<- ggplot(phantoms_new, aes(x=percent_occ)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=0.5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
plot2
ggsave(filename = "phantom_occupancy_in_unrarefied_density.tiff", plot = plot2,
       width = 17,
       height = 15, units = c("cm"),
       dpi = 300)



####################################################################################################
####################### Microbiome analysis of active community data using phyloseq ################
####################################################################################################



otu = read.csv("active_final_DNAabund_phantoms_threshold_detection5%.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 

any(rowSums(otu[])<1)
#FALSE
any(rowSums(otu[])<=1)
#TRUE, singletons present, KEEPing singletons


otu1 = otu[rowSums(otu[])>1, ,drop=FALSE]
otu1 ##no singletons
View(otu1) ##9962 OTUs in this dataset with no singletons (singletons defined here as OTUs present only in one sample in one copy, rowsums==1)
##We will keep the singletons

tax = read.csv("taxonomy-dn-99-edited.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("sample_metadata_dna_cdna.csv", sep=",", row.names=1)
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

any(taxa_sums(phyloseq_merged) <= 1) #true

any(taxa_sums(phyloseq_merged) < 1) #false


if (!requireNamespace("BiocManager", quietly = TRUE)) # Install BiocManager if not already installed. 
  install.packages("BiocManager")

BiocManager::install("microbiome")

# From GitHub
devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
microbiome::summarize_phyloseq(phyloseq_merged)

## Min. number of reads = 1915"

## Max. number of reads = 6024"

## Total number of reads = 910222"

## Average number of reads = 4293.5"

## Median number of reads = 4294"

## Sparsity = 0.903891699819058"

## Any OTU sum to 1 or less? YES"

## Number of singletons = 2426"

## Percent of OTUs that are singletons (i.e. exactly one read detected across all samples) 19.5834678721343"


## Number of sample variables are: 7"

## "treatment"     "crop"          "harvest"       "drought"       "planted"       "summarise"     "summarise_new"


# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(phyloseq_merged))
sum(sample_sum_df)
##all 212 samples "active" from dna abundance table-- 910222  (with singletons)
# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

smin <- min(sample_sums(phyloseq_merged))
smin # 1915
smean <- mean(sample_sums(phyloseq_merged))
smean # 4293.5
smax <- max(sample_sums(phyloseq_merged))
smax # 6024



#### Rarefaction curves

library(scales)
library(vegan)
library(dplyr)
rarecurve(t(otu_table(phyloseq_merged)), step=20, label = FALSE, col="cyan4") ##NOT USED IN FINAL MANUSCRIPT

###################################################
## Bacterial active community microbiome analysis
###################################################


########## Stacked barplots ##############

metadata<- read.csv("sample_metadata_dna_cdna.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones, this file includes all samples, bean and switchgrass
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final

phyloseq_class <- phyloseq_merged_final %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_class)

phyloseq_class$drought_new = factor(phyloseq_class$drought, levels=c("pre-drought","drought","well-watered"))
View(phyloseq_class)
devtools::install_github("doehm/evoPalette")
library(evoPalette)
launch_evo_palette()

palette_new41 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#fcfc00","#ffff9e","#ffb7ef","#fa7efc","#ae09ea","#521899",
                  "#a0fffc","#1e0047", "#CBD588", "#599861", "#508578","#FFD700","#FFA500",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD",
                  "#D14285")


library(evoPalette)

install.packages("ggh4x")

  
# install.packages("devtools")
devtools::install_github("teunbrand/ggh4x")

library(ggh4x)
plot<- ggplot(phyloseq_class, aes(x= treatment, y = Abundance, fill =Class)) + 
 # facet_grid(crop~planted+drought_new, scales="free") + theme(strip.text.x = element_text(size = 25))+ 
  facet_grid2(crop~planted+drought_new, scales="free", independent = "x")+ theme(strip.text.x = element_text(size = 35))+ 
  theme(strip.text.y = element_text(size = 35))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new41))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=25),
                                                     axis.title.y=element_text(size=35)) + 
  theme(legend.title = element_text(size=35))+ theme(legend.text = element_text(size=25))
plot 

###Figure S7

ggsave(filename = "barplot_class_decontam_all_removeempties2_newlabel.tiff", plot = plot,
       width = 75,
       height = 40, units = c("cm"),
       dpi = 300)


########
###Kruskall Wallis Test
########

metadata<- read.csv("sample_metadata_dna_cdna.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
#subset.metadata<- subset(metadata, crop=="switchgrass")
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final

phyloseq_merged_final


phyloseq_class <- phyloseq_merged_final %>% ##all 212 samples 
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} )%>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."


library(devtools)
install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq) ##under developement so didnt work
library(devtools)
install.packages("remotes")
remotes::install_github("TBrach/MicrobiomeX2") ##does not work

kruskal.test(phyloseq_class$Abundance ~ phyloseq_class$crop)


library(tidyverse)
physeq_bean<- phyloseq_class%>% filter(crop=="bean")


kruskal.test(physeq_bean$Abundance ~ physeq_bean$drought)


pairwise.wilcox.test(physeq_bean$Abundance, physeq_bean$drought, p.adjust.method = "fdr")



kruskal.test(physeq_bean$Abundance ~ physeq_bean$planted)

kruskal.test(physeq_bean$Abundance ~ physeq_bean$harvest)


physeq_sw<- phyloseq_class%>% filter(crop=="switchgrass")


kruskal.test(physeq_sw$Abundance ~ physeq_sw$drought)


kruskal.test(physeq_sw$Abundance ~ physeq_sw$planted)



kruskal.test(physeq_sw$Abundance ~ physeq_sw$harvest)



##subset to individual classes otherwise it is not making sense

##use simper to see which OTUs differ by crop, within crop which classes differ by planted , drought treatments.

##*********SIMPER RESULTS NOT REPORTED IN FINAL MANUSCRIPT
## Simper analysis, extract abundance matrix from a phyloseq object

##keep the abundances as raw abundances, not relative abundances. data already normalized to equal reads across samples

physeq_class_all <- phyloseq_merged_final %>%
  tax_glom(taxrank = "Class")   %>% 
  #transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()
# agglomerate at family/phylum level

View(physeq_class_all)

physeq_class<- physeq_class_all %>% select("Abundance", "Class", "Sample")

#write.csv(physeq_class, "physeq_class_simper.csv")

##Class"unclutured" has 2 values for each sample, from Desulfobacterota and Armatimonadota, so removed uncultured row for simper analysis 
##to do this I first had to rename uncultured in the csv file to _uncultured other wise it wasnt working ##i made an error with phyloseq object names so was not working, no need to rename to _uncultured, saves a step

#physeq_class<- read.csv("physeq_class_simper.csv")

physeq_class <- physeq_class %>% filter(Class!=" uncultured")

View(physeq_class)
physeq_wide_class<- pivot_wider(physeq_class, names_from=Sample, values_from=Abundance)

View(physeq_wide_class)

##identify duplicates (this code showed that the uncultured row was the problem)

physeq_wide_duplicates<- physeq_class %>%
  dplyr::group_by(Class, Sample) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 
View(physeq_wide_duplicates)

write.csv(physeq_wide_class, "physeq_wide_class_simper.csv") ##did sum calculations to see read counts per sample


physeq_wide_class_mat <- as.data.frame(physeq_wide_class)
rnames <-physeq_wide_class_mat[,1] 
rownames(physeq_wide_class_mat) <- rnames  # assign row names
View(physeq_wide_class_mat)
physeq_wide_class_new <- physeq_wide_class_mat %>% select(!Class)
View(physeq_wide_class_new)

physeq_wide_class_trans <- t(physeq_wide_class_new)
View(physeq_wide_class_trans)

#physeq_wide_class_trans <- cbind(SampleID = rownames(physeq_wide_class_trans), physeq_wide_class_trans, row.names=FALSE)
#View(physeq_wide_class_trans)

#rownames(physeq_wide_class_trans) <- NULL

#metadata<- read.csv("metadata_simper.csv", row.names=1) ##cannot use this as sample numbers dont match has to be 212 in metadata also
#View(metadata)
library(vegan)

sample_data = data.frame(sample_data(phyloseq_merged_final))

View(sample_data)
# running the simper analysis on the dataframe and the variable of interest 
crop_simper <- simper(physeq_wide_class_trans, sample_data$crop, permutations = 100)

# printing the top OTUs
print(crop_simper)


phyloseq_class <- phyloseq_merged_final %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()                                        # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  #arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

##filter taxa to Gamma, Actino, Alpha as per simper 

phyloseq_class_gamma<- phyloseq_class %>% filter(Class==" Gammaproteobacteria")

kruskal.test(phyloseq_class_gamma$Abundance ~ phyloseq_class_gamma$crop)


phyloseq_class_actino<- phyloseq_class %>% filter(Class==" Actinobacteria")

kruskal.test(phyloseq_class_actino$Abundance ~ phyloseq_class_actino$crop)


phyloseq_class_alpha<- phyloseq_class %>% filter(Class==" Alphaproteobacteria")

kruskal.test(phyloseq_class_alpha$Abundance ~ phyloseq_class_alpha$crop)


##SIMPER switchgrass
##subset to switchgrass drought and ww, this HAS to be done for both metadata and OTU file separately

##switchgrass drought
meta.sw <- sample_data[which(sample_data$crop=="switchgrass"),]

meta.sw.dr <- meta.sw[which(meta.sw$drought=="drought"),]


OTU.sw.dr <- physeq_wide_class_trans[which(sample_data$crop=="switchgrass" & sample_data$drought=="drought"),]



View(meta.sw.dr)

# running the simper analysis on the dataframe and the variable of interest 
planted.sw.simper <- simper(OTU.sw.dr, meta.sw.dr$planted, permutations = 100)

# printing the top OTUs
print(planted.sw.simper)


phyloseq_class_sw_dr_gamma<- phyloseq_class %>% filter(crop=="switchgrass" & drought=="drought" & Class==" Gammaproteobacteria")

kruskal.test(phyloseq_class_sw_dr_gamma$Abundance ~ phyloseq_class_sw_dr_gamma$planted)

View(phyloseq_class_sw_dr_gamma)


phyloseq_class_sw_dr_actino<- phyloseq_class %>% filter(crop=="switchgrass" & drought=="drought" & Class==" Actinobacteria")

kruskal.test(phyloseq_class_sw_dr_actino$Abundance ~ phyloseq_class_sw_dr_actino$planted)




phyloseq_class_sw_dr_alpha<- phyloseq_class %>% filter(crop=="switchgrass" & drought=="drought" & Class==" Alphaproteobacteria")

kruskal.test(phyloseq_class_sw_dr_alpha$Abundance ~ phyloseq_class_sw_dr_alpha$planted)


##switchgrass well-watered

meta.sw <- sample_data[which(sample_data$crop=="switchgrass"),]

meta.sw.ww <- meta.sw[which(meta.sw$drought=="well-watered"),]


OTU.sw.ww <- physeq_wide_class_trans[which(sample_data$crop=="switchgrass" & sample_data$drought=="well-watered"),]



View(meta.sw.ww)

# running the simper analysis on the dataframe and the variable of interest 
planted.sw.simper <- simper(OTU.sw.ww, meta.sw.ww$planted, permutations = 100)

# printing the top OTUs
print(planted.sw.simper)

phyloseq_class_sw_ww_gamma<- phyloseq_class %>% filter(crop=="switchgrass" & drought=="well-watered" & Class==" Gammaproteobacteria")

kruskal.test(phyloseq_class_sw_ww_gamma$Abundance ~ phyloseq_class_sw_ww_gamma$planted)


phyloseq_class_sw_ww_actino<- phyloseq_class %>% filter(crop=="switchgrass" & drought=="well-watered" & Class==" Actinobacteria")

kruskal.test(phyloseq_class_sw_ww_actino$Abundance ~ phyloseq_class_sw_ww_actino$planted)

phyloseq_class_sw_ww_alpha<- phyloseq_class %>% filter(crop=="switchgrass" & drought=="well-watered" & Class==" Alphaproteobacteria")

kruskal.test(phyloseq_class_sw_ww_alpha$Abundance ~ phyloseq_class_sw_ww_alpha$planted)





#### SIMPER bean 
##bean drought

meta.bean <- sample_data[which(sample_data$crop=="bean"),]

meta.bean.dr <- meta.sw[which(meta.sw$drought=="drought"),]


OTU.bean.dr <- physeq_wide_class_trans[which(sample_data$crop=="bean" & sample_data$drought=="drought"),]


View(meta.bean.dr)

# running the simper analysis on the dataframe and the variable of interest 
planted.bean.simper <- simper(OTU.bean.dr, meta.bean.dr$planted, permutations = 100)

# printing the top OTUs
print(planted.bean.simper)

  
phyloseq_class_bean_dr_gamma<- phyloseq_class %>% filter(crop=="bean" & drought=="drought" & Class==" Gammaproteobacteria")

kruskal.test(phyloseq_class_bean_dr_gamma$Abundance ~ phyloseq_class_bean_dr_gamma$planted)


phyloseq_class_bean_dr_actino<- phyloseq_class %>% filter(crop=="bean" & drought=="drought" & Class==" Actinobacteria")

kruskal.test(phyloseq_class_bean_dr_actino$Abundance ~ phyloseq_class_bean_dr_actino$planted)


phyloseq_class_bean_dr_alpha<- phyloseq_class %>% filter(crop=="bean" & drought=="drought" & Class==" Alphaproteobacteria")

kruskal.test(phyloseq_class_bean_dr_alpha$Abundance ~ phyloseq_class_bean_dr_alpha$planted)



##bean well-watered
meta.bean <- sample_data[which(sample_data$crop=="bean"),]

meta.bean.ww <- meta.bean[which(meta.bean$drought=="well-watered"),]


OTU.bean.ww <- physeq_wide_class_trans[which(sample_data$crop=="bean" & sample_data$drought=="well-watered"),]



View(meta.bean.ww)

# running the simper analysis on the dataframe and the variable of interest 
planted.bean.simper <- simper(OTU.bean.ww, meta.bean.ww$planted, permutations = 100)

# printing the top OTUs
print(planted.bean.simper)


phyloseq_class_bean_ww_gamma<- phyloseq_class %>% filter(crop=="bean" & drought=="well-watered" & Class==" Gammaproteobacteria")

kruskal.test(phyloseq_class_bean_ww_gamma$Abundance ~ phyloseq_class_bean_ww_gamma$planted)


phyloseq_class_bean_ww_actino<- phyloseq_class %>% filter(crop=="bean" & drought=="well-watered" & Class==" Actinobacteria")

kruskal.test(phyloseq_class_bean_ww_actino$Abundance ~ phyloseq_class_bean_ww_actino$planted)


phyloseq_class_bean_ww_alpha<- phyloseq_class %>% filter(crop=="bean" & drought=="well-watered" & Class==" Alphaproteobacteria")

kruskal.test(phyloseq_class_bean_ww_alpha$Abundance ~ phyloseq_class_bean_ww_alpha$planted)

 

######### Principal Coordinate analysis #########


## PCoaA analysis all samples bean and switchgrass

metadata<- read.csv("sample_metadata_dna_cdna.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final


##pcoa
phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged_final, 
  method = "PCoA", 
  distance = "bray"
)
##plot by crop
pcoa_crop<-plot_ordination(
  physeq = phyloseq_merged_final,##changed from phyloseq_rarefy
  ordination = phyloseq_pcoa,
  color = "crop",
  shape = "crop",
  title = "PCoA Switchgrass and Bean") + 
  #scale_fill_manual(name="drought", limits = c("pre-drought", "drought", "well-watered"),
  #values = c("#a65628", "red", "#ffae19",
  #"#4daf4a", "#1919ff", "darkorchid3", "magenta", "black" ))+
  scale_color_manual(values = c("#CC79A7", "#0072B2"))+
  geom_point(aes(color = crop, size=harvest)) +scale_shape_manual(values=c(16,17))+
  scale_size_manual(values = c(2,3,4,5,6,7))+
  geom_point() + theme(plot.title = element_text(hjust = 0.5))
library(ggplot2)
pcoa_crop 

ggsave(filename="pcoa-sw-bean-new.TIFF", plot=pcoa_crop, width=6.8, height=6, units="in", dpi=720)


## PCoA analysis within crop

## Plot by treatment (drought, planted levels) Bean

metadata<- read.csv("sample-metadata-switchgrass-cdna.csv")
metadata<- read.csv("sample-metadata-bean-cdna.csv")  ##change to bean when looking at bean samples 
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final


#pcoa

phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged_final, 
  method = "PCoA", 
  distance = "bray"
)

pcoa_bean<-plot_ordination(
  physeq = phyloseq_merged_final,##changed from phyloseq_rarefy
  ordination = phyloseq_pcoa,
  color = "drought",
  shape = "planted",
  title = "PCoA Bean"
) +  #scale_fill_manual(name="drought", limits = c("pre-drought", "drought", "well-watered"),
                                                    #values = c("#a65628", "red", "#ffae19",
                                                  #"#4daf4a", "#1919ff", "darkorchid3", "magenta", "black" ))+
  scale_color_manual(name="drought", limits = c("pre-drought", "drought", "watered"), values = c("#D55E00","#009E73","#56B4E9")) +
  geom_point(aes(color = drought, size=harvest)) +
  scale_shape_manual(values=c(16,17))+
  scale_size_manual(values = c(2,3,4,5,6,7))+
  theme(plot.title = element_text(hjust = 0.5))

pcoa_bean

ggsave(filename="pcoa-bean-new.TIFF", plot=pcoa_bean, width=6.8, height=6, units="in", dpi=720)

# Plot by treatment (drought, planted levels) Switchgrass

#pcoa

phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged_final, 
  method = "PCoA", 
  distance = "bray"
)
pcoa_sw<-plot_ordination(
  physeq = phyloseq_merged_final,##changed from phyloseq_rarefy
  ordination = phyloseq_pcoa,
  color = "drought",
  shape = "planted",
  title = "PCoA Switchgrass"
) +  #scale_fill_manual(name="drought", limits = c("pre-drought", "drought", "well-watered"),
                                                    #values = c("#a65628", "red", "#ffae19",
                                                  #"#4daf4a", "#1919ff", "darkorchid3", "magenta", "black" ))+
  scale_color_manual(name="drought", limits = c("pre-drought", "drought", "watered"), values = c("#D55E00","#009E73","#56B4E9")) +
  geom_point(aes(color = drought, size=harvest)) +
  scale_shape_manual(values=c(16,17))+
  scale_size_manual(values = c(2,3,4,5,6,7))+
  theme(plot.title = element_text(hjust = 0.5))

pcoa_sw

ggsave(filename="pcoa-sw-new.TIFF", plot=pcoa_sw, width=6.8, height=6, units="in", dpi=720)

##new publication figure

my_plot_list <- list(pcoa_crop, pcoa_bean, pcoa_sw)
plot_combined <- ggarrange(plotlist = my_plot_list, labels=c("A", "B", "C"), nrow=1, ncol=3,
                           common.legend = FALSE, heights = c(20, 20, 20)) + bgcolor("white") 

plot_combined

##Figure 4
ggsave(filename = "bean_sw_pcoa_combined.tiff", plot = plot_combined,
       width = 40,
       height = 11, units = c("cm"),
       dpi = 300)



#############################################
###Abundance dotplots over time##############
#############################################

metadata<- read.csv("bean-dr-ww-planted-unplanted.csv")

keep.samples <- as.vector(metadata$sample_id)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final

phyloseq_class <- phyloseq_merged_final %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_class)


##for dotplots 
phyloseq_class1<- phyloseq_class %>% group_by(summarise, Class, crop) %>% mutate(mean=mean(Abundance))


#Filter out less than 1% abundance 

phyloseq_class2<- phyloseq_class1 %>% filter(Class!="< 1% abund.")

View(phyloseq_class2) 

write.csv(phyloseq_class2, "bean_dr_ww_phyloseq_class2_dotplot_no1percent.csv") ##edited columns (drought and summarise new) so that pre-drought component exists for both drought and watered
phyloseq_class2<- read.csv("bean_dr_ww_phyloseq_class2_dotplot_no1percent.csv")


library(viridis)
bb<-c(0.300,0.250,0.200,0.150,0.100,0.050,0.001) # define breaks.
ll<-c("0.300","0.250","0.200","0.150","0.100","0.050","0.001") # labels. ##all samples representation checked for bean planted drought and watered such that the range does not miss any samples
plot_dotplot_beandrww<-
  ggplot(phyloseq_class2,aes(x=harvest,y=Class,size=mean,color=mean))+
  geom_point()+
  scale_x_discrete(limits = c("Day 0","Day 2","Day 3", "Day 4","Day 5", "Day 6"))+
  facet_grid(~summarise_new,scales="free")+ 
  theme(strip.text.x=element_text(size=10))+
  #scale_colour_gradient(low="#132B43",
                        #high="#56B1F7",limits=c(0.001,0.3),breaks=bb)+
  scale_colour_viridis(limits=c(0.001,0.3),breaks=bb)+
  scale_size_continuous(limits=c(0.001,0.3),labels=ll,breaks=bb, name="Mean Relative\nAbundance")+
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=9, angle = 90, vjust = 0.5, hjust=1))+
  #ggtitle("Bean Drought")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(color="Mean Relative\nAbundance")+
  ylab("Class") + 
  theme(axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10)) + 
  theme(legend.title = element_text(size=10))+ theme(legend.text = element_text(size=9))
  
plot_dotplot_beandrww


#Figure 6
ggsave(filename = "Fig6.tiff", plot = plot_dotplot_beandrww,
       width = 22 ,
       height = 16, units = c("cm"),
       dpi = 200)




#######################################
##### Statistical analyses ############
#######################################

##Community differences by crop type

metadata<- read.csv("sample_metadata_dna_cdna.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final


set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_final, method = "bray") 
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged_final)) 
View(sample_df)
set.seed(1)
adonis(phyloseq_bray ~ crop, data = sample_df)

##beta dispersion

metadata<- read.csv("sample_metadata_dna_cdna.csv")
subset.metadata_bean<- subset(metadata, crop=="bean") ##change to bean and sw as needed
subset.metadata_1<- subset.metadata %>% filter(drought!="pre-drought")##96 samples ##after knowing the posthoc tests I am removing one vraiable at a time to test the betadispersion with permutest to get DF, Pseudo F and dispersion, because post hoc Tukey does not give F statistic
subset.metadata_2<- subset.metadata %>% filter(drought!="watered")##56 samples
subset.metadata_3<- subset.metadata %>% filter(drought!="drought")##58 samples
subset.metadata_bean<- subset.metadata_bean %>% filter(drought!="pre-drought") ##remove pre-drought for time series testing
keep.samples <- as.vector(subset.metadata_bean$SampleID)
keep.samples

metadata<- read.csv("sample_metadata_dna_cdna.csv")
subset.metadata_sw<- subset(metadata, crop=="switchgrass") ##change to bean and sw as needed
subset.metadata_sw<- subset.metadata_sw %>% filter(drought!="pre-drought") ##remove pre-drought for time series dispersion
keep.samples <- as.vector(subset.metadata_sw$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
phyloseq_merged_final

set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_final, method = "bray") ##changed from _rarefy to _merged
phyloseq_bray


sample_df <- data.frame(sample_data(phyloseq_merged_final)) 
View(sample_df)


##betadispersion: Homogeneity of dispersion test by crop
set.seed(1)
beta_dispersion_crop <- betadisper(phyloseq_bray, sample_df$crop)
permutest(beta_dispersion_crop)


plot<- plot(beta_dispersion_crop, hull=FALSE, ellipse=TRUE, label=F, col=c("#CC79A7", "#0072B2"))
legend(x = "topleft", lty = c(4,6), text.font = 4, pch = c(16,17),
      col = c("#CC79A7", "#0072B2"),
       pt.bg = c("#CC79A7", "#0072B2"),text.col = "black", legend=c("bean", "switchgrass"))

##beta dispersion by drought ##subset to bean or switchgrass
set.seed(1)
beta_dispersion_drought <- betadisper(phyloseq_bray, sample_df$drought)
permutest(beta_dispersion_drought)

##Tukey post hoc when significant

drought_HSD <- TukeyHSD(beta_dispersion_drought)
drought_HSD

#plot
plot<- plot(beta_dispersion_drought, hull=FALSE, ellipse=TRUE, label=F, col=c("#D55E00","#009E73","#56B4E9"), pch=c(15,16,17))
legend(x = "topleft", lty = c(4,6), text.font = 4, pch = c(15,16,17),
       col = c("#D55E00","#009E73","#56B4E9"),
       pt.bg = c("#D55E00","#009E73","#56B4E9"),text.col = "black", legend=c("pre-drought", "drought", "watered"))                                                                         

##beta dispersion by planted ##subset to bean or switchgrass
set.seed(1)
beta_dispersion_planted <- betadisper(phyloseq_bray, sample_df$planted)
permutest(beta_dispersion_planted)

plot<- plot(beta_dispersion_planted, hull=FALSE, ellipse=TRUE, label=F, col=c("green", "blue"), pch=c(16,17))
legend(x = "topleft", lty = c(4,6), text.font = 4, pch = c(16,17),
       col = c("green", "blue"),
       pt.bg = c("green", "blue"),text.col = "black", legend=c("planted", "unplanted"))   

##by harvest ##subset to bean or switchgrass

set.seed(1)
beta_dispersion_harvest <- betadisper(phyloseq_bray, sample_df$harvest)
permutest(beta_dispersion_harvest)

##by summarise_new for each planted/drought/time condition #subset to bean or switchgrass
set.seed(1)
beta_dispersion_timeseries <- betadisper(phyloseq_bray, sample_df$summarise_new)
permutest(beta_dispersion_timeseries)
View(sample_df)

timeseries_HSD <- TukeyHSD(beta_dispersion_timeseries)
timeseries_HSD  


##betadispersion: Homogeneity of dispersion test
##beta <- betadisper(phyloseq_bray, sample_df$crop)
##beta
##set.seed(1) ##if set.seed is not set then result changes, due to permutations 
##beta_crop<- permutest(beta, pairwise=TRUE, permutations=999)
##beta_crop
##beta_crop_HSD<- TukeyHSD(beta) ##pairwise tukey test
##beta_crop_HSD

##plot(beta)
##boxplot(beta)



####### Bean  and switchgrass samples combined mantel correlations###########
  
  
# Calculate bray curtis distance matrix

phyloseq_merged ##212 samples and 12388 OTUs

set.seed(1)
phyloseq_bray <- phyloseq::distance(phyloseq_merged, method = "bray")
phyloseq_bray

library("ape")
set.seed(1)
random_tree = rtree(ntaxa(phyloseq_merged), rooted=TRUE, tip.label=taxa_names(phyloseq_merged))
plot(random_tree)
phyloseq_tree = merge_phyloseq(phyloseq_merged, random_tree)
phyloseq_tree


#calculate weighted unifrac distance matrix
phyloseq_wunifrac <- phyloseq::distance(phyloseq_tree, method = "wunifrac" )

##mantel tests to test correlation between wunifrac and BC (https://jkzorz.github.io/2019/07/08/mantel-test.html)
library(scales)
library(reshape2)
library(grid)
set.seed(1)
mantel_cor = mantel(phyloseq_bray, phyloseq_wunifrac, method = "spearman", permutations = 9999, na.rm = TRUE)
mantel_cor

mantel_cor = mantel(phyloseq_bray, phyloseq_wunifrac, method = "pearson", permutations = 9999, na.rm = TRUE)
mantel_cor

##permanova with wunifrac, comparing bean and switchgrass

sample_df <- data.frame(sample_data(phyloseq_merged)) 
View(sample_df)
set.seed(1)
adonis(phyloseq_wunifrac ~ crop, data = sample_df)


##plot correlation (NOT REPORTED IN FINAL MANUSCRIPT)



plot<- plot(phyloseq_bray,phyloseq_wunifrac,pch=16,cex=0.5,col="black",bty="l")



##prune to SW samples only
  
metadata<- read.csv("sample-metadata-switchgrass-cdna.csv")
keep.samples <- as.vector(metadata$sample_id)
keep.samples

phyloseq_merged_sw <- prune_samples(keep.samples, phyloseq_merged) ##changed from rarefy to merged
phyloseq_merged_sw


set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_sw, method = "bray")
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged_sw))
View(sample_df)

set.seed(1)
adonis(phyloseq_bray ~ planted, data = sample_df)

set.seed(1)
adonis(phyloseq_bray ~ drought, data = sample_df)

library(pairwiseAdonis)
set.seed(1)
pairwise.adonis2(phyloseq_bray ~ drought, data = sample_df)


  
##prune to bean samples only
  
metadata<- read.csv("sample-metadata-bean-cdna.csv")
keep.samples <- as.vector(metadata$sample_id)
keep.samples

phyloseq_merged_bean <- prune_samples(keep.samples, phyloseq_merged) #changed from rarefy to merged
phyloseq_merged_bean

set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_bean, method = "bray")
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged_bean))
View(sample_df)

set.seed(1)
adonis(phyloseq_bray ~ planted, data = sample_df) 
adonis(phyloseq_bray ~ drought, data = sample_df)

set.seed(1)
pairwise.adonis2(phyloseq_bray ~ drought, data = sample_df)


####################################################
########## Interaction Permanovas ##################
####################################################

##Switchgrass time versus drought treatment versus planted treatment interaction (using permanova)
  
metadata<- read.csv("SW-time-treat-nopredr.csv") ## SW planted and unplanted samples
keep.samples <- as.vector(metadata$sample_id)
keep.samples

phyloseq_sw_nopredr <- prune_samples(keep.samples, phyloseq_merged) 
phyloseq_sw_nopredr


set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_sw_nopredr, method = "bray")
phyloseq_bray


# Adonis test
sample_df <- data.frame(sample_data(phyloseq_sw_nopredr))
View(sample_df)
set.seed(1)
adonis(phyloseq_bray ~ planted*drought*harvest, by="margin", data = sample_df)



##Bean time versus drought treatment versus planted treatment interaction (using permanova)

metadata<- read.csv("bean-time-treat-nopredr.csv")
keep.samples <- as.vector(metadata$sample_id)
keep.samples

phyloseq_bean_nopredr <- prune_samples(keep.samples, phyloseq_merged) ##changed to _merged from _rarefy
phyloseq_bean_nopredr



set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_bean_nopredr, method = "bray")
phyloseq_bray


# Adonis test
sample_df <- data.frame(sample_data(phyloseq_bean_nopredr))
View(sample_df)

set.seed(1)
adonis(phyloseq_bray ~ planted*drought*harvest, by="margin", data = sample_df)




##linked resource webpage: https://stat.ethz.ch/pipermail/r-sig-ecology/2018-November/005830.html

############################################################################
###############Bray Curtis Similarity to pre-drought condition #############
############################################################################


##BC similarity to pre-drought over harvest time, bean drought and well watered samples (planted and unplanted)

metadata<- read.csv("bean-ww-planted-unplanted.csv") ##replace with "bean-drought-planted-unplanted.csv" when looking at drought samples

keep.samples <- as.vector(metadata$sample_id)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) 
phyloseq_merged_final 


sample_data(phyloseq_merged_final)

set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_final, method = "bray")
phyloseq_bray

install.packages("usedist")
library(usedist)## can be useful in this case for functions like dist_subset and dist_groups

install.packages("dendextend")
library(dendextend)
dist_long<-dist_long(phyloseq_bray)
dist_long
write.csv(dist_long, "dist_long_ww_bean_withsingletons_nosecondrarefaction.csv")
dist_long<-read.csv("dist_long_ww_bean_withsingletons_nosecondrarefaction.csv")##change dist_long column names to iso1 and iso2 for formatting, then read back in
metadata_new<-data.frame(sample_data(phyloseq_merged_final))
View(metadata_new)
write.csv(metadata_new, "metadata_dist_long_ww_bean_withsingletons_nosecondrarefaction.csv") 
metadata_new<- read.csv("metadata_dist_long_ww_bean_withsingletons_nosecondrarefaction_edited.csv")#change column names in metadata to group and sample_id, then read back in
install.packages("harrietr", dependencies = TRUE)
install.packages('devtools')
setRepositories(ind = 1:2)
devtools::install_github("andersgs/harrietr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")


library(harrietr)
library(ggtree)
update.packages("ggtree")
dist.df<- join_metadata(dist_long, metadata_new, isolate = 'sample_id', group = 'group', remove_ind=FALSE)
write.csv(dist.df,"dist_df_ww_bean_withsingletons_nosecondrarefaction.csv")##change planted.x and planted.y to planted_x and planted_y

##bean planted BC versus time 

distdfbean<- read.csv("dist_df_ww_bean_withsingletons_nosecondrarefaction.csv")
distdfbean.planted <- distdfbean  %>%
  filter(planted_x=="planted"&
           planted_y== "planted")

View(distdfbean.planted)

write.csv(distdfbean.planted, "distdfbean_planted_ww_withsingletons_nosecondrarefaction.csv")
beanpredr.planted<- distdfbean.planted  %>%
  filter(group_2=="pre-droughtplanted" &
           group_1!="pre-droughtplanted")

View(beanpredr.planted)

write.csv(beanpredr.planted, "beanpredr_planted_ww_withsingletons_nosecondrarefaction.csv") ##for well-watered treatment samples

write.csv(beanpredr.planted, "beanpredr_planted_dr_withsingletons_nosecondrarefaction.csv") ##for drought treatment samples

#bean unplanted BC versus time

distdfbean<- read.csv("dist_df_ww_bean_withsingletons_nosecondrarefaction.csv")
distdfbean.unplanted <- distdfbean  %>%
  filter(planted_x=="unplanted"&
           planted_y== "unplanted")

View(distdfbean.unplanted)

write.csv(distdfbean.unplanted, "distdfbean_unplanted_ww_withsingletons_nosecondrarefaction.csv")
beanpredr.unplanted<- distdfbean.unplanted  %>%
  filter(group_2=="pre-droughtunplanted" &
           group_1!="pre-droughtunplanted")

View(beanpredr.unplanted)



write.csv(beanpredr.unplanted, "beanpredr_unplanted_ww_withsingletons_nosecondrarefaction.csv") ##for well-watered treatment samples


write.csv(beanpredr.unplanted, "beanpredr_unplanted_dr_withsingletons_nosecondrarefaction.csv") ##for drought treatment samples




##diagram for BC similarity change by time, SW drought and well watered samples (planted and unplanted)


metadata<- read.csv("sw-drought-planted-unplanted.csv") ##change to "sw-ww-planted-unplanted.csv" for well-watered samples

keep.samples <- as.vector(metadata$sample_id)
keep.samples


phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) 
phyloseq_merged_final 

sample_data(phyloseq_merged_final)


set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_final, method = "bray")
phyloseq_bray

install.packages("usedist")
library(usedist)## can be useful in this case for functions like dist_subset and dist_groups

install.packages("dendextend")
library(dendextend)
dist_long<-dist_long(phyloseq_bray)
dist_long
write.csv(dist_long, "dist_long_dr_sw_withsingletons_nosecondrarefaction.csv")
dist_long<-read.csv("dist_long_dr_sw_withsingletons_nosecondrarefaction.csv")##change dist_long column names to iso1 and iso2 for formatting, then read back in
metadata_new<-data.frame(sample_data(phyloseq_merged_final))
View(metadata_new)
write.csv(metadata_new, "metadata_dist_long_dr_sw_withsingletons_nosecondrarefaction.csv") 
metadata_new<- read.csv("metadata_dist_long_dr_sw_withsingletons_nosecondrarefaction_edited.csv")#change column names in metadata to group and sample_id, then read back in
install.packages("harrietr", dependencies = TRUE)
install.packages('devtools')
setRepositories(ind = 1:2)
devtools::install_github("andersgs/harrietr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

library(harrietr)
library(ggtree)
update.packages("ggtree")
dist.df<- join_metadata(dist_long, metadata_new, isolate = 'sample_id', group = 'group', remove_ind=FALSE)
write.csv(dist.df,"dist_df_dr_sw_withsingletons_nosecondrarefaction.csv")##change planted.x and planted.y to planted_x and planted_y

##SW planted BC versus time 

distdfsw<- read.csv("dist_df_dr_sw_withsingletons_nosecondrarefaction.csv")
distdfsw.planted <- distdfsw  %>%
  filter(planted_x=="planted"&
           planted_y== "planted")

View(distdfsw.planted)

write.csv(distdfsw.planted, "distdfsw_planted_dr_withsingletons_nosecondrarefaction.csv")
swpredr.planted<- distdfsw.planted  %>%
  filter(group_2=="pre-droughtplanted" &
           group_1!="pre-droughtplanted") 

swpredr.planted1<- distdfsw.planted  %>%
  filter(group_1=="pre-droughtplanted" &
           group_2!="pre-droughtplanted") 

##combine swpredr.planted and swpredr.planted1 

swpredr.planted<-merge(swpredr.planted,swpredr.planted1, all.x = TRUE, all.y = TRUE)

View(swpredr.planted)

write.csv(swpredr.planted, "swpredr_planted_dr_withsingletons_nosecondrarefaction.csv") ##for drought treatment samples

write.csv(swpredr.planted, "swpredr_planted_ww_withsingletons_nosecondrarefaction.csv") ##for well-watered treatment samples


#SW unplanted BC versus time

distdfsw<- read.csv("dist_df_dr_sw_withsingletons_nosecondrarefaction.csv")
distdfsw.unplanted <- distdfsw  %>%
  filter(planted_x=="unplanted"&
           planted_y== "unplanted")

View(distdfsw.unplanted)

write.csv(distdfsw.unplanted, "distdfsw_unplanted_dr_withsingletons_nosecondrarefaction.csv")
swpredr.unplanted<- distdfsw.unplanted  %>%
  filter(group_2=="pre-droughtunplanted" &
           group_1!="pre-droughtunplanted")

swpredr.unplanted1<- distdfsw.unplanted  %>%
  filter(group_1=="pre-droughtunplanted" &
           group_2!="pre-droughtunplanted") 

##combine swpredr.planted and swpredr.planted1 

swpredr.unplanted<-merge(swpredr.unplanted,swpredr.unplanted1, all.x = TRUE, all.y = TRUE)

View(swpredr.unplanted)

write.csv(swpredr.unplanted, "swpredr_unplanted_dr_withsingletons_nosecondrarefaction.csv") ## for drought treatment samples

write.csv(swpredr.unplanted, "swpredr_unplanted_ww_withsingletons_nosecondrarefaction.csv") ## for well-watered treatment samples

##### Note: Merge the planted and unplanted datasets together - resulting in one dataset for bean drought (planted and unplanted) #######
#### and bean well-watered (planted and unplanted). Total 4 dataframes for bean and switchgrass #########################################


##### The merged file should be edited to indicate the replicates for each sample compared in the similarity matrix #######################################
##### and then only identical replicate comparisons should be retained, such that 1-1 rep (pre-drought - harvest 1), 2-2 rep, (predrought - harvest 1) ####
##### etc are retained. 1-2, 1-3 and other comparisons should not be included ##############################################################################
##### Add a column to this merged file to also indicate time (days) of harvest sampling under column "group" ###############################################
##### Also add a column to indicate similarity index as Bray curtis gives dissimilarity values, to do this simply subtract the "distance" value from 1 #######
##### Add a concat column to merge the identical replicate comparisons #####################################################################################

###### See example merge and merge-edited files on HPCC

##### Use the merged file to select the 1-1, 2-2 comparisons

data<-read.csv("bean-ww-BC-time-planted-unplanted-merge-withsingletons_nosecondrarefaction.csv") ##merged bean ww planted and unplanted dataset, similarly do for other spreadsheets

target<- c("11","22","33","44","55")

data1<-data %>% filter(concat %in% target)

View(data1)

write.csv(data1, "bean-ww-BC-time-planted-unplanted-merge-edited-withsingletons_nosecondrarefaction.csv") ###see sheets used for analysis here in the folder specified on HPCC

##plotting BC sim to pre-drought over harvest time

library(viridis)


###new code used for final plot 

library(viridis)
library(ggplot2)
sw_dr_BC_time<- read.csv("bean-sw-dr-ww-BC-time-merge-edited.csv") ##combine all datasheets obtained from code in line 2315
View(sw_dr_BC_time)
sw_dr_BC_time<- rename(sw_dr_BC_time, harvest = group)

sw_dr_BC_time$harvest <- recode(sw_dr_BC_time$harvest, 
                                `2` = 'Day 2',
                                `3` = 'Day 3',
                                `4` = 'Day 4',
                                `5` ='Day 5',
                                `6` = 'Day 6')


sw_dr_BC_time$drought <- recode(sw_dr_BC_time$drought, 'well-watered' = 'watered')

plot_bean_sw<- ggplot(sw_dr_BC_time, aes(x= harvest,y=similarity, linetype=drought, shape=drought,color=planted_x)) + theme_classic()+ scale_color_manual(values=c("#39568CFF", "#55C667FF")) +
  scale_shape_manual(values=c(1,16))+ #scale_size(range = c(1,10), # point size range
  #label = scales::dollar)+
  #expand_limits(shape = seq(2, 6, by = 1))+
  #scale_alpha_manual(values = c(2,3,4,5,6))+
  geom_point(aes(size=harvest)) +  scale_size_manual(values = c(2,3,4,5,6))+ ylim(0.1, 0.6)+ geom_smooth(method = lm, se=FALSE, aes(group=summarise), inherit.aes=TRUE) + #ggtitle("Bean")
  facet_grid(~crop)+theme(
    strip.text.x = element_text(
      size = 16
    ))+
  theme(plot.title = element_text(size= 18, hjust = 0.5))+
  labs(y = "Bray-Curtis similarity to pre-drought", color='planted')+ xlab("harvest")+theme(axis.text.x = element_text(colour = "black", size = 17,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(size = 16, colour = "black"), 
         axis.title.y= element_text(size=16, color="black"),
         axis.title.x= element_text(size=16, color="black"),
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =14, colour = "black"),
         legend.text = element_text(size = 12, colour = "black")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")), color= guide_legend(order = 1)) 


plot_bean_sw


##Figure S8
ggsave(filename = "bean_sw_drought_watered_compiled_new2.tiff", plot = plot_bean_sw,
       width = 25,
       height = 15, units = c("cm"),
       dpi = 300)



############****note: ONLY LINEAR MODEL R2 AND P VALUE REPORTED FOR EACH TREATMENT*TIME SERIES 
##ANCOVA for BC similarity versud time, catergorical variables- planted, covariate: time, categorical variable 2: crop, cat var 3: drought, dependent variable similarity


ancova<- read.csv("ancova_metadata_withsingletons_nosecondrarefaction.csv")

##type 3:
options(contrasts = c("contr.sum", "contr.poly"))
lm.crop<- lm(similarity~crop*planted*drought*days, data=ancova)
mod1.type3 <- car::Anova(lm.crop, type=3)
mod1.type3

mod1 <- aov(similarity~crop*planted*drought*days, data=ancova) ##test interaction
summary(mod1)

##if crop level is included minor interactions exist between planted:drought and crop:planted 

##separate bean and SW

###with_singletons_no_second_rarefaction

ancova_bean <-subset(ancova, crop=="bean")
mod1 <- aov(similarity~planted*drought*days, data=ancova_bean) ##test interaction
summary(mod1)

ancova_bean_dr <-subset(ancova, crop=="bean" & drought=="drought")
mod1 <- aov(similarity~planted*days, data=ancova_bean_dr) ##test interaction
summary(mod1)

##type 3:
options(contrasts = c("contr.sum", "contr.poly"))
lm.beandr<- lm(similarity~planted*days, data=ancova_bean_dr)
mod1.type3 <- car::Anova(lm.beandr, type=3)
mod1.type3

##bean dr type 3 : produces same interaction without type 3 so just stick to without type 3 and check for all conditions

##bean drought results run

##bean well-watered results run

##bean all sample results run

##no interaction for bean , slope is not different between treatments
#no interaction between time and planted , so no difference in slope

##withsingletons no second rarefaction

##bean all

mod2 <- aov(similarity~planted+drought+days, data=ancova_bean) ##if no interaction, test differences in intercept within categorical variable eg. sex, if significant then subset by sex and look for effects of covariate (snout)

summary(mod2)


## no difference in intercept: planted levels, difference seen in intercept by drought levels

anova(mod1,mod2)


##removing interaction term for bean doesnt significantly alter the fit of the model.
##therefore we can conclude that most parsimonious model is mod2.


##bean drought

mod2 <- aov(similarity~planted+days, data=ancova_bean_dr)
summary(mod2)


## no difference in intercept

anova(mod1,mod2)

##no difference when interaction is removed


options(contrasts = c("contr.sum", "contr.poly"))
lm.bean<- lm(similarity~planted+days, data=ancova_bean_dr)
mod2 <- car::Anova(lm.bean, type=2)
mod2 


anova (mod1, mod2) ##same warning as with switchgrass even though anova type 2 values are the same as without type 2.

##bean ww
mod2 <- aov(similarity~planted+days, data=ancova_bean_ww)

summary(mod2)


##no difference in intercept

anova(mod1,mod2)



##no difference when interaction is removed

options(contrasts = c("contr.sum", "contr.poly"))
lm.bean<- lm(similarity~planted+days, data=ancova_bean_ww)
mod2 <- car::Anova(lm.bean, type=2)
mod2

  

##Bean subset to drought and well watered ##with singletons no second rarefaction

dr_bean <- subset(ancova_bean, drought=="drought" )


reg1 <- lm(similarity~planted_y/days -1, data=dr_bean); summary(reg1) ##ONLY THESE RESULTS REPORTED IN FINAL MANUSCRIPT

##well-watered results

##drought results



##with singletons no second rarefaction

##switchgrass all

ancova_switchgrass<- subset(ancova, crop=="switchgrass")
mod1 <- aov(similarity~planted*drought*days, data=ancova_switchgrass) ##test interaction
summary(mod1)

options(contrasts = c("contr.sum", "contr.poly"))
lm.sw<- lm(similarity~planted*drought*days, data=ancova_switchgrass)
mod1.type3 <- car::Anova(lm.sw, type=3)
mod1.type3

##type 3 sw all results

##type 1 default sw all results


###no interaction for SW, slope is not different between treatments
##for SW
#no interaction between time and planted , so no difference in slope
  
#withsingletons no second rarefaction
  
mod2 <- aov(similarity~planted+drought+days, data=ancova_switchgrass) ##if no interaction, test differences in intercept within categorical variable eg. sex, if significant then subset by sex and look for effects of covariate (snout)

summary(mod2)



##significant differences in intercept for SW data **** (planted factor considered here)

anova(mod1,mod2) ##does removing the interaction affect the fit of the model? if not then proceed with model 2 and subset the data based on sex
##sepaarte regression for categirical variables testing effects of snout on y var, due to effects of sex on intercept sepaarte sexes and assess effect of snout on dependent variable (pelvic)


##with singletons no second rarefaction

  
##removing interaction term for SW doesnt significantly alter the fit of the model.
##therefore we can conclude that most parsimonious model is mod2.
  
##switchgrass drought
  
ancova_switchgrass_dr<- subset(ancova, crop=="switchgrass" & drought=="drought")
mod1 <- aov(similarity~planted*days, data=ancova_switchgrass_dr) ##test interaction
summary(mod1)




#spot checking type 3 anova :
options(contrasts = c("contr.sum", "contr.poly"))
lm.sw<- lm(similarity~planted*days, data=ancova_switchgrass_dr)
mod1.type3 <- car::Anova(lm.sw, type=3)
mod1.type3  




mod2 <- aov(similarity~planted+days, data=ancova_switchgrass_dr)
summary(mod2)  

 


anova (mod1, mod2)



options(contrasts = c("contr.sum", "contr.poly"))
lm.sw<- lm(similarity~planted+days, data=ancova_switchgrass_dr)
mod2 <- car::Anova(lm.sw, type=2)
mod2 


  

anova(mod1, mod2)




##switchgrass well-watered

ancova_switchgrass_ww<- subset(ancova, crop=="switchgrass" & drought=="well-watered")
mod1 <- aov(similarity~planted*days, data=ancova_switchgrass_ww) ##test interaction
summary(mod1)



mod2 <- aov(similarity~planted+days, data=ancova_switchgrass_ww)
summary(mod2)  


options(contrasts = c("contr.sum", "contr.poly"))
lm.sw<- lm(similarity~planted+days, data=ancova_switchgrass_ww)
mod2 <- car::Anova(lm.sw, type=2)
mod2 



anova(mod1, mod2) ##results without type 2, as both type 1 and type 2 gave same results
##however, this is because data is balanced here which was not the case for sw dr. 
##even though anova results were same here with or without type 2 , i still got the warning message as shown above saying : In anova.lmlist(object, ...) :
#models with response ‘"NULL"’ removed because response differs from model 1



##SW subset to drought and well watered ##with singletons no second rarefaction
  
##should I use these values for the slope values ("Estimate" values from table)?
##or should i not nest days within planted and just do a summary call on lm of individual subsets of data and get slopes
##DISCLAIMER:: The nested slopes calculated as nested within days give the exact same slope values as the one done separately.
##So just work with the subsetted data.

ww_sw <- subset(ancova_switchgrass, drought=="well-watered" )

##ONLY THESE RESULTS REPORTED IN FINAL MANUSCRIPT
reg1 <- lm(similarity~planted_y/days -1, data=ww_sw); summary(reg1) ##subset data, nest days within planted

##drought SW results

##well-watered SW results




###########################################################
###### Shoot biomass and gravimetric soil moisture #######
###########################################################


##SW biomass
sw_shootmass<-read.csv("Switchgrass_biomass.csv", sep=',', header=T, row.names=1)
View(sw_shootmass)
#rename column to remove units
library(tidyverse)
sw_shootmass <- sw_shootmass %>% filter(Treatment != "Pre-Drought") 
colnames(sw_shootmass)<-c("Harvest","Treatment","ShootMass")

## Bean biomass
bean_shootmass<-read.csv("Bean_biomass.csv", sep=',', header=T, row.names=1)
View(bean_shootmass)
#rename column to remove units
library(tidyverse)
bean_shootmass <- bean_shootmass %>% filter(Treatment != "Pre-Drought") 
colnames(bean_shootmass)<-c("Harvest","Treatment","ShootMass")


##3 way anova

View(sw_shootmass)
View(bean_shootmass)

library(ggplot2)
library(grid)
library(lattice)
install.packages("multcompView")
library(multcompView)
library(lsmeans)

options(contrasts = c("contr.sum", "contr.poly"))

model = lm(ShootMass~ Treatment + Harvest + Treatment:Harvest,
           data=bean_shootmass) ##replace with sw_shootmass when doing switchgrass analysis
shapiro.test(resid(model)) #W statistic > 0.9 means normality okay #check normality to flag outliers 
ee= as.matrix(resid(model))
qqnorm (ee) #look for outliers (far away from curve) and test for normality 
ee
install.packages("car")
library(car)
#options(contrasts = c("contr.sum", "contr.poly"))
Anova(model, type="III")



##########Plotting ########
########################


##SW biomass
sw_shootmass<-read.csv("Switchgrass_biomass.csv", sep=',', header=T, row.names=1)
View(sw_shootmass)
#rename column to remove units
library(tidyverse)
#sw_shootmass <- sw_shootmass %>% filter(Treatment != "Pre-Drought") #remove this code when plotting
colnames(sw_shootmass)<-c("Harvest","Treatment","ShootMass")

## Bean biomass
bean_shootmass<-read.csv("Bean_biomass.csv", sep=',', header=T, row.names=1)
View(bean_shootmass)
#rename column to remove units
library(tidyverse)
#bean_shootmass <- bean_shootmass %>% filter(Treatment != "Pre-Drought") #remove this code when plotting
colnames(bean_shootmass)<-c("Harvest","Treatment","ShootMass")

library(Rmisc)
sw_data1 <- summarySE(sw_shootmass, measurevar="ShootMass", groupvars=c("Harvest","Treatment"))
View(sw_data1)
#reshape data so Pre-treatment comes first in legend
sw_data1$Treatment <- factor(sw_data1$Treatment,levels = c("pre-drought","drought","watered"), ordered=TRUE)
library(ggplot2)
sw_data1$Treatment

bean_data1 <- summarySE(bean_shootmass, measurevar="ShootMass", groupvars=c("Harvest","Treatment"))
View(bean_data1)
#reshape data so Pre-treatment comes first in legend
bean_data1$Treatment <- factor(bean_data1$Treatment,levels = c("pre-drought","drought","watered"), ordered=TRUE)
library(ggplot2)
bean_data1$Treatment


## Gravimetric soil moisture

sw_soilmoisture<-read.csv("Switchgrass_soil.csv", sep=',', header=T, row.names=1)
View(sw_soilmoisture)
library(tidyverse)
sw_soilmoisture <- sw_soilmoisture %>% filter(Treatment != "Pre-Drought") 
View(sw_soilmoisture)
#rename column to remove units
colnames(sw_soilmoisture)<-c("Harvest","Treatment","Type","EnvMass","EnvWetMass","EnvDryMass")

sw_soilmoisture$percentmoisture<-((sw_soilmoisture$EnvWetMass-sw_soilmoisture$EnvMass)-(sw_soilmoisture$EnvDryMass-sw_soilmoisture$EnvMass))/(sw_soilmoisture$EnvDryMass-sw_soilmoisture$EnvMass)*100 

bean_soilmoisture<-read.csv("Bean_soil.csv", sep=',', header=T, row.names=1)
View(bean_soilmoisture)
library(tidyverse)
bean_soilmoisture <- bean_soilmoisture %>% filter(Treatment != "Pre-Drought") 
View(bean_soilmoisture)
#rename column to remove units
colnames(bean_soilmoisture)<-c("Harvest","Treatment","Type","EnvMass","EnvWetMass","EnvDryMass")

bean_soilmoisture$percentmoisture<-((bean_soilmoisture$EnvWetMass-bean_soilmoisture$EnvMass)-(bean_soilmoisture$EnvDryMass-bean_soilmoisture$EnvMass))/(bean_soilmoisture$EnvDryMass-bean_soilmoisture$EnvMass)*100 


##3 way anova

library(ggplot2)
library(grid)
library(lattice)
install.packages("multcompView")
library(multcompView)
options(contrasts = c("contr.sum", "contr.poly"))

model = lm(percentmoisture~ Type + Treatment + Harvest + Type:Treatment + Type:Harvest + Treatment:Harvest + Type:Treatment:Harvest,
           data=bean_soilmoisture) ## replace with sw_soilmoisture when doing switchgrass analysis
shapiro.test(resid(model)) #W statistic > 0.9 means normality okay #check normality to flag outliers 
ee= as.matrix(resid(model))
qqnorm (ee) #look for outliers (far away from curve) and test for normality 
ee
install.packages("car")
library(car)
#options(contrasts = c("contr.sum", "contr.poly"))
Anova(model, type="III")


########Plotting
##################

sw_soilmoisture<-read.csv("Switchgrass_soil.csv", sep=',', header=T, row.names=1)
View(sw_soilmoisture)
library(tidyverse)
#sw_soilmoisture <- sw_soilmoisture %>% filter(Treatment != "Pre-Drought") 
View(sw_soilmoisture)
#rename column to remove units
colnames(sw_soilmoisture)<-c("Harvest","Treatment","Type","EnvMass","EnvWetMass","EnvDryMass")

sw_soilmoisture$percentmoisture<-((sw_soilmoisture$EnvWetMass-sw_soilmoisture$EnvMass)-(sw_soilmoisture$EnvDryMass-sw_soilmoisture$EnvMass))/(sw_soilmoisture$EnvDryMass-sw_soilmoisture$EnvMass)*100 

bean_soilmoisture<-read.csv("Bean_soil.csv", sep=',', header=T, row.names=1)
View(bean_soilmoisture)
library(tidyverse)
#bean_soilmoisture <- bean_soilmoisture %>% filter(Treatment != "Pre-Drought") 
View(bean_soilmoisture)
#rename column to remove units
colnames(bean_soilmoisture)<-c("Harvest","Treatment","Type","EnvMass","EnvWetMass","EnvDryMass")

bean_soilmoisture$percentmoisture<-((bean_soilmoisture$EnvWetMass-bean_soilmoisture$EnvMass)-(bean_soilmoisture$EnvDryMass-bean_soilmoisture$EnvMass))/(bean_soilmoisture$EnvDryMass-bean_soilmoisture$EnvMass)*100 



library(Rmisc)
sw_data2 <- summarySE(sw_soilmoisture, measurevar="percentmoisture", groupvars=c("Harvest","Treatment","Type"),na.rm=TRUE)

#reshape data so Pre-treatment comes first in legend
sw_data2$Treatment <- factor(sw_data2$Treatment,levels = c("pre-drought","drought","watered"), ordered=TRUE)

sw_data2$Treatment
View(sw_data2)


bean_data2 <- summarySE(bean_soilmoisture, measurevar="percentmoisture", groupvars=c("Harvest","Treatment","Type"),na.rm=TRUE)

#reshape data so Pre-treatment comes first in legend
bean_data2$Treatment <- factor(bean_data2$Treatment,levels = c("pre-drought","drought","watered"), ordered=TRUE)

bean_data2$Treatment
View(bean_data2)


##compiled graph 1

##add a type column to the shootmass data for bean and switchgrass

bean_data1$Type<- c("planted")
bean_data1$crop<- c("Bean")
View(bean_data1)

sw_data1$Type<- c("planted")
sw_data1$crop<- c("Switchgrass")
View(sw_data1)

sw_data2$crop<- c("Switchgrass")
bean_data2$crop<- c("Bean")

View(sw_data2)
View(bean_data2)


data_compiled_percentmoisture<- full_join(bean_data2, sw_data2)

View(data_compiled_percentmoisture)

data_compiled_percentmoisture$measure<- c("Gravimetric Soil Moisture (%)")


data_compiled_shootbiomass<- full_join(sw_data1, bean_data1)

View(data_compiled_shootbiomass)
data_compiled_shootbiomass$measure<- c("Shoot Biomass (g)")

data_final<- full_join(data_compiled_percentmoisture, data_compiled_shootbiomass)
View(data_final)
data_final_1 = data_final %>% mutate(value=ifelse(is.finite(percentmoisture), percentmoisture, ShootMass))

View(data_final_1)

write.csv(data_final_1, "data_final_compiled.csv")

limits <- aes(ymax = value+se, ymin =value-se)
dodge <- position_dodge(width = 0.9)

data_final_moisture<- data_final_1 %>% filter (measure=="Gravimetric Soil Moisture (%)")
View(data_final_moisture)

data_final_biomass<- data_final_1 %>% filter (measure=="Shoot Biomass (g)")
View(data_final_biomass)

p1<-ggplot(data_final_moisture,aes(Harvest,value,fill=Treatment)) + geom_errorbar(limits, position = dodge, width = 0.25)+
  geom_bar(aes(fill=Treatment),position="dodge",stat="identity")+ 
  theme_classic()+
  facet_grid( crop~ measure+Type, scales="free_y" ) +
  theme(panel.background = element_blank())+ scale_x_discrete(limits = c("Day 0","Day 2","Day 3","Day 4","Day 5", "Day 6"))+
  scale_fill_manual(values=c("#D55E00","#009E73","#56B4E9")) + theme(axis.title.x = element_text(size=15)) + theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) + ylab ("")+ xlab("") + theme(axis.text.y=element_text(size=12)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=12))
p1



p2<-ggplot(data_final_biomass,aes(Harvest,value,fill=Treatment)) + geom_errorbar(limits, position = dodge, width = 0.25)+
  geom_bar(aes(fill=Treatment),position="dodge",stat="identity")+ 
  theme_classic()+
  facet_grid( crop~ measure+Type, scales="free_y" ) +
  theme(panel.background = element_blank())+ scale_x_discrete(limits = c("Day 0","Day 2","Day 3","Day 4","Day 5", "Day 6"))+
  scale_fill_manual(values=c("#D55E00","#009E73","#56B4E9")) + theme(axis.title.x = element_text(size=15)) + theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) + ylab ("")+ xlab("") + theme(axis.text.y=element_text(size=12)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=12))
p2




library(ggpubr)

plot_combined<- ggarrange(p1, p2,  labels = c("A", "B"),
                          common.legend = TRUE, legend = "bottom")

plot_combined


########## Figure 2#######

ggsave(filename = "bean_sw_moisture_biomass_new_colorblindfriendly_final_newlabel.tiff", plot = plot_combined,
       width = 22,
       height = 13, units = c("cm"),
       dpi = 300)


#####################################
########Activity dynamics############
#####################################

##load indicator species analysis package

install.packages("indicspecies")
library(indicspecies)

###########

##phyloseq object
library(phyloseq)
library(tidyverse)
otu = read.csv("active_final_DNAabund_phantoms_threshold_detection5%.csv", sep=",", row.names=1)

tax = read.csv("taxonomy-dn-99-edited.csv", sep=",", row.names=1)
tax = as.matrix(tax)

metadata = read.csv("sample_metadata_dna_cdna_1.csv", sep=",", row.names=1) # use for crop indicators 
metadata<- metadata %>% filter(drought!= "pre-drought", crop=="switchgrass") #use for within crop indicators without pre-drought

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

##subsetting options, not needed initially to get indicators
metadata<- read.csv("bean-ww-nopredr.csv")##change metadata file as needed to analyze for plots

keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) 
phyloseq_merged_final

##indicator species analysis

library(indicspecies)

GetIndicators <-function(dataframe, var){
  require(phyloseq); require(indicspecies); require(dplyr)
  #dataframe <- tax_glom(dataframe, taxrank = "Family") ##can be adjusted to indicate at family or class level, else default is OTU level
  otu <- as.data.frame(otu_table(dataframe, taxa_are_rows = TRUE))
  metadata = as(sample_data(dataframe), "data.frame")
  taxa <- as.data.frame(as.matrix(tax_table(dataframe)))
  multipatt <- multipatt(t(otu), metadata[,var], func = "r.g",
                         control=how(nperm=999), duleg=FALSE)
  multipatt -> multipatt_fdr
  multipatt_fdr$sign$p.value <- p.adjust(multipatt_fdr$sign$p.value, "fdr")
  #print(multipatt_fdr)
  multipatt_fdr$sign[which(
    multipatt_fdr$sign$p.value <= 0.05), ] -> indicator_taxa
  taxa$OTU <- rownames(taxa)
  data.frame(OTU = as.factor(row.names(indicator_taxa)), indicator_taxa) %>%
    dplyr::left_join(taxa, by="OTU") -> indicator_taxa
  rownames(indicator_taxa) <- indicator_taxa$OTU
  indicator_taxa <- arrange(indicator_taxa, desc(stat))
  return(indicator_taxa)
}


#GetIndicators(datafraem=phyloseq_object, var="Soil_location")


ind_bean <-
  GetIndicators(dataframe=phyloseq_merged, "summarise_new")

ind_sw <-
  GetIndicators(dataframe=phyloseq_merged, "summarise_new")

ind_crop <-
  GetIndicators(dataframe=phyloseq_merged, "crop")

ind_crop <-
  GetIndicators(dataframe=phyloseq_merged, "summarise_new1") ##using metadata file sample_metadata_dna_cdna_1.csv, this looks at all conditions across crops

##checking

View(ind_bean)
View(ind_sw)
View(ind_crop)
ind_bean %>% subset(x = ., s.drought_planted==1 & s.drought_unplanted==1)
ind_crop %>% subset(s.bean==1 & s.switchgrass==1) ##no common indicators across two plants when including ALL 212 samples in analysis, even at family level

##checking indicators across similar conditions between two crops

ind_crop %>% subset(s.bean.drought_planted==1 & s.switchgrass_drought_planted==1)
ind_crop %>% subset(s.bean.drought_unplanted==1 & s.switchgrass_drought_unplanted==1)
ind_crop %>% subset(s.bean.watered_planted==1 & s.switchgrass_watered_planted==1)
ind_crop %>% subset(s.bean.watered_unplanted==1 & s.switchgrass_watered_unplanted==1)


library(tidyverse)

write.csv(ind_bean, "bean_indicators.csv")

write.csv(ind_sw, "sw_indicators.csv")

write.csv(ind_crop, "crop_indicators.csv")


##plot indicator species
##first subset to taxa that are only associated with drought planted/drought unplanted/ww planted/ ww uplanted conditions in bean and sw, total 8 datasets. by subsetting s.planted.drought ==1 option etc (condition specific taxa) , ARE ANY TAXA SHARED ACROSS CONDITIONS BETWEEN CROPS? JOIN DATASETS AND SEE OVERLAP (EXAMPLE : COMPARE SW DROUGHT/PLANTED WITH BEAN DR/PLANTED?)
##similarly pick the taxa associated only with bean, bean ==1 from crop dataset, repeat for switchgrass (crop specific taxa), NOTE: no taxa are shared 
##merge or left join with active DNA abundance table by rownames and keep columns from abundance dataset
##pick top 20-30 OTUs by sorting in decreasing abundances
##max standardize abundances across rows
##plot bubble plot (using psmelt to long format - phyloseq or pivot longer - tidyverse)
##or plot heatmap


library(tidyverse)

#read in bean data from crop indicator table
crop_ind<- read.csv("bean_indicators_editedCol.csv") ##changed crop_indicators to bean_indicators and then selected planted drought, unplanted drought etc.
View(crop_ind)
crop_ind_sw <- crop_ind %>% filter(drought_planted==0 & drought_unplanted==0 & well_watered_planted==1 & well_watered_unplanted==0) ##when looking at unique associations with bean or switchgrass the condition should be switchgrass==1, && bean==0, and vice versa. select drought, planted conditions when looking at those conditions within plant, such as within bean
View(crop_ind_sw)

Active_DNA<- read.csv("active_final_DNAabund_phantoms_threshold_detection5%_edited.csv", row.names=1)
View(Active_DNA)
Active_new <- cbind(OTU = rownames(Active_DNA), Active_DNA)

ind_sw_abund<- left_join(crop_ind_sw, Active_new, by="OTU")
View(ind_sw_abund)
rownames(ind_sw_abund) <- ind_sw_abund[,1]
ind_sw_mat = ind_sw_abund[,16:ncol(ind_sw_abund)] ##chnage the column number as appropriate
View(ind_sw_mat)
#ind_sw_50 <- ind_sw_mat[order(rowSums(ind_sw_mat), decreasing = TRUE),]
#View(ind_sw_50)

###moved this to later so that decreasing order is set after subsetting to switchgrass/bean samples
#ind_bean_50 <- ind_bean_mat[order(rowSums(ind_bean_mat), decreasing = TRUE),]
#ind_bean_50 <- ind_bean_50[c(1:50),] ##select top 50 most abundant Classes
#View(ind_bean_50)

#ind_bean_20 <- ind_bean_mat[order(rowSums(ind_bean_mat), decreasing = TRUE),]
#ind_bean_20 <- ind_bean_20[c(1:20),] ##select top 50 most abundant Classes
#View(ind_bean_20)

####start again from here
##subset to correct samples before computing max standardization
ind_sw_long<- ind_sw_mat %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="Abundance", -ASV)#gather columns into key value pairs
View(ind_sw_long)

metadata<-read.csv("metadata_dna_cdna_new.csv")

newdf = left_join(ind_sw_long, metadata, by=c("SampleID"))
View(newdf)
new_df_sw<- newdf %>% filter(crop=="bean", drought=="well-watered", planted=="planted") ##change to new_df_bean or keep using the same object name, but change the parameters as needed for the crop and within crop factors, should match the condition in line 3679 

##pivot wider to plot averages, do the rowsums sorting after filtering to switchgrass? if i do it before then bean abundances may obscure the sw rowsum abudances (this is unlikely given that the OTUs picked are indicators for a specific crop and likely will have more abundances in that crop.)
new_df_sw_wide <- new_df_sw %>% select(ASV, Abundance, SampleID) %>%pivot_wider(names_from="SampleID", values_from="Abundance")  
View(new_df_sw_wide) 
rnames <-new_df_sw_wide$ASV

new_df_sw_wide<- new_df_sw_wide %>% select (!ASV)
rownames(new_df_sw_wide) <- rnames
View(new_df_sw_wide)

new_df_sw_wide<- as.matrix(new_df_sw_wide)
new_df_sw_wide <- new_df_sw_wide[order(rowSums(new_df_sw_wide), decreasing = TRUE),]
new_df_sw_wide <- new_df_sw_wide[c(1:50),] ##select top 50 most abundant Classes
View(new_df_sw_wide)



##merging with dna cdna data to get NAs inserted in dataframe
##using rarefied DNA and cDNA copies with 15k reads.
dna<- read.csv("DNAcopy.csv", header=TRUE, row.names = 1)
cdna <- read.csv("cDNAcopy.csv", header=TRUE, row.names = 1)

View(dna)
View(cdna)

dna_20taxa <- dna[rownames(dna)%in%rownames(new_df_sw_wide),]
View(dna_20taxa)

cdna_20taxa <- cdna[rownames(cdna)%in%rownames(new_df_sw_wide),]
View(cdna_20taxa)

dna1<- dna_20taxa %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="dna", -ASV)#gather columns into key value pairs
View(dna1)
cdna1<- cdna_20taxa %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="cdna", -ASV)
View(cdna1)
merge.df = full_join(dna1, cdna1, by=c("ASV", "SampleID"))
View(merge.df)

merge.df[is.na(merge.df)] <-0
View(merge.df)
metadata<-read.csv("metadata_dna_cdna_new.csv")

newdf = full_join(merge.df, metadata, by=c("SampleID"))
View(newdf)


finaldf = newdf %>%  mutate(ratioMethod1= ifelse(dna == 0 & cdna > 0, 100, cdna/dna)) %>%  
  mutate(dna2= ifelse(dna == 0, 1,dna))  %>%  mutate(ratioMethod2= cdna/dna2) %>% mutate(ratio_nophantoms = cdna/dna)
View(finaldf)
finaldf[is.na(finaldf)] <-0
finaldf$ratioMethod1[!is.finite(finaldf$ratioMethod1)] <- 0 ##using the code above for removing nas also works here as it changed the 0 denominator ratios from NA to 0, no infinite values only NaN
finaldf$ratioMethod2[!is.finite(finaldf$ratioMethod2)] <- 0 ##ratios computed by methods 1 and 2 will always be finite, because of the zero denominator being accounted for phantom taxa?
finaldf$ratio_nophantoms[!is.finite(finaldf$ratio_nophantoms)] <- 0
View(finaldf) 


##code taxa as inactive (0) or not detected in dna (NA)


finaldf_coded <- finaldf %>% #mutate(dna_code=ifelse(dna2==0, "not_detected", dna2)) 
  mutate(activity= ifelse(ratioMethod2 < 1, "inactive", "active"))
View(finaldf_coded) 

##Abundant_active dataset abundance combined with finaldf dataset 
Abundant_Active_new <- cbind(ASV = rownames(new_df_sw_wide), new_df_sw_wide)
View(Abundant_Active_new)
Abundant_Active_new<- as.data.frame(Abundant_Active_new)
Abundant_Active_long <-pivot_longer(Abundant_Active_new, !ASV, names_to="SampleID", values_to="AbundanceDNA_Active")
View(Abundant_Active_long)
merge_df_new<- merge(Abundant_Active_long, finaldf_coded, by=c("ASV", "SampleID"))
View(merge_df_new)
merge_df_new_NA<- merge_df_new %>% mutate(AbundanceDNA_Active_NA= ifelse(dna==0 & cdna==0, NA, AbundanceDNA_Active))
View(merge_df_new_NA)

##pivot wider

merge_df_new_NA_wide <- merge_df_new_NA %>% select(ASV, AbundanceDNA_Active_NA, treatment) %>%pivot_wider(names_from="treatment", values_from="AbundanceDNA_Active_NA")  
View(merge_df_new_NA_wide) 
data<- as.data.frame(merge_df_new_NA_wide)

##create matrix


#option 2, both options work
rnames <-data[,1] 
rownames(data) <- rnames  # assign row names
View(data)
mat_data<- data %>% select(!ASV)
View(mat_data)


#REMOVE OTUS that cant be plotted in heatmap, note which ones are removed, if there are only numbers nothing is removed. if there are NAs some may be removed.

#mat_data<- data.matrix(mat_data)


sum(is.na(as.matrix(dist(mat_data))))

giveNAs = which(is.na(as.matrix(dist(mat_data))),arr.ind=TRUE)
head(giveNAs)
mat_data[c(1,17),]

tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
  sum(is.na(as.matrix(dist(mat_data[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))] ##https://stackoverflow.com/questions/61469201/pheatmap-won-t-cluster-rows-na-nan-inf-in-foreign-function-call-arg-10 
rmv ## the intention with this code is not to remove rows but to see for which ASVs are we not able to calculate a distance for the heatmap
mat_data = mat_data[-as.numeric(rmv),] ##use this line only is something needs to be removed
View(mat_data)


####merge taxonomy with OTU hash id

taxonomy<- read.csv("taxonomy-dn-99-edited.csv", row.names=1)
View(taxonomy)
taxonomy_family<- taxonomy %>% select(c("Family", "Class"))
View(taxonomy_family)

mat_data_merge <- merge(mat_data, taxonomy_family, by="row.names")
View(mat_data_merge)
mat_data_merge$OTUclass <- paste(mat_data_merge$Class, mat_data_merge$Row.names)
View(mat_data_merge)

rnames <-mat_data_merge$OTUclass

rownames(mat_data_merge) <- rnames  # assign row names
View(mat_data_merge)
mat_data_merge_otuclass<- mat_data_merge %>% select(!Row.names)%>% select (!Family) %>% select(!Class) %>% select(!OTUclass)
View(mat_data_merge_otuclass)

write.csv(mat_data_merge_otuclass, "mat_data_merge_otuclass_bean_dr_planted.csv")
write.csv(mat_data_merge_otuclass, "mat_data_merge_otuclass_bean_dr_unplanted_ind.csv")

write.csv(mat_data_merge_otuclass, "mat_data_merge_otuclass_bean_ww_planted_ind.csv")
write.csv(mat_data_merge_otuclass, "mat_data_merge_otuclass_bean_ww_unplanted_ind.csv")

if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
install_github("raivokolde/pheatmap")
library(pheatmap)

library(vegan)

mat_data_merge_otuclass<- read.csv("mat_data_merge_otuclass_bean_dr_planted_ind.csv", row.names=1)

mat_data_merge_otuclass<- as.data.frame(mat_data_merge_otuclass)
#mat_data_merge_otuclass<- as_tibble(mat_data_merge_otuclass)

class(mat_data_merge_otuclass)
##replace NAs with zero
matrix<- as.matrix(mat_data_merge_otuclass)
matrix
View(matrix)
matrix[is.na(matrix)] <-0
View(matrix)

#matrix_new<-data.matrix(matrix)
#View(matrix_new)
##compute max standardization

matrix_max_standardize_sw <-decostand(matrix, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(matrix_max_standardize_sw)

##put back NAs where they were before
matrix1<- as.matrix(mat_data_merge_otuclass)
matrix1
View(matrix1)
matrix1[is.na(matrix1)] <- -100
View(matrix1)
mask <- matrix1 == -100
mask

matrix_max_standardize_sw[mask] <- matrix1[mask]
View(matrix_max_standardize_sw)

matrix_max_standardize_sw[matrix_max_standardize_sw == -100]<- NA
View(matrix_max_standardize_sw) ##gives the standardized values with NA


###

library(viridis)

##pivot longer, merge with metadata
##convert to dataframe if doing this part

#sw_long<- matrix_max_standardize_sw %>%
  #tibble::rownames_to_column(var="ASV") %>%
 # tidyr::gather(key="SampleID", value="RA", -ASV)#gather columns into key value pairs
#View(sw_long)

#metadata<-read.csv("sample_metadata_dna_cdna_edited.csv")

#new_df_sw = left_join(sw_long, metadata, by=c("SampleID"))
#View(new_df_sw)

#compute averages by harvest and summarise new, might not work if there are NAs, might need the zeros . In that case dont put back the NAs after max standardization
#newdf_mean_sw<- new_df_sw %>% group_by(summarise, ASV) %>% summarise(mean=mean(RA))
#View(newdf_mean_sw)

##pivot wider to plot averages
#sw_ind_wide <- newdf_mean_sw %>% select(ASV, mean, summarise) %>%pivot_wider(names_from="summarise", values_from="mean")  
#View(sw_ind_wide) 

#class(sw_ind_wide)

#sw_ind_mat<- as.matrix(sw_ind_wide)
#class(sw_ind_mat)
#View(sw_ind_mat)

#rnames <-sw_ind_wide$ASV

#sw_ind_wide<- sw_ind_wide %>% select (!ASV)
#rownames(sw_ind_wide) <- rnames
#View(sw_ind_wide)


#plot

##bean planted drought, using default clustering "complete" linkage

#create the breaks
bk2 = unique(c(0.00000000, seq(0.00000001, 1.00000000, length=20)))

col1 = colorRampPalette(c("gray"))(1)

col2 =viridis(20)

colors2 <- c(col1, col2)



plot_heatmap_bean_planted_drought <- pheatmap::pheatmap(matrix_max_standardize_sw, ##if plotting sample by sample and not the means then use matrix_max_standardize_sw
                                        cluster_cols = F, border_color = "white", na_col = "white", clustering_distance_rows="euclidean",clustering_method = "complete",
                                        #main="Bean indicator species (Unique to planted+drought condition)",
                                        color = colors2,
                                        breaks = bk2
                                        )
                                                
                                        
library(pvclust)

mat_t<- t(matrix_max_standardize_sw)
View(mat_t)
result_bean_drought <- pvclust(mat_t, method.dist="euclidian", method.hclust="complete", nboot=10000, parallel=TRUE)
result_bean_drought


plot(result_bean_drought) #NOT USED IN FINAL MANUSCRIPT
pvrect(result_bean_drought, alpha=0.95)
seplot(result_bean_drought, identify=TRUE) #Values on the edges of the clustering are p-values (%). Red values are AU p-values, and green values are BP values. Clusters with AU larger than 95% are highlighted by rectangles, which are strongly supported by data.

##trying with Ward.D

mat_t<- t(matrix_max_standardize_sw)
View(mat_t)
result_bean_drought <- pvclust(mat_t, method.dist="euclidian", method.hclust="ward.D", nboot=10000, parallel=TRUE)
result_bean_drought

plot(result_bean_drought)
pvrect(result_bean_drought, alpha=0.95)
seplot(result_bean_drought, identify=TRUE) ##ward gives one au.se error value of 0.5 which is quite large and it has lower au. p value. clustering remains same for Ward and complete clustering approaches so sticking to complete

##bean dr unplanted (NOT USED IN FINAL MANUSCRIPT)

#create the breaks
bk2 = unique(c(0.00000000, seq(0.00000001, 1.00000000, length=20)))

#set different color vectors for each interval
col1 = colorRampPalette(c("gray"))(1)

col2 =viridis(20)

colors2 <- c(col1, col2)

plot_heatmap_bean_unplanted_drought <- pheatmap::pheatmap(matrix_max_standardize_sw, ##if plotting sample by sample and not the means then use matrix_max_standardize_sw
                                                        cluster_cols = F, border_color = "white", na_col = "white", clustering_distance_rows="euclidean",clustering_method = "complete",
                                                        #main="Bean indicator species (Unique to planted+drought condition)",
                                                        color = colors2,
                                                        breaks = bk2
)



##Figure 5B

ggsave(filename = "bean_ind_planted_drought_newXlabel_ind1.tiff", plot = plot_heatmap_bean_planted_drought,
       width = 30,
       height = 25, units = c("cm"),
       dpi = 300)


#ggsave(filename = "bean_ind_unplanted_drought_newXlabel_ind1.tiff", plot = plot_heatmap_bean_unplanted_drought,
      # width = 30,
      # height = 25, units = c("cm"),
      # dpi = 300)

#ggsave(filename = "bean_ind_planted_ww_newXlabe_ind.tiff", plot = plot_heatmap_bean_planted_ww,
      # width = 30,
      # height = 25, units = c("cm"),
       #dpi = 300)


#ggsave(filename = "bean_ind_unplanted_ww_newXlabel_ind.tiff", plot = plot_heatmap_bean_unplanted_ww,
       #width = 30,
       #height = 25, units = c("cm"),
       #dpi = 300)
## Note: for the heatmaps within crop (for each condition) we only plotted for bean. because for sw there were no indicator species for any condition.





##Venn diagrams for overlap and exclusion of OTUs across different conditions in bean and switchgrass
##venn diagran should have 4 conditions or circles (planted_dr, planted_ww, unplanted_dr, unplanted_ww)

if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")


install.packages("VennDiagram")
library(ggVennDiagram)

library(VennDiagram)

bean_ind<- read.csv("bean_indicators_editedCol.csv", header = TRUE)
View(bean_ind)
class(bean_ind)
library(tibble)
bean_ind <- as_tibble(bean_ind)
library(tidyverse)


drought_planted<- bean_ind %>% filter(drought_planted==1)
drought_unplanted= bean_ind %>% filter(drought_unplanted==1)
watered_planted= bean_ind %>% filter(well_watered_planted==1)
watered_unplanted=bean_ind %>% filter(well_watered_unplanted==1)



#Define sets for diagram
SET1 <- drought_planted$OTU
SET2 <- drought_unplanted$OTU
SET3 <- watered_planted$OTU
SET4 <- watered_unplanted$OTU


#Draw the diagram 

##Figure 5A (REMADE IN POWERPOINT)
venn.diagram(list(drought_planted=SET1, drought_unplanted=SET2, watered_planted=SET3, watered_unplanted=SET4),main.cex=4, 
             sub.cex=3,
                   fill = c("red", "green", "blue", "yellow"),
                   alpha = c(0.5, 0.5, 0.5, 0.5),
                   filename='venn_bean_ind.tiff', height = 6000 , 
                   width = 6000 , resolution = 700, units="px", 
                   imagetype="tiff", compression = "lzw", output=TRUE)


###checking the OTUs for each set for abundance plotting. the only plot most important here is the drought planted set with unique OTUs not shared in any other set.

dr_planted_1= bean_ind %>% filter(drought_planted==1 & drought_unplanted==0 & well_watered_planted==0 & well_watered_unplanted==0)
View(dr_planted_1)
dr_unplanted_1 = bean_ind %>% filter(drought_planted==0  & drought_unplanted==1 & well_watered_planted==0 & well_watered_unplanted==0)
ww_planted_1 = bean_ind %>% filter(drought_planted==0  & drought_unplanted==0 & well_watered_planted==1 & well_watered_unplanted==0)
ww_unplanted_1 =bean_ind %>% filter(drought_planted==0  & drought_unplanted==0 & well_watered_planted==0 & well_watered_unplanted==1)



#########################################################
################Alpha Diversity #########################
#########################################################


##richness (observed OTUs) estimates, inverse simpson

library(ggplot2)
library(dplyr)
library(vegan)
library(reshape2)
library(scales)
library(grid)
library(readxl)


# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(phyloseq_merged)
trials = 100
richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(phyloseq_merged)
evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(phyloseq_merged)
# It is always important to set a seed when you subsample so your result is replicable, note of caution that this will yield different results of alpha diversity depending on the version of R being used. This is because set.seed function can vary across R versions. Here I am reporting results from R v.3.4.0 
set.seed(3)
for (i in 1:100) {
  # Subsample
  #r <- rarefy_even_depth(phyloseq_merged, sample.size = 1936, verbose = FALSE, replace = TRUE) #no rarefaction a second time.
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(phyloseq_merged, measures = "Observed"))) ##changed r to phyloseq_merged, no second rarefaction
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(phyloseq_merged, measures = "InvSimpson"))) ##changed r to phyloseq_merged, no second rarefaction
  evenness[ ,i] <- even
}
# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)
# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd) 
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)
alpha <- rbind(rich_stats, even_stats)

s <- read.csv("sample_metadata_dna_cdna.csv")
alphadiv <- merge(alpha, s, by = "SampleID") 
write.csv(alphadiv,file = "alphadiv_withsingleton_nosecondrarefaction.csv")

alphadiv_sw<-alphadiv[alphadiv$crop=="switchgrass",]
dim(alphadiv_sw)
alphadiv_bean<-alphadiv[alphadiv$crop=="bean",]
dim(alphadiv_bean)
View(alphadiv_sw)
View(alphadiv_bean)
write.csv(alphadiv_bean, file="alphadiv_bean_withsingleton_nosecondrarefaction.csv")
write.csv(alphadiv_sw, file="alphadiv_switchgrass_withsingleton_nosecondrarefaction.csv")

group_sw<-alphadiv_sw%>%
  group_by(summarise,drought,measure) %>%
  summarise(mean_update = mean(mean))
write.csv(group_sw, file="alphadivmean_sw_withsingleton_nosecondrarefaction.csv")
group_bean<-alphadiv_bean%>%
  group_by(summarise, drought,measure) %>%
  summarise(mean_update = mean(mean))
write.csv(group_bean, file="alphadivmean_bean_withsingleton_nosecondrarefaction.csv")

### Inverse Simpson bean
alphadiv_bean <- read.csv("alphadiv_bean_withsingleton_nosecondrarefaction.csv")
pd<-position_dodge(0.7)
alphadiv_bean_InvS =  alphadiv_bean %>% filter(measure=="Inverse Simpson")
alphadiv_bean_InvS$drought <- factor(alphadiv_bean_InvS$drought,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)
View(alphadiv_bean_InvS)


####Inverse Simpson Switchgrass

alphadiv_sw <- read.csv("alphadiv_switchgrass_withsingleton_nosecondrarefaction.csv")
pd<-position_dodge(0.7)
alphadiv_sw_InvS =  alphadiv_sw %>% filter(measure=="Inverse Simpson")
alphadiv_sw_InvS$drought <- factor(alphadiv_sw_InvS$drought,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)
View(alphadiv_sw_InvS)


####Bean Richness
alphadiv_bean <- read.csv("alphadiv_bean_withsingleton_nosecondrarefaction.csv")
pd<-position_dodge(0.7)
alphadiv_bean_rich =  alphadiv_bean %>% filter(measure=="Richness")
alphadiv_bean_rich$drought <- factor(alphadiv_bean_rich$drought,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)
View(alphadiv_bean_rich)



###Switchgrass richness
alphadiv_sw <- read.csv("alphadiv_switchgrass_withsingleton_nosecondrarefaction.csv")
pd<-position_dodge(0.7)
alphadiv_sw_rich =  alphadiv_sw %>% filter(measure=="Richness")
alphadiv_sw_rich$drought <- factor(alphadiv_sw_rich$drought,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)
View(alphadiv_sw_rich)


####
####Alpha diversity ANOVA statistics
####


library(ggplot2)
library(grid)
library(lattice)
library(multcompView)
library(tidyverse)
library(dplyr)
alphadiv_rich_old <- read.csv("bean_sw_rich_even_all.csv") ## combined richness and evenness estimates of all samples
alphadiv_rich <- alphadiv_rich_old %>% filter(drought!="pre-drought") ## if pre-drought is not removed then error in anova saying "Error in Anova.III.lm(mod, error, singular.ok = singular.ok, ...) : "
#"there are aliased coefficients in the model", possibly unbalanced data
pd<-position_dodge(0.7)

alphadiv_rich =  alphadiv_rich %>% filter(measure=="Inverse Simpson" & crop=="bean")
alphadiv_rich =  alphadiv_rich %>% filter(measure=="Inverse Simpson" & crop=="switchgrass")
alphadiv_rich =  alphadiv_rich %>% filter(measure=="Richness" & crop=="bean")
alphadiv_rich =  alphadiv_rich %>% filter(measure=="Richness" & crop=="switchgrass")

##Conduct the following code for each of the datasets selected above
View(alphadiv_rich)
options(contrasts = c("contr.sum", "contr.poly")) ### must set this before running the model for type III tests
model = lm(mean~ planted+drought + harvest 
           + planted:drought + planted:harvest + drought:harvest + planted:drought:harvest,
           data=alphadiv_rich)

shapiro.test(resid(model)) #W statistic > 0.9 means normality okay #check normality to flag outliers 
ee= as.matrix(resid(model))
qqnorm (ee) #look for outliers (far away from curve) and test for normality 
ee

library(car)
#options(contrasts = c("contr.sum", "contr.poly"))
Anova(model, type=c("III"))




#########################
####Metabolomics analysis 
#########################

##For PCA plots 
library(MetabolAnalyze)
PCA<- read.csv("PCA and PLSDA_all features_POS_transpose.csv", row.names=1)
View(PCA)
class(PCA)
pca_mat<- as.matrix(PCA)
View(pca_mat)
pca_mat_one<- pca_mat+1
View(pca_mat_one)
log_pca<- log10(pca_mat_one) ##log normalize
View(log_pca)
pca_scale<- scaling(log_pca, type = "pareto") ##pareto scale
View(pca_scale)
install.packages("NormalizeMets")
#ppca.metabol(pca_scale, minq=1, maxq=2, scale = "none", epsilon = 0.1, 
#plot.BIC = FALSE, printout=TRUE)
pca<-prcomp(pca_scale, center=F, scale=F)
pcaresults<- summary(pca)

summary(pca)
pcaresults
scree.data<- as.data.frame(pcaresults$importance)
score.data<- as.data.frame(pcaresults$x)
scree.data
score.data
loadings.data<- as.data.frame(pcaresults$rotation)
loadings.data
write.csv(scree.data, "pca_scree.csv")
write.csv(score.data, "pca_scores_prcomp.csv")
write.csv(loadings.data, "pca_loadings_prcomp.csv")
data<- read.csv("pca_scores_prcomp.csv", header=T)
View(data)
data<- data[, c(1:4)]
metadata<- read.csv("metadata_PCA_PLSDA.csv", row.names=1)
View(metadata)
planted<- as.vector(metadata$planted)
drought<-as.vector(metadata$drought)
harvest<- as.vector(metadata$harvest)
data_metadata<- cbind (data, planted, drought, harvest)
View(data_metadata)

##ggplot
library(viridis)

plot2<-ggplot(data_metadata, aes(PC1, PC3, shape=planted, color=drought))+  
  #scale_fill_manual(name="drought", limits = c("pre-drought", "drought", "well-watered"),
  #values = c("#a65628", "red", "#ffae19",
  #"#4daf4a", "#1919ff", "darkorchid3", "magenta", "black" ))+
  scale_color_manual(name="drought", limits = c("pre-drought", "drought", "watered"), values = c("#D55E00","#009E73","#56B4E9")) +
  geom_point(aes(color = drought, size=harvest)) +
  scale_shape_manual(values=c(16,17))+
  scale_size_manual(values = c(2,3,4,5,6,7))+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Scores Plot")+
  ylab("PC3 0.4%")+ ##add the variance from summary(pca) function
  xlab("PC1 96.5%")+ ##add the variance from summary(pca) function
  stat_ellipse(level=0.65, type="t", linetype=2)
plot2

##Figure S4B
ggsave(filename="pca-metabolite.TIFF", plot=plot2, width=6.8, height=6, units="in", dpi=300)




##calculating eucidean distances to supplement PCA with PERMANOVA

metadata<- read.csv("metadata_PCA_PLSDA.csv")##change as needed between with and without pre drought samples
#metadata_planted<- metadata%>% filter(planted=="planted")
distm<- dist(pca_scale, method="euclidean")

##permanova using euclidean distance (above)

library(vegan)
adonis <- adonis2(distm ~metadata[,colnames(metadata)=='planted'])
adonis



adonis <- adonis2(distm ~metadata[,colnames(metadata)=='drought']) 
adonis
##differences between pre-dr drought and watered samples


adonis <- adonis2(distm ~metadata[,colnames(metadata)=='harvest'])
adonis




##checking differences between drought and watered samples within planted soil samples

library(MetabolAnalyze)
PCA<- read.csv("PCA and PLSDA_all features_POS_transpose_nopredr_onlyplanted.csv", row.names=1)
View(PCA)
class(PCA)
pca_mat<- as.matrix(PCA)
View(pca_mat)
pca_mat_one<- pca_mat+1
View(pca_mat_one)
log_pca<- log10(pca_mat_one)
View(log_pca)
pca_scale<- scaling(log_pca, type = "pareto")
View(pca_scale)


metadata<- read.csv("metadata_PCA_PLSDA_nopredr_onlyplanted.csv")##change as needed between with and without pre drought samples

distm<- dist(pca_scale, method="euclidean")
adonis <- adonis2(distm ~metadata[,colnames(metadata)=='drought']) 
adonis

## without pre-drought and unplanted samples, checking differences between drought and watered planted samples --not significant





##HEATMAP 
##USING THE TOP 50 LOADINGS FROM XINGXING 


heatmap<- read.csv("Heatmap_TOP50PLSDAcoef_POS_version2.csv", row.names=1)##ADD LOG TRANSFORM
View(heatmap)
heat_mat<- as.matrix(heatmap)
View(heat_mat)
heat_mat_one<- heat_mat+1
View(heat_mat_one)
log_heat<- log10(heat_mat_one)
View(log_heat)



# We will scale the data by each row from 0 (lowest) to 1 (highest)

##install packages
library("circlize")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library("ComplexHeatmap")#lib.loc="~/R/win-library/3.5"
metadata<- read.csv("metadata_PCA_PLSDA.csv")

##use log_heat data
for(i in 1:nrow(log_heat)){
  min_value = min(log_heat[i,])
  max_value = max(log_heat[i,])
  log_heat[i,] <- (log_heat[i,]-min_value)/(max_value-min_value)
}

log_heat_scale <- as.matrix(log_heat)
View(log_heat_scale)

fontsize<-1.1

ha_column = HeatmapAnnotation(Treatment = metadata$Label,
                              col = list(Treatment = c("pre-drought_planted" =  "coral", "pre-drought_unplanted" = "gold", 
                                                       "drought_unplanted"="cyan","watered_unplanted"="darkblue",
                                                       "drought_planted"="green","watered_planted"="darkgreen"), labels=c("pre-drought_planted", "pre-drought_unplanted",
                                                                                                                          "drought_planted", "drought_unplanted", "watered_planted", "watered_unplanted")),
                              annotation_legend_param = list(title = " ",title_gp = gpar(fontsize = 14),
                                                             labels_gp = gpar(fontsize = 14)),show_legend = c(TRUE))

map_colors<-colorRampPalette(viridis(12))
ht=Heatmap(log_heat_scale,show_row_names = T,show_column_names = F,cluster_columns = TRUE,
           col = colorRamp2(breaks=c(0,0.2,0.4,0.6,0.8,1), hcl_palette = "viridis") , #row_dend_width = unit(6,"cm"),
           clustering_method_rows = "ward.D",
           row_names_gp = gpar(fontsize = 8),
           clustering_distance_rows ="euclidean",row_dend_width=unit(2,"cm"),
           heatmap_legend_param = list(at = c(1, 0), 
                                       labels = c("high", "low"),
                                       title = " ", title_gp = gpar(fontsize = 6),
                                       labels_gp = gpar(fontsize = 8),
                                       legend_direction = "horizontal"),
           top_annotation = ha_column,rect_gp = gpar(col = "black", lty = 1, lwd = 1),km=3)

##cluster cols by ward (have not done this yet)
##ht=Heatmap(log_heat_scale,show_row_names = T,show_column_names = F,cluster_columns = TRUE,
           ##clustering_method_columns = "ward.D",
           ##clustering_distance_columns = "euclidean",
          ## col = colorRamp2(breaks=c(0,0.2,0.4,0.6,0.8,1), hcl_palette = "viridis") , #row_dend_width = unit(6,"cm"),
           ##clustering_method_rows = "ward.D",
           ##row_names_gp = gpar(fontsize = 8),
           ##clustering_distance_rows ="euclidean",row_dend_width=unit(2,"cm"),
          ## heatmap_legend_param = list(at = c(1, 0), 
                                       ##labels = c("high", "low"),
                                       ##title = " ", title_gp = gpar(fontsize = 6),
                                       ##labels_gp = gpar(fontsize = 8),
                                      ## legend_direction = "horizontal"),
           ##top_annotation = ha_column,rect_gp = gpar(col = "black", lty = 1, lwd = 1),km=3)

##Figure S4D
draw(ht, heatmap_legend_side = "bottom") ##export heatmap as pdf

library(pvclust)
##pvclust to check confidece in cluster of features
mat_t<- t(log_heat_scale)
View(mat_t)
result_metab <- pvclust(mat_t, method.dist="euclidian", method.hclust="ward.D", nboot=10000, parallel=TRUE)
result_metab

plot(result_metab) ##Figure NOT USED IN FINAL MANUSCRIPT
pvrect(result_metab, alpha=0.80) ##this cluster does not exactly match the heatmap but seems like cluster 1 has 85% confidence and the rest clustered together, cluster 1 has two or three features placed differently in the heatmap, cluster 3 one part is at 100% other at 76%, cluster 2 is at 94% cluster so pretty good
seplot(result_metab, identify=TRUE) 



###pvclust to check confidence in cluster of samples

result_metab_sample <- pvclust(log_heat_scale, method.dist="euclidian", method.hclust="complete", nboot=10000, parallel=TRUE)
result_metab_sample

plot(result_metab_sample) #Figure NOT USED IN FINAL MANUSCRIPT
pvrect(result_metab_sample, alpha=0.95) ##planted samples cluster is at 95% confidence with method  "complete", also clustering of all predrought samples are at 95% confidence, only one drought sample in that cluster as shown in the heatmap
seplot(result_metab_sample, identify=TRUE) 



##3D PCA
                    
##SCATTERPLOT3D
install.packages("scatterplot3d") # Install
library("scatterplot3d") # load
                    
scatterplot3d(data_metadata[,2:4], angle=55)
                    
write.csv(data_metadata, "data_metadata.csv") #sorted the drought column alphabetically but did not help with anything later so nothing was changed here basically 
                    
data_metadata<- read.csv("data_metadata.csv")
                    
colors <- c("#009E73","#D55E00","#56B4E9")
colors <- colors[as.numeric(as.factor(data_metadata$drought)) ]
colors
                    
shapes = c(20,19)
shapes <- shapes[as.numeric(as.factor(data_metadata$harvest))]
shapes
                    
shapes1 = c(16,17)
shapes1 <- shapes1[as.numeric(as.factor(data_metadata$planted))]
shapes1
                    
data_metadata$drought <- factor(data_metadata$drought, levels = c("pre-drought", "drought", "watered"))
                    
data_metadata$harvest <- factor(data_metadata$harvest, levels = c("Day 0", "Day 6"))
                    
data_metadata$planted <- factor(data_metadata$planted, levels = c("planted", "unplanted"))


##Figure S4C, save as PDF                  
plot<- scatterplot3d(data_metadata[,2:4], pch = shapes1, color=colors, cex.symbols=2)
#legend("bottomright", legend = levels(as.factor(data_metadata$harvest)),
# pch = c(20,19), cex=1)
legend("topleft", legend = levels(as.factor(data_metadata$planted)),
pch =  c(16,17),inset=c(0.001, 0), cex=1, box.lty = 0, title="planted")
legend("topleft", legend = levels(as.factor(data_metadata$drought)),
col =c("#D55E00","#009E73","#56B4E9"), pch=16, inset=c(0.002, 0.15), cex=1, box.lty = 0, title="drought")
                    

########################                    
##qubit dna concentration 
########################


##add geom_smooth for stats, fit a linear model and trend over time- V3 edit pending

data<- read.csv("qubit_plot.csv")
data<- as.data.frame(data)
class(data)

library(tidyverse)
View(data)
data
library(ggplot2)


plot_bean_sw<- ggplot(data, aes(x= harvest,y=final_conc, linetype=drought, shape=drought,color=planted)) + theme_classic()+ scale_color_manual(values=c("#39568CFF", "#55C667FF")) +
  scale_shape_manual(values=c(1,16))+ #scale_size(range = c(1,10), # point size range
                                                 #label = scales::dollar)+
  #expand_limits(shape = seq(2, 6, by = 1))+
  #scale_alpha_manual(values = c(2,3,4,5,6))+
  geom_point(size=2) +  #scale_size_manual(values = c(6))+ 
  geom_smooth(method = lm, se=FALSE, aes(group=summarise), inherit.aes=TRUE) + #ggtitle("Bean")
  facet_grid(~crop)+theme(
    strip.text.x = element_text(
      size = 16
    ))+
  theme(plot.title = element_text(size= 18, hjust = 0.5))+
  labs(y=expression(paste(''*mu~'g DNA ml'^-1*'nucleic acid')), color='planted')+ xlab("harvest")+theme(axis.text.x = element_text(colour = "black", size = 17,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(size = 16, colour = "black"), 
         axis.title.y= element_text(size=16, color="black"),
         axis.title.x= element_text(size=16, color="black"),
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =14, colour = "black"),
         legend.text = element_text(size = 12, colour = "black")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")), color= guide_legend(order = 1)) 
                                                                                       

plot_bean_sw

##Figure S6
ggsave(filename = "DNA_conc_revised2.tiff", plot = plot_bean_sw,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)


#calculating slopes, intercepts, and R2


library(tidyverse)
bean_dr_planted<- data %>% filter(drought=="drought" & planted=="planted" & crop=="Bean")

summary(lm(bean_dr_planted, formula = final_conc ~ harvest_1))

bean_dr_unplanted<- data %>% filter(drought=="drought" & planted=="unplanted" & crop=="Bean")

summary(lm(bean_dr_unplanted, formula = final_conc ~ harvest_1))

bean_ww_planted<- data %>% filter(drought=="watered" & planted=="planted" & crop=="Bean")

summary(lm(bean_ww_planted, formula = final_conc ~ harvest_1))

bean_ww_unplanted<- data %>% filter(drought=="watered" & planted=="unplanted" & crop=="Bean")

summary(lm(bean_ww_unplanted, formula = final_conc ~ harvest_1))


sw_dr_planted<- data %>% filter(drought=="drought"  & planted=="planted" & crop=="Switchgrass")
summary(lm(sw_dr_planted, formula = final_conc ~ harvest_1))

sw_dr_unplanted<- data %>% filter(drought=="drought"  & planted=="unplanted" & crop=="Switchgrass")
summary(lm(sw_dr_unplanted, formula = final_conc ~ harvest_1))

sw_ww_planted<- data %>% filter(drought=="watered"  & planted=="planted" & crop=="Switchgrass")
summary(lm(sw_ww_planted, formula = final_conc ~ harvest_1))

sw_ww_unplanted<- data %>% filter(drought=="watered"  & planted=="unplanted" & crop=="Switchgrass")
summary(lm(sw_ww_unplanted, formula = final_conc ~ harvest_1))

##Main manuscript Figure notes
##Figure 1 was created by GLBRC communications team
##Figure 3 was generated by collaborator





##Supplemental manuscript figure notes
##Figure S1 was experiment images 
##Figure S2 generated by GLBRC communications team
##Figure S4A was generated by collaborator using software noted in methods section of MS, Inkscape was used to combine plots
##Figure S5 generated by collaborator using software noted in methods section of MS

################################################
#######Figures and stats not reported in final MS
###############################################

####################################################################################################
####################### Percent active OTUs detected at threshold >=1 and >1 in rarefied dataset ###
######################  Percent singletons and singleton phantoms in rarefied dataset ##############
####################################################################################################


##plotting treatment differences with % phantom taxa, %  singleton phantom taxa and % active otus greater than equal to 1 ratio
##percent active otus greater than 1 using 2 methods 

dna<- read.csv("DNAcopy.csv", header=TRUE, row.names = 1)
cdna <- read.csv("cDNAcopy.csv", header=TRUE, row.names = 1)
View(dna)
View(cdna)

dna1<- dna %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="dna", -ASV)#gather columns into key value pairs
View(dna1)
cdna1<- cdna %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="cdna", -ASV)
View(cdna1)
merge.df = full_join(dna1, cdna1, by=c("ASV", "SampleID"))
View(merge.df)##10313164 taxa which is 212 samples * 48647 taxa

merge.df[is.na(merge.df)] <-0
View(merge.df)
metadata<-read.csv("metadata_dna_cdna_new.csv")

newdf = full_join(merge.df, metadata, by=c("SampleID"))
View(newdf)


finaldf = newdf %>% filter(newdf$crop == "switchgrass") %>%  mutate(ratioMethod1= ifelse(dna == 0 & cdna > 0, 100, cdna/dna)) %>%  
  mutate(dna2= ifelse(dna == 0, 1,dna))  %>%  mutate(ratioMethod2= cdna/dna2) %>% mutate(ratio_nophantoms = cdna/dna)
View(finaldf)
finaldf[is.na(finaldf)] <-0
finaldf$ratioMethod1[!is.finite(finaldf$ratioMethod1)] <- 0 ##using the code above for removing nas also works here as it changed the 0 denominator ratios from NA to 0, no infinite values only NaN
finaldf$ratioMethod2[!is.finite(finaldf$ratioMethod2)] <- 0
finaldf$ratio_nophantoms[!is.finite(finaldf$ratio_nophantoms)] <- 0
View(finaldf) #5205229 rows: 48647*107 for SW, 5107935 for bean 48647*105

library(dplyr) ##use package version 1.0.7 , this is important to have the code below work


##normalize calculations sample by sample
##Switchgrass
finaldf2 <- finaldf %>% group_by(SampleID) %>% 
  mutate(percent_gr_1_method1 = (sum(ratioMethod1 > 1)/sum(dna+cdna!=0)) *100)  %>% 
  group_by(SampleID) %>% 
  mutate(percent_gr_1_method2 = (sum(ratioMethod2 > 1)/sum(dna+cdna!=0)) *100) %>% 
  group_by(SampleID) %>%
  mutate(percent_gr_1_nophantoms = (sum(ratio_nophantoms > 1)/sum(dna+cdna!=0)) *100) %>% 
  group_by(SampleID) %>% 
  mutate(percent_phantoms = (sum(dna == 0 & cdna>0))/sum(dna+cdna!=0) *100) %>% 
  group_by(SampleID) %>% 
  mutate(percent_phantom_singleton = (sum(dna == 0 & cdna==1))/sum(dna+cdna!=0) *100)

View(finaldf2)
write.csv(finaldf2, "SW_phantomstats_allfinite_greaterthan1_normalizedtotalotus_samplebysample_new.csv")



##check for ratio>1
df_test<- finaldf2 %>% filter(SampleID=="X_125")
View(df_test)
df_test1<- df_test %>% filter(dna+cdna!=0)
View(df_test1)#2137 otus
df_test2<- df_test %>% filter(ratio_nophantoms > 1)
View(df_test2)## 135  otus
##comes out as 6.3172672 % for percent_gr_eq_1_nophantoms column- so calculation is normalized by sample now

finaldf2 <- finaldf %>% group_by(SampleID) %>% 
  mutate(percent_gr_eq_1_method1 = (sum(ratioMethod1 >= 1)/sum(dna+cdna!=0)) *100)  %>% 
  group_by(SampleID) %>% 
  mutate(percent_gr_eq_1_method2 = (sum(ratioMethod2 >= 1)/sum(dna+cdna!=0)) *100) %>% 
  group_by(SampleID) %>%
  mutate(percent_gr_eq_1_nophantoms = (sum(ratio_nophantoms >= 1)/sum(dna+cdna!=0)) *100) %>% 
  group_by(SampleID) %>% 
  mutate(percent_phantoms = (sum(dna == 0 & cdna>0))/sum(dna+cdna!=0) *100) %>% 
  group_by(SampleID) %>% 
  mutate(percent_phantom_singleton = (sum(dna == 0 & cdna==1))/sum(dna+cdna!=0) *100)

View(finaldf2)
write.csv(finaldf2, "SW_phantomstats_allfinite_greaterthaneq1_normalizedtotalotus_samplebysample_new.csv")

##check for ratio>=1
df_test<- finaldf2 %>% filter(SampleID=="X_125")
View(df_test)
df_test1<- df_test %>% filter(dna+cdna!=0)
View(df_test1)#2137 otus
df_test2<- df_test %>% filter(ratio_nophantoms >= 1)
View(df_test2)##138 otus
##comes out as 6.45765 % for percent_gr_eq_1_nophantoms column- so calculation is normalized by sample now

#bean 
finaldf4 <- finaldf %>% group_by(SampleID) %>% 
  mutate(percent_gr_1_method1 = (sum(ratioMethod1 > 1)/sum(dna+cdna!=0)) *100)  %>% 
  group_by(SampleID) %>% 
  mutate(percent_gr_1_method2 = (sum(ratioMethod2 > 1)/sum(dna+cdna!=0)) *100) %>% 
  group_by(SampleID) %>%
  mutate(percent_gr_1_nophantoms = (sum(ratio_nophantoms > 1)/sum(dna+cdna!=0)) *100) %>% 
  group_by(SampleID) %>% 
  mutate(percent_phantoms = (sum(dna == 0 & cdna>0))/sum(dna+cdna!=0) *100) %>% 
  group_by(SampleID) %>% 
  mutate(percent_phantom_singleton = (sum(dna == 0 & cdna==1))/sum(dna+cdna!=0) *100)
View(finaldf4)
write.csv(finaldf4, "Bean_phantomstats_allfinite_greaterthan1_normalizedtotalotus_samplebysample_new.csv")


##check for >1

df_test3<- finaldf4 %>% filter(SampleID=="X_127")
View(df_test3)
df_test4<- df_test3 %>% filter(dna+cdna!=0)
View(df_test4)#  2363 otus
df_test5<- df_test3 %>% filter(ratio_nophantoms > 1)
View(df_test5)## 245 otus
##comes out as 10.368176  % for percent_gr_1_nophantoms column- so calculation is normalized by sample now


finaldf4 <- finaldf %>% group_by(SampleID) %>% 
  mutate(percent_gr_eq_1_method1 = (sum(ratioMethod1 >= 1)/sum(dna+cdna!=0)) *100)  %>% 
  group_by(SampleID) %>% 
  mutate(percent_gr_eq_1_method2 = (sum(ratioMethod2 >= 1)/sum(dna+cdna!=0)) *100) %>% 
  group_by(SampleID) %>%
  mutate(percent_gr_eq_1_nophantoms = (sum(ratio_nophantoms >= 1)/sum(dna+cdna!=0)) *100) %>% 
  group_by(SampleID) %>% 
  mutate(percent_phantoms = (sum(dna == 0 & cdna>0))/sum(dna+cdna!=0) *100) %>% 
  group_by(SampleID) %>% 
  mutate(percent_phantom_singleton = (sum(dna == 0 & cdna==1))/sum(dna+cdna!=0) *100)
View(finaldf4)
write.csv(finaldf4, "Bean_phantomstats_allfinite_greaterthaneq1_normalizedtotalotus_samplebysample_new.csv")

##check for >=1

df_test3<- finaldf4 %>% filter(SampleID=="X_127")
View(df_test3)
df_test4<- df_test3 %>% filter(dna+cdna!=0)
View(df_test4)#  2363 otus
df_test5<- df_test3 %>% filter(ratio_nophantoms >= 1)
View(df_test5)## 302 otus
##comes out as 12.78036  % for percent_gr_eq_1_nophantoms column- so calculation is normalized by sample now

###boxplots for ratio threshold >1, SW

#long format to combine method 1 and 2 together for facetting

sw_phantoms<- read.csv("SW_phantomstats_allfinite_greaterthan1_normalizedtotalotus_samplebysample_new.csv")
library(tidyverse)

library(dplyr)

sw_phantoms$harvest <- recode(sw_phantoms$harvest, predrought = 'Day 0', 
                        harvest1 = 'Day 2',
                        harvest2 = 'Day 3',
                        harvest3 = 'Day 4',
                        harvest4 ='Day 5',
                        harvest5 = 'Day 6')


sw_phantoms$drought <- recode(sw_phantoms$drought, 'well-watered' = 'watered',predrought = 'pre-drought')

View(sw_phantoms)


sw_phantoms_longer<- sw_phantoms %>% pivot_longer(cols=starts_with("percent_gr"),
                                                  names_to = "method",
                                                  values_to="percent_active")
View(sw_phantoms_longer)
df_test<- sw_phantoms_longer%>% filter(treatment=="PreDrP1")
View(df_test)

#library(ggplot2)
#library(car)
#packageDescription("car")
#facetlabel <- c("method1 (threshold >= 1)", "method2 (threshold >=1)")

##Make a modified copy of the original data to change the facet names in method column
#sw_phantoms_longer<- as.tibble(sw_phantoms_longer)
sw_phantoms_longer <- sw_phantoms_longer %>%
  mutate(method = recode(method,
                         "percent_gr_1_method1" = "method1 (threshold > 1)",
                         "percent_gr_1_method2" = "method2 (threshold > 1)",
                         "percent_gr_1_nophantoms" = "no_phantoms_included(threshold > 1)"
  ))


View(sw_phantoms_longer)

##boxplot of ratio threshold >1 , bean


bean_phantoms<- read.csv("Bean_phantomstats_allfinite_greaterthan1_normalizedtotalotus_samplebysample_new.csv")
library(tidyverse)
bean_phantoms$harvest <- recode(bean_phantoms$harvest, predrought = 'Day 0', 
                             harvest1 = 'Day 2',
                             harvest2 = 'Day 3',
                             harvest3 = 'Day 4',
                             harvest4 ='Day 5',
                             harvest5 = 'Day 6')


bean_phantoms$drought <- recode(bean_phantoms$drought, 'well-watered' = 'watered',predrought = 'pre-drought')


bean_phantoms_longer<- bean_phantoms %>% pivot_longer(cols=starts_with("percent_gr"),
                                                      names_to = "method",
                                                      values_to="percent_active")
View(bean_phantoms_longer)
df_test<- bean_phantoms_longer%>% filter(treatment=="PreDrP1")
View(df_test)

library(ggplot2)

#facetlabel <- c("method1 (threshold >= 1)", "method2 (threshold >=1)")

##Make a modified copy of the original data to change the facet names in method column

bean_phantoms_longer <- bean_phantoms_longer %>%
  mutate(method = recode(method,
                         "percent_gr_1_method1" = "method1 (threshold > 1)",
                         "percent_gr_1_method2" = "method2 (threshold > 1)",
                         "percent_gr_1_nophantoms" = "no_phantoms_included(threshold > 1)"
  ))

View(bean_phantoms_longer)

##reorder levels of harvest
level_order <- c('Day 0', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6')
#plot
plot1 <- ggplot(sw_phantoms_longer, aes(x=factor(harvest, level=level_order), y=percent_active,fill=drought)) + facet_grid(planted~method) +
  theme(strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent active OTUs") + ggtitle("Switchgrass")+ theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y=element_text(size=17),
                                      axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))  + scale_fill_manual(limits = c("pre-drought", "drought", "watered"), values=c("#D55E00","#009E73", "#56B4E9"))
plot1
ggsave(filename = "SW_percent_active_threshold>1_colorblindfriendly_newlabel.tiff", plot = plot1,
       width = 24,
       height = 18, units = c("cm"),
       dpi = 300)

##reorder levels of harvest
level_order <- c('Day 0', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6')
#plot
plot2 <- ggplot(bean_phantoms_longer, aes(x=factor(harvest, level=level_order), y=percent_active,fill=drought)) + facet_grid(planted~method) +
  theme(strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 12))+ geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent active OTUs") + ggtitle("Bean")+theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y=element_text(size=17),
                                      axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))  + scale_fill_manual(limits = c("pre-drought", "drought", "watered"), values=c("#D55E00","#009E73", "#56B4E9"))
plot2
ggsave(filename = "bean_percent_active_threshold>1_colorblindfriendly_newlabel.tiff", plot = plot2,
       width = 24,
       height = 18, units = c("cm"),
       dpi = 300)


##generate plots 1 and 2 again with threshold >=1 , SW

sw_phantoms<- read.csv("SW_phantomstats_allfinite_greaterthaneq1_normalizedtotalotus_samplebysample_new.csv")
library(tidyverse)

library(dplyr)

sw_phantoms$harvest <- recode(sw_phantoms$harvest, predrought = 'Day 0', 
                        harvest1 = 'Day 2',
                        harvest2 = 'Day 3',
                        harvest3 = 'Day 4',
                        harvest4 ='Day 5',
                        harvest5 = 'Day 6')


sw_phantoms$drought <- recode(sw_phantoms$drought, 'well-watered' = 'watered',predrought = 'pre-drought')

View(sw_phantoms)


sw_phantoms_longer<- sw_phantoms %>% pivot_longer(cols=starts_with("percent_gr"),
                                                  names_to = "method",
                                                  values_to="percent_active")
View(sw_phantoms_longer)
df_test<- sw_phantoms_longer%>% filter(treatment=="PreDrP1")
View(df_test)

#library(ggplot2)
#library(car)
#packageDescription("car")
#facetlabel <- c("method1 (threshold >= 1)", "method2 (threshold >=1)")

##Make a modified copy of the original data to change the facet names in method column
#sw_phantoms_longer<- as.tibble(sw_phantoms_longer)
sw_phantoms_longer <- sw_phantoms_longer %>%
  mutate(method = recode(method,
                         "percent_gr_eq_1_method1" = "method1 (threshold >= 1)",
                         "percent_gr_eq_1_method2" = "method2 (threshold >= 1)",
                         "percent_gr_eq_1_nophantoms" = "no_phantoms_included(threshold >= 1)"
  ))


View(sw_phantoms_longer)

##boxplot of ratio threshold >=1 , bean


bean_phantoms<- read.csv("Bean_phantomstats_allfinite_greaterthaneq1_normalizedtotalotus_samplebysample_new.csv")
library(tidyverse)
bean_phantoms$harvest <- recode(bean_phantoms$harvest, predrought = 'Day 0', 
                             harvest1 = 'Day 2',
                             harvest2 = 'Day 3',
                             harvest3 = 'Day 4',
                             harvest4 ='Day 5',
                             harvest5 = 'Day 6')


bean_phantoms$drought <- recode(bean_phantoms$drought, 'well-watered' = 'watered',predrought = 'pre-drought')


bean_phantoms_longer<- bean_phantoms %>% pivot_longer(cols=starts_with("percent_gr"),
                                                      names_to = "method",
                                                      values_to="percent_active")
View(bean_phantoms_longer)
df_test<- bean_phantoms_longer%>% filter(treatment=="PreDrP1")
View(df_test)

library(ggplot2)

#facetlabel <- c("method1 (threshold >= 1)", "method2 (threshold >=1)")

##Make a modified copy of the original data to change the facet names in method column

bean_phantoms_longer <- bean_phantoms_longer %>%
  mutate(method = recode(method,
                         "percent_gr_eq_1_method1" = "method1 (threshold >= 1)",
                         "percent_gr_eq_1_method2" = "method2 (threshold >= 1)",
                         "percent_gr_eq_1_nophantoms" = "no_phantoms_included(threshold >= 1)"
  ))

View(bean_phantoms_longer)

##reorder levels of harvest
level_order <- c('Day 0', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6')
#plot
plot1 <- ggplot(sw_phantoms_longer, aes(x=factor(harvest, level=level_order), y=percent_active,fill=drought)) + facet_grid(planted~method) +
  theme(strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent active OTUs") + ggtitle("Switchgrass")+ theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y=element_text(size=17),
                                      axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))  + scale_fill_manual(limits = c("pre-drought", "drought", "watered"), values=c("#D55E00","#009E73", "#56B4E9"))
plot1
ggsave(filename = "SW_percent_active_threshold>=1_colorblindfriendly_newlabel.tiff", plot = plot1,
       width = 24,
       height = 18, units = c("cm"),
       dpi = 300)

##reorder levels of harvest
level_order <- c('Day 0', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6')
#plot
plot2 <- ggplot(bean_phantoms_longer, aes(x=factor(harvest, level=level_order), y=percent_active,fill=drought)) + facet_grid(planted~method) +
  theme(strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 12))+ geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent active OTUs") + ggtitle("Bean")+theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y=element_text(size=17),
                                      axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))  + scale_fill_manual(limits = c("pre-drought", "drought", "watered"), values=c("#D55E00","#009E73", "#56B4E9"))
plot2
ggsave(filename = "bean_percent_active_threshold>=1_colorblindfriendly_newlabel.tiff", plot = plot2,
       width = 24,
       height = 18, units = c("cm"),
       dpi = 300)



##boxplots of percent phantoms and singleton phantoms

sw_phantoms<- read.csv("SW_phantomstats_allfinite_greaterthaneq1_normalizedtotalotus_samplebysample_new.csv")

sw_phantoms$harvest <- recode(sw_phantoms$harvest, predrought = 'Day 0', 
                              harvest1 = 'Day 2',
                              harvest2 = 'Day 3',
                              harvest3 = 'Day 4',
                              harvest4 ='Day 5',
                              harvest5 = 'Day 6')


sw_phantoms$drought <- recode(sw_phantoms$drought, 'well-watered' = 'watered',predrought = 'pre-drought')



level_order <- c('Day 0', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6')
#plot
plot3 <- ggplot(sw_phantoms, aes(x=factor(harvest, level=level_order), y=percent_phantoms,fill=drought)) + facet_grid(~planted) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent phantom OTUs") + ggtitle("Switchgrass Percent Phantom OTUs")+ theme(plot.title = element_text(hjust = 0.5))+theme(axis.text.y=element_text(size=10),
                                                 axis.title.y=element_text(size=15)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=13))+ scale_fill_manual(limits = c("pre-drought", "drought", "watered"), values=c("#D55E00","#009E73", "#56B4E9"))
  
plot3 
ggsave(filename = "percent_phantoms_SW_colorblindfriendly_newlabel.tiff", plot = plot3,
       width = 18,
       height = 10, units = c("cm"),
       dpi = 300)

plot4 <- ggplot(sw_phantoms, aes(x=factor(harvest, level=level_order), y=percent_phantom_singleton,fill=drought)) + facet_grid(~planted) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent phantom singleton OTUs") + ggtitle("Switchgrass Percent Phantom Singleton OTUs")+ theme(plot.title = element_text(hjust = 0.5))+theme(axis.text.y=element_text(size=10),
                                       axis.title.y=element_text(size=15)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=13))+ scale_fill_manual(limits = c("pre-drought", "drought", "watered"), values=c("#D55E00","#009E73", "#56B4E9"))
  
plot4
ggsave(filename = "percent_phantom_singleton_SW_colorblindfriendly_newlabel.tiff", plot = plot4,
       width = 18,
       height = 10, units = c("cm"),
       dpi = 300)

#ylim(0,80) ##percent phantoms
#ylim(0,20) ##percent phantom singleton

bean_phantoms<- read.csv("Bean_phantomstats_allfinite_greaterthaneq1_normalizedtotalotus_samplebysample_new.csv")


bean_phantoms$harvest <- recode(bean_phantoms$harvest, predrought = 'Day 0', 
                                harvest1 = 'Day 2',
                                harvest2 = 'Day 3',
                                harvest3 = 'Day 4',
                                harvest4 ='Day 5',
                                harvest5 = 'Day 6')


bean_phantoms$drought <- recode(bean_phantoms$drought, 'well-watered' = 'watered',predrought = 'pre-drought')

level_order <- c('Day 0', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6')
#plot
plot5 <- ggplot(bean_phantoms, aes(x=factor(harvest, level=level_order), y=percent_phantoms,fill=drought)) + facet_grid(~planted) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent phantom OTUs") + ggtitle("Bean Percent Phantom OTUs")+ theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y=element_text(size=10),
                                       axis.title.y=element_text(size=15)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=13)) + scale_fill_manual(limits = c("pre-drought", "drought", "watered"), values=c("#D55E00","#009E73", "#56B4E9"))
  
plot5 
ggsave(filename = "percent_phantoms_bean_colorblindfriendly_newlabel.tiff", plot = plot5,
       width = 18,
       height = 10, units = c("cm"),
       dpi = 300)

plot6 <- ggplot(bean_phantoms, aes(x=factor(harvest, level=level_order), y=percent_phantom_singleton,fill=drought)) + facet_grid(~planted) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent phantom singleton OTUs") + ggtitle("Bean Percent Phantom Singleton OTUs")+ theme(plot.title = element_text(hjust = 0.5))+theme(axis.text.y=element_text(size=10),
                                       axis.title.y=element_text(size=15)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=13)) + scale_fill_manual(limits = c("pre-drought", "drought", "watered"), values=c("#D55E00","#009E73", "#56B4E9"))
  
plot6
ggsave(filename = "percent_phantom_singleton_bean_colorblindfriendly_newlabel.tiff", plot = plot6,
       width = 18,
       height = 10, units = c("cm"),
       dpi = 300)

#ylim(0,80) ##percent phantoms
#ylim(0,20) ##percent phantom singleton

##for arranging in grid for publication
library(ggpubr)
percent_active_bythreshold<- ggarrange(plot2, plot1, ncol=1, nrow=2, widths=c(1,1),heights=c(1,1), common.legend = TRUE, legend="bottom", labels= c("A", "B"))
View(percent_active_bythreshold)
ggsave("combined_bean_sw_percentactive>1_bythreshold_newlabel.tiff",percent_active_bythreshold,width=22,height=26, units="cm")

##repeat for threshold >=1 with new plots 1 and 2
percent_active_bythreshold<- ggarrange(plot2, plot1, ncol=1, nrow=2, widths=c(1,1),heights=c(1,1), common.legend = TRUE, legend="bottom", labels= c("A", "B"))
View(percent_active_bythreshold)
ggsave("combined_bean_sw_percentactive>=1_bythreshold_newlabel.tiff",percent_active_bythreshold,width=22,height=26, units="cm")


#ggsave("Figure3.eps",Figure3,width=6,height=4, units="in")

percent_phantoms_combined<- ggarrange(plot5, plot6, plot3, plot4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom", labels= c("A", "B", "C", "D"))
View(percent_phantoms_combined)
ggsave("combined_bean_sw_percentphantoms_newlabel.tiff",percent_phantoms_combined,width=27,height=27, units="cm")



#######
##ANCOVA by drought treatments (drought versus watered)
######

ancova<- read.csv("ancova_metadata_withsingletons_nosecondrarefaction.csv")

##slope significances
ancova_bean <-subset(ancova, crop=="switchgrass" & planted_y=="unplanted")
mod1 <- aov(similarity~drought*days, data=ancova_bean) ##test interaction
summary(mod1)

##intercept significances
mod2 <- aov(similarity~drought+days, data=ancova_bean)
summary(mod2)


##slopes and intercept significances (regression)

reg1 <- lm(similarity~drought/days -1, data=ancova_bean); summary(reg1) 


##ANCOVA by crop treatments (bean vs switchgrass)

ancova<- read.csv("ancova_metadata_withsingletons_nosecondrarefaction.csv")

##slope significances
ancova_bean <-subset(ancova, drought=="well-watered" & planted_y=="unplanted")
mod1 <- aov(similarity~crop*days, data=ancova_bean) ##test interaction
summary(mod1)


##intercept significances
mod2 <- aov(similarity~crop+days, data=ancova_bean)
summary(mod2)

##slopes and intercept significances (regression)

reg1 <- lm(similarity~crop/days -1, data=ancova_bean); summary(reg1) 




#####Alpha diversity unused code in final manuscript


##combine bean and switchgrass richness and evenness together in panel

#alphadiv_bean_InvS, alphadiv_bean_rich, alphadiv_sw_InvS, alphadiv_sw_rich
View(alphadiv_bean_InvS)
alphadiv_bean_IS_new<- alphadiv_bean_InvS %>% mutate(Crop=ifelse(crop=="bean", "Bean", NA)) 
View(alphadiv_bean_IS_new)

View(alphadiv_bean_rich)
alphadiv_bean_rich_new<- alphadiv_bean_rich %>% mutate(Crop=ifelse(crop=="bean", "Bean", NA)) 
View(alphadiv_bean_rich_new)

View(alphadiv_sw_InvS)
alphadiv_sw_IS_new<- alphadiv_sw_InvS %>% mutate(Crop=ifelse(crop=="switchgrass", "Switchgrass", NA))
View(alphadiv_sw_IS_new)

View(alphadiv_sw_rich)
alphadiv_sw_rich_new<- alphadiv_sw_rich %>% mutate(Crop=ifelse(crop=="switchgrass", "Switchgrass", NA))
View(alphadiv_sw_rich_new)

alphadiv_combined= rbind(alphadiv_sw_rich_new, alphadiv_bean_rich_new, alphadiv_sw_IS_new, alphadiv_bean_IS_new)
View(alphadiv_combined)

alphadiv_combined$harvest <- factor(alphadiv_combined$harvest,levels = c("pre-drought","harvest1","harvest2","harvest3", "harvest4", "harvest5"), ordered=TRUE)

write.csv(alphadiv_combined, "alphadiv_combined.csv") ##edit to have pre drought in drought and well watered levels
alphadiv_combined<- read.csv("alphadiv_combined.csv")
alphadiv_combined$harvest <- factor(alphadiv_combined$harvest,levels = c("pre-drought","harvest1","harvest2","harvest3", "harvest4", "harvest5"), ordered=TRUE)

alphadiv_combined<- read.csv("alphadiv_combined_LabelChange.csv") ##rename harvest to days in csv file 
alphadiv_combined$harvest <- factor(alphadiv_combined$harvest,levels = c("Day 0","Day 2","Day 3","Day 4", "Day 5", "Day 6"), ordered=TRUE)

View(alphadiv_combined)

alphadiv_rich<- alphadiv_combined %>% filter (measure=="Richness (Number of observed OTUs)")
View(alphadiv_rich)
#sw_dr_BC_time<- rename(sw_dr_BC_time, Day = group)
#alphadiv_combined$harvest <- recode(alphadiv_combined$harvest, 
# `2` = 'Day 2',
#  `3` = 'Day 3',
#  `4` = 'Day 4',
# `5` ='Day 5',
# `6` = 'Day 6')



#sw_dr_BC_time$drought <- recode(sw_dr_BC_time$drought, 'well-watered' = 'watered')

alphadiv_rich$crop<- recode(alphadiv_rich$crop, 'switchgrass' = 'Switchgrass')

alphadiv_rich$crop<- recode(alphadiv_rich$crop, 'bean' = 'Bean')

#sw_dr_BC_time$Day<-as.factor(as.numeric(sw_dr_BC_time$Day))

plot_bean_sw<- ggplot(alphadiv_rich, aes(x= harvest,y=mean, linetype=drought, shape=drought,color=planted)) + theme_classic()+ scale_color_manual(values=c("#39568CFF", "#55C667FF")) +
  scale_shape_manual(values=c(1,16))+ #scale_size(range = c(1,10), # point size range
  #label = scales::dollar)+
  #expand_limits(shape = seq(2, 6, by = 1))+
  #scale_alpha_manual(values = c(2,3,4,5,6))+
  geom_point(aes(size=harvest)) +  scale_size_manual(values = c(1,2,3,4,5,6))+ geom_smooth(method = lm, se=FALSE, aes(group=summarise_new), inherit.aes=TRUE) + #ggtitle("Bean")
  facet_grid(~crop)+theme(
    strip.text.x = element_text(
      size = 16
    ))+
  theme(plot.title = element_text(size= 18, hjust = 0.5))+
  labs(y = "Richness (Number of observed OTUs)", color='planted')+ xlab("harvest")+theme(axis.text.x = element_text(colour = "black", size = 17,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(size = 16, colour = "black"), 
         axis.title.y= element_text(size=16, color="black"),
         axis.title.x= element_text(size=16, color="black"),
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =14, colour = "black"),
         legend.text = element_text(size = 12, colour = "black")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")), color= guide_legend(order = 1)) 


plot_bean_sw


##Figure S14
ggsave(filename="alphadiv_combined_newlabel_sameasBCsim_FINAL.TIFF", plot=plot_bean_sw, width=25, height=15, units="cm", dpi=300)






############
##ANCOVA for alpha diversity over time, catergorical variables- planted, covariate: time, categorical variable 2: crop, cat var 3: drought, dependent variable richness/evenness

ancova<- read.csv("ancova_alphadiv.csv")
View(ancova)
ancova<- ancova %>% filter(measure=="Richness") ##change to inverse simpson when testing that variable
#ancova<- ancova %>% filter(measure=="Inverse Simpson")
mod1 <- aov(mean~crop*planted*drought*days, data=ancova) ##test interaction
summary(mod1)

###subset to bean and sw
ancova_bean <-subset(ancova, crop=="bean")

##if crop level is included, interactions exist between planted:drought and crop:planted (same as beta diversity similarity changes)

##separate bean and SW

###with_singletons_no_second_rarefaction

ancova_bean <-subset(ancova, crop=="bean")
mod1 <- aov(mean~planted*drought*days, data=ancova_bean) ##test interaction
summary(mod1)

##no interaction for bean , slope is not different between treatments (in this case drought and also planted levels)
#no interaction between time and planted , and time and drought, so no difference in slope


##with singletons no second rarefaction

mod2 <- aov(mean~planted+drought+days, data=ancova_bean) ##if no interaction, test differences in intercept within categorical variable eg. sex, if significant then subset by sex and look for effects of covariate (snout)

summary(mod2)

## difference in intercept: planted levels


anova(mod1,mod2)

##removing interaction term for bean doesnt significantly alter the fit of the model.
##therefore we can conclude that most parsimonious model is mod2.


##bean subset to drought 
ancova_bean_dr <-subset(ancova, crop=="bean" & drought=="drought")
mod1 <- aov(mean~planted*days, data=ancova_bean_dr) ##test interaction
summary(mod1)


mod2 <- aov(mean~planted+days, data=ancova_bean_dr)
summary(mod2)  


options(contrasts = c("contr.sum", "contr.poly"))
lm.bean<- lm(mean~planted+days, data=ancova_bean_dr)
mod2 <- car::Anova(lm.bean, type=2)
mod2 

##since same results for both type 2 and type 1, use the results with type 1

anova(mod1, mod2)


##bean well-watered

ancova_bean_ww <-subset(ancova, crop=="bean" & drought=="well-watered")
mod1 <- aov(mean~planted*days, data=ancova_bean_ww) ##test interaction
summary(mod1)

mod2 <- aov(mean~planted+days, data=ancova_bean_ww)
summary(mod2)  


options(contrasts = c("contr.sum", "contr.poly"))
lm.bean<- lm(mean~planted+days, data=ancova_bean_ww)
mod2 <- car::Anova(lm.bean, type=2)
mod2 


anova(mod1, mod2) ##since same results for type 1 and type 2 used the type 1


##Bean subset to drought and well watered ##with singletons no second rarefaction

dr_bean <- subset(ancova_bean, drought=="drought")


reg1 <- lm(mean~planted/days -1, data=dr_bean); summary(reg1) ##subset data, nest days within planted

##drought bean results
##wellwatered bean results


##with singletons no second rarefaction
ancova_switchgrass<- subset(ancova, crop=="switchgrass")
mod1 <- aov(mean~planted*drought*days, data=ancova_switchgrass) ##test interaction
summary(mod1)


###no interaction for SW, slope is not different between treatments
##for SW
#no interaction between time and planted , or time and drought so no difference in slope, minor interaction between planted and drought but not significant.

#with singletons no second rarefaction

mod2 <- aov(mean~planted+drought+days, data=ancova_switchgrass) ##if no interaction, test differences in intercept within categorical variable eg. sex, if significant then subset by sex and look for effects of covariate (snout)

summary(mod2)

##no significant differences in intercept for SW data **** (planted factor and drought factor considered here)

anova(mod1,mod2) ##does removing the interaction affect the fit of the model? if not then proceed with model 2 and subset the data based on sex
##sepaarte regression for categirical variables testing effects of snout on y var, due to effects of sex on intercept sepaarte sexes and assess effect of snout on dependent variable (pelvic)


##with singletons no second rarefaction


##removing interaction term for SW doesnt significantly alter the fit of the model.
##therefore we can conclude that most parsimonious model is mod2.


##SW subset to drought 
ancova_sw_dr <-subset(ancova, crop=="switchgrass" & drought=="drought")
mod1 <- aov(mean~planted*days, data=ancova_sw_dr) ##test interaction
summary(mod1)

mod2 <- aov(mean~planted+days, data=ancova_sw_dr)
summary(mod2)  


options(contrasts = c("contr.sum", "contr.poly"))
lm.sw<- lm(mean~planted+days, data=ancova_sw_dr)
mod2 <- car::Anova(lm.sw, type=2)
mod2 


##since same results for both type 2 and type 1, use the results with type 1

anova(mod1, mod2)

##SW  well-watered

ancova_sw_ww <-subset(ancova, crop=="switchgrass" & drought=="well-watered")
mod1 <- aov(mean~planted*days, data=ancova_sw_ww) ##test interaction
summary(mod1)


mod2 <- aov(mean~planted+days, data=ancova_sw_ww)
summary(mod2)  


options(contrasts = c("contr.sum", "contr.poly"))
lm.sw<- lm(mean~planted+days, data=ancova_sw_ww)
mod2 <- car::Anova(lm.sw, type=2)
mod2 

anova(mod1, mod2) ##since same results for type 1 and type 2 used the type 1

##SW subset to drought and well watered ##with singletons no second rarefaction

dr_sw <- subset(ancova_switchgrass, drought=="drought" )

reg1 <- lm(mean~planted/days -1, data=dr_sw); summary(reg1) ##subset data, nest days within planted

##drought SW results
##wellwatered SW results


### Resistance Calculation from Jackson's code (Dormancy dispersal paper)
##Calculate resistance based on BC similarity at harvest timepoint 5 (this was the final time point of stress)


T5_map_DNA <- read.csv("resistance.csv")

T5_map_DNA_Disturbed <- T5_map_DNA[T5_map_DNA$drought=="Y",]

View(T5_map_DNA_Disturbed)

#y0_T5_SW_planted <-  mean(T5_map_DNA$ActiveSimilarity_to_predr[T5_map_DNA$drought=="N"& T5_map_DNA$planted=="planted" & T5_map_DNA$crop=="switchgrass"])

#y0_T5_SW_unplanted <-  mean(T5_map_DNA$ActiveSimilarity_to_predr[T5_map_DNA$drought=="N"& T5_map_DNA$planted=="unplanted" & T5_map_DNA$crop=="switchgrass"])

#y0_T5_bean_planted <-  mean(T5_map_DNA$ActiveSimilarity_to_predr[T5_map_DNA$drought=="N"& T5_map_DNA$planted=="planted" & T5_map_DNA$crop=="bean"])

#y0_T5_bean_unplanted <-  mean(T5_map_DNA$ActiveSimilarity_to_predr[T5_map_DNA$drought=="N"& T5_map_DNA$planted=="unplanted" & T5_map_DNA$crop=="bean"])

#yl_T5_SW_planted <- T5_map_DNA_Disturbed$ActiveSimilarity_to_predr[T5_map_DNA_Disturbed$planted=="planted" & T5_map_DNA_Disturbed$crop=="switchgrass"]

#yl_T5_SW_unplanted <-  T5_map_DNA_Disturbed$ActiveSimilarity_to_predr[T5_map_DNA_Disturbed$planted=="unplanted" & T5_map_DNA_Disturbed$crop=="switchgrass"]

#yl_T5_bean_planted <-  T5_map_DNA_Disturbed$ActiveSimilarity_to_predr[T5_map_DNA_Disturbed$planted=="planted" & T5_map_DNA_Disturbed$crop=="bean"]

#yl_T5_bean_unplanted <-  T5_map_DNA_Disturbed$ActiveSimilarity_to_predr[T5_map_DNA_Disturbed$planted=="unplanted" & T5_map_DNA_Disturbed$crop=="bean"]


##when comparing planted drought (treatment) and unplanted drought (control) treatments for resistance calculations
yl_T5_SW_planted <- T5_map_DNA_Disturbed$ActiveSimilarity_to_predr[T5_map_DNA_Disturbed$planted=="planted" & T5_map_DNA_Disturbed$crop=="switchgrass"]

y0_T5_SW_unplanted <-  mean(T5_map_DNA_Disturbed$ActiveSimilarity_to_predr[T5_map_DNA_Disturbed$planted=="unplanted" & T5_map_DNA_Disturbed$crop=="switchgrass"])

yl_T5_bean_planted <-  T5_map_DNA_Disturbed$ActiveSimilarity_to_predr[T5_map_DNA_Disturbed$planted=="planted" & T5_map_DNA_Disturbed$crop=="bean"]

y0_T5_bean_unplanted <-  mean(T5_map_DNA_Disturbed$ActiveSimilarity_to_predr[T5_map_DNA_Disturbed$planted=="unplanted" & T5_map_DNA_Disturbed$crop=="bean"])


##subset data

SW_planted <- T5_map_DNA_Disturbed[T5_map_DNA_Disturbed$planted=="planted" & T5_map_DNA_Disturbed$crop=="switchgrass",]

#SW_unplanted <-  T5_map_DNA_Disturbed[T5_map_DNA_Disturbed$planted=="unplanted" & T5_map_DNA_Disturbed$crop=="switchgrass",]

bean_planted <-  T5_map_DNA_Disturbed[T5_map_DNA_Disturbed$planted=="planted" & T5_map_DNA_Disturbed$crop=="bean",]

#bean_unplanted <-  T5_map_DNA_Disturbed[T5_map_DNA_Disturbed$planted=="unplanted" & T5_map_DNA_Disturbed$crop=="bean",]



library(tidyverse)

#put all data frames into list

##resistance for planted bean and sw treatments during drought compared to unplanted controls 

SW_planted$ActiveResistance <- 1- ((2*abs(y0_T5_SW_unplanted-yl_T5_SW_planted))/(y0_T5_SW_unplanted+abs(y0_T5_SW_unplanted - yl_T5_SW_planted))) ## should the last part of the equation have + or -? Does not match with the formula on the paper

#SW_unplanted$ActiveResistance <- 1- ((2*abs(y0_T5_SW_unplanted-yl_T5_SW_unplanted))/(y0_T5_SW_unplanted+abs(y0_T5_SW_unplanted - yl_T5_SW_unplanted))) ## should the last part of the equation have + or -? Does not match with the formula on the paper

bean_planted$ActiveResistance <- 1- ((2*abs(y0_T5_bean_unplanted-yl_T5_bean_planted))/(y0_T5_bean_unplanted+abs(y0_T5_bean_unplanted - yl_T5_bean_planted))) ## should the last part of the equation have + or -? Does not match with the formula on the paper

#bean_unplanted$ActiveResistance <- 1- ((2*abs(y0_T5_bean_unplanted-yl_T5_bean_unplanted))/(y0_T5_bean_unplanted+abs(y0_T5_bean_unplanted - yl_T5_bean_unplanted))) ## should the last part of the equation have + or -? Does not match with the formula on the paper

df_all <- rbind(SW_planted, bean_planted)

View(df_all)

write.csv(df_all, "resistance_plot.csv")

RS_plot <- ggplot(df_all, aes(x=group_1, y=ActiveResistance))+ facet_grid(~crop) +
  geom_boxplot()+ ylim(0,1)+
  geom_point()+
  ylab(label = "Resistance")+
  xlab(label = NULL)+
  theme_bw()+
  scale_x_discrete(labels=c("Day 6"))+
  theme(axis.text.x = element_text(size=15))
RS_plot 

### Figure S19

ggsave(filename = "resistance_newlabel.tiff", plot = RS_plot,
       width = 8,
       height = 8, units = c("cm"),
       dpi = 300)

#anova

anova<- aov(ActiveResistance ~ crop, data=df_all )

summary(anova)

###########################
###########################
## Pie charts (did not use for the final version)
###########################
###########################


metadata<- read.csv("bean_sw_dr_withpredr.csv")

View(metadata)
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) 
phyloseq_merged_final


phyloseq_class <- phyloseq_merged_final %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_class)


phyloseqclass1<- phyloseq_class %>% group_by(summarise, Class, planted, crop, harvest) %>% summarise(mean=mean(Abundance)) 

test<- phyloseqclass1 %>% filter(summarise=="harvest5planted" & planted == "planted" & crop=="bean")
View(test)

total<-sum(test$mean)
total
phyloseqclass2<- test %>% group_by(summarise, planted) %>% mutate(percent=(mean/(sum(test$mean)))*100) 

library(dplyr)
detach_package("plyr")
library(dplyr)


phyloseqclass2<- phyloseqclass1 %>% dplyr::group_by(summarise, planted, crop, harvest) %>% dplyr::mutate(percent=(mean/(sum(mean)))*100)

data<- phyloseqclass2 %>% filter(summarise=="harvest1planted" & planted=="planted" & crop=="bean") 

View(data)
sum(data$percent)
sum(phyloseqclass2$percent)

View(phyloseqclass1)

View(phyloseqclass2)
View(phyloseqclass_1)


sum(phyloseq_class$Abundance)
unique(phyloseq_class$OTU)



##pie chart
library(ggplot2)
library(scales)
library(RColorBrewer)
library(viridis)
library(rcartocolor)

library(Polychrome)
set.seed(935234)
P42 <- createPalette(42, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P42)
P42 <- sortByHue(P42)
P42 <- as.vector(t(matrix(P42, ncol=3)))
swatch(P42)
names(P42) <- NULL



phyloseqclass2$harvest <- factor(phyloseqclass2$harvest, levels=c("pre-drought", "harvest1", "harvest2", "harvest3", "harvest4", "harvest5"))

bp<- ggplot(phyloseqclass2, aes(x="", y=percent, fill=Class))+geom_bar(width = 1, stat = "identity")
bp

pie <- bp + facet_grid(crop+planted~harvest, labeller = label_wrap_gen(width=25))+ theme(strip.background =element_rect(fill="white"), strip.text.x = element_text(size = 20, colour = "black"))+geom_col(position = "fill") +coord_polar("y", start=0)
pie

pie_bean_sw_dr_withPredr_planted_unplanted <-pie  + theme_void()+ theme(strip.background =element_blank(), strip.text.x = element_text(size = 13, colour = "black", vjust=1), strip.placement = "inside")+theme(axis.title.x=element_blank(),
                                                                                                                                                                                                                axis.text.x=element_blank(),
                                                                                                                                                                                                                axis.ticks.x=element_blank())+ theme(axis.title.y=element_blank(),
                                                                                                                                                                                                                                                     axis.text.y=element_blank(),
                                                                                                                                                                                                                                                     axis.ticks.y=element_blank(), panel.grid=element_blank())+ scale_fill_manual(values=P42)+guides(fill=guide_legend(title="Class"))+
  theme(legend.text=element_text(size=12),
        legend.title=element_text(size=12))

pie_bean_sw_dr_withPredr_planted_unplanted



ggsave(pie_bean_sw_dr_withPredr_planted_unplanted, file="piepanel_bean_sw_dr_withPredr_planted_unplanted_neworder.TIFF", unit=c("in"), height=20, width=30, dpi=300)



###############                    
##Bubble plot for actinos and alphaproteobacterias by OTU using max standardization
##############

library(tidyverse)

#read in OTU table

Active_DNA<- read.csv("active_final_DNAabund_phantoms_threshold_detection5%.csv", row.names=1)
View(Active_DNA)
#Active_new <- cbind(OTU = rownames(Active_DNA), Active_DNA)

active_long<-Active_DNA %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="Abundance", -ASV)#gather columns into key value pairs
View(active_long)

metadata<-read.csv("sample_metadata_dna_cdna.csv")

newdf = left_join(active_long, metadata, by=c("SampleID"))
View(newdf)
new_df_bean<- newdf %>% filter(crop=="bean" & drought=="drought" & planted=="planted")
View(new_df_bean)


##pivot wider
new_df_bean_wide <- new_df_bean %>% select(ASV, Abundance, SampleID) %>%pivot_wider(names_from="SampleID", values_from="Abundance")  
View(new_df_bean_wide) 

rnames <-new_df_bean_wide$ASV

new_df_bean_wide<- new_df_bean_wide %>% select (!ASV)
rownames(new_df_bean_wide) <- rnames

write.csv(new_df_bean_wide,"new_df_bean_wide.csv")

new_df_bean_wide<- read.csv("new_df_bean_wide.csv", row.names=1)
new_df_bean_wide = new_df_bean_wide[rowSums(new_df_bean_wide[])>0, ,drop=FALSE]
View(new_df_bean_wide)
#new_df_bean_wide<- as.data.frame(new_df_bean_wide)
#rownames(new_df_bean_wide) <- rnames


##subset OTUs to Actinos and alphas

##merge taxonomy with OTU hash id

taxonomy<- read.csv("taxonomy-dn-99-edited.csv", row.names=1)
View(taxonomy)
taxonomy_family<- taxonomy %>% select(c("Family", "Class"))
View(taxonomy_family)

merge_tax <- merge(new_df_bean_wide, taxonomy_family, by="row.names")
View(merge_tax)
merge_tax$OTUclass <- paste(merge_tax$Class, merge_tax$Row.names)
View(merge_tax)

merge_tax_edited<- merge_tax %>% filter(Class==" Alphaproteobacteria") ##Change to Actino as needed
View(merge_tax_edited)
rnames <-merge_tax_edited$OTUclass

rownames(merge_tax_edited) <- rnames  # assign row names
View(merge_tax_edited)
merge_tax_edited<- merge_tax_edited %>% select(!Row.names)%>% select (!Family) %>% select(!Class) %>% select(!OTUclass)
View(merge_tax_edited)


merge_tax_edited1 <- merge_tax_edited[order(rowSums(merge_tax_edited), decreasing = TRUE),]
merge_tax_20 <- merge_tax_edited1[c(1:20),] ##select top 50 most abundant Classes
View(merge_tax_20)

#max_standardization

matrix<- as.matrix(merge_tax_20)
matrix
View(matrix)
library(vegan)
merge_tax_maxstand <-decostand(matrix, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(merge_tax_maxstand)
library(pheatmap)
library(viridis)

merge_tax_maxstand <- as.data.frame(merge_tax_maxstand)
##pivot longer, merge with metadata

bean_long<- merge_tax_maxstand %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="RA", -ASV)#gather columns into key value pairs
View(bean_long)

metadata<-read.csv("sample_metadata_dna_cdna.csv")
View(metadata)
new_df_bean = left_join(bean_long, metadata, by=c("SampleID"))
View(new_df_bean)
##new_df_bean<- newdf %>% filter(crop=="bean")


#compute averages by harvest and summarise new
newdf_mean_bean<- new_df_bean %>% group_by(summarise, ASV) %>% summarise(mean=mean(RA))
View(newdf_mean_bean)

newdf_mean_bean$summarise <- recode(newdf_mean_bean$summarise, 'Day 2planted' = 'Day2',
                                    'Day 3planted' = 'Day3',
                                    'Day 4planted' = 'Day4',
                                    'Day 5planted' = 'Day5',
                                    'Day 6planted' = 'Day6'
) 

##plot

##dotplot code  

bb<-c(0.85, 0.70, 0.60,0.50,0.40,0.30,0.20,0.10,0.01) # define breaks. change for Alpha and Actino
ll<-c(0.85, 0.70, 0.60,0.50,0.40,0.30,0.20,0.10,0.01) # labels. change for Alpha and Actino
plot_dotplot_beandrought<-
  ggplot(newdf_mean_bean,aes(x=summarise,y=ASV,size=mean,color=mean))+
  geom_point()+
  #scale_x_discrete(limits = c("pre-drought","harvest1","harvest2", "harvest3","harvest4", "harvest5"))+
  #facet_grid(~planted,scales="free")+ 
  theme(strip.text.x=element_text(size=10))+
  #scale_colour_gradient(low="#132B43",
  #high="#56B1F7",limits=c(0.001,0.3),breaks=bb)+
  scale_colour_viridis(limits=c(0.01,0.85),breaks=bb)+
  scale_size_continuous(limits=c(0.01,0.85),labels=ll,breaks=bb, name="Mean Relative\nAbundance")+
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=9, angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Bean Drought (Alphaproteobacteria)")+ ##change title for alpha and actino
  theme(plot.title = element_text(hjust = 0.5))+
  labs(color="Mean Relative\nAbundance")+
  ylab("Class-OTU") + xlab("harvest")+
  theme(axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10)) + 
  theme(legend.title = element_text(size=10))+ theme(legend.text = element_text(size=9))





plot_dotplot_beandrought

##Supplementary Figure S17A and B
ggsave(filename = "dotplotOTU_Alpha_bean_dr_planted_newlabel_FINAL_1.tiff", plot = plot_dotplot_beandrought, ##change filename for actino
       width = 17 ,
       height = 16, units = c("cm"),
       dpi = 200)

##checked these plots so they include all abundance data points except 0 

                    






