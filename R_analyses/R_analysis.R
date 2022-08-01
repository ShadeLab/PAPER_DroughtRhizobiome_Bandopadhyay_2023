#######################################################################################
####Bacterial community Analysis for Drought Experiment using Bean and Switchgrass#####
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
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point() ##Supplemental Figure S2

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
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") ##Supplemental Figure S3
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
rarecurve<-rarecurve(t(otu_table(phyloseq_merged1)), step=20, label = FALSE, col="cyan4") ## Supplemental Figure S7

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


finaldf = newdf %>% filter(newdf$crop == "switchgrass") %>%  mutate(ratioMethod1= ifelse(dna == 0 & cdna > 0, 100, cdna/dna)) %>%  #replace "switchgrass" with "bean" when filtering to bean samples
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
##comes out as 6.3172672 % for percent_gr_1_nophantoms column- so calculation is normalized by sample now

##Switchgrass
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


#bean 
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


###boxplots for ratio threshold >1 , SW

#long format to combine method 1 and 2 together for facetting

sw_phantoms<- read.csv("SW_phantomstats_allfinite_greaterthan1_normalizedtotalotus_samplebysample_new.csv")
library(tidyverse)
sw_phantoms_longer<- sw_phantoms %>% pivot_longer(cols=starts_with("percent_gr"),
                                                  names_to = "method",
                                                  values_to="percent_active")
View(sw_phantoms_longer)
df_test<- sw_phantoms_longer%>% filter(treatment=="PreDrP1")
View(df_test)

library(ggplot2)
library(car)
packageDescription("car")


##Make a modified copy of the original data to change the facet names in method column

sw_phantoms_longer <- sw_phantoms_longer %>%
  mutate(method = recode(method,
                         "percent_gr_1_method1" = "method1 (threshold > 1)",
                         "percent_gr_1_method2" = "method2 (threshold > 1)",
                         "percent_gr_1_nophantoms" = "no_phantoms_included(threshold > 1)"
  ))

View(sw_phantoms_longer)

##boxplot of ratio threshold > 1, bean

bean_phantoms<- read.csv("Bean_phantomstats_allfinite_greaterthan1_normalizedtotalotus_samplebysample_new.csv")
library(tidyverse)
bean_phantoms_longer<- bean_phantoms %>% pivot_longer(cols=starts_with("percent_gr"),
                                                      names_to = "method",
                                                      values_to="percent_active")
View(bean_phantoms_longer)
df_test<- bean_phantoms_longer%>% filter(treatment=="PreDrP1")
View(df_test)

library(ggplot2)



##Make a modified copy of the original data to change the facet names in method column

bean_phantoms_longer <- bean_phantoms_longer %>%
  mutate(method = recode(method,
                         "percent_gr_1_method1" = "method1 (threshold > 1)",
                         "percent_gr_1_method2" = "method2 (threshold > 1)",
                         "percent_gr_1_nophantoms" = "no_phantoms_included(threshold > 1)"
  ))

View(bean_phantoms_longer)

##reorder levels of harvest
level_order <- c('predrought', 'harvest1', 'harvest2', 'harvest3', 'harvest4', 'harvest5')
#plot
plot1 <- ggplot(sw_phantoms_longer, aes(x=factor(harvest, level=level_order), y=percent_active,fill=drought)) + facet_grid(planted~method) +
  theme(strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent active OTUs") + ggtitle("Switchgrass")+ theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y=element_text(size=17),
                                                                                                             axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))  + scale_fill_manual(limits = c("predrought", "drought", "well-watered"), values=c("#D55E00","#009E73", "#56B4E9"))
plot1
ggsave(filename = "SW_percent_active_threshold>1_colorblindfriendly.tiff", plot = plot1,
       width = 24,
       height = 18, units = c("cm"),
       dpi = 300)

##reorder levels of harvest
level_order <- c('predrought', 'harvest1', 'harvest2', 'harvest3', 'harvest4', 'harvest5')
#plot
plot2 <- ggplot(bean_phantoms_longer, aes(x=factor(harvest, level=level_order), y=percent_active,fill=drought)) + facet_grid(planted~method) +
  theme(strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 12))+ geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent active OTUs") + ggtitle("Bean")+theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y=element_text(size=17),
                                                                                                     axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))  + scale_fill_manual(limits = c("predrought", "drought", "well-watered"), values=c("#D55E00","#009E73", "#56B4E9"))
plot2
ggsave(filename = "bean_percent_active_threshold>1_colorblindfriendly.tiff", plot = plot2,
       width = 24,
       height = 18, units = c("cm"),
       dpi = 300)



###boxplots for ratio threshold >= 1 , SW

#long format to combine method 1 and 2 together for facetting

sw_phantoms<- read.csv("SW_phantomstats_allfinite_greaterthan1_normalizedtotalotus_samplebysample_new.csv")
library(tidyverse)
sw_phantoms_longer<- sw_phantoms %>% pivot_longer(cols=starts_with("percent_gr"),
                                                  names_to = "method",
                                                  values_to="percent_active")
View(sw_phantoms_longer)
df_test<- sw_phantoms_longer%>% filter(treatment=="PreDrP1")
View(df_test)

library(ggplot2)
library(car)
packageDescription("car")


##Make a modified copy of the original data to change the facet names in method column

sw_phantoms_longer <- sw_phantoms_longer %>%
  mutate(method = recode(method,
                         "percent_gr_eq_1_method1" = "method1 (threshold >= 1)",
                         "percent_gr_eq_1_method2" = "method2 (threshold >= 1)",
                         "percent_gr_eq_1_nophantoms" = "no_phantoms_included(threshold >= 1)"
  ))

View(sw_phantoms_longer)

##boxplot of ratio threshold >= 1, bean

bean_phantoms<- read.csv("Bean_phantomstats_allfinite_greaterthan1_normalizedtotalotus_samplebysample_new.csv")
library(tidyverse)
bean_phantoms_longer<- bean_phantoms %>% pivot_longer(cols=starts_with("percent_gr"),
                                                      names_to = "method",
                                                      values_to="percent_active")
View(bean_phantoms_longer)
df_test<- bean_phantoms_longer%>% filter(treatment=="PreDrP1")
View(df_test)

library(ggplot2)



##Make a modified copy of the original data to change the facet names in method column

bean_phantoms_longer <- bean_phantoms_longer %>%
  mutate(method = recode(method,
                         "percent_gr_eq_1_method1" = "method1 (threshold >= 1)",
                         "percent_gr_eq_1_method2" = "method2 (threshold >= 1)",
                         "percent_gr_eq_1_nophantoms" = "no_phantoms_included(threshold >= 1)"
  ))

View(bean_phantoms_longer)

##reorder levels of harvest
level_order <- c('predrought', 'harvest1', 'harvest2', 'harvest3', 'harvest4', 'harvest5')
#plot
plot1 <- ggplot(sw_phantoms_longer, aes(x=factor(harvest, level=level_order), y=percent_active,fill=drought)) + facet_grid(planted~method) +
  theme(strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent active OTUs") + ggtitle("Switchgrass")+ theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y=element_text(size=17),
                                                                                                             axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))  + scale_fill_manual(limits = c("predrought", "drought", "well-watered"), values=c("#D55E00","#009E73", "#56B4E9"))
plot1
ggsave(filename = "SW_percent_active_threshold>=1_colorblindfriendly.tiff", plot = plot1,
       width = 24,
       height = 18, units = c("cm"),
       dpi = 300)

##reorder levels of harvest
level_order <- c('predrought', 'harvest1', 'harvest2', 'harvest3', 'harvest4', 'harvest5')
#plot
plot2 <- ggplot(bean_phantoms_longer, aes(x=factor(harvest, level=level_order), y=percent_active,fill=drought)) + facet_grid(planted~method) +
  theme(strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 12))+ geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent active OTUs") + ggtitle("Bean")+theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y=element_text(size=17),
                                                                                                     axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))  + scale_fill_manual(limits = c("predrought", "drought", "well-watered"), values=c("#D55E00","#009E73", "#56B4E9"))
plot2
ggsave(filename = "bean_percent_active_threshold>=1_colorblindfriendly.tiff", plot = plot2,
       width = 24,
       height = 18, units = c("cm"),
       dpi = 300)



##boxplots of percent phantoms and singleton phantoms

sw_phantoms<- read.csv("SW_phantomstats_allfinite_greaterthaneq1_normalizedtotalotus_samplebysample_new.csv")

level_order <- c('predrought', 'harvest1', 'harvest2', 'harvest3', 'harvest4', 'harvest5')
#plot
plot3 <- ggplot(sw_phantoms, aes(x=factor(harvest, level=level_order), y=percent_phantoms,fill=drought)) + facet_grid(~planted) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent phantom OTUs") + ggtitle("Switchgrass Percent Phantom OTUs")+ theme(plot.title = element_text(hjust = 0.5))+theme(axis.text.y=element_text(size=10),
                                                                                                                                  axis.title.y=element_text(size=15)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=13))+ scale_fill_manual(limits = c("predrought", "drought", "well-watered"), values=c("#D55E00","#009E73", "#56B4E9"))

plot3 
ggsave(filename = "percent_phantoms_SW_colorblindfriendly.tiff", plot = plot3,
       width = 18,
       height = 10, units = c("cm"),
       dpi = 300)

plot4 <- ggplot(sw_phantoms, aes(x=factor(harvest, level=level_order), y=percent_phantom_singleton,fill=drought)) + facet_grid(~planted) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent phantom singleton OTUs") + ggtitle("Switchgrass Percent Phantom Singleton OTUs")+ theme(plot.title = element_text(hjust = 0.5))+theme(axis.text.y=element_text(size=10),
                                                                                                                                                      axis.title.y=element_text(size=15)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=13))+ scale_fill_manual(limits = c("predrought", "drought", "well-watered"), values=c("#D55E00","#009E73", "#56B4E9"))

plot4
ggsave(filename = "percent_phantom_singleton_SW_colorblindfriendly.tiff", plot = plot4,
       width = 18,
       height = 10, units = c("cm"),
       dpi = 300)



bean_phantoms<- read.csv("Bean_phantomstats_allfinite_greaterthaneq1_normalizedtotalotus_samplebysample_new.csv")

level_order <- c('predrought', 'harvest1', 'harvest2', 'harvest3', 'harvest4', 'harvest5')
#plot
plot5 <- ggplot(bean_phantoms, aes(x=factor(harvest, level=level_order), y=percent_phantoms,fill=drought)) + facet_grid(~planted) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent phantom OTUs") + ggtitle("Bean Percent Phantom OTUs")+ theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.y=element_text(size=10),
                                                                                                                            axis.title.y=element_text(size=15)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=13)) + scale_fill_manual(limits = c("predrought", "drought", "well-watered"), values=c("#D55E00","#009E73", "#56B4E9"))

plot5 
ggsave(filename = "percent_phantoms_bean_colorblindfriendly.tiff", plot = plot5,
       width = 18,
       height = 10, units = c("cm"),
       dpi = 300)

plot6 <- ggplot(bean_phantoms, aes(x=factor(harvest, level=level_order), y=percent_phantom_singleton,fill=drought)) + facet_grid(~planted) +
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+geom_boxplot() + ylim(0,80)+ theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=10, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Percent phantom singleton OTUs") + ggtitle("Bean Percent Phantom Singleton OTUs")+ theme(plot.title = element_text(hjust = 0.5))+theme(axis.text.y=element_text(size=10),
                                                                                                                                               axis.title.y=element_text(size=15)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=13)) + scale_fill_manual(limits = c("predrought", "drought", "well-watered"), values=c("#D55E00","#009E73", "#56B4E9"))

plot6
ggsave(filename = "percent_phantom_singleton_bean_colorblindfriendly.tiff", plot = plot6,
       width = 18,
       height = 10, units = c("cm"),
       dpi = 300)



##for arranging in grid for publication

##Supplemental Figure 5
library(ggpubr)
percent_active_bythreshold_gr_eq_1<- ggarrange(plot2, plot1, ncol=1, nrow=2, widths=c(1,1),heights=c(1,1), common.legend = TRUE, legend="bottom", labels= c("A", "B"))
View(percent_active_bythreshold)
ggsave("combined_bean_sw_percentactive>=1_bythreshold.tiff",percent_active_bythreshold,width=22,height=26, units="cm")

##Supplemental Figure 6
percent_active_bythreshold_gr_1<- ggarrange(plot2, plot1, ncol=1, nrow=2, widths=c(1,1),heights=c(1,1), common.legend = TRUE, legend="bottom", labels= c("A", "B"))
View(percent_active_bythreshold)
ggsave("combined_bean_sw_percentactive>1_bythreshold.tiff",percent_active_bythreshold,width=22,height=26, units="cm")

#Supplemental Figure 4
percent_phantoms_combined<- ggarrange(plot5, plot6, plot3, plot4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom", labels= c("A", "B", "C", "D"))
View(percent_phantoms_combined)
ggsave("combined_bean_sw_percentphantoms.tiff",percent_phantoms_combined,width=27,height=27, units="cm")






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
rarecurve(t(otu_table(phyloseq_merged)), step=20, label = FALSE, col="cyan4") ##Supplemental Figure S10

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


plot<- ggplot(phyloseq_class, aes(x= treatment, y = Abundance, fill =Class)) + 
  facet_grid(crop~planted+drought_new, scales="free") + theme(strip.text.x = element_text(size = 25))+ 
  theme(strip.text.y = element_text(size = 25))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new41))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=17),
                                     axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))
plot 

 
### Supplemental Figure 12
ggsave(filename = "barplot_class_decontam_bean+sw+dr+planted_with_singleton_nosecondrarefaction.tiff", plot = plot,
       width = 70,
       height = 30, units = c("cm"),
       dpi = 300)


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
  physeq = phyloseq_merged_final,
  ordination = phyloseq_pcoa,
  color = "crop",
  shape = "crop",
  title = "PCoA Switchgrass and Bean") + 
  scale_color_manual(values = c("#CC79A7", "#0072B2"))+
  geom_point(aes(color = crop), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))
library(ggplot2)
pcoa_crop 


### For main Figure 4 A
ggsave(filename="pcoa-sw-bean-new-with-singleton-nosecondrarefaction-final-colorblindfriendly.TIFF", plot=plot2, width=6.8, height=6, units="in", dpi=720)

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

#plot
pcoa_bean<-plot_ordination(
  physeq = phyloseq_merged_final,
  ordination = phyloseq_pcoa,
  color = "drought",
  shape = "planted",
  title = "PCoA Bean"
) +
  scale_color_manual(name="drought", limits = c("pre-drought", "drought", "well-watered"), values = c("#D55E00","#009E73","#56B4E9")) +
  geom_point(aes(color = drought), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))

pcoa_bean

### For main Figure 4 B
ggsave(filename="pcoa-bean-new-with-singleton-nosecondrarefaction-final-colorblindfriendly.TIFF", plot=plot2, width=6.8, height=6, units="in", dpi=720)



# Plot by treatment (drought, planted levels) Switchgrass

#pcoa

phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged_final, 
  method = "PCoA", 
  distance = "bray"
)
#plot
pcoa_sw<-plot_ordination(
  physeq = phyloseq_merged_final,
  ordination = phyloseq_pcoa,
  color = "drought",
  shape = "planted",
  title = "PCoA Switchgrass"
) + 
  scale_color_manual(name="drought", limits = c("pre-drought", "drought", "well-watered"), values = c("#D55E00","#009E73","#56B4E9")) +
  geom_point(aes(color = drought), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))

pcoa_sw

### For main Figure 4 C
ggsave(filename="pcoa-sw-new-with-singleton-nosecondrarefaction-final-colorblindfriendly.TIFF", plot=plot2, width=6.8, height=6, units="in", dpi=720)



#############################################
###Abundance dotplots over time##############
#############################################

metadata<- read.csv("sw-drought-planted-unplanted.csv")
metadata<- read.csv("bean-drought-planted-unplanted.csv") ##change metadata file as needed between bean and switchgrass
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
phyloseq_class1<- phyloseq_class %>% group_by(summarise, Class) %>% mutate(mean=mean(Abundance))

write.csv(phyloseq_class1, "bean_phyloseq_class1.csv")

write.csv(phyloseq_class1, "sw_phyloseq_class1.csv")

##combine and read back in to plot together or modify the "group_by" function when making phyloseq_class1 object(see code above when making phyloseq_class1 object for dotplot)

##final dotplot code
library(viridis)
plot_dotplot_swdrought<- ggplot(phyloseq_class1, aes(x= harvest, y = Class, size=mean, color=mean)) + 
  facet_grid(~planted, scales="free") + theme(strip.text.x = element_text(size = 10))+ 
  geom_point()+ 
  scale_x_discrete(limits = c("pre-drought","harvest1","harvest2", "harvest3","harvest4", "harvest5"))+
  scale_colour_viridis(limits=c(0.00, 0.35),  #limits=c(0.00, 0.25)
                       oob = scales::squish, direction=1)+ scale_size(range= c(1, 6), limits= c(0.00, 0.35),name="Mean Relative Abundance") +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=9, angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Switchgrass Drought")+ # change to "Bean Drought" when doing bean samples
  theme(plot.title = element_text(hjust = 0.5))+
  labs(color="Mean Relative Abundance")+
  ylab("Class") + 
  theme(axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10)) + 
  theme(legend.title = element_text(size=10))+ theme(legend.text = element_text(size=9))

plot_dotplot_beandrought ##For Figure 4 D
plot_dotplot_swdrought ## For Figure 4 E

ggsave(filename = "dotplotclass_sw_dr_with_singletons_nosecondrarefaction_FINAL.tiff", plot = plot_dotplot_swdrought,
       width = 17 ,
       height = 16, units = c("cm"),
       dpi = 200)

ggsave(filename = "dotplotclass_bean_dr_with_singletons_nosecondrarefaction_FINAL.tiff", plot = plot_dotplot_beandrought,
       width = 17 ,
       height = 16, units = c("cm"),
       dpi = 200)

permanova_stats <- read.csv("crop_permanova_stats.csv", header = TRUE) ##table for Figure 4 F
table1<- ggtexttable(permanova_stats, rows = NULL, theme = ttheme(base_size = 16))
table1


## Main manuscript figure 4

library(ggpubr)
my_plot_list <- list(pcoa_crop, pcoa_bean, pcoa_sw, plot_dotplot_beandrought, plot_dotplot_swdrought, table1)
plot_combined <- ggarrange(plotlist = my_plot_list, labels=c("A", "B", "C", "D", "E", "F"),
                           common.legend = FALSE, heights = c(10,15)) + bgcolor("white") 

plot_combined


## Combined Figure 4
ggsave(filename = "bean_sw_pcoa_dotplots_new_colorblindfriendly_final.tiff", plot = plot_combined,
       width = 44,
       height = 25, units = c("cm"),
       dpi = 300)


#######################################
##### Statistical anslyses ############
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



#######bean  and switchgrass samples combined mantel correlations###########
  
  
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


##plot correlation

##Supplemental Figure S11

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

## For Figure 5 A, B and Supplemental Figure 13.
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

sw_dr_BC_time<- read.csv("sw-dr-BC-time-planted-unplanted-merge-edited-withsingletons_nosecondrarefaction.csv")
plot_sw_dr_BC<- ggplot(sw_dr_BC_time, aes(x= group,y=similarity, color=planted_x)) + theme_classic()+ scale_color_manual(values=c("#39568CFF", "#55C667FF")) +
  geom_point(aes(), size=2) + ylim(0.1, 0.6)+geom_smooth(method = lm, se=FALSE) + ggtitle("Switchgrass Drought") + theme(plot.title = element_text(size= 18, hjust = 0.5))+
  labs(y = "Bray-Curtis similarity to pre-drought", x = "Harvest time point (days)", color='planted') + theme(axis.text.x = element_text(face = "bold",colour = "black", size = 17,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(face = "bold", size = 17, colour = "black"), 
         axis.title= element_text(face = "bold", size = 18, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =14, face = "bold", colour = "black"),
         legend.text = element_text(size = 12, face = "bold", colour = "black")) 
plot_bean_ww_BC ##plot for bean well-watered condition
plot_sw_ww_BC ## plot for sw well-watered condition 
plot_sw_dr_BC ## plot for sw drought condition
plot_bean_dr_BC ## plot for bean drought condition

##Combine for Figure 5 and Supplemental Figure 13
library(ggpubr)

plot_combined_dr<- ggarrange(plot_bean_dr_BC, plot_sw_dr_BC,labels = c("A", "B"))

plot_combined_ww<- ggarrange(plot_bean_ww_BC, plot_sw_ww_BC,labels = c("A", "B"))


### Figure 5 A, B

ggsave(filename = "bean_sw_drought_BCsim_FINAL.tiff", plot = plot_combined_dr,
       width = 30,
       height = 15, units = c("cm"),
       dpi = 300)


##Supplemental Figure 13 A, B

ggsave(filename = "bean_sw_ww_BCsim_FINAL.tiff", plot = plot_combined_ww,
       width = 30,
       height = 15, units = c("cm"),
       dpi = 300)






############
##ANCOVA for BC similarity versud time, catergorical variables- planted, covariate: time, categorical variable 2: crop, cat var 3: drought, dependent variable similarity

ancova<- read.csv("ancova_metadata_withsingletons_nosecondrarefaction.csv")
mod1 <- aov(similarity~crop*planted*drought*days, data=ancova) ##test interaction
summary(mod1)



##if crop level is included minor interactions exist between planted:drought and crop:planted 

##separate bean and SW

ancova_bean <-subset(ancova, crop=="bean")
mod1 <- aov(similarity~planted*drought*days, data=ancova_bean) ##test interaction
summary(mod1)



##no interaction for bean , slope is not different between treatments
#no interaction between time and planted , so no difference in slope



mod2 <- aov(similarity~planted+drought+days, data=ancova_bean) ##if no interaction seen above, test differences in intercept within categorical variable eg. planted, if significant then subset by planted treatment and look for effects of covariate (harvest time)

summary(mod2)


## no difference in intercept: planted levels


ancova_switchgrass<- subset(ancova, crop=="switchgrass")
mod1 <- aov(similarity~planted*drought*days, data=ancova_switchgrass) ##test interaction
summary(mod1)



###no interaction for SW, slope is not different between treatments for SW
###no interaction between time and planted , so no difference in slope
  
  
mod2 <- aov(similarity~planted+drought+days, data=ancova_switchgrass) ##if no interaction seen above, test differences in intercept within categorical variable eg. planted, if significant then subset by planted treatment and look for effects of covariate (harvest time)

summary(mod2)


##significant differences in intercept for SW data **** : planted factor 

anova(mod1,mod2) ##does removing the interaction affect the fit of the model? if not then proceed with model 2 and subset the data based on the two planted levels: planted and unplanted
##separate regression for categorical variable: planted factor, testing effects of time on y var, due to effects of planted factor on intercept 
## separate planted and unplanted datasets and and assess effect of time on dependent variable (BC similarity)

  
##removing interaction term for SW doesn't significantly alter the fit of the model.
##therefore we can conclude that most parsimonious model is mod2.
  
  
##SW subset to drought and well watered
  
ww_sw <- subset(ancova_switchgrass, drought=="well-watered" )
dr_sw <- subset(ancova_switchgrass, drought=="drought" )

reg1 <- lm(similarity~planted_y/days -1, data=ww_sw); summary(reg1) ##subset data, nest days within planted
reg2 <- lm(similarity~planted_y/days -1, data=dr_sw); summary(reg2) ##subset data, nest days within planted


#### planted treatment in switchgrass is significantly affected by the harvest time for well watered condition
## neither planted nor unplanted treatment is significantly affected by the harvest time for drought condition


#################################
######## DESeq analysis #########


#####Bean Drought DeSeq analysis 

metadata<- read.csv("bean-dr-nopredr.csv") 
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) 
phyloseq_merged_final

head(sample_data(phyloseq_merged_final)$planted, 25)

if (!requireNamespace("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
packageVersion("DESeq2")

##convert phyloseq format to DESEq dataset with dispersions estimated using the experimental design formula
diagdds = phyloseq_to_deseq2(phyloseq_merged_final, ~ planted)
diagdds$planted
diagdds

diagdds$planted <-factor(diagdds$planted, levels = c("unplanted", "planted"))
diagdds$planted #make sure that Control is the first level in the treatment factor, so that the
#default log2 fold changes are calculated as treatment over control and not the other way around.

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

##The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

#The following results function call creates a table of the results of the tests. Very fast. 
#The hard work was already stored with the rest of the DESeq2-related data in our latest version 
#of the diagdds object (see above). I then order by the adjusted p-value, removing the entries 
#with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display.


res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_merged_final)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

sigtab_positives<- sigtab %>% filter(!log2FoldChange<0)
View(sigtab_positives) 


#Let's look at the OTUs that were positively enriched in the planted bean soil samples compared to the unplanted ones. The following makes a nice ggplot2 summary of the results.

library("ggplot2")
install.packages("Polychrome")
library(Polychrome)
P52<- createPalette(52,  c("#010101", "#ff0000"))
names(P52) <- NULL
library(viridis)


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "P52", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab_positives$log2FoldChange, sigtab_positives$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_positives$Phylum = factor(as.character(sigtab_positives$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_positives$log2FoldChange, sigtab_positives$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_positives$Genus = factor(as.character(sigtab_positives$Genus), levels=names(x))

View(sigtab_positives) 


## Figure 5 E - Bean drought positively enriched OTUs in planted treatments compared to unplanted

plot_bean_dr<- ggplot(sigtab_positives, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  #scale_x_reordered() +
  geom_point(size=3) + scale_color_manual(values=as.vector(P52))+ 
  ggtitle("Bean Drought")+
  theme(plot.title = element_text(hjust = 0.5, size=18))+
  theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_text(size=18))+
  theme(axis.text.y=element_text(size=13), axis.title.y=element_text(size=18)) + 
  theme(legend.title = element_text(size=18))+ theme(legend.text = element_text(size=13)) +
  guides(col=guide_legend(ncol=1)) +theme(strip.text.x = element_text(size = 30)) +geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  ylab("Log2FoldChange")+coord_flip() 

ggsave(filename = "Deseq2-bean-dr-nopredr-unplanted-planted-newordering-withsingletons-nosecondrarefaction.tiff", plot = plot_bean_dr,
       width = 90,
       height = 70, units = c("cm"),
       dpi = 300)


##### Switchgrass Drought DeSeq analysis 

metadata<- read.csv("sw-dr-nopredr.csv") 
keep.samples <- as.vector(metadata$SampleID)
keep.samples

phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) 
phyloseq_merged_final

head(sample_data(phyloseq_merged_final)$planted, 25)

##convert phyloseq format to DESEq dataset with dispersions estimated using the experimental design formula
diagdds = phyloseq_to_deseq2(phyloseq_merged_final, ~ planted)
diagdds$planted
diagdds

diagdds$planted <-factor(diagdds$planted, levels = c("unplanted", "planted"))
diagdds$planted #make sure that Control is the first level in the treatment factor, so that the
#default log2 fold changes are calculated as treatment over control and not the other way around.

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

##The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

#The following results function call creates a table of the results of the tests. Very fast. 
#The hard work was already stored with the rest of the DESeq2-related data in our latest version 
#of the diagdds object (see above). I then order by the adjusted p-value, removing the entries 
#with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display.


res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_merged_final)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

sigtab_positives<- sigtab %>% filter(!log2FoldChange<0)
View(sigtab_positives) 



#Let's look at the OTUs that were positively enriched in the planted switchgrass soil samples compared to the unplanted ones. The following makes a nice ggplot2 summary of the results.

library("ggplot2")
install.packages("Polychrome")
library(Polychrome)
P52<- createPalette(52,  c("#010101", "#ff0000"))
names(P52) <- NULL
library(viridis)


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "P52", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab_positives$log2FoldChange, sigtab_positives$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_positives$Phylum = factor(as.character(sigtab_positives$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab_positives$log2FoldChange, sigtab_positives$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_positives$Genus = factor(as.character(sigtab_positives$Genus), levels=names(x))

View(sigtab_positives) 


## Figure 5 F - Switchgrass drought positively enriched OTUs in planted treatments compared to unplanted

plot_sw_dr<- ggplot(sigtab_positives, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  #scale_x_reordered() +
  geom_point(size=3) + scale_color_manual(values=as.vector(P52))+ 
  ggtitle("Switchgrass Drought")+
  theme(plot.title = element_text(hjust = 0.5, size=18))+
  theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_text(size=18))+
  theme(axis.text.y=element_text(size=13), axis.title.y=element_text(size=18)) + 
  theme(legend.title = element_text(size=18))+ theme(legend.text = element_text(size=13)) +
  guides(col=guide_legend(ncol=1)) +theme(strip.text.x = element_text(size = 30)) +geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  ylab("Log2FoldChange")+coord_flip() 

ggsave(filename = "Deseq2-sw-dr-nopredr-unplanted-planted-newordering-withsingletons-nosecondrarefaction.tiff", plot = plot_sw_dr,
       width = 90,
       height = 70, units = c("cm"),
       dpi = 300)


plot_sw_dr
plot_bean_dr


## Figure 5 E, F (edited in inkscape)
plot_combined_deseq_positives<- ggarrange(plot_bean_dr, plot_sw_dr,labels = c("E", "F"))
ggsave(filename = "bean_sw_dr_deseqPositivies_FINAL.tiff", plot = plot_combined_deseq_positives,
       width = 45,
       height = 15, units = c("cm"),
       dpi = 300)


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

###Plot Inverse Simpson bean
alphadiv_bean <- read.csv("alphadiv_bean_withsingleton_nosecondrarefaction.csv")
pd<-position_dodge(0.7)
alphadiv_bean_InvS =  alphadiv_bean %>% filter(measure=="Inverse Simpson")
alphadiv_bean_InvS$drought <- factor(alphadiv_bean_InvS$drought,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)
View(alphadiv_bean_InvS)

IS_bean<-ggplot(alphadiv_bean_InvS, aes(y=mean, x= harvest, fill=drought)) + #color = Treatment, group = Treatment
  facet_grid(~planted+drought, scales="free", space="free") + geom_boxplot()+
  ylab("Inverse Simpson") +
  ggtitle("Inverse Simpson - Bean")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("#D55E00","#009E73", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=12)) +
  theme(axis.text=element_text(size=12))+
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=12)) + theme_classic()
IS_bean

ggsave(filename="IS_bean_300dpi_newcolor_withsingleton_nosecondrarefaction_colorblindfriendly.TIFF", plot=IS_bean, width=15, height=6.8, units="in", dpi=300)

####Plot Inverse Simpson Switchgrass

alphadiv_sw <- read.csv("alphadiv_switchgrass_withsingleton_nosecondrarefaction.csv")
pd<-position_dodge(0.7)
alphadiv_sw_InvS =  alphadiv_sw %>% filter(measure=="Inverse Simpson")
alphadiv_sw_InvS$drought <- factor(alphadiv_sw_InvS$drought,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)
View(alphadiv_sw_InvS)

cbPalette <- c("#D55E00","#009E73", "#56B4E9")

IS_sw<-ggplot(alphadiv_sw_InvS, aes(y=mean, x= harvest, fill=drought)) + #color = Treatment, group = Treatment
  facet_grid(~planted+drought, scales="free", space="free") + geom_boxplot()+
  ylab("Inverse Simpson") +
  ggtitle("Inverse Simpson - Switchgrass")+ 
  theme(plot.title = element_text(hjust = 0.5))+
scale_fill_manual(values=c("#D55E00","#009E73", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=12)) +
  theme(axis.text=element_text(size=12))+
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=12)) + theme_classic()

IS_sw
ggsave(filename="IS_sw_300dpi_newcolor_withsingleton_nosecondrarefaction_colorblindfriendly.TIFF", plot=IS_sw, width=15, height=6.8, units="in", dpi=300)

####Bean Richness
alphadiv_bean <- read.csv("alphadiv_bean_withsingleton_nosecondrarefaction.csv")
pd<-position_dodge(0.7)
alphadiv_bean_rich =  alphadiv_bean %>% filter(measure=="Richness")
alphadiv_bean_rich$drought <- factor(alphadiv_bean_rich$drought,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)
View(alphadiv_bean_rich)


rich_bean<-ggplot(alphadiv_bean_rich, aes(y=mean, x= harvest, fill=drought)) +
  facet_grid(~planted+drought, scales="free", space="free") + geom_boxplot()+
  ylab("Richness(observed OTUs)") +
  ggtitle("Richness - Bean")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("#D55E00","#009E73", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=12)) +
  theme(axis.text=element_text(size=12))+
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=12)) +theme_classic()
rich_bean

ggsave(filename="rich_bean_300dpi_newcolor_withsingleton_nosecondrarefaction_colorblindfriendly.TIFF", plot=rich_bean, width=15, height=6.8, units="in", dpi=300)

###Switchgrass richness
alphadiv_sw <- read.csv("alphadiv_switchgrass_withsingleton_nosecondrarefaction.csv")
pd<-position_dodge(0.7)
alphadiv_sw_rich =  alphadiv_sw %>% filter(measure=="Richness")
alphadiv_sw_rich$drought <- factor(alphadiv_sw_rich$drought,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)
View(alphadiv_sw_rich)


rich_sw<-ggplot(alphadiv_sw_rich, aes(y=mean, x= harvest, fill=drought)) +
  facet_grid(~planted+drought, scales="free", space="free") + geom_boxplot()+
  ylab("Richness(observed OTUs)") +
  ggtitle("Richness - Switchgrass")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("#D55E00","#009E73", "#56B4E9"))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=12)) +
  theme(axis.text=element_text(size=12))+
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=12)) +theme_classic()

rich_sw
ggsave(filename="rich_sw_300dpi_newcolor_withsingleton_nosecondrarefaction_colorblindfriendly.TIFF", plot=rich_sw, width=15, height=6.8, units="in", dpi=300)

## Figure 6 main manuscript

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

alphadiv_combined_plot<-ggplot(alphadiv_combined, aes(y=mean, x= harvest, fill=drought)) +
  facet_grid(planted~Crop+measure, scales="free", space="free") + geom_boxplot()+
scale_fill_manual(values=c("#D55E00","#009E73", "#56B4E9"))+ theme_classic()+
  ylab("Alpha Diversity")+
  xlab("")+
  ggtitle("Alpha Diversity Estimates")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=12)) +
  theme(axis.text=element_text(size=12))+
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=12)) 

alphadiv_combined_plot

ggsave(filename="alphadiv_combined_withsingleton_nosecondrarefaction_colorblindfriendly.TIFF", plot=alphadiv_combined_plot, width=10, height=6.8, units="in", dpi=300)



####Alpha diversity statistics

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
RS_plot <- ggplot(df_all, aes(x=group_1, y=ActiveResistance))+ facet_grid(~crop) +
  geom_boxplot()+ ylim(0,1)+
  geom_point()+
  ylab(label = "Resistance")+
  xlab(label = NULL)+
  theme_bw()+
  scale_x_discrete(labels=c("Harvest 5"))+
  theme(axis.text.x = element_text(size=15))
RS_plot 

### Supplemental Figure 14

ggsave(filename = "resistance.tiff", plot = RS_plot,
       width = 8,
       height = 8, units = c("cm"),
       dpi = 300)



## Pie charts for Figure 5 


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


### Figure 5 C and D, edited in Inkscape for combined final plot 

ggsave(pie_bean_sw_dr_withPredr_planted_unplanted, file="piepanel_bean_sw_dr_withPredr_planted_unplanted_neworder.TIFF", unit=c("in"), height=20, width=30, dpi=300)


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
sw_data1$Treatment <- factor(sw_data1$Treatment,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)
library(ggplot2)
sw_data1$Treatment

bean_data1 <- summarySE(bean_shootmass, measurevar="ShootMass", groupvars=c("Harvest","Treatment"))
View(bean_data1)
#reshape data so Pre-treatment comes first in legend
bean_data1$Treatment <- factor(bean_data1$Treatment,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)
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
sw_data2$Treatment <- factor(sw_data2$Treatment,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)

sw_data2$Treatment
View(sw_data2)


bean_data2 <- summarySE(bean_soilmoisture, measurevar="percentmoisture", groupvars=c("Harvest","Treatment","Type"),na.rm=TRUE)

#reshape data so Pre-treatment comes first in legend
bean_data2$Treatment <- factor(bean_data2$Treatment,levels = c("pre-drought","drought","well-watered"), ordered=TRUE)

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
  theme(panel.background = element_blank())+ scale_x_discrete(limits = c("pre-drought","harvest1","harvest2","harvest3","harvest4", "harvest5"))+
  scale_fill_manual(values=c("#D55E00","#009E73","#56B4E9")) + theme(axis.title.x = element_text(size=15)) + theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) + ylab ("")+ xlab("") + theme(axis.text.y=element_text(size=12)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=12))
p1



p2<-ggplot(data_final_biomass,aes(Harvest,value,fill=Treatment)) + geom_errorbar(limits, position = dodge, width = 0.25)+
  geom_bar(aes(fill=Treatment),position="dodge",stat="identity")+ 
  theme_classic()+
  facet_grid( crop~ measure+Type, scales="free_y" ) +
  theme(panel.background = element_blank())+ scale_x_discrete(limits = c("pre-drought","harvest1","harvest2","harvest3","harvest4", "harvest5"))+
  scale_fill_manual(values=c("#D55E00","#009E73","#56B4E9")) + theme(axis.title.x = element_text(size=15)) + theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, ncol=1)) + ylab ("")+ xlab("") + theme(axis.text.y=element_text(size=12)) + 
  theme(legend.title = element_text(size=15))+ theme(legend.text = element_text(size=12))
p2




library(ggpubr)

plot_combined<- ggarrange(p1, p2,  labels = c("A", "B"),
                          common.legend = TRUE, legend = "bottom")

plot_combined


########## Main Figure 2#######

ggsave(filename = "bean_sw_moisture_biomass_new_colorblindfriendly_final.tiff", plot = plot_combined,
       width = 22,
       height = 13, units = c("cm"),
       dpi = 300)


#####################################
########Activity dynamics############
#####################################


library(tidyverse)

##using rarefied DNA and cDNA copies with 15k reads.
dna<- read.csv("DNAcopy.csv", header=TRUE, row.names = 1)
cdna <- read.csv("cDNAcopy.csv", header=TRUE, row.names = 1)

Active_DNA<- read.csv("active_final_DNAabund_phantoms_threshold_detection5%.csv", row.names=1)
View(Active_DNA)
library(gplots) ##for heatmap2 function
Abundant_Active <- Active_DNA[order(rowSums(Active_DNA), decreasing = TRUE),]
Abundant_Active <- Abundant_Active[c(1:50),] ##select top 50 most abundant Classes
View(Abundant_Active)


View(dna)
View(cdna)

dna_50taxa <- dna[rownames(dna)%in%rownames(Abundant_Active),]
View(dna_50taxa)

cdna_50taxa <- cdna[rownames(cdna)%in%rownames(Abundant_Active),]
View(cdna_50taxa)

dna1<- dna_50taxa %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="dna", -ASV)#gather columns into key value pairs
View(dna1)
cdna1<- cdna_50taxa %>%
  tibble::rownames_to_column(var="ASV") %>%
  tidyr::gather(key="SampleID", value="cdna", -ASV)
View(cdna1)
merge.df = full_join(dna1, cdna1, by=c("ASV", "SampleID"))
View(merge.df)## 50*212= 10,600 entries

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

##code taxa as inactive (0) or not detected in dna (NA)

finaldf_coded <- finaldf %>% #mutate(dna_code=ifelse(dna2==0, "not_detected", dna2)) 
  mutate(activity= ifelse(ratioMethod2 < 1, "inactive", "active"))
View(finaldf_coded) 

##Abundant_active dataset abundance combined with finaldf dataset 
Abundant_Active_new <- cbind(ASV = rownames(Abundant_Active), Abundant_Active)
View(Abundant_Active_new)
Abundant_Active_long <-pivot_longer(Abundant_Active_new, !ASV, names_to="SampleID", values_to="AbundanceDNA_Active")
View(Abundant_Active_long)
merge_df_new<- merge(Abundant_Active_long, finaldf_coded, by=c("ASV", "SampleID"))
View(merge_df_new)
merge_df_new_NA<- merge_df_new %>% mutate(AbundanceDNA_Active_NA= ifelse(dna==0 & cdna==0, NA, AbundanceDNA_Active))
View(merge_df_new_NA)


##heatmap coding when including phantoms, becomes 0 for inactive and the usual DNA abundance for active taxa is used for everything else thats active

##when phantoms are accounted , there is no undetected taxa, so all i need is to code inactive as 0, and the rest should just be the active abundance filtered to DNA data



##plot heatmap for drought samples in switchgrass only
drought_switchgrass_planted<- merge_df_new_NA %>% filter (crop=="switchgrass", drought == "drought", planted=="planted") ##select data as needed for drought, well-watered, planted or unplanted samples
View(drought_switchgrass_planted)

drought_sw_planted_heat <- drought_switchgrass_planted %>% select(ASV, AbundanceDNA_Active_NA, treatment) %>%pivot_wider(names_from="treatment", values_from="AbundanceDNA_Active_NA")  
View(drought_sw_planted_heat) 
data<- as.data.frame(drought_sw_planted_heat)


##try this option 1
rnames <-data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names
View(mat_data)

#OR 
#option 2, both options work
rownames(data) <- rnames  # assign row names
View(data)
mat_data<- data %>% select(!ASV)
View(mat_data)

mat_data<- data.matrix(mat_data)

sum(is.na(as.matrix(dist(mat_data))))

giveNAs = which(is.na(as.matrix(dist(mat_data))),arr.ind=TRUE)
head(giveNAs)
mat_data[c(1,17),]

tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
  sum(is.na(as.matrix(dist(mat_data[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))]
rmv
#[1] "31" "72" "74" "84" "23" "64" (for 100 taxa)

#[1] "17" "45" (for 50 top taxa) ##45 seems to NOT be rowsums==0, then why is this being removed? ## the intention here is not to remove 
##rows but to see for which ASVs are we not able to calculate a distance for the heatmap


mat_data = mat_data[-as.numeric(rmv),]
View(mat_data)

##merge taxonomy with OTU hash id

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

if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
install_github("raivokolde/pheatmap")
library(pheatmap)





##max standardization, margin should be 1 not 2, 1 means calculation across rows, 2 means by column.
##since ASVs are rows for us, we should calculate max standardization by row
library(vegan)

##replace NAs with zero
##then compute max standardization
## put back NAs where they were before
##compute heatmap

##replace NAs with zero
matrix<- as.matrix(mat_data_merge_otuclass)
matrix
View(matrix)
matrix[is.na(matrix)] <-0
View(matrix)

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

##plot heatmap
##use when needing png image

png("heatmap_active_with_phantoms_sw_drought_planted_viridis_with_title.png",    # create PNG for the heat map        
    width = 8*300,        # 5 x 300 pixels
    height = 7*300,
    res = 300,            # 300 pixels per inch
    pointsize = 100)        # smaller font size


if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
install_github("raivokolde/pheatmap")

library(viridis)
#pheatmap(matrix_max_standardize, color=scale_fill_viridis_d(option="viridis"), cluster_cols=FALSE) ## including all top 100 taxa (or less due to removed rows, see code in lines 250-268), coz subsetting is not working with NA in the matrix, all NAs are coded as grey only, color change not working
#ggsave(filename="heatmap_class_sw_drought_planted.TIFF", plot= plot, width=15, height=8, unit=c("cm"), dpi=300)


##dev.off() ##use when needing png image 
pheatmap::pheatmap(matrix_max_standardize_sw,
                   cluster_cols = F,
                   na_col = "white",
                   border_color = "white",
                   main="Bean drought",
                   color = viridis(n = 256, alpha = 1, 
                                   begin = 0, end = 1, option = "viridis"
                   ))

dev.off()


##use for combining plots (drought planted samples)
plot_heatmap_bean <- pheatmap::pheatmap(matrix_max_standardize_bean, ##DO THE ABOVE ANALYSIS BY SUBSETTING TO BEAN DROUGHT PLANTED SAMPLES THEN COMBINE AS BELOW
                                        cluster_cols = F,
                                        na_col = "white",
                                        border_color = "white",
                                        main="Bean drought",
                                        color = viridis(n = 256, alpha = 1, 
                                                        begin = 0, end = 1, option = "viridis"
                                        ))
plot_heatmap_sw<- pheatmap::pheatmap(matrix_max_standardize_sw, #DO THE ABOVE ANALYSIS BY SUBSETTING TO SWITCHGRASS DROUGHT PLANTED SAMPLES THEN COMBINE AS BELOW
                                     cluster_cols = F,
                                     na_col = "white",
                                     border_color = "white",
                                     main="Switchgrass drought",
                                     color = viridis(n = 256, alpha = 1, 
                                                     begin = 0, end = 1, option = "viridis"
                                     ))


a <- list(plot_heatmap_bean[[4]])
a[[2]] <- plot_heatmap_sw[[4]]
z <- do.call(grid.arrange,a)
plot(z)


##Main manuscript Figure 7
ggsave(filename = "combined_heatmap_colorblindfriendly_final_planted_drought_bean_sw.tiff", plot = z,
       width = 25,
       height = 35, units = c("cm"),
       dpi = 300)


## Supplemental Figure S15, use for combining plots (drought unplanted samples)
plot_heatmap_bean <- pheatmap::pheatmap(matrix_max_standardize_bean, ##DO THE ABOVE ANALYSIS BY SUBSETTING TO BEAN DROUGHT UNPLANTED SAMPLES 
                                        cluster_cols = F,
                                        na_col = "white",
                                        border_color = "white",
                                        main="Bean drought - unplanted",
                                        color = viridis(n = 256, alpha = 1, 
                                                        begin = 0, end = 1, option = "viridis"
                                        ))
plot_heatmap_sw<- pheatmap::pheatmap(matrix_max_standardize_sw, ##DO THE ABOVE ANALYSIS BY SUBSETTING TO SWITCHGRASS DROUGHT UNPLANTED SAMPLES 
                                     cluster_cols = F,
                                     na_col = "white",
                                     border_color = "white",
                                     main="Switchgrass drought - unplanted",
                                     color = viridis(n = 256, alpha = 1, 
                                                     begin = 0, end = 1, option = "viridis"
                                     ))


a <- list(plot_heatmap_bean[[4]])
a[[2]] <- plot_heatmap_sw[[4]]
z <- do.call(grid.arrange,a)
plot(z)


ggsave(filename = "combined_heatmap_colorblindfriendly_final_unplanted_drought_bean_sw.tiff", plot = z,
       width = 25,
       height = 35, units = c("cm"),
       dpi = 300)







