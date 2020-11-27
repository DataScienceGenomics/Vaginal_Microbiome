#Phyloseq Obj Creation
#V34

library(CountClust)
library(dplyr)
library(readxl)
library(metagenomeSeq)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(pheatmap)

#Import OTU file
otu_mat_col_types <- c("text", rep("numeric", 111))##last 112 are int
otu_mat<- read_excel("~/Maternal_V34_100samples_OTU_table_mc2_w_tax_no_pynast_failures_12212018.xlsx", col_types=otu_mat_col_types )

#Import Taxa file
tax_mat_col_types<- rep("text", 13)##all 13 are int
tax_mat<- read_excel("~/Maternal_V34_100samples_Tax_table_mc2_w_tax_no_pynast_failures_12212018.xlsx", col_types=tax_mat_col_types)

tax_mat$Species <- gsub(" ", "", tax_mat$Species)
idx1 <- which(is.na(tax_mat$Species))
tax_mat$Species[idx1] <- tax_mat$Genus[idx1]
tax_mat<-tax_mat %>% subset(!is.na(Genus))
tax_mat<-tax_mat %>% subset(Genus != "g__")
tax_mat$SpeciesM<-tax_mat$Species
tax_mat$SpeciesM<-gsub("g__", "", tax_mat$SpeciesM)
tax_mat$SpeciesM<-gsub("s__", "", tax_mat$SpeciesM)
tax_mat$SpeciesM<-gsub(" ", "", tax_mat$SpeciesM)
tax_mat$Species <- gsub("^[gs]__", "", tax_mat$Species)

#Import Sample_data file
samples_df <- read.table("~/samples_df_V34.txt", sep ="\t",header = TRUE)

samples_df<-dplyr::filter(samples_df, !grepl("Water_040518", ID))
samples_df<-dplyr::filter(samples_df, !grepl("Water", ID))

#Get rid of duplicate IDs
samples_df<-dplyr::filter(samples_df, !grepl("IPM1150", IPM_ID))
samples_df<-dplyr::filter(samples_df, !grepl("IPM1151", IPM_ID))
samples_df<-dplyr::filter(samples_df, !grepl("IPM1203", IPM_ID))
samples_df<-dplyr::filter(samples_df, !grepl("IPM1152", IPM_ID))
samples_df<-dplyr::filter(samples_df, !grepl("IPM1204", IPM_ID))
samples_df<-dplyr::filter(samples_df, !grepl("IPM1205", IPM_ID))
samples_df<-dplyr::filter(samples_df, !grepl("IPM1153", IPM_ID))
samples_df<-dplyr::filter(samples_df, !grepl("IPM1162", IPM_ID))
samples_df<-dplyr::filter(samples_df, !grepl("IPM1215", IPM_ID))
samples_df<-dplyr::filter(samples_df, !grepl("IPM1163", IPM_ID))
samples_df<-dplyr::filter(samples_df, !grepl("IPM1219", IPM_ID))


#Define row names
names(otu_mat)<-str_replace_all(names(otu_mat), c(" " = "_" , "," = "" ))
row.names(otu_mat) <- otu_mat$OTU_ID
otu_mat1 <- otu_mat
otu_mat <- otu_mat %>% dplyr::select (-OTU_ID)
rownames(otu_mat)<-rownames(otu_mat1)
otu_mat<-as.data.frame(otu_mat)
otu_mat <- otu_mat[, colSums(otu_mat) > 0]

#Idem for two other matrixes
names(tax_mat)<-str_replace_all(names(tax_mat), c(" " = "_" , "," = "" ))
tax_mat<-as.data.frame(tax_mat)
row.names(tax_mat) <- tax_mat$OTU_ID
#tax_mat1<-tax_mat
tax_mat <- tax_mat %>% dplyr::select (-OTU_ID)
tax_mat<-tax_mat %>% dplyr::select(-c(NCBI_ID,NCBI_taxonomy,DB_ID,Center,Center_SEQ))
tax_mat = tax_mat[,c("Kingdom","Phylum","Class","Order","Family","Genus","Species","SpeciesM")]


otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)


#Transform to phyloseq item
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

TAX = tax_table(tax_mat)
rownames(samples_df)<-samples_df$IPM_ID
samples_df <- samples_df %>% dplyr::select (-IPM_ID)

samples1 = sample_data(samples_df)

carbom2020 <- phyloseq(OTU, TAX, samples1)

carbom2020

#Get rid of contam taxa
carbomVtaxa = subset_taxa(
  carbom2020,
  Species!="Staphylococcus_epidermis")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="Staphylococcus_hominis")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="Prevotella_intermedia")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Pediococcus")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Catonella")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Dermabacter")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="Porphyromonas_endodontalis")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="s__viridiflava")
carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Jeotgalicoccus")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="viridiflava")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Erythromicrobium")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Dermacoccus")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="Planomicrobium")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Calothrix")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Phormidium")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Rhodobacter")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="Aurantimonas")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="Phormidium")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="epidermidis")

carbom1 = subset_taxa(
  carbomVtaxa,
  Species!="Staphylococcus_epidermidis")

#Filter taxa
carbom1 = filter_taxa(carbom1, function(x) sum(x > 2) > (0.10*length(x)), TRUE)

carbom1


#V12

otu_mat <- read.table("~/M_generated_new_OTU_File19thApr2020.txt", sep ="\t",header = TRUE)
table(is.na(otu_mat))
otu_mat[is.na(otu_mat)] <- 0
tax_mat_col_types<- rep("text", 13)##all 13 are int
tax_mat<- read_excel("~/Maternal_V12_417Samples_TAXNOMY_assignment_no_pynast_failures_1_1_2019.xlsx", col_types=tax_mat_col_types)

idx <- which(tax_mat$Species == "s__")
tax_mat$Species[idx] <- tax_mat$Genus[idx]

idx1 <- which(is.na(tax_mat$Species))
tax_mat$Species[idx1] <- tax_mat$Genus[idx1]


samples_df <- read.table("~//samples_df_V12.txt", sep ="\t",header = TRUE)
sample_IPM<-samples_df%>%select(IPM_ID)
rownames(samples_df)<-samples_df$IPM_ID
samples_df<-samples_df%>%select(-IPM_ID)
rownames(samples_df)<-sample_IPM$IPM_ID

#Define row names
otu_mat <- otu_mat %>% remove_rownames
row.names(otu_mat) <- otu_mat$OTU.ID
k<-as.character(row.names(otu_mat))
otu_mat<-otu_mat[-1]
row.names(otu_mat) <- k
otu_mat<-as.data.frame(otu_mat)
table(is.na(otu_mat))


#Idem for two other matrixes
names(tax_mat)<-str_replace_all(names(tax_mat), c(" " = "_" , "," = "" ))
tax_mat<-as.data.frame(tax_mat)
row.names(tax_mat) <- tax_mat$OTU_ID
k1<-as.character(row.names(tax_mat))
tax_mat<-tax_mat[-1]
row.names(tax_mat) <- k1
tax_mat<-as.data.frame(tax_mat)
tax_mat<-tax_mat %>% subset(!is.na(Genus))
tax_mat<-tax_mat %>% subset(Genus != "g__")
tax_mat<-tax_mat %>% select(-c(NCBI_ID,NCBI_taxonomy,DB_ID,Center,Center_SEQ))
tax_mat = tax_mat[,c("Kingdom","Phylum","Class","Order","Family","Genus","Species")]
table(is.na(tax_mat))


#Idem for two other matrixes

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)


#Transform to phyloseq item
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

TAX = tax_table(tax_mat)
samples1 = sample_data(samples_df)

carbom2019 <- phyloseq(OTU, TAX, samples1)

carbom2019

#Get rid of contam taxa
carbomVtaxa = subset_taxa(
  carbom2019,
  Species!="s__Staphylococcus_epidermis")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="s__ Staphylococcus_hominis")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="s__Prevotella_intermedia")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Pediococcus")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Catonella")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Dermabacter")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="s__Porphyromonas_endodontalis")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="s__viridiflava")
carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Jeotgalicoccus")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Species!="s__viridiflava")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Erythromicrobium")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Dermacoccus")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="Planomicrobium")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Calothrix")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Phormidium")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="g__Rhodobacter")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="Aurantimonas")

carbomVtaxa = subset_taxa(
  carbomVtaxa,
  Genus!="Phormidium")

carbom = subset_taxa(
  carbomVtaxa,
  Species!="s__epidermidis")

carbom

#filtering criteria
carbom = filter_taxa(carbom2019, function(x) sum(x >2 ) > (0.10*length(x)), TRUE)
carbom
