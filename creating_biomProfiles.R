suppressMessages(library(tidyverse))

####This Rscript is to arrange salmon's quantification results into CAMI's biom format#### 

####Reading the metadata with results from taxonkit####
####Make sure this metadata file has a common column with the quantification results####
ncbi <- read.delim("/vol/projects/salampal/data/microbiome_assembly/ncbiNomenclature_gtdbReference.tsv",
                   header=T,sep="\t",stringsAsFactors = F)

####Splitting the taxonomy based on a delimiter. Taxonkit is ";"####
ncbiSplit <- separate(ncbi, ncbiReformed, into=c("kingdom","phylum","class","order","family","genus","species"), 
                      sep = ";", remove = FALSE,convert = FALSE, fill = "right")
ncbiSplit <- separate(ncbiSplit, taxIDReformed, into=c("kingdom_taxID","phylum_taxID","class_taxID","order_taxID","family_taxID","genus_taxID","species_taxID"), 
                      sep = ";", remove = FALSE,convert = FALSE, fill = "right")

####Function to create the biom format####
####salmonQuant <- Salmon's quantification file
####salmonMeta <- to pick the total number of reads to later estimate the percentage of reads
####outFileName <- prefix name of the biom file
####loopNo <- biom files have a pattern with respect to their goldstandards. 
####The initial headers in the biom format should match for OPAL to compare.####

creatingProfiles <- function(salmonQuant,salmonMeta,outFileName,loopNo){
  #Read the quant file.
  salmon_quant <- read.delim(salmonQuant,header=T,sep="\t",stringsAsFactors=F)
  #merge the metadata and drop the lines not present in the quant file.
  salmon_quant_joining <- left_join(salmon_quant,ncbiSplit,by=c("Name" = "accession"))
  #Read all the lines from the meta file
  noOfReads <- readLines(salmonMeta)
  #Pick the line with total number of reads
  reads <- as.numeric(gsub(",","",str_split(noOfReads[36],pattern=": ",simplify = T)[[2]]))
  
  #Warning message arise to different R versions (I have worked with only 3.5.1/2) and tidyverse library
  #Now on at each taxa level I summarise the percentage of reads and combine the required columns as per biom format
  #Example:
  #@@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE
  #2	superkingdom	2	Bacteria	98.4757
  k <- salmon_quant_joining %>% group_by(kingdom,kingdom_taxID) %>%
    summarise(salmon_count=sum(NumReads)) %>%
    mutate(salmon_count_percentage=round((100*(salmon_count/reads)),digits = 4))
  k <- subset(k,k$kingdom!="")
  df <- cbind("@@TAXID"=k$kingdom_taxID,"RANK"=rep("superkingdom",nrow(k)),"TAXPATH"=k$kingdom_taxID,"TAXPATHSN"=k$kingdom,"PERCENTAGE"=k$salmon_count_percentage)
  
  #From phylum on the summary changes. Removing any rows that do not have a taxa name or ID at that level (blanks and na's).
  
  p <- salmon_quant_joining %>% group_by(kingdom,phylum,kingdom_taxID,phylum_taxID) %>%
    summarise(salmon_count=sum(NumReads)) %>%
    mutate(salmon_count_percentage=round((100*(salmon_count/reads)),digits = 4))
  p <- subset(p,p$phylum!="")
  p <- subset(p,!is.na(p$phylum_taxID))
  p <- unite(p, TAXPATHSN,1:2,sep = "|", remove = TRUE, na.rm = FALSE)
  p <- unite(p, TAXPATH,2:3,sep = "|", remove = F, na.rm = FALSE)
  p_df <- cbind("@@TAXID"=p$phylum_taxID,"RANK"=rep("phylum",nrow(p)),"TAXPATH"=p$TAXPATH,"TAXPATHSN"=p$TAXPATHSN,"PERCENTAGE"=p$salmon_count_percentage)
  df <- rbind(df,p_df)
  
  c <- salmon_quant_joining %>% group_by(kingdom,phylum,class,kingdom_taxID,phylum_taxID,class_taxID) %>%
    summarise(salmon_count=sum(NumReads)) %>%
    mutate(salmon_count_percentage=round((100*(salmon_count/reads)),digits = 4))
  c <- subset(c,c$class!="")
  c <- subset(c,!is.na(c$class_taxID))
  c <- unite(c, TAXPATHSN,1:3,sep = "|", remove = TRUE, na.rm = FALSE)
  c <- unite(c, TAXPATH,2:4,sep = "|", remove = F, na.rm = FALSE)
  c_df <- cbind("@@TAXID"=c$class_taxID,"RANK"=rep("class",nrow(c)),"TAXPATH"=c$TAXPATH,"TAXPATHSN"=c$TAXPATHSN,"PERCENTAGE"=c$salmon_count_percentage)
  df <- rbind(df,c_df)
  
  o <- salmon_quant_joining %>% group_by(kingdom,phylum,class,order,kingdom_taxID,phylum_taxID,class_taxID,order_taxID) %>%
    summarise(salmon_count=sum(NumReads)) %>%
    mutate(salmon_count_percentage=round((100*(salmon_count/reads)),digits = 4))
  o <- subset(o,o$order!="")
  o <- subset(o,!is.na(o$order_taxID))
  o <- unite(o, TAXPATHSN,1:4,sep = "|", remove = TRUE, na.rm = FALSE)
  o <- unite(o, TAXPATH,2:5,sep = "|", remove = F, na.rm = FALSE)
  o_df <- cbind("@@TAXID"=o$order_taxID,"RANK"=rep("order",nrow(o)),"TAXPATH"=o$TAXPATH,"TAXPATHSN"=o$TAXPATHSN,"PERCENTAGE"=o$salmon_count_percentage)
  df <- rbind(df,o_df)
  
  f <- salmon_quant_joining %>% group_by(kingdom,phylum,class,order,family,kingdom_taxID,phylum_taxID,class_taxID,order_taxID,family_taxID) %>%
    summarise(salmon_count=sum(NumReads)) %>%
    mutate(salmon_count_percentage=round((100*(salmon_count/reads)),digits = 4))
  f <- subset(f,f$family!="")
  f <- subset(f,!is.na(f$family_taxID))
  f <- unite(f, TAXPATHSN,1:5,sep = "|", remove = TRUE, na.rm = FALSE)
  f <- unite(f, TAXPATH,2:6,sep = "|", remove = F, na.rm = FALSE)
  f_df <- cbind("@@TAXID"=f$family_taxID,"RANK"=rep("family",nrow(f)),"TAXPATH"=f$TAXPATH,"TAXPATHSN"=f$TAXPATHSN,"PERCENTAGE"=f$salmon_count_percentage)
  df <- rbind(df,f_df)
  
  g <- salmon_quant_joining %>% group_by(kingdom,phylum,class,order,family,genus,kingdom_taxID,phylum_taxID,class_taxID,order_taxID,family_taxID,genus_taxID) %>%
    summarise(salmon_count=sum(NumReads)) %>%
    mutate(salmon_count_percentage=round((100*(salmon_count/reads)),digits = 4))
  g <- subset(g,g$genus!="")
  g <- subset(g,!is.na(g$genus_taxID))
  g <- unite(g, TAXPATHSN,1:6,sep = "|", remove = TRUE, na.rm = FALSE)
  g <- unite(g, TAXPATH,2:7,sep = "|", remove = F, na.rm = FALSE)
  g_df <- cbind("@@TAXID"=g$genus_taxID,"RANK"=rep("genus",nrow(g)),"TAXPATH"=g$TAXPATH,"TAXPATHSN"=g$TAXPATHSN,"PERCENTAGE"=g$salmon_count_percentage)
  df <- rbind(df,g_df)
  
  s <- salmon_quant_joining %>% group_by(kingdom,phylum,class,order,family,genus,species,kingdom_taxID,phylum_taxID,class_taxID,order_taxID,family_taxID,genus_taxID,species_taxID) %>%
    summarise(salmon_count=sum(NumReads)) %>%
    mutate(salmon_count_percentage=round((100*(salmon_count/reads)),digits = 4))
  s <- subset(s,s$species!="")
  s <- subset(s,!is.na(s$species_taxID))
  s <- unite(s, TAXPATHSN,1:7,sep = "|", remove = TRUE, na.rm = FALSE)
  s <- unite(s, TAXPATH,2:8,sep = "|", remove = F, na.rm = FALSE)
  s_df <- cbind("@@TAXID"=s$species_taxID,"RANK"=rep("species",nrow(s)),"TAXPATH"=s$TAXPATH,"TAXPATHSN"=s$TAXPATHSN,"PERCENTAGE"=s$salmon_count_percentage)
  df <- rbind(df,s_df)
  
  #Convering a tibble to a dataframe
  df_df <- as.data.frame(df,stringsAsFactors = F)
  #Include headers as per biom format and keep the SampleID same as the one in goldstandards.
  headers <- c("#CAMI Taxonomic Profiling Output",paste("@SampleID:",loopNo,sep=""),"@Version:0.9.3","@Ranks:superkingdom|phylum|class|order|family|genus|species|strain",
               "@TaxonomyID:GTDB_r89","@__program__:Salmon")
  cat(headers,sep="\n",file = outFileName)
  write.table(df_df,file = outFileName,append=T,row.names = F,quote = F,sep="\t")
}

#Looping the function to create biom format profiles for all the samples in the toy datasets

for (i in "High") {
  whichComp <- paste("dataset_",i[1],"_Complexity",sep="")
  pathList <-  paste("/vol/projects/salampal/data/microbiome_assembly/cami_toy_datasets/",whichComp,"/",sep="")
  for (sketch in c(10)) {
    sp <- paste(pathList,whichComp,"_genomeSequences_mashScreen_GTDB_s",sketch,"K_k21_mm",sep="")
    for (mm in c(1)) {
      mmp <- paste(sp,mm,sep="")
      for (j in 1:5) {
        qfileList <- paste(mmp,"_sample",j,"_salmonQuant/quant.sf",sep="")
        mfileList <- paste(mmp,"_sample",j,"_salmonQuant/aux_info/meta_info.json",sep="")
        outName <- paste(mmp,"_sample",j,"_tpmCutoff_gt65_profile.tsv",sep="")
        creatingProfiles(qfileList,mfileList,outName,paste("HC_Sample",j,sep=""))
      }
    }
  }
}


for (i in "Low") {
  whichComp <- paste("dataset_",i,"_Complexity",sep="")
  pathList <-  paste("/vol/projects/salampal/data/microbiome_assembly/cami_toy_datasets/",whichComp,"/",sep="")
  for (sketch in c(10,4)) {
    sp <- paste(pathList,whichComp,"_genomeSequences_mashScreen_GTDB_s",sketch,"K_k21_mm",sep="")
    for (mm in c(1:7)) {
      mmp <- paste(sp,mm,sep="")
      qfileList <- paste(mmp,"_salmonQuant/quant.sf",sep="")
      mfileList <- paste(mmp,"_salmonQuant/aux_info/meta_info.json",sep="")
      outName <- paste(mmp,"_profile.tsv",sep="")
    }
  }
}

for (i in "Medium") {
  whichComp <- paste("dataset_",i,"_Complexity",sep="")
  pathList <-  paste("/vol/projects/salampal/data/microbiome_assembly/cami_toy_datasets/",whichComp,"/",sep="")
  for (sketch in c(10,4)) {
    sp <- paste(pathList,whichComp,"_genomeSequences_mashScreen_GTDB_s",sketch,"K_k21_mm",sep="")
    for (mm in c(1:7)) {
      mmp <- paste(sp,mm,sep="")
      for (j in 1:2) {
        qfileList <- paste(mmp,"_sample",j,"_salmonQuant/quant.sf",sep="")
        mfileList <- paste(mmp,"_sample",j,"_salmonQuant/aux_info/meta_info.json",sep="")
        outName <- paste(mmp,"_sample",j,"_profile.tsv",sep="")
        creatingProfiles(qfileList,mfileList,outName,paste("MC_Sample",j,"_180bp",sep=""))
      }
    }
  }
}
