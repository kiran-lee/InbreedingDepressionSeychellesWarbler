# This script investigates pedigree and genomic relatedness in ~1900 island-confined Seychelles warblers.
R.Version()
library(readr)
library(tidyverse)
library(dplyr)
library(readxl)
library(ggplot2)
library("data.table") 
setwd("~/Documents/GitHub/InProgressGenomicsInbreedingSeychellesWarblers/Data")


#Create the dataframe for the analysis
## We want to create dataframe of sample filenames (.bam files), their BirdID,and Coverage.
CoverageExtraSamples <- read.delim("coveragefilenameallsortedextrasamples.txt",sep=" ",header=F)
names(CoverageExtraSamples)[names(CoverageExtraSamples) == 'V1'] <- 'SeqID'
names(CoverageExtraSamples)[names(CoverageExtraSamples) == 'V2'] <- 'Coverage'
CoverageExtraSamples$Filepath<-CoverageExtraSamples$SeqID
CoverageExtraSamples$Filepath <- gsub("^.{0,4}", "", CoverageExtraSamples$Filepath)
CoverageExtraSamples$SeqID<-sub('.', '', CoverageExtraSamples$SeqID)
CoverageExtraSamples<-CoverageExtraSamples %>% separate(SeqID, c('Filepath1', 'Filepath2', 'Filepath3', 'Filepath4','Plate','Filepath5','SeqID'), sep = '/', convert = TRUE)
CoverageExtraSamples = subset(CoverageExtraSamples, select = c(Plate, SeqID, Coverage, Filepath))
CoverageExtraSamples$ID<-CoverageExtraSamples$SeqID

##Files to match BirdIDs
Identifiers<-read_excel("SheffieldSubmissions.xlsx")
Identifiers26076<-read_excel("ID 26076_Sample information table.xlsx")
PilotIdentifiers<-read_excel("SamplesForPilotTargetCapture_290119_sortBTN_Qubit.xlsx")
LIMS26629renamed<-read.table("lims26629renamed.txt",sep=" ",header=F)
LIMS26757p1raw4<-read.csv("Samples for Sequencing 25072023.csv",colClasses=c("NULL",NA,NA,NA,NA,NA))
MissingLIMS26757p1raw4 <- data.frame(BirdID=c(6572,6145,6373,5904,6144,6651), 
                                     FieldRing=c(NA,NA,NA,NA,NA,NA),
                                     BTO=c(NA,NA,NA,NA,NA,NA), 
                                     BloodID=c(8178, 5883,6287,5620,5880,7239), 
                                     BloodTubeNumber=c(NA,NA,NA,NA,NA,NA), 
                                     stringsAsFactors=FALSE)
LIMS26757p1raw4<-rbind(LIMS26757p1raw4,MissingLIMS26757p1raw4)
BloodID<-read.csv("BloodID.csv")

##Clean ID numbers
LIMS26629renamed<-LIMS26629renamed %>% arrange(V1)
LIMS26629renamed$SeqID<-CoverageExtraSamples$SeqID[CoverageExtraSamples$Plate=='LIMS26629']
LIMS26629renamed<-LIMS26629renamed %>%
  mutate(SeqID=sort(SeqID))
LIMS26629renamed$V1<-sub(".*\\-", "", LIMS26629renamed$V1)
colnames(LIMS26629renamed)[colnames(LIMS26629renamed) == 'V1'] <- 'ID'

#Fix 6 samples that were mislabeled, as identified by Rowan
LIMS26629renamed$ID[which(LIMS26629renamed$SeqID == "100_ACAAGAACCT-CGATACTGAA_L002__all_mapped_rehead.bam")] <- 2544
LIMS26629renamed$ID[which(LIMS26629renamed$SeqID == "101_AGAGTATGTG-AGATGGCTTC_L002__all_mapped_rehead.bam")] <- 3227
LIMS26629renamed$ID[which(LIMS26629renamed$SeqID == "96_CAACCATACA-ACCGGTTATA_L002__all_mapped_rehead.bam")] <- 764
LIMS26629renamed$ID[which(LIMS26629renamed$SeqID == "97_GTAGGCCGTT-GCCACTGTCT_L002__all_mapped_rehead.bam")] <- 2688
LIMS26629renamed$ID[which(LIMS26629renamed$SeqID == "98_CGGATTGATC-AGTCACAACA_L002__all_mapped_rehead.bam")] <- 2493
LIMS26629renamed$ID[which(LIMS26629renamed$SeqID == "99_ACTGGCAAGA-TGTTGTCCAT_L002__all_mapped_rehead.bam")] <- 2522

CoverageExtraSamples$ID[CoverageExtraSamples$Plate=="LIMS26629"]<-NA
CoverageExtraSamples <- merge(CoverageExtraSamples,LIMS26629renamed,by="SeqID", all = TRUE) %>%
  mutate(ID = coalesce(ID.x, ID.y)) %>%
  select (-c(ID.x, ID.y))



PilotSequences<-subset(CoverageExtraSamples, CoverageExtraSamples$Plate=="LIMS24675"|CoverageExtraSamples$Plate=="LIMS25133")
PilotSequences$ID<-trimws(sapply(strsplit(PilotSequences$ID, "_"), `[[`, 2))
MissingBloodIDSequences<-subset(CoverageExtraSamples, CoverageExtraSamples$Plate == "LIMS26076p4"|CoverageExtraSamples$Plate == "LIMS26076raw")
MissingBloodIDSequences$ID<-as.numeric(sapply(strsplit(MissingBloodIDSequences$ID, "_"), "[[", 1))
BloodIDSequences<-subset(CoverageExtraSamples, CoverageExtraSamples$Plate!="LIMS24675"& CoverageExtraSamples$Plate!="LIMS25133"&CoverageExtraSamples$Plate != "LIMS26076p4"&CoverageExtraSamples$Plate != "LIMS26076raw" )
BloodIDSequences$ID <- sub('_repeat','',BloodIDSequences$ID)
BloodIDSequences$ID <- sub('-repeat','',BloodIDSequences$ID)
BloodIDSequences$ID<-trimws(sapply(strsplit(BloodIDSequences$ID, "_"), `[[`, 1))
BloodIDSequences$ID<-trimws(sapply(strsplit(BloodIDSequences$ID, "-"), `[[`, 2))
BloodIDSequences$ID<-sub('.+-(.+)', '\\1', BloodIDSequences$ID)

##Pilot sequences use Blood Tube Number, the rest use BloodID
PilotSequences$Identifier<-paste("BloodTubeNumber")
MissingBloodIDSequences$Identifier<-paste("BloodID")
BloodIDSequences$Identifier<-paste("BloodID")

##Make IDs numeric
PilotSequences$ID<-as.numeric(PilotSequences$ID)
BloodIDSequences$ID<-as.numeric(BloodIDSequences$ID)
MissingBloodIDSequences$ID<-as.numeric(MissingBloodIDSequences$ID)
Identifiers26076$Sample_number<-as.numeric(Identifiers26076$`Sample number`)
Identifiers26076$BloodID<-as.numeric(Identifiers26076$`Sample name`)

#Join BirdIDs
PilotSequences <- PilotSequences %>% 
  left_join(select(PilotIdentifiers, BirdID, BloodTubeNumber), by = c("ID" = "BloodTubeNumber"))

BloodIDSequences$BloodID<-BloodIDSequences$ID
BloodIDSequences <- BloodIDSequences %>% 
  left_join(select(Identifiers, BirdID, BloodID), by = c("ID" = "BloodID"))
names(BloodIDSequences)[names(BloodIDSequences) == 'BirdID.x'] <- 'BirdID'
BloodIDSequences<-BloodIDSequences[!duplicated(BloodIDSequences), ]

MissingBloodIDSequences <- MissingBloodIDSequences %>% 
  left_join(select(Identifiers26076, Sample_number, BloodID), by = c("ID" = "Sample_number"))
MissingBloodIDSequences <- MissingBloodIDSequences %>% 
  left_join(select(Identifiers, BirdID, BloodID), by = c("BloodID" = "BloodID"))

##Concatenate into one file (2018 individuals including duplicates)
MissingBloodIDSequencesFormatted= subset(MissingBloodIDSequences, select = c(Plate, SeqID, Coverage, Filepath, BloodID, Identifier, BirdID))
colnames(MissingBloodIDSequencesFormatted)[colnames(MissingBloodIDSequencesFormatted) == 'BloodID'] <- 'ID'
BloodIDSequencesFormatted= subset(BloodIDSequences, select = c(Plate, SeqID, Coverage, Filepath, ID, Identifier, BirdID))
SequencedIndividualsBirdIDs=rbind(PilotSequences,MissingBloodIDSequencesFormatted,BloodIDSequencesFormatted)

##Add in missing BirdIDs from LIMS26757p1raw4
SequencedIndividualsBirdIDsExtra<- SequencedIndividualsBirdIDs %>% 
  left_join(select(LIMS26757p1raw4, BirdID, BloodID), by = c("ID" = "BloodID")) %>%
  mutate(BirdID = coalesce(BirdID.x, BirdID.y)) %>%
  select (-c(BirdID.x, BirdID.y))

##Add in inbreeding coefficients
##ROH calculated using PLINK:  --allow-extra-chr --bfile mergedimputedchromosomes --homozyg-density 200 --homozyg-gap 300 --homozyg-het 2 --homozyg-kb 1300 --homozyg-snp 50 --homozyg-window-het 2 --homozyg-window-missing 4 --homozyg-window-snp 50 --out mergedimputedchromosomes
MySequencedIndividuals<-read_excel("mergedimputedchromosomes.hom.indiv.xlsx")
MySequencedIndividuals$IID<-sub('.', '', MySequencedIndividuals$IID)
MySequencedIndividuals<-MySequencedIndividuals %>% separate(IID, c('Filepath1', 'Filepath2', 'Filepath3', 'Plate','Filepath4','SeqID'), sep = '/', convert = TRUE)
MySequencedIndividuals = subset(MySequencedIndividuals, select = c(Plate, SeqID, NSEG, KB, KBAVG))
MySequencedIndividuals$ID<-MySequencedIndividuals$SeqID
MySequencedIndividuals$FROH<-MySequencedIndividuals$KB/1091184475*1000
SequencedIndividualsBirdIDsExtraROH<- SequencedIndividualsBirdIDsExtra %>% 
  full_join(select(MySequencedIndividuals, FROH, SeqID), by = c("SeqID" = "SeqID"))

##Add in life-history traits
#Lifespan ----
##Read files
BirthDate <- read_csv("BirthDate27032023.csv", col_types = cols(BirthDate = col_date(format = "%d/%m/%Y"))) #In query table, this is BirdID
LastSeenYear <- read_csv("CurrentBTOextended27032023.csv") #In query table, this is CurrentBTOextended

##Make terms
###BirthYear from BirthDate
BirthDate <- BirthDate %>% 
  mutate(BirthYear = format(BirthDate, "%Y")) %>%
  mutate(BirthYear = as.numeric(BirthYear))

###Lifespan
Lifespan <- merge(BirthDate,LastSeenYear,by="BirdID", all = TRUE) %>% 
  mutate(LastSeenYea = as.numeric(LastSeenYea)) %>%
  mutate(BirthYear = as.numeric(BirthYear)) %>%
  mutate(Lifespan = LastSeenYea - BirthYear) %>%
  filter(LastSeenYea < 2022)  %>%
  select(BirdID,Lifespan,BirthYear,LastSeenYea)
colnames(Lifespan)[colnames(Lifespan) == 'LastSeenYea'] <- 'LastSeenYear'

##Link lifespan to dataframe
SequencedIndividualsBirdIDsExtraROH<- merge(SequencedIndividualsBirdIDsExtraROH, Lifespan, by="BirdID", all = TRUE)
SurvivingBirdsLifespan<-subset(SequencedLifespan, SequencedLifespan$Lifespan>0)

##Link n offspring to dataframe
Offspring <- read_csv("Offspring27032023.csv",col_types = cols(BirthDate = col_date(format = "%d/%m/%Y")))
ROcount <- Offspring %>% filter(Confidence > 80) %>% count(Parent) %>% rename(BirdID = Parent,ReproductiveOutput = n)
SequencedIndividualsBirdIDsExtraROH<- merge(SequencedIndividualsBirdIDsExtraROH, ROcount, by="BirdID", all = TRUE)

##Deduplicated dataset. For duplicated samples, pick the one with best coverage
SequencedIndividualsBirdIDsExtraDeduplicated<- SequencedIndividualsBirdIDsExtraROH %>% 
  group_by(BirdID) %>%
  top_n(1, abs(Coverage))

SequencedIndividualsBirdIDsExtraDuplicates <-SequencedIndividualsBirdIDsExtraROH[duplicated(SequencedIndividualsBirdIDsExtraROH$BirdID)|duplicated(SequencedIndividualsBirdIDsExtraROH$BirdID, fromLast=TRUE),]
SequencedIndividualsBirdIDsExtraDuplicates<-SequencedIndividualsBirdIDsExtraDuplicates %>% arrange((BirdID), Coverage)

##Quick visual inspection of inbreeding (FROH) depression (Lifespan and N offspring produced in lifetime)
ROHxLifespan<-ggplot(dplyr::filter(SequencedIndividualsBirdIDsExtraDeduplicated), aes(x=FROH, y=Lifespan)) + geom_point() + geom_smooth(method='lm') 
ROHxLifespan
ggsave("ROHxLifespan.png")

ROHxReproductiveOutput<-ggplot(SequencedIndividualsBirdIDsExtraDeduplicated , aes(x=FROH, y=ReproductiveOutput)) + geom_point() + geom_smooth(method='lm')
ROHxReproductiveOutput

##Does coverage affect inbreeding coefficient? Yes up to between 1.5x-2x.
sum((SequencedIndividualsBirdIDsExtra$Coverage>0 & SequencedIndividualsBirdIDsExtra$FROH>0), na.rm=TRUE)
ROHxCoverage<-ggplot(SequencedIndividualsBirdIDsExtra , aes(x=Coverage, y=FROH)) + geom_point() + geom_smooth(method='lm')+
  annotate("text",x=13, y=0.45, label= "All samples, n=1920")
ROHxCoverage

sum((SequencedIndividualsBirdIDsExtra$Coverage>0.5 & SequencedIndividualsBirdIDsExtra$FROH>0), na.rm=TRUE)
ROHxCoverage0.5x<-ggplot(dplyr::filter(SequencedIndividualsBirdIDsExtra, Coverage > 0.5) , aes(x=Coverage, y=FROH)) + geom_point() + geom_smooth(method='lm')+
  annotate("text",x=13, y=0.45, label= "Coverage>0.5x, n=1826")
ROHxCoverage0.5x

sum((SequencedIndividualsBirdIDsExtra$Coverage>1 & SequencedIndividualsBirdIDsExtra$FROH>0), na.rm=TRUE)
ROHxCoverage1x<-ggplot(dplyr::filter(SequencedIndividualsBirdIDsExtra, Coverage > 1) , aes(x=Coverage, y=FROH)) + geom_point() + geom_smooth(method='lm')+
  annotate("text",x=13, y=0.45, label= "Coverage>1x, n=1728")
ROHxCoverage1x

sum((SequencedIndividualsBirdIDsExtra$Coverage>1.5 & SequencedIndividualsBirdIDsExtra$FROH>0), na.rm=TRUE)
ROHxCoverage1.5x<-ggplot(dplyr::filter(SequencedIndividualsBirdIDsExtra, Coverage > 1.5) , aes(x=Coverage, y=FROH)) + geom_point() + geom_smooth(method='lm')+
  annotate("text",x=13, y=0.45, label= "Coverage>1.5x, n=1496")
ROHxCoverage1.5x

sum((SequencedIndividualsBirdIDsExtra$Coverage>2 & SequencedIndividualsBirdIDsExtra$FROH>0), na.rm=TRUE)
ROHxCoverage2x<-ggplot(dplyr::filter(SequencedIndividualsBirdIDsExtra, Coverage > 2) , aes(x=Coverage, y=FROH)) + geom_point() + geom_smooth(method='lm')+
  annotate("text",x=13, y=0.45, label= "Coverage>2x, n=1154")
ROHxCoverage2x

##Add in duplicates as found by KING
KINGduplicates<-read.delim("duplicates.txt",sep="",header=T,row.names=NULL)
DuplicateSamplesGTcheckSequencedCoverage$KINGduplicate = with(DuplicateSamplesGTcheckSequencedCoverage,
                                                              ifelse(Filepath.1 %in% KINGduplicates$IID1, "KINGduplicate",
                                                                     ifelse(Filepath.1 %in% KINGduplicates$IID2, "KINGduplicate",
                                                                            ifelse(Filepath.2 %in% KINGduplicates$IID1, "KINGduplicate",
                                                                                   ifelse(Filepath.2 %in% KINGduplicates$IID2, "KINGduplicate", NA)
                                                                            ))
                                                              )
)

DuplicateSamplesGTcheckSequencedCoverage$KINGkinship = with(DuplicateSamplesGTcheckSequencedCoverage,
                                                            ifelse(Filepath.1 %in% KINGduplicates$IID1, KINGduplicates$KINSHIP[match(Filepath.1,KINGduplicates$IID1)],
                                                                   ifelse(Filepath.1 %in% KINGduplicates$IID2, KINGduplicates$KINSHIP[match(Filepath.1,KINGduplicates$IID2)],
                                                                          ifelse(Filepath.2 %in% KINGduplicates$IID1, KINGduplicates$KINSHIP[match(Filepath.2,KINGduplicates$IID1)],
                                                                                 ifelse(Filepath.2 %in% KINGduplicates$IID2, KINGduplicates$KINSHIP[match(Filepath.2,KINGduplicates$IID1)], NA)
                                                                          ))
                                                            )
)

DuplicateSamplesGTcheckSequenced$KINGkinship = with(DuplicateSamplesGTcheckSequenced,
                                                    ifelse(Filepath.1 %in% KINGduplicates$IID1, KINGduplicates$KINSHIP[match(Filepath.1,KINGduplicates$IID1)],
                                                           ifelse(Filepath.1 %in% KINGduplicates$IID2, KINGduplicates$KINSHIP[match(Filepath.1,KINGduplicates$IID2)],
                                                                  ifelse(Filepath.2 %in% KINGduplicates$IID1, KINGduplicates$KINSHIP[match(Filepath.2,KINGduplicates$IID1)],
                                                                         ifelse(Filepath.2 %in% KINGduplicates$IID2, KINGduplicates$KINSHIP[match(Filepath.2,KINGduplicates$IID1)], NA)
                                                                  ))
                                                    )
)
##To do: add in other variables that might influence inbreeding depression
#Fixed effects: Territory
#Random effects: Year of sample


##Output the dataframe used
write.table(SequencedIndividualsBirdIDsExtra, file = "SequencedIndividualsBirdIDsExtra.txt", sep = "\t",
            col.names = T, row.names = F, quote = FALSE)

### Script to choose unrelated samples to resequence and create a reference panel ###
##Choose samples to resequence
BestCoverageUnrelated<-read.table("unrelatedimputedchromosomes250sampled2nddegree.king.cutoff.in.id")
BestCoverageUnrelated$SeqID<-basename(BestCoverageUnrelated$V2)
BestCoverageUnrelatedSampleID<- merge(BestCoverageUnrelated, SequencedIndividualsBirdIDsExtraDeduplicated, by="SeqID", all = FALSE)
BestCoverageUnrelatedToResequence<-BestCoverageUnrelatedSampleID %>%
  select(SeqID, Plate, ID, Identifier, BirdID, Coverage )
###Add in plate and well numbers to help find libraries
BestCoverageUnrelatedToResequenceWellID<- merge(BestCoverageUnrelatedToResequence, Identifiers, by="BirdID", all = FALSE)
Identifiers26076<- Identifiers26076 %>% 
  left_join(select(BloodIDSequences, BloodID, BirdID) , by = c("BloodID" = "BloodID"))
BestCoverageUnrelatedToResequenceWellIDPlate <- merge(BestCoverageUnrelatedToResequenceWellID, Identifiers26076, by= "BirdID", all= TRUE)
BestCoverageUnrelatedToResequenceWellIDPlate <- subset(BestCoverageUnrelatedToResequenceWellIDPlate, !is.na(BestCoverageUnrelatedToResequenceWellIDPlate$SeqID))
BestCoverageUnrelatedToResequenceWellIDPlate <- BestCoverageUnrelatedToResequenceWellIDPlate[-c(18, 20, 21), ]

BestCoverageUnrelatedToResequenceWellIDPlate <- subset(BestCoverageUnrelatedToResequenceWellIDPlate,BestCoverageUnrelatedToResequenceWellIDPlate$Identifier=="BloodID")

BestCoverageUnrelatedToResequenceWellIDPlate <- subset(BestCoverageUnrelatedToResequenceWellIDPlate,BestCoverageUnrelatedToResequenceWellIDPlate$Plate.x!="LIMS24675")
BestCoverageUnrelatedToResequenceWellIDPlate <- subset(BestCoverageUnrelatedToResequenceWellIDPlate,BestCoverageUnrelatedToResequenceWellIDPlate$Plate.x!="LIMS25133")

openxlsx::write.xlsx(BestCoverageUnrelatedToResequenceWellIDPlate, file ="BestCoverageUnrelatedToResequenceWellIDPlate.xlsx", quote=FALSE)
getwd()



###########################################################################################
#Sequioa
library(sequoia)  
library(usethis) 
usethis::edit_r_environ()

#First, for duplicate samples, pick the sample with highest coverage and then use PLINK to subset these and create a .raw file
DeduplicatedforSequioa<-subset(SequencedIndividualsBirdIDsExtraDeduplicated, Plate!="LIMS26757p1raw4" & Plate!="LIMS26757p4")
RenameforSequioa<-data.frame(rep("0",length(DeduplicatedforSequioa$Filepath)),DeduplicatedforSequioa$Filepath,rep("0",length(DeduplicatedforSequioa$Filepath)),DeduplicatedforSequioa$BirdID)
DeduplicatedforSequioa<-data.frame(rep("0",length(DeduplicatedforSequioa$Filepath)),DeduplicatedforSequioa$Filepath)

write.table(DeduplicatedforSequioa, file = "DeduplicatedforSequioa.txt", sep = "\t",
            col.names = F, row.names = F, quote = FALSE)
write.table(RenameforSequioa, file = "RenameforSequioa.txt", sep = "\t",
            col.names = F, row.names = F, quote = FALSE)

#Read in SNPs in .raw format
sw_GenoM_family <- GenoConvert(InFile = "imputedinputfile_for_sequoia.raw", InFormat="raw")

#Read in pedigree from masterbayes
masterped <- read.csv("pedFINAL_var124567_combined_unique_20230413.csv")
head(masterped)
names(masterped)[names(masterped) == 'dad'] <- 'sire'
names(masterped)[names(masterped) == 'mum'] <- 'dam'

#Create a more confident masterped
masterpedconfident<-subset(masterped, masterped$p>0.9)

#drop unnecessary column
masterped <- masterped[-c(4)]
masterpedconfident <- masterpedconfident[-c(4)]

#Read in life history data
sw_LifeHistData <- read.csv("LifeHistoryData.csv", header = TRUE)
sw_LifeHistData<-sw_LifeHistData[!duplicated(sw_LifeHistData),]
sw_LifeHistData$Sex <- case_when(sw_LifeHistData$Sex == "0" ~ 1, 
                                 sw_LifeHistData$Sex == "1" ~ 2, 
                                 TRUE ~ 3)
#convert years
for (n in seq(length(sw_LifeHistData$BirthDate))) {
  sw_LifeHistData$Year[n] <- as.numeric(substr(sw_LifeHistData$BirthDate[n], 7, 11))
}

#drop birds born before 1991 as they are not sampled
sw_LifeHistData <- sw_LifeHistData[sw_LifeHistData$Year > 1990,]

#drop birds not in cousin
sw_LifeHistData <- sw_LifeHistData[sw_LifeHistData$Island == "CN",]

#drop unnecessary column
sw_LifeHistData <- sw_LifeHistData[,-2]
sw_LifeHistData <- sw_LifeHistData[,-3]

sw_LifeHistData<-na.omit(sw_LifeHistData)


#Looks like lots of duplicate BirdIDs so create new life hist file
BirthDate <- read_csv("BirthDate27032023.csv", col_types = cols(BirthDate = col_date(format = "%d/%m/%Y"))) #In query table, this is BirdID
BirdIDSexYear <- read_excel("BirdIDSexYear.xlsx")
lifehist<- merge(BirthDate, BirdIDSexYear, by= "BirdID", all= TRUE)%>%
  select(BirdID,Sex,BirthDate)
lifehist<-lifehist[!duplicated(lifehist),]
lifehist$Sex <- case_when(lifehist$Sex == "0" ~ 1, 
                          lifehist$Sex == "1" ~ 2, 
                          TRUE ~ 3)
lifehist$BirthDate<-format(as.Date(lifehist$BirthDate, format="%Y/%m/%d"),"%Y")
colnames(lifehist)[colnames(lifehist) == 'BirthDate'] <- 'Year'
lifehist<-lifehist[-c(1, 2), ]  
lifehist$BY.min<-lifehist$Year
lifehist$BY.max<-lifehist$Year
lifehist<-merge(lifehist, LastSeenYear, by = "BirdID", all = TRUE)
lifehist<-subset(lifehist,select=-c(6,7,9))
colnames(lifehist)[colnames(lifehist) == 'LastSeenYea'] <- 'Year.last'

# And a quick check of how many parents we *should* be able to assign based on this test family...
sum(is.na(masterped$sire) == FALSE)
sum(is.na(masterped$dam) == FALSE)

# 1. Initial diagnostics --------------------------------------------------

# Check genotype calls across the dataset (-9 = missing calls)
summary.factor(as.factor(sw_GenoM_family))/(nrow(sw_GenoM_family)*ncol(sw_GenoM_family))

# Check for Mendelian errors in the data at the locus level
GenoM_checks <- SnpStats(sw_GenoM_family, Ped = masterped[,c("id", "dam", "sire")])

ggplot(GenoM_checks, aes(x = OHdam)) + geom_histogram() + theme_classic()
ggplot(GenoM_checks, aes(x = OHsire)) + geom_histogram() + theme_classic()
sw_OHLLR_DATABASE <- CalcOHLLR(Pedigree = masterped[,c("id", "dam", "sire")], GenoM = sw_GenoM_family)
ggplot(sw_OHLLR_DATABASE, aes(y = LLRdam, x = OHdam)) + geom_point(alpha = 0.6) + theme_classic()
ggplot(sw_OHLLR_DATABASE, aes(y = LLRsire, x = OHsire)) + geom_point(alpha = 0.6) + theme_classic()
SummarySeq(sw_OHLLR_DATABASE, Panels = "OH")

#Looks like 62/676 SNPs don't work well so remove
PoorSNPs<-subset(GenoM_checks, GenoM_checks$Err.hat>0.1 )
PoorSNPs<-rownames_to_column(PoorSNPs, var = "SNP")
PoorSNPs$SNP<- gsub("SNP","",PoorSNPs$SNP)
paste0( PoorSNPs$SNP, collapse=",")
dropSNP<-c(366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,486,487,488,489,490,491,492,493,494,495,496,497,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541)

sw_GenoM_family <- GenoConvert(InFile =sw_GenoM_family, InFormat="seq",
                               dropcol = c(366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,486,487,488,489,490,491,492,493,494,495,496,497,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541))

#Redo check. Looks much better without those 62 SNPs.
GenoM_checks <- SnpStats(sw_GenoM_family, Ped = masterped[,c("id", "dam", "sire")])
ggplot(GenoM_checks, aes(x = OHdam)) + geom_histogram() + theme_classic()
ggplot(GenoM_checks, aes(x = OHsire)) + geom_histogram() + theme_classic()
sw_OHLLR_DATABASE <- CalcOHLLR(Pedigree = masterped[,c("id", "dam", "sire")], GenoM = sw_GenoM_family)
ggplot(sw_OHLLR_DATABASE, aes(y = LLRdam, x = OHdam)) + geom_point(alpha = 0.6) + theme_classic()
ggplot(sw_OHLLR_DATABASE, aes(y = LLRsire, x = OHsire)) + geom_point(alpha = 0.6) + theme_classic()
SummarySeq(sw_OHLLR_DATABASE, Panels = "OH")


# 2. Run basic sequoia ------------------------------------

# Note modules used to build up from basic input check to analysis:
# i.   'pre': Input check
# ii.  'dup': check for duplicates
# iii. 'par': parentage assignment
# iv.  'ped': full pedigree reconstruction

# Just parentage assignment (Module = "par")
sw_family_sequoia_justparents <- sequoia(GenoM = sw_GenoM_family, LifeHistData = lifehist, Module = "par", quiet="verbose", args.AP= list(MaxAgeParent = c(12, 12)))
sw_family_pedigree_justparents_comparison <- PedCompare(masterped[,c("id", "dam", "sire")], sw_family_sequoia_justparents$Pedigree, Plot=TRUE, Symmetrical = TRUE)
sw_family_pedigree_justparents_comparison <- PedCompare(masterpedconfident[,c("id", "dam", "sire")], sw_family_sequoia_justparents$Pedigree, Plot=TRUE, Symmetrical = TRUE)

SummarySeq(SeqList = sw_family_sequoia_justparents)

#Investigate offspring to dam mismatches and birth year of offspring, dam1 and dam2, last seen year of offspring, dam1 and dam2, the territory ID (of offspring first seen date and last seen date, of dam first seen date and last seen date)
mismatch<-sw_family_pedigree_justparents_comparison[["Mismatch"]]
mismatch<-merge(mismatch, Lifespan , by.x= "id", by.y = "BirdID", all= FALSE)
colnames(mismatch)[colnames(mismatch) == 'BirthYear'] <- 'OffspringBirthYear'
colnames(mismatch)[colnames(mismatch) == 'LastSeenYear'] <- 'OffspringLastSeenYear'
mismatch<-merge(mismatch, Lifespan , by.x= "dam.1", by.y = "BirdID", all= FALSE)
colnames(mismatch)[colnames(mismatch) == 'BirthYear'] <- 'MasterbayesDamBirthYear'
colnames(mismatch)[colnames(mismatch) == 'LastSeenYear'] <- 'MasterbayesDamLastSeenYear'
mismatch<-merge(mismatch, Lifespan , by.x= "dam.2", by.y = "BirdID", all= FALSE)
colnames(mismatch)[colnames(mismatch) == 'BirthYear'] <- 'SequoiaDamBirthYear'
colnames(mismatch)[colnames(mismatch) == 'LastSeenYear'] <- 'SequoiaDamLastSeenYear'
mismatch$MasterbayesDamOffspringYearClash<-mismatch$MasterbayesDamLastSeenYear-mismatch$OffspringBirthYear
mismatch$SequoiaDamOffspringYearClash<-mismatch$SequoiaDamLastSeenYear-mismatch$OffspringBirthYear
OffspringTerritories<- BirdIDSexYear %>% 
  group_by(BirdID) %>%
  top_n(-1, abs(Year)) %>%
  select(BirdID,TerritoryID)
DamTerritories<- BirdIDSexYear %>% 
  group_by(BirdID) %>%
  top_n(1, abs(Year))%>%
  select(BirdID,TerritoryID)
mismatch<-merge(mismatch, OffspringTerritories , by.x= "id", by.y = "BirdID", all= FALSE)
colnames(mismatch)[colnames(mismatch) == 'TerritoryID'] <- 'OffspringTerritoryID'
mismatch<-merge(mismatch, DamTerritories , by.x= "dam.1", by.y = "BirdID", all= FALSE)
colnames(mismatch)[colnames(mismatch) == 'TerritoryID'] <- 'MasterbayesDamTerritoryID'
mismatch<-merge(mismatch, DamTerritories , by.x= "dam.2", by.y = "BirdID", all= FALSE)
colnames(mismatch)[colnames(mismatch) == 'TerritoryID'] <- 'SequoiaDamTerritoryID'

# Full pedigree reconstruction (Module = "ped")
sw_family_sequoia <- sequoia(GenoM = sw_GenoM_family, LifeHistData = lifehist, Module = "ped", quiet="verbose", args.AP= list(MaxAgeParent = c(12, 12)))
sw_family_pedigree_comparison <- PedCompare(masterped[,c("id", "dam", "sire")], sw_family_sequoia$Pedigree, Plot=TRUE)$MergedPed

# 3. Sequoia parameters: improving performance ---------------------------- HAVEN'T TESTED THIS YET- KL.

## Sequoia has a number of parameters that you can play with, but three most important are probably:
# 1. Err = estimated genotyping error rate
# 2. Tfilter = threshold LLR between proposed relationship vs unrelated, to select candidate relatives
# 3. Tassign = minimum LLR required for acceptance of proposed relationship, relative to next most likely relationship.

# I found in my data that playing with Err rates immediately solved my assignment issues, as part of the problem
# was that there were clear genotyping errors in the data due to low coverage. To directly assess this, look at 
# assignment rates in this family across different error rates. 

# The purpose of this section is to assess how sequoia responds to different error parameters. 

seq_error <- c(0.001, 0.01, 0.05, 0.1, 0.2)
seq_loopinfo <- data.frame("Error_rate"= NA_integer_,
                           "DAM - Match - assigned correctly" = NA_integer_, "DAM - Mismatch - assigned incorrectly" = NA_integer_,
                           "DAM - Mismatch - null assignment incorrect" = NA_integer_, "DAM - Mismatch - assigned new incorrect" = NA_integer_,
                           "DAM - Match - null assignment correct" = NA_integer_,
                           "SIRE - Match - assigned correctly" = NA_integer_, "SIRE - Mismatch - assigned incorrectly" = NA_integer_,
                           "SIRE - Mismatch - null assignment incorrect" = NA_integer_, "SIRE - Mismatch - assigned new incorrect" = NA_integer_,
                           "SIRE - Match - null assignment correct" = NA_integer_)

for (error in 1:length(seq_error)) {
  # Run sequoia
  seq_run <- sequoia(GenoM = sw_GenoM_family, LifeHistData = lifehist, Module = "par", Err = seq_error[error], quiet = TRUE)
  seq_comp <- PedCompare(masterped[,c("id", "dam", "sire")], seq_run$Pedigree, Plot=FALSE)$MergedPed
  
  seq_comp <- seq_comp[!(grepl("F00", seq_comp$id) | grepl("M00", seq_comp$id)),]
  
  # Dam summary stats
  dam_seq_match <- sum(seq_comp$dam.class == "Match")
  dam_seq_mismatch <- sum(seq_comp$dam.class == "Mismatch")
  dam_seq_notassigned <- sum(seq_comp$dam.class == "P1only")
  dam_seq_assigneddummy <- sum(seq_comp$dam.class == "P2only")
  dam_seq_assignednull <- sum(seq_comp$dam.class == "_")
  
  # Sire summary stats
  sire_seq_match <- sum(seq_comp$sire.class == "Match")
  sire_seq_mismatch <- sum(seq_comp$sire.class == "Mismatch")
  sire_seq_notassigned <- sum(seq_comp$sire.class == "P1only")
  sire_seq_assigneddummy <- sum(seq_comp$sire.class == "P2only")
  sire_seq_assignednull <- sum(seq_comp$sire.class == "_")
  
  seq_loopinfo <- rbind(seq_loopinfo, c(seq_error[error],
                                        dam_seq_match, dam_seq_mismatch, dam_seq_notassigned, dam_seq_assigneddummy, dam_seq_assignednull,
                                        sire_seq_match, sire_seq_mismatch, sire_seq_notassigned, sire_seq_assigneddummy, sire_seq_assignednull))
}

rm(dam_seq_assigneddummy, dam_seq_assignednull, dam_seq_match, dam_seq_mismatch, dam_seq_notassigned, 
   sire_seq_assigneddummy, sire_seq_assignednull, sire_seq_match, sire_seq_mismatch, sire_seq_notassigned, 
   error, seq_comp, seq_run, seq_error)

seq_loopinfo <- seq_loopinfo[!is.na(seq_loopinfo$Error_rate),]
seq_loopinfo_long <- melt(seq_loopinfo, id = c("Error_rate"))
seq_loopinfo_long$value_prop <- as.numeric(seq_loopinfo_long$value)

seq_loopinfo_long$Subset <- gsub("\\..*", "", seq_loopinfo_long$variable)
seq_loopinfo_long$Assignment_category <- gsub("^[^.]+.", "", seq_loopinfo_long$variable)

seq_loopinfo_long$Assignment_category <- factor(seq_loopinfo_long$Assignment_category, levels = c("..Match...assigned.correctly", "..Match...null.assignment.correct", "..Mismatch...assigned.incorrectly", 
                                                                                                  "..Mismatch...null.assignment.incorrect", "..Mismatch...assigned.new.incorrect"))

seq_loopinfo_long_SUMMARY <- seq_loopinfo_long %>% group_by(Error_rate, Assignment_category) %>% summarize(Subset = "Total", total_value = sum(value), total_value_prop = sum(value_prop))

# Graphing the proportion of match/mismatch assignments with both dams and sires
ggplot(seq_loopinfo_long_SUMMARY, aes(x = Error_rate, y = total_value_prop, col = Assignment_category, group = Assignment_category)) +
  geom_point() + geom_line(size = 1) + 
  geom_abline(linetype = "dashed", slope = 0, size = 1, 
              intercept = sum(is.na(sw.extended.family$dam), is.na(sw.extended.family$sire)), 
              color = "#4247a9") +
  geom_abline(linetype = "dashed", slope = 0, size = 1, 
              intercept = sum(is.na(sw.extended.family$dam) == FALSE, is.na(sw.extended.family$sire) == FALSE), 
              color = "#008CCC") +
  scale_color_manual(values = c("#008CCC", "#4247a9", "#800020", "#C31F48", "#d9b3bc")) +
  theme_classic() 

# The same but dams only
ggplot(subset(seq_loopinfo_long, seq_loopinfo_long$Subset == "DAM"), aes(x = Error_rate, y = value, col = Assignment_category, group = variable)) +
  geom_point() + geom_line(size = 1) + 
  scale_color_manual(values = c("#008CCC", "#4247a9", "#800020", "#C31F48", "#d9b3bc")) +
  geom_abline(linetype = "dashed", slope = 0, size = 1, 
              intercept = sum(is.na(sw.extended.family$dam)), 
              color = "#4247a9") +
  geom_abline(linetype = "dashed", slope = 0, size = 1,  
              intercept = sum(is.na(sw.extended.family$dam) == FALSE), 
              color = "#008CCC") +
  theme_classic() +
  ggtitle("Dam assignments only")


# Note that this is set up to *only* do parentage assignment, not full reconstruction. 
# With this dataset, I get the same conclusion using both module = "par" and module = "ped", 
# but the parentage assignment is substantially faster. Plus due to the small number of samples,
# full reconstruction tends to be overenthusiastic with filling null parents. This tends to resolve
# when doing full reconstruction where all possible parents are genotyped, as sequoia can then
# actually find the correct parents. 

# Conclusion: running sequoia with Err = 0.03 seems to optimise assignments, by maximising the
# number of correct matches, while minimising mismatches and missed assignments. 
# Let's have a closer look:

# Just parentage assignment (Module = "par")
sw_family_sequoia_justparents_err <- sequoia(GenoM = sw_GenoM_family, LifeHistData = lifehist, Module = "par", Err = 0.03)
sw_family_pedigree_justparents_err_comparison <- PedCompare(masterped[,c("id", "dam", "sire")], sw_family_sequoia_justparents_err$Pedigree, Plot=TRUE)$MergedPed

# Full pedigree reconstruction (Module = "ped")
sw_family_sequoia_err <- sequoia(GenoM = sw_GenoM_family, LifeHistData = lifehist, Module = "ped", Err = 0.03)
sw_family_pedigree_err_comparison <- PedCompare(masterped[,c("id", "dam", "sire")], sw_family_sequoia_err$Pedigree, Plot=TRUE)$MergedPed

### Extra imputation accuracy checks. Don't bother. Because we are comparing against samples which did not sequence well in the first place. ###

SequencedIndividualsBirdIDsExtraDuplicates$Repeat<-c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,3,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2)
library("reshape2")
SequencedIndividualsBirdIDsExtraDuplicatesWide <- reshape(SequencedIndividualsBirdIDsExtraDuplicates, idvar = "BirdID", timevar = "Repeat", direction = "wide")
SequencedIndividualsBirdIDsExtraDuplicatesWide$Coverage.1<-as.numeric(SequencedIndividualsBirdIDsExtraDuplicatesWide$Coverage.1)
SequencedIndividualsBirdIDsExtraDuplicatesWideCoverage <- subset(SequencedIndividualsBirdIDsExtraDuplicatesWide,SequencedIndividualsBirdIDsExtraDuplicatesWide$Coverage.1>0.1)

#FROH of duplicate samples shows some positive correlation
ROHxROH<-ggplot(SequencedIndividualsBirdIDsExtraDuplicatesWideCoverage , aes(x=FROH.1, y=FROH.2)) + geom_point() + geom_smooth(method='lm') + expand_limits(x = 0, y = 0)
ROHxROH

#If if I were to re-run STITCH using these samples I missed, then we would get on average 0.945x extra coverage per sample.
mean(SequencedIndividualsBirdIDsExtraDuplicatesWide$Coverage.2-SequencedIndividualsBirdIDsExtraDuplicatesWide$Coverage.1[SequencedIndividualsBirdIDsExtraDuplicatesWide$Plate.1=="LIMS26757p1raw4"])
CovxCov<-ggplot(dplyr::filter(SequencedIndividualsBirdIDsExtraDuplicatesWide, Plate.1=="LIMS26757p1raw4") , aes(x=Coverage.1, y=Coverage.2)) + geom_point() + geom_smooth(method='lm')
CovxCov

#Create file of duplicates for gtcheck to check concordance
DuplicateSamplesGTcheckSequenced<-subset(SequencedIndividualsBirdIDsExtraDuplicatesWide, SequencedIndividualsBirdIDsExtraDuplicatesWide$Plate.1!="LIMS26757p1raw4")
DuplicateSamplesGTcheckSequenced<-subset(DuplicateSamplesGTcheckSequenced, DuplicateSamplesGTcheckSequenced$Plate.2!="LIMS26757p1raw4")
DuplicateSamplesGTcheck<-subset(DuplicateSamplesGTcheckSequenced, select = c(Filepath.1, Filepath.2))
write.table(DuplicateSamplesGTcheck, file = "DuplicateSamplesGTcheck.txt", sep = "\t",
            col.names = F, row.names = F, quote = FALSE)

#Add in concordance for duplicates
ConcordanceDuplicates <- read.delim("ndiscordanceduplicatesproportion.txt",sep=" ",header=T,col.names=c("nDiscordant","n","Concordance"),row.names=NULL)
ConcordanceDuplicates$Filepath.1<-DuplicateSamplesGTcheck$Filepath.1
ConcordanceDuplicates$Filepath.2<-DuplicateSamplesGTcheck$Filepath.2
SequencedIndividualsBirdIDsExtraDuplicatesWide<- SequencedIndividualsBirdIDsExtraDuplicatesWide %>% 
  left_join(select(ConcordanceDuplicates, nDiscordant, n, Concordance), by = c("Filepath.1" = "Filepath.1"))

DuplicateSamplesGTcheckSequenced$nDiscordant<-ConcordanceDuplicates$nDiscordant
DuplicateSamplesGTcheckSequenced$n<-ConcordanceDuplicates$n
DuplicateSamplesGTcheckSequenced$Concordance<-ConcordanceDuplicates$Concordance

#Impact of duplicate concordance on FROH
hist(DuplicateSamplesGTcheckSequenced$Concordance, breaks=40)
ROHxROH<-ggplot(dplyr::filter(DuplicateSamplesGTcheckSequenced, Coverage.1>0) , aes(x=FROH.1, y=FROH.2, color=Concordance)) + geom_point() + geom_smooth(method='lm')
ROHxROH

CovxCon<-ggplot(dplyr::filter(DuplicateSamplesGTcheckSequenced) , aes(x=(Coverage.2-Coverage.1), y=Concordance)) + geom_point() + geom_smooth(method='lm')
CovxCon

#Randomly take 29 non-identical samples to compare genotypes
Random29Samples<-subset(SequencedIndividualsBirdIDsExtraDeduplicated,SequencedIndividualsBirdIDsExtraDeduplicated$Plate!="LIMS26757p1raw4")
Random29Samples<-Random29Samples[sample(nrow(Random29Samples), 29),]
Random29Samples$NonDuplicate<-DuplicateSamplesGTcheck$Filepath.1
Random29Samples<-subset(Random29Samples, select = c(Filepath, NonDuplicate))
duplicated(Random29Samples) #check these are non-duplicate pairs
write.table(Random29Samples, file = "Random29Samples.txt", sep = "\t",
            col.names = F, row.names = F, quote = FALSE)

#Add in concordance of non-duplicates
ConcordanceNonDuplicates <- read.delim("ndiscordancenonduplicatesproportion.txt",sep=" ",header=T,col.names=c("nDiscordant","n","Concordance"),row.names=NULL)
Random29Samples$nDiscordant<-ConcordanceNonDuplicates$nDiscordant
Random29Samples$n<-ConcordanceNonDuplicates$n
Random29Samples$Concordance<-ConcordanceNonDuplicates$Concordance
hist(Random29Samples$Concordance, breaks=40)
Random29Samples<-Random29Samples %>%
  mutate(Comparison="NonDuplicate")
DuplicateSamplesGTcheckSequenced<-DuplicateSamplesGTcheckSequenced %>%
  mutate(Comparison="Duplicate")

#Now compare concordance of duplicates and non-duplicates in histogram
DuplicatesConcordance<-data.frame(DuplicateSamplesGTcheckSequenced$Concordance,DuplicateSamplesGTcheckSequenced$Comparison)
NonDuplicatesConcordance<-data.frame(Random29Samples$Concordance,Random29Samples$Comparison)
colnames(DuplicatesConcordance)[colnames(DuplicatesConcordance) == 'DuplicateSamplesGTcheckSequenced.Concordance'] <- 'Concordance'
colnames(DuplicatesConcordance)[colnames(DuplicatesConcordance) == 'DuplicateSamplesGTcheckSequenced.Comparison'] <- 'Comparison'
colnames(NonDuplicatesConcordance)[colnames(NonDuplicatesConcordance) == 'Random29Samples.Concordance'] <- 'Concordance'
colnames(NonDuplicatesConcordance)[colnames(NonDuplicatesConcordance) == 'Random29Samples.Comparison'] <- 'Comparison'
NonDupDup<-rbind(DuplicatesConcordance,NonDuplicatesConcordance)
hist(NonDupDup$Concordance)
ggplot(NonDupDup, aes(x = Concordance, fill = Comparison)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 54)

##Compare duplicate pairs with non-duplicate pairs but call rate >0.99 and MAF>0.01
###Add in concordance for duplicates MAF>0.01
ConcordanceDuplicatesMAF <- read.delim("ndiscordanceduplicatesproportionMAF.txt",sep=" ",header=T,col.names=c("nDiscordant","n","Concordance"),row.names=NULL)
ConcordanceDuplicatesMAF$Filepath.1<-DuplicateSamplesGTcheck$Filepath.1
ConcordanceDuplicatesMAF$Filepath.2<-DuplicateSamplesGTcheck$Filepath.2

DuplicateSamplesGTcheckSequenced$nDiscordant<-ConcordanceDuplicatesMAF$nDiscordant
DuplicateSamplesGTcheckSequenced$n<-ConcordanceDuplicatesMAF$n
DuplicateSamplesGTcheckSequenced$Concordance<-ConcordanceDuplicatesMAF$Concordance
hist(DuplicateSamplesGTcheckSequenced$Concordance, breaks=40)


##Add in concordance of non-duplicates MAF>0.01
ConcordanceNonDuplicatesMAF <- read.delim("ndiscordancenonduplicatesproportionMAF.txt",sep=" ",header=T,col.names=c("nDiscordant","n","Concordance"),row.names=NULL)
Random29Samples$nDiscordant<-ConcordanceNonDuplicatesMAF$nDiscordant
Random29Samples$n<-ConcordanceNonDuplicatesMAF$n
Random29Samples$Concordance<-ConcordanceNonDuplicatesMAF$Concordance
hist(Random29Samples$Concordance, breaks=40)
Random29Samples<-Random29Samples %>%
  mutate(Comparison="NonDuplicate")
DuplicateSamplesGTcheckSequenced<-DuplicateSamplesGTcheckSequenced %>%
  mutate(Comparison="Duplicate")

#Now compare concordance of duplicates and non-duplicates in histogram
DuplicatesConcordanceMAF<-data.frame(DuplicateSamplesGTcheckSequenced$Concordance,DuplicateSamplesGTcheckSequenced$Comparison)
NonDuplicatesConcordanceMAF<-data.frame(Random29Samples$Concordance,Random29Samples$Comparison)
colnames(DuplicatesConcordanceMAF)[colnames(DuplicatesConcordanceMAF) == 'DuplicateSamplesGTcheckSequenced.Concordance'] <- 'Concordance'
colnames(DuplicatesConcordanceMAF)[colnames(DuplicatesConcordanceMAF) == 'DuplicateSamplesGTcheckSequenced.Comparison'] <- 'Comparison'
colnames(NonDuplicatesConcordanceMAF)[colnames(NonDuplicatesConcordanceMAF) == 'Random29Samples.Concordance'] <- 'Concordance'
colnames(NonDuplicatesConcordanceMAF)[colnames(NonDuplicatesConcordanceMAF) == 'Random29Samples.Comparison'] <- 'Comparison'
NonDupDupMAF<-rbind(DuplicatesConcordanceMAF,NonDuplicatesConcordanceMAF)
hist(NonDupDupMAF$Concordance)
ggplot(dplyr::filter(NonDupDupMAF, Coverage.1>0.1), aes(x = Concordance, fill = Comparison)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 54)

### LOOKS LIKE I NEED TO DOUBLE CHECK COMPARISON FILE  ###

DuplicateSamplesGTcheckSequencedCoverage<-subset(DuplicateSamplesGTcheckSequenced, DuplicateSamplesGTcheckSequenced$Coverage.1 > 0.1)
write.table(DuplicateSamplesGTcheckSequencedCoverage, file = "DuplicateSamplesGTcheck.txt", sep = "\t",
            col.names = F, row.names = F, quote = FALSE)

###How does coverage affect concordance?
DuplicateSamplesGTcheckSequenced$MeanCoverage<-mean(DuplicateSamplesGTcheckSequenced$Coverage.1,DuplicateSamplesGTcheckSequenced$Coverage.2)
CovxCon<-ggplot(dplyr::filter(DuplicateSamplesGTcheckSequenced) , aes(x=Coverage.1, y=Concordance)) + geom_point() + geom_smooth(method='lm')
CovxCon

##How does the plate it is sequenced on affect concordance?
sum(DuplicateSamplesGTcheckSequenced$Plate.1==DuplicateSamplesGTcheckSequenced$Plate.2)
sum(DuplicateSamplesGTcheckSequenced$Plate.1!=DuplicateSamplesGTcheckSequenced$Plate.2)
DuplicateSamplesGTcheckSequenced$SamePlate<-ifelse(DuplicateSamplesGTcheckSequenced$Plate.1==DuplicateSamplesGTcheckSequenced$Plate.2, "Same", "Different")
PlatexCon<-ggplot(dplyr::filter(DuplicateSamplesGTcheckSequenced) , aes(x=SamePlate, y=Concordance)) + geom_point() + geom_smooth(method='lm')
PlatexCon

#How does concordance affect FROH?
ConxFROH<-ggplot(dplyr::filter(DuplicateSamplesGTcheckSequenced) , aes(x=FROH.1, y=FROH.2, color=Concordance)) + geom_point() + geom_smooth(method='lm')
ConxFROH

#Create file of duplicates for gtcheck to check concordance
DuplicateExtraSamplesGTcheck<-subset(SequencedIndividualsBirdIDsExtraDuplicatesWide, select = c(Filepath.1, Filepath.2))
write.table(DuplicateExtraSamplesGTcheck, file = "DuplicateExtraSamplesGTcheck.txt", sep = "\t",
            col.names = F, row.names = F, quote = FALSE)

### Additional concordance tests using SnpSift- Looks great! ###
snpsiftconcordance<-read.table("concordance_773truth_773testreheadered.by_sample.txt",header=T,check.names = FALSE)
snpsiftconcordance$concordance<-(snpsiftconcordance$`REF/REF`+snpsiftconcordance$`ALT_1/ALT_1`+snpsiftconcordance$`ALT_2/ALT_2`)/(snpsiftconcordance$`REF/REF`+snpsiftconcordance$`ALT_1/ALT_1`+snpsiftconcordance$`ALT_2/ALT_2`+snpsiftconcordance$`REF/ALT_1`+snpsiftconcordance$`REF/ALT_2`+snpsiftconcordance$`ALT_1/REF`+snpsiftconcordance$`ALT_1/ALT_2`+snpsiftconcordance$`ALT_2/REF`+snpsiftconcordance$`ALT_2/ALT_1`)
snpsiftconcordance$nsites<-(snpsiftconcordance$`REF/REF`+snpsiftconcordance$`ALT_1/ALT_1`+snpsiftconcordance$`ALT_2/ALT_2`+snpsiftconcordance$`REF/ALT_1`+snpsiftconcordance$`REF/ALT_2`+snpsiftconcordance$`ALT_1/REF`+snpsiftconcordance$`ALT_1/ALT_2`+snpsiftconcordance$`ALT_2/REF`+snpsiftconcordance$`ALT_2/ALT_1`)
mean(snpsiftconcordance$concordance)

