# InbreedingDepressionSeychellesWarbler
 
These are the data and scripts used to investigate inbreeding depression post-bottlenek in Seychelles warblers using 1,936 low coverage, whole-genome sequences that have been imputed to recover missing genotypes and verified for sample identity as in https://github.com/kiran-lee/SeychellesWarblerGenomicToolkit. 

![Cousin 1900s](Other/File-photo-Cousins-coconut-plantations-were-begun-in-the-early-1900s-thegem-blog-default.jpg)

Photo from Nature Seychelles.

We trace the demographic history, showing a drastic bottleneck to Ne = 12.6, 9 generations ago.
We use model-based approach to calculate individual FROH using RZooRoH and investigate how this covaries with key fitness traits:  annual survival, annual fecundity, lifespan and lifetime fecundity, including relevant variables showing severe inbreeding depression.
We perform a GWAS-like study investigating effects of SNPs within ROH on fitness traits.
We evaluate if any inbreeding avoidance by mate choice takes place.
We show evidence this research can be applied in choosing which surviving, genotyped individuals to translocate.

Perhaps the most useful file to other researchers is a dataframe (frohdata.xlsx) of all birds in the Seychelles warbler database (BirdID), the plate number they were sequenced on (Plate), sequences sample name as named by Liverpool University (SeqID),  sequencing coverage (Coverage), tube number of sample used to link SeqID to BirdID  (ID), whether the ID used was BloodID or BloodTubeNumber (Identifier) FROH inbreeding coefficient estimated from PLINK (FROH_all), and inbreeding estimated by RZooRoH (FROH_2.5, FROH_5, FROH_25, _FROH_50, FROH_250, FROH_500, _FROH_2500), lifespan (Lifespan), year of birth (BirthYear), year it was last seen (LastSeenYear), lifetime offspring produced (ReproductiveOutput), natal year death rate as an environmental control of lifetime hardship faced (BirthYearDeathRate), age of mother at birth (MaternalAge).
