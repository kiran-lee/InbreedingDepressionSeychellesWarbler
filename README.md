# InbreedingDepressionSeychellesWarbler
 
These are the data and scripts used to investigate individual inbreeding depression using ~1900 low coverage, whole-genome Seychelles warbler sequences that have been imputed to recover missing genotypes. The pipeline to impute SNPs and verify accuracy can be found in my other repository: https://github.com/kiran-lee/SNPsSeychellesWarbler/tree/main and the imputed SNPs can be found as a compressed vcf file on Zenodo: https://zenodo.org/records/12570527. 

An explanation of the files in Data can be found in the InbreedingDepressionSeychellesWarblers.R script.

Perhaps the most useful file to other researchers is a dataframe (Seychelles_warbler_traits.txt) of all birds in the Seychelles warbler database (BirdID), the plate number they were sequenced on (Plate), sequences sample name as named by Liverpool University (SeqID),  sequencing coverage (Coverage), tube number of sample used to link SeqID to BirdID  (ID), whether the ID used was BloodID or BloodTubeNumber (Identifier) FROH inbreeding coefficient (FROH),  lifespan (Lifespan), year of birth (BirthYear), year it was last seen (LastSeenYear) and lifetime offspring produced (ReproductiveOutput).
