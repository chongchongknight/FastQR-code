
################################################################################
## important!!!!!!!
## For the use of data please see the instructions: https://www.cancerrxgene.org/, 
## All the user should have the permission before usage.
## Also see the paper: Iorio, Francesco, et al. "A landscape of pharmacogenomic interactions in cancer." Cell 166.3 (2016): 740-754.
################################################################################











library(data.table)




# Y -drug response
drug=fread("v17.3_fitted_dose_response.csv")
dim(drug)
# [1] 224202     13
dim(drug[which(drug$DRUG_NAME=="Gefitinib"),])

table(drug[which(drug$PUTATIVE_TARGET=="EGFR"), "DRUG_NAME"])
# Cetuximab Erlotinib Gefitinib Pelitinib 
# 899       392       878       976 

# X-RNA
RNA=fread("sanger1018_brainarray_ensemblgene_rma.txt", header=T)
dim(RNA)
# [1] 17737  1019 (n=1018)

# Z- tissue types
cov=fread("Cell_Lines_Details_CSV.csv")
dim(cov)
# [1] 1002   13
table(cov[, 8]) # tissue 
# aero_dig_tract              bone            breast  digestive_system 
# 1                79                40                52                52 
# kidney   large_intestine          leukemia              lung        lung_NSCLC 
# 34                51                85                22               111 
# lung_SCLC          lymphoma           myeloma    nervous_system     neuroblastoma 
# 66                70                18                57                32 
# pancreas              skin       soft_tissue           thyroid urogenital_system 
# 32                58                21                16               105 

cov=cov[which(cov[,8]!=""),]
# [1] 1001   13

# length(table(drug$DRUG_ID))
# [1] 266 # 226 different drugs; for each drug, run analysis
# range(table(drug$DRUG_ID))
# [1] 388 989 # the sample sizes for each drug. 



# Example - drug 1

# ID run for each drug in DRUG_ID and for each gene  - "I updated drugIDindex for y"
drugIDindex=1530

drugsub=drug[which(drug$DRUG_ID==drugIDindex),]
y_COSMIC_ID=data.frame(drugsub[, "COSMIC_ID"])[,1]
x_COSMIC_ID=colnames(RNA)[-1] 
z_COSMIC_ID=data.frame(cov[,2])[,1]
xyz_ID=intersect(intersect(y_COSMIC_ID, x_COSMIC_ID), z_COSMIC_ID)

# X, Y, Z data
y=data.frame(drugsub[match(xyz_ID, drugsub$COSMIC_ID), "LN_IC50"])[,1]
x=data.frame(RNA[,match(xyz_ID, colnames(RNA)), with=F])
z=data.frame(cov[match(xyz_ID, as.character(cov$'COSMIC identifier')), 8])[,1]

# For each row of x, run analysis to obtain p-value; Here is where you revise to LP procedure
i=1
summary(lm(y~as.numeric(x[i,])+as.factor(z)))$coef[2,4]

unique(z)

lm(y~as.numeric(x[i,])+as.factor(z))


unique(z)



### the above code show how we clean the data, we also give a pseudo data for the example

#data = cbind(y, z, t(x))

#dim(x) 17737   835, we choose 5000 column
#length(z) 835
#length(y) 835


data = read.csv("pseudodata.csv")[, -1]


## For the implementation of the data please choose one of drug and follow the 
## extraction process above to use the corresponding package
