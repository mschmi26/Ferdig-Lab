#install.packages("dplyr")
library(dplyr)
library(tidyverse)

#read in data
setwd('/Users/macke/Documents/RSA/MKK2835xNHP1337')
eastman  = read.csv("Malaria Parental Lines.csv", 
                          header = TRUE, 
                          #row.names = 1, 
                          stringsAsFactors = FALSE, 
                          fileEncoding="UTF-8-BOM")
head(eastman)
#eastman = tibble(eastman)
#subest eastman data to the important columns 
eastman_subset = eastman %>%
  select(Protocol.Name,
         Sample.Name,
         Sample.ID,
         Primary.MOA,
         IC50..uM.,
         AC50..uM.,
         AUC,
         Hill.Coef,
         CC.v2,
         Efficacy,
         Structure)

#replace blanks with NA
eastman_subset[eastman_subset == ""] = NA
#eastman_subset[eastman_subset == "null"] = NA
# Loop over rows of the rows in the eastman dataframe
for (row in 1:nrow(eastman_subset)) {
  #make the important varablies in the current row their own seperate value used in the if statements
  curve <- eastman_subset[row, "CC.v2"]
  efficacy  <- eastman_subset[row, "Efficacy"]
  para  <- eastman_subset[row, "Protocol.Name"]
  drug  <- eastman_subset[row, "Sample.Name"]
  NCATS_ID <- eastman_subset[row, "Sample.ID"]
  # if curve class of the current row is an NA, OR curve class is not -1.1, -1.2, -2.1, -2.2 OR effeicacy is greater than -80 then make ac50 ic50 and hill coef an NA in the dataframe
  if(is.na(curve) == TRUE || curve %in% c(-1.1, -1.2, -2.1, -2.2) == FALSE || efficacy > -80) {
    eastman_subset[row, "AC50..uM."] = NA
    eastman_subset[row, "IC50..uM."] = NA #value in eastman_subset dataframn at current row in column "IC50_uM_" is turned to NA
    eastman_subset[row, "Hill.Coef"] = NA
  } #if there isn't a normal drug name in Sample.ID, make sample.ID the NCATS identifier 
  if (drug == "null") { 
    eastman_subset[row, "Sample.Name"] = NCATS_ID
  }
}
#replace nulls with NA
eastman_subset[eastman_subset == "null"] = NA

#group data by parasite
eastman_grouped = eastman_subset %>% group_by(Protocol.Name) 
#split the dateframe into seperate dataframes by the groups just made
eastman_splits = group_split(eastman_grouped)




#make a list of parasite names
parasites = c("mal31", "KH004", "MKK2835", "NF54", "NHP1337", "NHP4026")
#loop through tibbles seperated by drug. add drug name to the end of AC_50_, IC_50_, AUC_, Hill_Coef_
drug_list = list()
for (tib in 1:length(eastman_splits)){
  ###add drug to mapping columns
  AC50_uM_tmp = (paste(parasites[tib], "_AC50_uM", sep = ''))
  IC50_uM_tmp = (paste(parasites[tib], "_IC50_uM", sep = ''))
  AUC_tmp = (paste(parasites[tib], "_AUC", sep = ''))
  Hill_Coef_tmp = (paste(parasites[tib], "_Hill_Coef", sep = ''))
  CC_v2_tmp = (paste(parasites[tib], "_CC_v2", sep = ''))
  Efficacy_tmp = (paste(parasites[tib], "_Efficacy", sep = ''))
  

  current = eastman_splits[[tib]] %>% rename(!! AC50_uM_tmp := AC50..uM.,
                                             !! IC50_uM_tmp := IC50..uM.,
                                             !! AUC_tmp := AUC,
                                             !! Hill_Coef_tmp := Hill.Coef,
                                             !! CC_v2_tmp := CC.v2,
                                             !! Efficacy_tmp := Efficacy)
  
  drug_list[[tib]] <- current;drug_list
}

#phenotypes_df = cbind(phenotypes_df, drug_list[[tibs]][,5:8])

# for each parasite dataframe remove rows with an "NA"
mal31_df = subset(drug_list[[1]], mal31_IC50_uM != "NA") %>% select(-1, -6, -7, -8)
KH004_df = subset(drug_list[[2]], KH004_IC50_uM != "NA") %>% select(-1, -6, -7, -8)
MKK2835_df = subset(drug_list[[3]], MKK2835_IC50_uM != "NA") %>% select(-1, -6, -7, -8)
NHP1337_df = subset(drug_list[[5]], NHP1337_IC50_uM != "NA") %>% select(-1, -6, -7, -8)
NF54GFP = subset(drug_list[[5]], NHP1337_IC50_uM != "NA") %>% select(-1, -6, -7, -8)


MKK2835_NHP1337_df = merge(MKK2835_df, NHP1337_df, by=c("Sample.Name","Sample.ID","Primary.MOA","Structure")) %>%
  add_column(fold_change=NA, .after = "Primary.MOA") %>%
  add_column(uM_difference=NA, .after = "fold_change")
MKK2835_NHP1337_df = MKK2835_NHP1337_df %>% mutate(fold_change = pmax(as.numeric(MKK2835_NHP1337_df$NHP1337_IC50_uM), as.numeric(MKK2835_NHP1337_df$MKK2835_IC50_uM))/pmin(as.numeric(MKK2835_NHP1337_df$NHP1337_IC50_uM), as.numeric(MKK2835_NHP1337_df$MKK2835_IC50_uM)))
MKK2835_NHP1337_df = MKK2835_NHP1337_df %>% mutate(uM_difference = pmax(as.numeric(MKK2835_NHP1337_df$NHP1337_IC50_uM), as.numeric(MKK2835_NHP1337_df$MKK2835_IC50_uM))-pmin(as.numeric(MKK2835_NHP1337_df$NHP1337_IC50_uM), as.numeric(MKK2835_NHP1337_df$MKK2835_IC50_uM)))




write.csv(MKK2835_NHP1337_df, "MKK2835_NHP1337_df.csv", row.names = F)






###add phenotype columns to new data frame
phenotypes_df <- data.frame(matrix(ncol = 0, nrow = 68))
names = drug_list[[1]][,1]
phenotypes_df = data.frame(names[,-1], row.names=names$Protocol.Name)

for (tibs in 1:length(drug_list)){
  if(nrow(drug_list[[tibs]]) != 68){  
    print(tibs)
    print(nrow(drug_list[[tibs]]))
  }  
  phenotypes_df = cbind(phenotypes_df, drug_list[[tibs]][,5:8])
}

write.csv(phenotypes_df, "all_eastman_phenotypes.csv")
