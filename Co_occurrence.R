library(maftools)
install.packages("readxl")
library(readxl)

#Input data
maf = read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_mut') 
clin = read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin') 

#Filtering
clin <- clin[clin$VHL == "1", ]

data <- maf[maf$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode, ]
read <- read.maf(maf = data, clinicalData = clin)
plot <- somaticInteractions(
  maf = read,
  countStats = TRUE,
  top = 100, 
  pvalue = c(0.05, 0.1),
  fontSize = 0.7
)


