# Install and load packages
install.packages("readxl")
install.packages("maftools")
library(readxl)
library(maftools)


# Input data
maf = read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_mut') 
clin = read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin') 

#Filtering by VHL Status
clin <- clin[clin$VHL == "1", ]

# Read maf file
data <- maf[maf$Tumor_Sample_Barcode %in% clin$Tumor_Sample_Barcode, ]
read <- read.maf(maf = data, clinicalData = clin)

#Generate the plot
plot <- somaticInteractions(
  maf = read,
  countStats = TRUE,
  top = 100, 
  pvalue = c(0.05, 0.1),
  fontSize = 0.7
)