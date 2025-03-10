if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("maftools")
library(maftools)
install.packages("readxl")
library(readxl)

cagekid.maf = read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_mut') 
cagekid.clin = read_excel('/Users/jessie/Desktop/Cagekid/Cagekid_All.xlsx', sheet = 'Cagekid_clin') 
all = read.maf(maf = cagekid.maf, clinicalData = cagekid.clin)

# oncoplots for all samples
oncoplot(maf = all, top = 18)

#VHLwt
VHLwt_clin <- cagekid.clin[cagekid.clin$VHL == "0", ]
VHLwt_maf_data <- cagekid.maf[cagekid.maf$Tumor_Sample_Barcode %in% VHLwt_clin$Tumor_Sample_Barcode, ]
VHLwt <- read.maf(maf = VHLwt_maf_data, clinicalData = VHLwt_clin)
oncoplot(maf = VHLwt, top = 11)

#VHLmut
VHLmut_clin <- cagekid.clin[cagekid.clin$VHL == "1", ]
VHLmut_maf_data <- cagekid.maf[cagekid.maf$Tumor_Sample_Barcode %in% VHLmut_clin$Tumor_Sample_Barcode, ]
VHLmut <- read.maf(maf = VHLmut_maf_data, clinicalData = VHLmut_clin)
oncoplot(maf = VHLmut, top = 15)


# Get the top 15 mutated genes excluding VHL
top_genes <- getGeneSummary(VHLmut)$Hugo_Symbol[1:16]  # Get top 16 to adjust for exclusion
top_genes <- top_genes[top_genes != "VHL"]  # Remove VHL
# Generate oncoplot without displaying VHL
oncoplot(maf = VHLmut, genes = top_genes)




#VHL+1
VHL_1_clin <- cagekid.clin[cagekid.clin$VHLCat == "VHLmut+3", ]
VHL_1_maf_data <- cagekid.maf[cagekid.maf$Tumor_Sample_Barcode %in% VHL_1_clin$Tumor_Sample_Barcode, ]
VHL_1 <- read.maf(maf = VHL_1_maf_data, clinicalData = VHL_1_clin)
oncoplot(maf = VHL_1, top = 18)




oncoplot(maf = entire, top = 15)
oncoplot(maf = VHLwt, top = 15)
oncoplot(maf = VHLmut, top = 15)
oncoplot(maf = VHL0, top = 15)
oncoplot(maf = VHL1, top = 15)
oncoplot(maf = VHL2, top = 15)
oncoplot(maf = VHL3, top = 15)
















