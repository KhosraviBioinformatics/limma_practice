"""
@author: AmirMohammad Khosravi
@Date: 1404\03\16 
@description: practicing limma
"""

setwd("C:/Users/khosravi/OneDrive/Desktop/DEGproject")

#Packages
install.packages("BiocManager")
BiocManager::install("limma")

#read data
data_file <- read.csv("T2vsT1finalRsyn1.csv")
View(data_file)
colnames(data_file)
print(which(grepl("Tt", names(data_file))))
print(which(grepl("To", names(data_file))))

##############################
#Data manipulation
##############################

# Reorder the data base on control & case(controls come first):
# Find column indices
genes <- which(grepl("X", names(data_file)))
to_cols <- which(grepl("To", names(data_file)))
tt_cols <- which(grepl("Tt", names(data_file)))

# Combine the indices in the desired order
new_order <- c(genes,to_cols, tt_cols)

# Reorder the dataframe
#data_reordered <- cbind(gene = genes, data_file)
data_reordered <- data_file[, new_order]
View(data_reordered)
print(which(grepl("Tt", names(data_reordered))))
print(which(grepl("To", names(data_reordered))))

#check if having a duplicate of gene symbols or not:
anyDuplicated(data_reordered[[1]])
View(data_reordered[[1]])
# Save the ordered data:
write.csv(data_reordered, "data_reordered.csv", row.names = FALSE)

##############################
# Limma
##############################
# Load the package
library(limma)

# Extract expression matrix (remove the first column which contains gene names)
expr_data <- as.matrix(data_reordered[, -1])
rownames(expr_data) <- data_reordered[[1]]
# Check for 0 data
any(expr_data == 0) 
# Define sample groups
# Suppose first 'n' columns are control (To), remaining are case (Tt)
num_controls <- length(which(grepl("To", names(data_reordered))))
num_cases <- length(which(grepl("Tt", names(data_reordered))))

group <- factor(c(rep("Control", num_controls), rep("case", num_cases)))
group
design <- model.matrix(~0+group)
View(design)
colnames(design) <- c("case" , "Control")

# Fit linear model
fit <- lmFit(expr_data, design)
cont <- makeContrasts(contrasts = "case - Control" , levels = design)
fit2 <- contrasts.fit(fit = fit , contrasts = cont)
fit3 <- eBayes(fit2)

# Check top DEGs
top_table <- topTable(fit3, number=Inf, adjust.method="fdr")
View(top_table)

# Save results
write.csv(top_table, "limma_results.csv")

