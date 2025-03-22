###### loading required libraries 
library(data.table)
library(MOFA2)

######## Load data from files
file1 <- "GSE168404_miRNA_C-vs-P.all.txt"
file2 <- "GSE168404_mRNA_C-vs-P.all.txt"
file3 <- "GSE168404_lncRNA_C-vs-P.all.txt"

######### Read the data from the files
data1 <- fread(file1)
data2 <- fread(file2)
data3 <- fread(file3)

######## Remove the second column from data2 as two identifiers are there
data2 <- data2[, -2, with = FALSE]

######## Function to set row names and convert to matrix
set_rownames <- function(df) {
  df <- as.data.frame(df)  # Ensure it's a data frame
  rownames(df) <- df[, 1]  # Set first column as row names
  df <- df[, -1]  # Remove the first column
  return(as.matrix(df))  # Convert to matrix
}

######### Apply function to each dataset
data1 <- set_rownames(data1)
data2 <- set_rownames(data2)
data3 <- set_rownames(data3)

####### checking for duplicated in identifiers 

check_duplicates <- function(mat) {
  duplicates <- duplicated(rownames(mat))
  if (any(duplicates)) {
    cat("Duplicates found in row names:\n")
    print(rownames(mat)[duplicates])
  } else {
    cat("No duplicates found in row names.\n")
  }
}

######### Check for duplicates in each dataset
cat("Checking data1:\n")
check_duplicates(data1)

cat("\nChecking data2:\n")
check_duplicates(data2)

cat("\nChecking data3:\n")
check_duplicates(data3)

###### Normalization 

data_logged1.logged =log2(data1 + 1)
data_logged2.logged =log2(data2 + 1)
data_logged3.logged =log2(data3 + 1)


################ scaling 
data1_scaled=scale(data_logged1.logged, center=TRUE , scale=TRUE)
data2_scaled=scale(data_logged2.logged, center=TRUE , scale=TRUE)
data3_scaled=scale(data_logged3.logged, center=TRUE , scale=TRUE)

########## Boxplots for Raw, Logged, and Scaled Data

par(mfrow=c(1,3)) # Arrange plots in a row

boxplot(data1, main="Raw miRNA", outline=FALSE)
boxplot(data_logged1.logged, main="Log2 Normalized miRNA", outline=FALSE)
boxplot(data1_scaled, main="Scaled & normalized miRNA", outline=FALSE)

par(mfrow=c(1,3)) 

boxplot(data2, main="Raw mRNA", outline=FALSE)
boxplot(data_logged2.logged, main="Log2 Normalized mRNA", outline=FALSE)
boxplot(data2_scaled, main="Scaled & normalized mRNA", outline=FALSE)

par(mfrow=c(1,3))

boxplot(data3, main="Raw lncRNA", outline=FALSE)
boxplot(data_logged3.logged, main="Log2 Normalized lncRNA", outline=FALSE)
boxplot(data3_scaled, main="Scaled & normalized lncRNA", outline=FALSE)

########## Histograms for Data Distribution

par(mfrow=c(1,3)) 

hist(as.vector(data1), breaks=50, main="Raw miRNA", xlab="Expression", col="lightblue")
hist(as.vector(data_logged1.logged), breaks=50, main="Log2 Transformed miRNA", xlab="Expression", col="lightgreen")
hist(as.vector(data1_scaled), breaks=50, main="Scaled & normalized miRNA", xlab="Expression", col="lightcoral")

par(mfrow=c(1,3)) 

hist(as.vector(data2), breaks=50, main="Raw mRNA", xlab="Expression", col="lightblue")
hist(as.vector(data_logged2.logged), breaks=50, main="Log2 Transformed mRNA", xlab="Expression", col="lightgreen")
hist(as.vector(data2_scaled), breaks=50, main="Scaled & normalized mRNA", xlab="Expression", col="lightcoral")

par(mfrow=c(1,3)) 

hist(as.vector(data3), breaks=50, main="Raw lncRNA", xlab="Expression", col="lightblue")
hist(as.vector(data_logged3.logged), breaks=50, main="Log2 Transformed lncRNA", xlab="Expression", col="lightgreen")
hist(as.vector(data3_scaled), breaks=50, main="Scaled & normalized lncRNA", xlab="Expression", col="lightcoral")


###### Adding the data sets to one list

data_list <- list(
  miRNA = data1_scaled,  
  mRNA = data2_scaled,  
  lncRNA = data3_scaled  
)

######## selecting top variable features foe each data set

highly_variable_features <- apply(data_list$mRNA, 1, var)
selected_features <- names(sort(highly_variable_features, decreasing = TRUE))[1:10000]  # Keep top 
data_list$mRNA <- data_list$mRNA[selected_features, ]

highly_variable_features <- apply(data_list$lncRNA, 1, var)
selected_features <- names(sort(highly_variable_features, decreasing = TRUE))[1:10000]  # Keep top 
data_list$lncRNA <- data_list$lncRNA[selected_features, ]

highly_variable_features <- apply(data_list$miRNA, 1, var)
selected_features <- names(sort(highly_variable_features, decreasing = TRUE))[1:2000]  # Keep top 
data_list$miRNA <- data_list$miRNA[selected_features, ]

######## MOFA object with filtered data

MOFAobject <- create_mofa(data_list)

####### Visualize data structure
plot_data_overview(MOFAobject)

########## Data options
data_opts <- get_default_data_options(MOFAobject)
#data_opts$scale_views <- FALSE  # Data is already Z-transformed

######### Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

######### Training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$maxiter <- 10000  # More iterations for better convergence
train_opts$convergence_mode <- "slow"  # More precise optimization
train_opts$seed <- 42

######## Prepare MOFA object
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

######### Train the model
MOFAobject <- run_mofa(MOFAobject, use_basilisk = TRUE)

######## Save the model
outfile <- paste0(getwd(), "/model.hdf5")
saveRDS(MOFAobject, outfile)

######## Plot the variance explained per factor
plot_variance_explained(MOFAobject)

######## Plot the variance explained per factor
plot_variance_explained(MOFAobject, plot_total = TRUE)
plot_factor_cor(MOFAobject)


#####ANOVA
factors <- get_factors(MOFAobject)$group1  # Extract factors
factors_df <- as.data.frame(factors)  # Convert to data frame
factors_df$Sample <- rownames(factors_df)  # Add sample IDs as a column

metadata <- read.csv("MetaDatagrad.csv", row.names = 1)
head(metadata)

# Ensure row names are properly aligned
merged_data <- merge(metadata, factors_df, by = "row.names")  

# Rename the new column created from row names
colnames(merged_data)[1] <- "Sample"

factors_df <- na.omit(factors_df)
table(merged_data$group)  # Check unique group values
merged_data$group <- as.factor(merged_data$group)  # Ensure it's a factor
merged_data <- merged_data[!is.na(merged_data$group), ]
zero_var_factors <- apply(factors_df, 2, var) == 0
factors_df <- factors_df[, !zero_var_factors]  # Remove zero-variance factors


anova_results <- lapply(colnames(factors_df), function(factor) {
  model <- aov(as.formula(paste(factor, "~ group")), data = merged_data)
  summary(model)
})



