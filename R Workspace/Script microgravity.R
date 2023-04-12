#########
# Notes #
#########

# In this experiment, we have the following statistical parameters of interest:
# - Explanatory variable (variable that influences gene expression. We can have
#   explanatory variables that are covariates - they are numerical and
#   continuous- or factors - which are discrete, as in this case. The different
#   values that a factor can take are called 'levels'): Our explanatory variable
#   is the the gravity condition, and it has three levels - `Earth`,
#   `Space_0g` and `Space_1g`. For a specific library, these levels are either
#   0 or 1 depending on whether that library comes from a cell culture in the
#   respective gravity condition.
# - Design matrix (used to define the structure of the statistical model. It
#   stores the value of the explanatory variables for each library): Named
#   `design` in this code, it has the shape:
#                Earth   Space_0g Space_1g
#            1     0        1        0
#            2     0        1        0
#            3     0        1        0
#            4     0        0        1
#            5     0        0        1
#            6     0        0        1
#            7     1        0        0
#            8     1        0        0
#            9     1        0        0
#   This is the matrix that is used in order to find the parameters of the linear
#   model (in this case, the 'means model') that fit to our data. Fitting the model
#   means solving the equation:
#                      DESIGN_MATRIX * PARAMETER_MATRIX = GENE_DATA,
#   Where the PARAMETER_MATRIX is the matrix that stores in its single column
#   (for this case, as we only have one explanatory variable to study) the values
#   the model parameters that we want to find.
#
#   Note that in this case we are using a design matrix without an intercept term
#   (the model uses a pure combination of the values of the explanatory variable
#   with no constant additive term). However, as our variable is a factor, the
#   resultant model is equivalent to the model we would obtain if we were to use
#   a design matrix with an intercept term but using a 'means-reference model'
#   instead of the 'means model', there the reference would be one of the levels
#   of the explanatory variable (`Earth`, for example).
# - Contrast matrix: once we have fitted our model, we obtain the model parameters.
#   These are like 'weights' can give us an idea of how the gene expression is in
#   each of the levels of the explanatory variable (actually, the model parameters
#   in this case are exactly the average expression of each gene for each of the
#   levels of the explanatory variable). However, sometimes it is more useful to
#   study the difference between the gene expression in the different levels, more
#   than the absolute values of each level separately. The contrast matrix specifies
#   linear operations that should be performed between the model parameters to
#   calculate specific values of interest. A single contrast matrix can contain
#   many of these linear combinations of interest.
#   In our case, the contrast matrix contains one contrast (named 'microgravity')
#   which gives the difference `Space_0g - Earth`:
#                     Levels     microgravity
#                       Earth              -1
#                       Space_0g            1
#                       Space_1g            0
#
# After computing the contrast, we find the following data:
#                              microgravity
#                       Down           1487
#                       NotSig        13585
#                       Up             1053
# It means that the contrast defined in the contrast matrix is significantly
# negative in 1487 genes (these genes are down regulated in Space_0g),
# significantly positive in 1053 genes (these genes are up regulated in
# Space_0g) and not significantly different in 13585 genes.

#########
# Setup #
#########

# install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", repos=NULL, type="source")
BiocManager::install("edgeR")
BiocManager::install("Glimma")
install.packages("gplots")
BiocManager::install("org.Hs.eg.db")
install.packages("RColorBrewer")

# load installed packages
library(limma)
library(edgeR)
library(Glimma)
library(gplots)
library(org.Hs.eg.db)
library(RColorBrewer)

# Use the setwd() function to set your working directory.
getwd() # function to find out what your current working directory is.
setwd("/home/bernat/Desktop/Uni/Third\ Year/Biophysics\ II/Project/Response\ to\ Microgravity/R\ Workspace") # function to set your working directory (path).
getwd() # after you've set the working directory, verify it by calling the getwd() function.

# Prepare count matrix
seqdata <- read.delim("/home/bernat/Desktop/Uni/Third\ Year/Biophysics\ II/Project/Response\ to\ Microgravity/R\ Workspace/GSE188793_counts.txt", stringsAsFactors = FALSE)
sampleinfo <- read.delim("/home/bernat/Desktop/Uni/Third\ Year/Biophysics\ II/Project/Response\ to\ Microgravity/R\ Workspace/SampleInfo.txt")
head(seqdata)
sampleinfo
# Remove first column from seqdata
countdata <- seqdata[, -(1:1)]
# Store GeneID as rownames
rownames(countdata) <- seqdata[, 1]
# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata), start = 1, stop = 7)
head(countdata)

##################
# Data Filtering #
##################
# Remove low expression genes

# Obtain CPMs (counts per million: https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/cpm)
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)
# Summary of how many TRUEs there are in each row
# There are 12787 genes that have TRUEs in all samples.
sum(rowSums(thresh) == ncol(thresh))
table(rowSums(thresh))
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 3
# Make a subset of the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep, ]
summary(keep)
dim(counts.keep)
# Let's have a look and see whether our threshold of 0.5 does indeed correspond
# to a count of about 5-10 (We will look at the first sample)
plot(myCPM[, 1], countdata[, 1], ylim = c(0, 50), xlim = c(0, 3))
abline(v = 0.5)

###################
# Quality Control #
###################

# Convert counts to DGEList object
y <- DGEList(counts.keep)
# Library size information is stored in the samples slot
y$samples
y$samples$lib.size
# Represent plot of logCPM
barplot(y$samples$lib.size, names = colnames(y), las = 2)
# Add a title to the plot
title("Barplot of library sizes")

# Get log2 counts per million
logcounts <- cpm(y, log = TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab = "", ylab = "Log2 counts per million", las = 2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h = median(logcounts), col = "blue")
title("Boxplots of logCPMs (unnormalised)")

# Principal Component Analysis
# Interactive Glimma PCA Plot
labels <- paste(sampleinfo$SampleName, sampleinfo$Gravity)
group <- paste(sampleinfo$Gravity, sep = ".")
group <- factor(group)
glMDSPlot(y, labels = labels, groups = group, folder = "mds")

# Hierarchical clustering
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing = TRUE))[1:500]
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var, ]
## Get some nicer colours
mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple", "orange")[sampleinfo$Gravity]
# Plot the heatmap
heatmap.2(highly_variable_lcpm, col = rev(morecols(50)), trace = "none", main = "Top 500 most variable genes across samples", ColSideColors = col.cell, scale = "row")
# Save the heatmap
png(file = "High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm, col = rev(morecols(50)), trace = "none", main = "Top 500 most variable genes across samples", ColSideColors = col.cell, scale = "row")
dev.off()

########################
# Sample Normalization #
########################
# Apply normalization to DGEList object
y <- calcNormFactors(y)
y$samples
# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group) # See notes section
# Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
design
# Transform data
par(mfrow = c(1, 1))
v <- voom(y, design, plot = TRUE)
v

# represent boxplot normalized vs non normalized
par(mfrow = c(1, 2))
boxplot(logcounts, xlab = "", ylab = "Log2 counts per million", las = 2, main = "Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h = median(logcounts), col = "blue")
boxplot(v$E, xlab = "", ylab = "Log2 counts per million", las = 2, main = "Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h = median(v$E), col = "blue")

####################################
# Differential Expression Analysis #
####################################
# In this section, we study the variation in the gene expression of the different
# libraries that we have (3 for ISS 1G, 3 for ISS 1G and 3 for ground 1G). We want to
# find genes that are over-expressed or under-expressed in microgravity.

# Fit the linear model
fit <- lmFit(v)
names(fit)
cont.matrix <- makeContrasts(microgravity = Space_0g - Earth, levels = design)
cont.matrix

# Compute Contrast (find up/down/insignificantly-regulated genes)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont) # Compute t-statistics (among others)
dim(fit.cont)
summa.fit <- decideTests(fit.cont) # Classify as up, down or not significant according to t-statistics.
summary(summa.fit)

# Sort genes
topTable(fit.cont, coef = "microgravity", sort.by = "p")

# Annotation
columns(org.Hs.eg.db)
ann <- select(org.Hs.eg.db, keys = rownames(fit.cont), columns = c("ENTREZID", "SYMBOL", "GENENAME"), keytype = "ENSEMBL")
ann <- ann[!duplicated(ann$ENSEMBL), ] # Make sure there are no duplicate annotations
head(ann)
dim(ann)
table(ann$ENSEMBL == rownames(fit.cont))
fit.cont$genes <- ann
topTable(fit.cont, coef = "microgravity", sort.by = "p")
limma.res <- topTable(fit.cont, coef = "microgravity", sort.by = "p", n = "Inf")

# Export to CSV
write.csv(limma.res, file = "microgravity.csv", row.names = FALSE)

# Interactive Glimma Volcano Plot
group2 <- group
levels(group2)
glXYPlot(
  x = fit.cont$coefficients[, 1], y = fit.cont$lods[, 1],
  xlab = "logFC", ylab = "B", main = "microgravity",
  counts = y$counts, groups = group2, status = summa.fit[, 1],
  side.main = "ENSEMBL", folder = "volcano"
)
