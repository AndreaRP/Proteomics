# Heatmap of Ad-hoc Pathways:
# Get the list of DE genes of those pathways
base_dir <- "./ANALYSIS/03.CustomPathways/"
# Then we tell the program about the pathways that we have.
list <- list.files(base_dir, pattern="*list.txt")
list
pathways <- c("Calcium", "ECM", "Hemo", "Mitochondrial", "ROS and DNA Binding", "Sarcomeric Cytoskeleton")
#pathways.1 <- c("Calcium", "ECM", "Hemo", "ROS and DNA Binding", "Sarcomeric Cytoskeleton")
prots <- data.frame(pathway=pathways, path=file.path(base_dir, list), stringsAsFactors = FALSE)
prots
# Search in each list the top 3 significant fdr in D21/D0 contrast
# Select the genes with fdr <0.05 in all the dataset
completeTableD21_WT <- topTable(fit2,coef=2,number=nrow(exprs), adjust="fdr")

# Define function to extract top 5 proteins by absolute logFC of the String pathway
# it is passed (eg.: Mitochondria)
getTopProts <- function(pathway, n){
  # Read the pathway list of proteins
  pathway.prot.names <- scan(prots$path[prots$pathway==pathway], what="", sep="\n")
  # Select the proteins from the DE subset. This way we get the intersection between DE and pathway prots.
  pathway.prots <- completeTableD21_WT[pathway.prot.names,]
  # Sort the pathway proteins by asbolute logFC:
  #   Get the order
  order <- order(abs(pathway.prots$logFC), decreasing = TRUE)
  #   Sort by that order and select first n
  #pathway.prots[order,][1:5,]
  # Add info about the pathway
  df <- data.frame(pathway.prots[order,][1:n,], pathway=pathway)
  return(df)
}

# Create a new data frame to store the info of all the pathways
custom.pathways <- data.frame()
#pathways[c(1:3,5:6)]
for (path in pathways[c(1:3,5:6)]){
  # Add to the merged data.frame
  custom.pathways <- rbind(custom.pathways,getTopProts(path, 5))
}
# Add the mitochondrial prots, we get the first 6 because we are going to eliminate 1
custom.pathways <- rbind(custom.pathways,getTopProts(pathways[4],6))
custom.pathways
# As we want to skip Q61941, we get the first 6, and then eliminate the one with that ID
which(rownames(custom.pathways)=="Q61941")
custom.pathways <- custom.pathways[-26,]

# Get Zq values (normalized) of those proteins
custom.subset <- merge(custom.pathways, exprs, by="row.names", all.x = TRUE)

# Discard unwanted columns. (We keep logFC to reorder)
custom.subset <- custom.subset[,c(-3:-7)]
# Reorder by logFC so that the ones increasing the most appear first
custom.subset <- custom.subset[order(custom.subset$pathway, custom.subset$logFC, decreasing = TRUE),]
# Save the pathway vector to use in the heatmap
pathways.ordered <- custom.subset$pathway
# Assign row names
rownames(custom.subset) <- custom.subset[,1]

# Delete unnecessary columns (names, logFC and pathway)
custom.subset <- custom.subset[,-1:-3]
custom.subset

# Finally, plot
# Transform to matrix so it can be plotted
custom.subset <- data.matrix(custom.subset, rownames.force = TRUE)
# Define which color is assigned to what pathway (3 proteins each pathway, so each color is assigned 3 times)
groups.color <- c(rep(palette[1],5), rep(palette[2],5), rep(palette[3],5), rep(palette[4],5), 
                  rep(palette[5], 5), rep(palette[6], 5))
groups.color

############################### ANNOTATION ################################ 
# Biomart
library('biomaRt')
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")
# listAttributes(ensembl)
# listFilters(ensembl)
biomart.prots <- getBM(attributes=c("uniprotswissprot", "external_gene_name"), # Info I want to get
                       filters="uniprotswissprot", # Type of info I'm providing
                       values=rownames(custom.subset), # Info I'm providing
                       mart=ensembl)
biomart.prots
# Eliminate duplicate entries
biomart.prots <- biomart.prots[!duplicated(biomart.prots$uniprotswissprot),]
# Change rownames to SwissProtID
rownames(biomart.prots) <- biomart.prots[,1]
# Some don't have gene name
setdiff(rownames(custom.subset), biomart.prots$uniprotswissprot)

# UniProt
# dplyr causes conflict with UniProts "select" method
detach("package:dplyr", unload=TRUE)
library("UniProt.ws")
# Check m musculus taxId
#availableUniprotSpecies(pattern="musculus")
# Create conection with that Id
up <- UniProt.ws(taxId=10090)
# Type of data (id) I will provide and vector with the ids
kt <- "UNIPROTKB"
keys <- rownames(custom.subset)
# Info I want to retrieve
# columns(up)
columns <- c("GENES")
# Query
uniprot.prots <- select(up, keys, columns, kt)
# Change rownames to SwissProtID
rownames(uniprot.prots) <- uniprot.prots[,1]
uniprot.prots
# Merge results
library("dplyr")

colnames(uniprot.prots) <- c("UniProtId", "GeneName")
colnames(biomart.prots) <- c("UniProtId", "GeneName")
uniprot.prots
biomart.prots
all <- full_join(uniprot.prots, biomart.prots, by = "UniProtId")
all
# Populate Second column with first column names where NA
all$GeneName.y[is.na(all$GeneName.y)] <- as.character(all$GeneName.x[is.na(all$GeneName.y)])
# Get the first name for each gene
all$GeneName.y <- sapply(strsplit(all$GeneName.y, " "), '[', 1)
# write.table(all, file="./ANALYSIS/03.CustomPathways/annotation.txt", sep = "\t", quote=FALSE, row.names = F)

pdf(file="./RESULTS/heatmap_custom_de_prots_5_ordered.pdf")
#png("Plot3.png", width = 1200, height = 1200, units = "px")
hm2 <- heatmap.2(custom.subset, trace = "none", scale = "row", col=hmcol, 
                 dendrogram = "column", # Only calculate dendrogram for rows (samples)
                 Rowv = NULL, # Don't reorder
                 #Colv = FALSE,
                 rowsep=c(5, 10, 15, 20, 25, 30), # Where the row separator should be
                 sepcolor = "white", # row separator color
                 colRow = "black", # Row groups colors
                 RowSideColors = groups.color # Sidebar indicating groups
                 , labRow = all$GeneName.y # Row label to use
) 

#hm2$colDendrogram
#c(custom.subset$pathways)
dev.off()

# Print legend
pdf(file="./RESULTS/heatmap_custom_de_prots_legend_5_ordered.pdf")
plot.new()
plot(legend("topleft",legend=unique(pathways.ordered), fill=palette[1:6], 
            border=FALSE, bty="n", y.intersp = 0.8, cex=0.8, title = "Canonical Pathway"))
dev.off()

# Export zq values table with gene names:
merged <- data.frame(custom.subset, GeneName = all$GeneName.y)
merged <- as.matrix(merged)
merged
rm(columns) 
write.table(merged, file="./ANALYSIS/03.CustomPathways/custom_pathways_top_de.txt", sep = "\t", quote=FALSE, row.names = F)
columnas <- c("SwissProtId", colnames(merged))



