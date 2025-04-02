###This part is executed on the linux server.
###Execute the following 2 lines in the terminal, in the directory with the STAR count files.
# ls *.tab | cat > files.txt
# paste *.tab > merged.tsv

###Open R and execute the following block of code. I executed it on my personal computer.
dat <- read.table('rnaQuest/data/merged.tsv')
#Set ENSG column as rownames for dataset.
row.names(dat) <- dat$V1
head(dat)
#Each counts file from STAR has four columns. The first is the ENSG ID, the second is the unstranded counts that we are looking for, the
#third and fourth columns are the stranded counts (forward and reverse, respectively). We want to keep the second column of each sample.
dat <- dat[,c(FALSE,TRUE,FALSE,FALSE)]
head(dat)

#Add sample ID to column names of dat
files <- read.delim2('rnaQuest/data/files.txt', header = FALSE, sep = "")
files <- sub("\\..*", "", files$V1)
colnames(dat) <- files
head(dat)


###Replace ENSG ID with gene symbol with annotation library.
# install.packages("BiocManager")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")  # For human; use org.Mm.eg.db for mouse

library(AnnotationDbi)
library(org.Hs.eg.db)

##Annotation
# Example list of ENSG IDs
ensg_ids <- rownames(dat)

# Convert ENSG IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = ensg_ids, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

print(gene_symbols)
head(gene_symbols)
dat <- cbind(gene_symbols,dat)
dat$gene_symbols[is.na(dat$gene_symbols)] <- "NA"
head(dat)

#Write out merged.tsv
write.table(dat, file = "rnaQuest/data/merged.tsv", row.names=TRUE, sep="\t")
