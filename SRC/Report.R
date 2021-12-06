# Download libraries 
library(ASpli)
library(GenomicFeatures)

# GTF pre-processing 
genomeTxDb <- makeTxDbFromGFF(file = "/Users/anmol/Desktop/gencode.v34.annotation.gtf", format = "gtf")
features <- binGenome( genomeTxDb )
# BAM and target file 
BAMFiles <- "/Users/anmol/Desktop/file.bam"
targets <- data.frame( 
  row.names = paste0('Sample',c(1:12)),
  bam = BAMFiles,
  f1 = c( 'A','A','A','A','A','A',
          'B','B','B','B','B','B'),
  f2 = c( 'C','C','C','D','D','D',
          'C','C','C','D','D','D'),
  stringsAsFactors = FALSE)
getConditions(targets)
mBAMs <- data.frame(bam      = sub("[02]","",targets$bam[c(1,4,7,10)]),
                    condition= c("A_H","A_D","B_C","B_D"))

# Read counts against annotated features 
gbcounts  <- gbCounts( features = features, 
                       targets = targets, 
                       minReadLength = 100, maxISize = 50000,
                       libType="SE", 
                       strandMode=0)

# Junction based de-novo countring and splicing signal estimation
asd   <- jCounts(counts = gbcounts, 
                 features = features, 
                 minReadLength = 100,
                 libType="SE", 
                 strandMode=0)

# Differential gene expression and bin usage signal estimation 
gb      <- gbDUreport(counts=gbcounts,
                      contrast = c( 1, -1, -1, 1 ) )

# Differential junction usage analysis 
jdur    <- jDUreport(asd, 
                     contrast =  c( 1, -1, -1, 1 ) ,
                     mergedBams = mBAMs)

# Bin and junction signal integration 
sr      <- splicingReport(gb, jdur, counts =gbcounts    )

# Summary of integration of splicing signals along genomic-regions 
is      <- integrateSignals(sr,asd)

# Export results 
exportIntegratedSignals( is, output.dir="aspliExample",
                         sr, gbcounts, features, asd, mBAMs)

