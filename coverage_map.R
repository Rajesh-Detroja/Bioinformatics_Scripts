## Reference: https://blog.liang2.tw/posts/2016/01/plot-seq-depth-gviz/

# BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

library(data.table)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Gviz")

library(Gviz)

args <- commandArgs(trailingOnly = TRUE)
filename <- (args[1])

bedgraph_dt <- fread(filename,
  col.names = c('chromosome', 'start', 'end', 'value')
)

for (hg38_chr in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY")){
#for (hg38_chr in c("chr1")){
# Specifiy the range to plot
thechr <- hg38_chr
#st <- 41176e3
#en <- 41324e3

bedgraph_dt_one_chr <- bedgraph_dt[chromosome == thechr]
dtrack <- DataTrack(
  range = bedgraph_dt_one_chr,
  genome = 'hg38',
  type="histogram",
#  window="auto",
  col.histogram="darkblue",
  fill.histogram="darkblue",
  name="Coverage"
)

itrack <- IdeogramTrack(
  genome = "hg38", chromosome = thechr
)

gtrack <- GenomeAxisTrack()

png(filename=paste(hg38_chr,".png", sep = ""), 
    type="cairo",
    units="in", 
    width=8, 
    height=5, 
    pointsize=15, 
    res=300)

plotTracks(
  list(dtrack, gtrack, itrack),
  main=paste(gsub('.filtered.bedgraph','', filename)),
  cex.main=1,
  fontface.main=2,
  col.main="black",
  margin=6
  #from = st, to = en
)

dev.off()

}
