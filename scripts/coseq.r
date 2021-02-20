library(coseq)
option_list <- list(
    make_option(c("-i", "--count_file"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."   ),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                 )				
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)


setwd(opt$outdir)
counts <- read.table(opt$count_file,header=T,row.names=1,sep='\t')
##########################
##########################
run_kmeans <- coseq(object=counts, K=2:20, transformation="logclr",model="kmeans")
pdf("kmeans.pdf")
plot(run_kmeans)
dev.off()
write.csv(tcounts(run_kmeans),"kmeans_modules.csv")
write.csv(coseqFullResults(run_kmeans),"fullresults-kmeans.csv")
##########################
##########################
##########################
run_arcsin_normal <- coseq(object=counts, K=2:4, iter=5, transformation="arcsin",
model="Normal",GaussianModel = "Gaussian_pk_Lk_Bk")
pdf("run_arcsin_normal.pdf")
plot(run_arcsin_normal)
dev.off()
write.csv(tcounts(run_arcsin_normal),"run_arcsin_normal.csv")
write.csv(coseqFullResults(run_arcsin_normal),"fullresults-run_arcsin_normal.csv")
##########################
##########################
##########################
run_logit_normal <- coseq(object=counts, K=2:4, iter=5, transformation="logit",
model="Normal",GaussianModel = "Gaussian_pk_Lk_Bk")
pdf("run_logit_normal.pdf")
plot(run_logit_normal)
dev.off()
write.csv(tcounts(run_logit_normal),"run_logit_normal.csv")
write.csv(coseqFullResults(run_logit_normal),"fullresults-run_logit_normal.csv")