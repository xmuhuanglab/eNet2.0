# Rscript  args.R --stage='healthy disease' --metadata.file='0.metadata.csv' --atac.file='0.cre.csv' --atac.raw='YES' --rna.exists='YES' --rna.file='0.rna.csv' --rna.raw='YES' --genomeArg='hg19' --cutoff=90 --Hscore=75 --nCores=2 --task='111' --inputdir='/home/rstudio/kitematic/R/Comparison/' --dir='/home/rstudio/kitematic/R/Comparison/'

library("argparse")

# 创建参数解析对象
parser <- ArgumentParser()

# 设置参数
parser$add_argument("--stage", type="character", default='healthy disease',
                    help="Two names of cell stages for comparison, separated by ' ' ")

parser$add_argument("--metadata.file",
                    help="The metadata of cell names and stages")

parser$add_argument("--atac.file",
                    help="The file name of scATAC-seq matrix")

parser$add_argument("--atac.raw", default='YES',
                    help="Is the scATAC-seq matrix raw count?")

parser$add_argument("--rna.exists", type="character", default='YES', 
                    help="Does the scRNA-seq matrix exist?")

parser$add_argument("--rna.file", 
                    help="The file name of scRNA-seq matrix")

parser$add_argument("--rna.raw", default='YES',
                    help="Is the scRNA-seq matrix raw count?")

parser$add_argument("--genomeArg", 
                    help="The reference genome used. It can only be hg19/hg38 for human, and mm9/mm10 for mouse.")

parser$add_argument("--cutoff", type='integer', default=90, 
                    help="The co-accessibilty threshould value in the format of quantile, which will be used to determine significantly co-accessibile enhancer pairs.")

parser$add_argument("--Hscore", type='integer', default=75, 
                    help="The Hub score threshould value in the format of quantile, which will be used to define hub enhancers.")

parser$add_argument("--nCores", type='integer', default=2, 
                    help="The number of cores.")

parser$add_argument("--task", 
                    help="The task number.")

parser$add_argument("--dir", 
                    help="Initial working path.")

parser$add_argument("--gdir", default='/home/rstudio/kitematic/genome/',
                    help="The path to store the genome file.")

parser$add_argument("--inputdir", default='/home/rstudio/kitematic/',
                    help="The path to store the genome file.")

# 调用解析器，此时args就被赋值为命令行参数输入的相应值
argsAll <- parser$parse_args()

stage <- argsAll$stage
metadata.file <- argsAll$metadata.file
atac.file <- argsAll$atac.file
atac.raw <- argsAll$atac.raw
rna.exists <- argsAll$rna.exists
rna.file <- argsAll$rna.file
rna.raw <- argsAll$rna.raw
genomeArg <- argsAll$genomeArg
cutoff <- argsAll$cutoff
Hscore <- argsAll$Hscore
nCores <- argsAll$nCores
task <- argsAll$task
dir <- argsAll$dir
gdir <- argsAll$gdir
inputdir <- argsAll$inputdir
rm(parser)
rm(argsAll)


dir.create(paste0(dir,task,'/out/'), recursive = T)
save.image(file = paste0('args_',task,'.RData'))

