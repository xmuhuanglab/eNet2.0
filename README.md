# Introduction of eNet 2.0
eNet 2.0 is a comprehensive tool for enhancer network analysis using scATAC-seq data, making up of four main parts: 
1) enhancer network construction
2) network topology analysis
3) network comparison analysis
4) network dynamics analysis

![image](https://github.com/Abigail1l1l1l/eNet-2.0/assets/63504517/69587234-bf0e-44cb-8af7-a9f5666314bd)


# Tutorial
## Pull the docker image and run docker container
```
docker pull hlinhlin/enet2.0:latest
docker run -it hlinhlin/enet2.0 /bin/bash
```
## Network construction and topology analysis
```r
### Set working directory
cd /home/rstudio/kitematic/R/Topology/
### Prepare parameters
Rscript  args.R --atac.file='0.cre.csv' --atac.raw='YES' --rna.exists='YES' --rna.file='0.rna.csv' --rna.raw='YES' --genomeArg='hg19' --cutoff=90 --Hscore=75 --nCores=2 --task='test1' --inputdir='/home/rstudio/kitematic/R/Topology/' --dir='/home/rstudio/kitematic/R/Topology/'
### Step0: Preprare input files
Rscript Step0_Input.R --task='test1' 
### Step1: Calculate chromatin co-accessibility
Rscript Step1_Edge.R --task='test1' 
### Step2: Identify E-P links
Rscript Step2_Node.R --task='test1'
### Step3: Build enhancer networks
Rscript Step3_Network.R --task='test1' 
### Step4: Identify hub enhancers
Rscript Step4_Module.R --task='test1'
```
#### Check output
```shell
ls ./test1/out/
# 0.cre.mat.Normalized.RDS  0.rna.mat.Normalized.RDS  1.conns.Rdata         2.genePeakTabFilt.Rdata  3.NetworkInfo.Rdata  4.EnhPair.Rdata
# 0.cre.mat.RDS             0.rna.mat.RDS             2.GenePeakcorr.Rdata  2.GPPair.Rdata           3.NetworkList.Rdata  4.ModEnhInfo.Rdata
```

## Network comparison analysis
```r
### Set working directory
cd /home/rstudio/kitematic/R/Comparison/
### Prepare parameters
Rscript  args.R --stage='healthy disease' --metadata.file='0.metadata.csv' --atac.file='0.cre.csv' --atac.raw='YES' --rna.exists='YES' --rna.file='0.rna.csv' --rna.raw='YES' --genomeArg='hg19' --cutoff=90 --Hscore=75 --nCores=2 --task='test2' --inputdir='/home/rstudio/kitematic/R/Comparison/' --dir='/home/rstudio/kitematic/R/Comparison/'
### Step0: Preprare input files
Rscript  Step0_Input.R --task='test2' 
### Step1: Identify E-P links
Rscript  Step1_Edge.R --task='test2' 
### Step2: Identify E-P links
Rscript  Step2_Node.R --task='test2' 
### Step3: Build enhancer networks
Rscript  Step3_Network.R --task='test2' 
### Step4: Identify differential enhancer networks 
Rscript  Step4_DiffMod.R --task='test2'
```
#### Check output
```shell
ls ./test2/out/
# 0.cre.mat.Normalized.RDS  0.rna.mat.RDS          2.GenePeakcorr.Rdata     3.disease.NetworkInfo.Rdata  3.healthy.NetworkList.Rdata
# 0.cre.mat.RDS             1.disease.conns.Rdata  2.genePeakTabFilt.Rdata  3.disease.NetworkList.Rdata
# 0.rna.mat.Normalized.RDS  1.healthy.conns.Rdata  2.GPPair.Rdata           3.healthy.NetworkInfo.Rdata
```

## Network dynamics analysis
```r
### Set working directory
cd /home/rstudio/kitematic/R/Dynamics/
### Prepare parameters
Rscript  args.R --stage='stage1 stage2 stage3 stage4' --metadata.file='0.metadata.csv' --atac.file='0.cre.csv' --atac.raw='YES' --rna.exists='YES' --rna.file='0.rna.csv' --rna.raw='YES' --genomeArg='hg19' --cutoff=90 --Hscore=75 --nCores=2 --task='test3' --inputdir='/home/rstudio/kitematic/R/Dynamics/' --dir='/home/rstudio/kitematic/R/Dynamics/'
### Step0: Preprare input files
Rscript  Step0_Input.R --task='test3' 
### Step1: Identify E-P links
Rscript  Step1_Edge.R --task='test3' 
### Step2: Identify E-P links
Rscript  Step2_Node.R --task='test3' 
### Step3: Build enhancer networks
Rscript  Step3_Network.R --task='test3' 
### Step4: Identify differential enhancer networks 
Rscript  Step4_DyNet.R --task='test3' 
```
#### Check output
```shell
ls ./test3/out/
# 0.cre.mat.Normalized.RDS  0.rna.mat.RDS         1.stage3.conns.Rdata  2.genePeakTabFilt.Rdata     3.stage1.NetworkList.Rdata  3.stage3.NetworkInfo.Rdata  3.stage4.NetworkList.Rdata
# 0.cre.mat.RDS             1.stage1.conns.Rdata  1.stage4.conns.Rdata  2.GPPair.Rdata              3.stage2.NetworkInfo.Rdata  3.stage3.NetworkList.Rdata  4.DyNet.Rdata
# 0.rna.mat.Normalized.RDS  1.stage2.conns.Rdata  2.GenePeakcorr.Rdata  3.stage1.NetworkInfo.Rdata  3.stage2.NetworkList.Rdata  3.stage4.NetworkInfo.Rdata  4.NetworkDynamicsInfo.Rdata
```

