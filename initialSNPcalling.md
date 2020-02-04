# SNPcalling.md

This script process the data from raw to an initial set of SNPs using Stacks v2.41. All the related metadata can be requested from dutoit.ludovic@gmail.com.



### Demultiplexing (process_radtags)

```
#bash
cd raw 
ln -s ../source_files/Worms/*gz .
cd ..

#!/bin/sh
process_radtags -P -p raw/ -o ./samples/ -b ~/repos/scripts/polychaetes_GBS/metadata/barcodes.txt -e pstI -r -c -q --inline_inline -T 8
```

```
370408774 total sequences
 39463860 barcode not found drops (10.7%)
   269317 low quality read drops (0.1%)
  1229681 RAD cutsite not found drops (0.3%)
329445916 retained reads (88.9%)
```


### Concatenating the output sample by sample

```
#bash
cd samples
rm *rem*
rm Empty*
```

then I move them all into the concat folder concatenating forward and reverse together.


```
#bash
mkdir concat
cd samples
rm *rem*
rm Empty*
```

Then I concat samples by samples all passing reads wether reverse or forward (code not shown here but available on request)


### Run denovo_map.pl 

Preliminary analyses suggest parameters don't lead to massive changes in the number of SNPs. We therefore run denovo_map.pl with default settings.
```
#bash
cut -f 3 barcodes.txt | grep -v Empty  | awk '{print $0,"\tsinglepop"}' - > popmap_all.txt
mkdir output_M2
echo '#!/bin/sh' > cleanM2.sh
echo "denovo_map.pl --samples concat/ --popmap popmap_all.txt  -o output_M2  -M 2 -n 2 -m 3 -T 16"  >> cleanM2.sh
sbatch -A uoo00116 -t 10-00:00:00 --partition=long -J cleanM2 -c 16 --mem=64G cleanM2.sh # specific to mahuika and ludovic.dutoit

```

Two samples ,  1.1.2.16 and 3.2.2.27  have absolutely no data after processing, causing the end of the pipeline to crash. We therefore run the end of the pipeline manually after removing those 2 samples.

First, we create a new population map without those two samples:

```
grep 1.1.2.16 popmap_all.txt -v  >   popmap_all_nonull.txt
grep 3.2.2.27  popmap_all_nonull.txt -v  >   popmap_all_nonull2.txt
```
We then run the last steps of the pipeline.


```
tsv2bam -P output_M2/ -M popmap_all_nonull2.txt
gstacks -P output_M2 -M popmap_all_nonull2.txt

```

We then run populations, excluding no samples and only the SNPs with more than 65% heterozygosity as likely collapsed paralogs.

```
populations -P output_M2/ -M popmap_all_nonull2.txt  --vcf  --max-het 0.65 # then filter it with 
```
This creates a *source* unfiltered dataset for which we can try different type of filtering.  We are interested in removing low quality individuals ( very few SNPs sequenced) and very low quality SNPs ( very few individuals covered).

The unfiltered dataset is at 
The SNPs are in [output/unfiltered/](output/unfiltered/). 


### Filtering

**Exploration**

I then use a little bit of R code to investigate how much SNps we have across how many individuals?

```
library("VariantAnnotation")
data<-readVcf("populations.snps.vcf")


numbermissing<-function(x){
  return(length(grep("\\./\\.",x)))
}

countsofmissing<-apply(geno(data)$GT,2,numbermissing)
nonmissing<-dim( geno(data)$GT)[1]-countsofmissing
### get list of removals ad then re run as white list with  proportion of missing look
summary(nonmissing)
quantile(nonmissing,seq(0,1,0.05))
# 0%  5%  10% 15% 20% 25% 30% 35%
# 12.0  235.2 478.6 625.6 944.0 1179.0  1571.4  2132.0
# 40% 45% 50% 55% 60% 65% 70% 75%
# 3586.2  5286.2  7602.0  10359.8 13533.0 16869.2 19376.6 23608.0
# 80% 85% 90% 95% 100%
# 28335.2 33788.0 50114.0 71672.4 150823.0

```
Without strong knowledge of the systems and the study goal, it is hard to make a fair call on filtering criteria (i.e. which individuals we will exclude). We do know that we might not need very many SNPs but that large amount of missing data are problematic for many analyses that imputes genotypes in the presence of missing data. (i.e. structure).


I choose to make a dataset where I removed the worse 1/3 of samples, therefore excluding any sample with less than 1919 SNPs:

```
quantile(nonmissing,0.33)
lowest_samples <-names(which(nonmissing<1918.48))# remove the lowest 33%
```

I explore what is left and what do if we remove those in R still.
```
prefilter_gt<-geno(data)$GT # create a genotype matrix
sample_filtered<- prefilter_gt[,-which(colnames(prefilter_gt)%in%lowest_samples)] # remove those samples
propofmissingpersnp<-apply(sample_filtered,1,numbermissing)/dim(sample_filtered)[2] # what is the proportion of missing data snp by snp
length(which(propofmissingpersnp<0.4)) # 60% or more of individuals have data, how many SNP is that?
```

***Removing the 33% individuals with the lowest number of SNPs and then excluding SNP for which less than 60% of individuals have data seem to be a fair call for a filtering step ***

To do so, I create a final popmap without those individuals

```
popmap<-read.table("../popmap_all_nonull2.txt")
new_popmap<-popmap[-which(as.character(popmap[,1])%in%lowest_samples),]
write.table(new_popmap,"../popmap_filtered33percent.txt",sep="\t",row.names=F,col.names=F,quote=F)
```
And finally back in bash I create those SNP files ion the vcf, structure, plink and treemix formats

```
mkdir filtered33percent
populations -P output_M2/ -M popmap_filtered33percent.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -r 0.6 -O filtered33percent # then filter it with 
```

  
We obtained a final dataset of 319 SNPs for 1857 individuals identified in [popmap_filtered33percent.txt](output/filtered33percent/popmap_filtered33percent.txt).

The SNPs are in [output/filtered33percent/](output/filtered33percent/)

 Some loci have very many SNPs across their 76bp. I did not exclude them as I know nothing about the genetic diversity of polychaetes:

```
#number of SNPs per loci    1 2 3 4 5 6 7 8 9 10  11  12  13  14  15  16  17  18
#Numberof loci with n SNPs    284 165 78  53  28  18  9 10  7 8 4 2 3 3 1 4 1 1
```

**Removing 40% of individuals but getting more SNPs**

To show a contrasting pictures, I can retrieve more SNPs by removing 

back into R :
```
library("VariantAnnotation")
data<-readVcf("populations.snps.vcf")


numbermissing<-function(x){
  return(length(grep("\\./\\.",x)))
}

countsofmissing<-apply(geno(data)$GT,2,numbermissing)
nonmissing<-dim( geno(data)$GT)[1]-countsofmissing

hist
quantile(nonmissing,0.4)
    40%
3586.2
lowest_samples <-names(which(nonmissing<3586.2))# remove the lowest 40%

prefilter_gt<-geno(data)$GT # create a genotype matrix
sample_filtered<- prefilter_gt[,-which(colnames(prefilter_gt)%in%lowest_samples)] # remove those samples
propofmissingpersnp<-apply(sample_filtered,1,numbermissing)/dim(sample_filtered)[2] # what is the proportion of missing data snp by snp
length(which(propofmissingpersnp<0.4)) # 60% or more of individuals have data, how many SNP is that?
#[1] 2900
popmap<-read.table("../popmap_all_nonull2.txt")
new_popmap<-popmap[-which(as.character(popmap[,1])%in%lowest_samples),]
write.table(new_popmap,"../popmap_filtered40percent.txt",sep="\t",row.names=F,col.names=F,quote=F)  

### Add a plot for visualisation
png("numberofSNPsperindsbeforefiltering.png")
par(mfrow=c(1,2))
 hist(nonmissing, main =" Overall " ,xlab="Number of SNPs", ylab= "Number of individuals",breaks=1000,col="black")
 hist(nonmissing,  main ="Zoomed in" ,xlab="Number of SNPs", ylab= "Number of individuals",breaks=1000,xlim=c(0,10000),col="grey")
 dev.off()
```

Finally we create the dataset using:
```
mkdir filtered40percent
populations -P output_M2/ -M popmap_filtered40percent.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -r 0.6 -O filtered40percent # then filter it with 
```


