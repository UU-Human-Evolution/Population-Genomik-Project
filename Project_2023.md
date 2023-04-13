# Population Genomik Project [March 2023]

## Instructions for project work

The goal of this population genomics project is to familiarize you with the workflows, different types of commonly used analyses as well as a way of thinking about how one approaches questions such as trying to explain the underlying forces that affect populations in general. We know that people tend to migrate and throughout history they have done that a lot. Similar populations get isolated and with time may become more distant, distant peoples admix with time may become more closely related - and trying to decipher human population history is usually not an easy task.

The structure of the project work will be organized as follows: 

- you need an username & password to log into Uppmax;
- up to you to decide if you want to work on your own laptops or in Hubben;
- we have an introductory meeting where we discuss the framework of the exercize in depth;
- in the few weeks you have until the course ends feel free to send in emails/slack with specific questions if something is not working (or you're just curious);
- we can also set up meeting times that work for all of us;
- the end product should be a Scientific report of ~10 pages, including Abstract, Introduction, Materials & Methods, Results, Discussion, References, Appendices (length is not crucial here, substance is);
- A short presentation on the 17th of March;
- the deadline for the written report is the day of the presentation.

## Unmasking the mystery

The specific point of the exercise is for you to identify the populations you have been given in the **'Unknown.fam'**. In order to achieve this, you have been provided a bunch of references (three different datasets) that should come in handy as a way to compare the unknown to the known. Both reference and target populations can be found in the DATA directory:
```
/proj/uppmax2023-2-18/MYSTERY_QUEST/DATA/Reference_datasets
/proj/uppmax2023-2-18/MYSTERY_QUEST/DATA/Mystery_populations
```
Helpful scripts are there to aid you on your quest can be found:
```
/proj/uppmax2023-2-18/MYSTERY_QUEST/SCRIPTS
```
The unknown samples are a mysterious bunch. They look different from the rest don't they? (*hint - check the genotyping rate - i.e. missingness*).
Your first task is to figure out what types of samples you have stumbled upon and why they are (potentially) very valuable, yet they look the way they do (*the real Sherlocks among you would have noticed by now something odd by just glancing at their names!*) 

In the steps that will lead you towards unmasking the unknown populations you will need to: 

1. Check for related individuals in your **reference dataset**;
2. Filter the related individuals from the **reference dataset** - *we have nothing against family, but individuals that are in close kinship will affect our downstream analysis and we don't want that*;
3. Filter the **reference** individuals & SNPs (missingness, HWE, MAF) - *word of caution: do not apply the 'maf' filtering until you have your final references dataset!*;
   *Remember: The **unknown** samples **do not** go through any filtering!*
4. Merge the reference datasets to each other;
6. Unify SNP names across your dataset - change the SNP name into the names based on position. This is needed when merging datasets from different SNP chips (the same position can have different names on different chips);
7. After merging together your reference dataset and you've finished filtering it - merge it together with your **unknown** data;
8. Before running PCA/Admixture you need to prune your data for SNPs in LD (pruning is explained in the admixture manual https://dalexander.github.io/admixture/admixture-manual.pdf);
9. Run a PCA on just on your "reference" dataset as well as on the combined "reference + Unknown dataset";
10. Run a projected PCA;
11. Run Admixture (remember - Admixture takes a lot of time! Organize your time wisely);
12. Read up on some literature - it will help for writing the report! (Pay special attention to the key figures in them):
[Schlebusch 2012](https://pubmed.ncbi.nlm.nih.gov/22997136/)
[Gurdasani_2015](https://www.nature.com/articles/nature13997) 
[Patin 2017](https://science.sciencemag.org/content/356/6337/543)
[Vicente_2021](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-01193-z)
13. Don't panic - there's plenty of time. Or is there?

## Some technicalities

* You should have R studio installed on your computer

Before we start here are some basic Unix/Linux commands if you are not used to working in a Unix-style terminal:

### Logging into SNOWY & working on the server:

```
ssh -AX user@rackham.uppmax.uu.se
```
After which you type in your password. *And yes, it'll not show up when typing it, but if you type it in correctly & press enter, you're in.*

Note: When acessing SNOWY from Rackhams login nodes you must always use the flag -M for all SLURM commands.
Examples:

```
- squeue -M snowy
- jobinfo -M snowy
- sbatch -M snowy slurm_script_file
- scancel -u username -M snowy
- interactive -A projectname -M snowy -p node -n 32 -t 01:00:00
```

Note: It is recommended to load all your modules in your job script file. This is even more important when running on Snowy since the module environment is not the same on the Rackham login nodes as on Snowy compute nodes.
You can read up more on how to use SNOWY : [Snowy_User_guide](https://www.uppmax.uu.se/support/user-guides/snowy-user-guide/)

Generally, for any bigger job (bigger Plink jobs, PCA, Admixture) **always work in interactive mode**.

### How to run interactively on a compute node:

Ex.:
```
   interactive -A uppmax2023-2-18 -M snowy -p node -n 32 -t 02:00:00
```
### Moving about:
```
    cd – change directory
    pwd – display the name of your current directory
    ls – list names of files in a directory
```
### File/Directory manipulations:
```
    cp – copy a file
    mv – move or rename files or directories
    mkdir – make a directory
    rm – remove files or directories
    rmdir – remove a directory
```
### Display file content:
```
    cat - concatenate file contents to the screen or to a file
    less - open file for viewing
    more 
    tail
    head
    realpath - checking the path to a file 
```
### You can check your jobs with:

```
jobinfo -u YOUR_USERNAME -M snowy
```
When creating and editing a file in text editor you can use ```nano``` or ```vim``` or something else.

# PART 0 Knowing your way around & using Plink   	    

**Beware** - it is vital to not apply any filters to your unknown samples - we don't know anything about them at this point! 

FYI: Link to PLINK site:[https://www.cog-genomics.org/plink2](https://www.cog-genomics.org/plink2)

PLINK is a software for fast and efficient filtering, merging, editing of large SNP datasets, and exports to different outputs directly usable in other programs. Stores huge datasets in compact binary format from which it reads-in data very quickly.

### Running the program:

If working on Rackham type in the following:
```
module load bioinfo-tools
module load plink/1.90b4.9 

```
The software is already pre-installed on Rackham, you just have to load it

Try it:
```
plink
```

### This is what a basic command structure looks like in Plink: 
``` 
plink --filetypeflag filename --commandflag commandspecification --outputfilecommand --out outputfilename
```
For example to do filtering of missing markers at 10% frequency cutoff (reading-in a bed format file called file1, doing the filtering, writing to a ped format file called file2):
```
plink --bfile file1 --geno 0.1 --recode --out file2
```

### Input formats:
File format summary:

ped format: usual format (.ped and .fam)
.ped contains marker and genotype info and .fam files contain sample info

bed format: binary/compact ped format (.fam .bim and .bed)
(.fam - sample info  .bim - marker info  and .bed - genotype info in binary format

tped format: transposed ped format (.tfam and .tped files)

tfam sample info, .tped marker and genotype info in transposed format

lgen format: long format (see manual, not used that often)

### Setup
Navigate to the directory you want to work in.

```
cd /proj/uppmax2023-2-18 #Uppmax project for this course
```
If you don't already have a working directory for this course then create one now:
```
mkdir your_unique_team_name
```

### The course material can be found at the following path:
```
/proj/uppmax2023-2-18/MYSTERY_QUEST
```

Copy the datasets from the directory “DATA” to your working folder by **changing** the command below while you are in your working folder. The **dot** means copy the contents "here". Do this for all your datasets. 

*Advice - it might be easier for you if all your datasets (.bed +.bim +.fam files for all four, instead of having to go in and out of the different folders) are in the same directory.

```
cp /full_path_to_course_material/unk1.bed .
cp /full_path_to_course_material/unk1.bim .
cp /full_path_to_course_material/unk1.fam .
```

**The files above contain SNPs from a number of unknown individuals in bed file format. You are going to figure out the ancestry of these population groups during this practical.**
 
Look at the `.bim` (markers) and `.fam` (sample info) files by typing:

```
less Unknown.bim
``` 
do the same for the `.fam` file

```
less Unknown.fam 
```
As mentioned before the `.bim` file store the variant markers present in the data and `.fam` lists the samples. (you can try to look at the .bed as well, but this file is in binary format and cannot be visualized as text if you want to see the specific genotype info you must export the bed format to a ped format)

Read in a  bed file dataset and convert to ped format by typing/pasting in:

```
plink --bfile Unknown --recode --out Unknown_ped 
```

### Look at the info that Plink prints to the screen. 
How many SNPs are there in the data? How many individuals? 
How much missing data? What do their names suggest?

Look at the missingness information of each individual and SNP by typing:

```
plink  --bfile Unknown --missing --out test1miss

```
Look at the two generated files by using the `less` command (q to quit)

```
less test1miss.imiss
less test1miss.lmiss
```

The `.imiss` contains the individual missingness and the `.lmiss` the marker missingness
Do you understand the columns of the files? The last three columns are the number of missing, the total, and the fraction of missing markers and individuals for the two files respectively

Look at the first few lines of the newly generated .map (fileset variant information file) and .ped (Paternal ID Maternal ID Sex Phenotype) files using the more command as demonstrated above.

Read in bed/ped file and convert to tped:

```
plink --bfile Unknown --recode transpose --out Unknown_tped 
plink --file Unknown_ped --recode transpose --out Unknown_tped 
```
**Do you see the difference in the two commands above for reading from a bed (--bfile) and reading from a ped (--file) file. 
Which one takes longer to read-in?**

Look at the first few lines of the  `.tfam` and `.tped` files by using the `less` command

**Can you see what symbol is used to encode missing data?**

Note- try to always work with bed files, they are much smaller and takes less time to read in. 
See this for yourself by inspecting the difference in file sizes:

```
ls -lh * 
```

Plink can convert to other file formats as well, you can have a look in the manual for the different types of conversions

# PART 1: Merging, Filtering, QC steps

## Step 1 Filtering for missing data
You need to do all of the following steps for each of the three reference datasets separately.
*Remember, you are only filtering the reference datasets. Not the Unknown.*

First, we filter for marker missingness:

Paste in the command below to filter out markers with more than 10% missing data

```
plink --bfile YOUR_DATASET --geno 0.1 --make-bed --out YOUR_DATASET2 
```

Look at the screen output, how many SNPs were excluded?

## Step 2 Filter for individual missingness

Paste in the command below to filter out individual missingness

```
plink --bfile YOUR_DATASET2 --mind 0.15 --make-bed --out YOUR_DATASET3 
```

Look at the screen output, how many individuals were excluded?

## Step 3 MAF filtering

To filter for minimum allele frequency is not always optimal, especially if you are going to merge your data with other datasets in which the alleles might be present. But we will apply a MAF filter in this case
(remember not to filter the 'Unknown' samples)

Filter data for a minimum allele frequency of 1% by pasting in:

```
plink --bfile YOUR_DATASET3 --maf 0.01 --make-bed --out YOUR_DATASET4 
```

How many SNPs are left at this point?

## Step 4 Filtering for SNPs out of Hardy-Weinberg equilibrium

Most likely, SNPs out of HWE usually indicates problems with the genotyping. However, to avoid filtering out SNPs that are selected for/against in certain groups (especially when working with case/control data) filtering HWE per group is recommended. After, only exclude the common SNPs that fall out of the HWE in the different groups - (OPTIONAL). But for reasons of time, we will now just filter the entire dataset for SNPs that aren’t in HWE with a significance value of 0.001

```
plink --bfile YOUR_DATASET4 --hwe 0.001 --make-bed --out YOUR_DATASET5
```

Look at the screen. How many SNPs were excluded?

If you only what to look at the HWE stats you can do as follows. By doing this command you can also obtain the observed and expected heterozygosities. 

```
plink --bfile YOUR_DATASET5 --hardy --out hardy_YOUR_DATASET5
```
Look at file hardy_YOUR_DATASET5.hwe, see if you understand the output?

There are additional filtering steps that you can go further. PLINK site on the side lists all the cool commands that you can use to treat your data. Usually, we also filter for related individuals and do a sex-check on the X-chromosome to check for sample mix-ups. 

## Step 5 Filtering out related individuals

In the `SCRIPTS` folder, there is a script called `sbatch_KING.sh` that can be used to run [KING](http://people.virginia.edu/~wc9c/KING/manual.html) You can have a look inside for instructions on how to run the script. Check out the manual and try and figure out how the software works.
After you have run the script have a look at the produced output files and figure out how to remove the related individuals. (*Hint - keep via Plink might be a good option) 
*P.S. You don't need to do run KING for the Unknown dataset.*

First load ```bioinfo-tools``` and then ```module load KING```.
The script for KING is as follows:

```
#!/bin/bash -l
#
#
#SBATCH -J king
#SBATCH -t 12:00:00
#SBATCH -A uppmax2023-2-18
#SBATCH -n 8
king -b $1  --unrelated
```
You can submit it as a job by doing:
```
sbatch -M snowy THIS_SCRIPT.sh YOUR_DATASET5.bed
```
Look at the output from KING & keep the unrelated individuals.
Checkpoint *At this point, when you exclude the related individuals you should be at YOUR_DATASET_6*


## Step 6 Data Merging & strand flipping

The next step would be to start merging your data with comparative datasets. But first you need to merge the three reference datasets to eachother, filter and clean-up, before you finally merge the **combined** to the **Unknown** dataset. 

Usually, when you merge your data with another dataset there are strand issues. The SNPs in the other dataset might be typed on the reverse DNA strand and yours on the forward, or vice versa. Therefore you need to flip the strand to the other orientation for all the SNPs where there is a strand mismatch. One should not flip C/G and A/T SNPs because one cannot distinguish reverse and forward orientation (i.e. C/G becomes G/C unlike other SNPs i.e. G/T which become C/A). Therefore before merging and flipping all A/T and C/G SNPs must be excluded. However, this can be a problem since some of your SNPs in your dataset may be monomorphic when you don't apply the MAF filter. I.E in the bim file they will appear as C 0 (with 0 meaning missing). So you don't know what kind of SNP it is, it can be C G or C T for instance if it is C G it needs to be filtered out but not if it is C T.
Therefore, before merging our data to other datasets it is important to first merge your data with a fake / reference_individual, that you prepare, which is heterozygous at every SNP position. This “fake” reference individual you can easily prepare from the SNP info file you get from the genotyping company or your own genomics processing software (such as Genome Studio from Illumina). You can also prepare it from data downloaded for each SNP from a web-database such as dbSNP. 

### A Proper merging path
First we prepare for merging the datasets.

### A.1 Rename SNPs 
When merging datasets from different chips, the same position can have different names on different chips. The rename_SNP.py script can be used on the .bim files to change the SNP name into the names based on position. Luckily for you, this has already been done for the three reference datasets. Otherwise there would have been just the additional step of renaming the SNPs for each dataset prior to merging. But since we don't know where the Unknown dataset came from, we want to conduct this step for it. Check the bim file for the Unknown dataset and see the SNP names before you rename them.

In order to run rename_SNP.py first load ```module load python/2.7.11``` (Issues might occur if not using the right version of Python).

```python rename_SNP.py Unknown.bim```
After which you should get that the output with replaced SNP names has been written to Unknown.bim_pos_names. You can delete the original bim (or rename it if you're feeling unsure) & edit the ```Unknown.bim_pos_names``` to be ```Unknown.bim```.

Check to see what happened after the script did its job.

### A.2 Merging the reference datasets to the RefInd

We want to merge the our data with a set of reference populations that we get from an already published study such as HapMap data or HGDP population data. Many of the sites archiving the data provided them in PLINK format as well. For this practical, we selected a few populations from HapMap and HGDP to compare your Unknown populations to.

*Have you noticed that PLINK sometimes generates a .nof file and in the log file output the following is mentioned 902 SNPs with no founder genotypes observed Warning, MAF set to 0 for these SNPs (see --nonfounders) This is the monomorphic SNPs in your data.*

So our first step will be merging with a fake reference (RefInd) that we prepared from SNP info data beforehand. 
* Copy the RefInd1 files (bim, bed, and fam) from the DATA folder to your working folder. 
* Copy the three reference datasets (bim, bed and fam) to the working folder. Look at the .fam files, do you recognize some of these populations? There are two HapMap and three HGDP populations (You can read up on what they are).
* **Just as you did for the Unknowns, rename the SNPs for the RefInd as well, using the same script as above**

The next step is to extract your SNPs of interest from the RefInd. Use the Schlebusch_2012 dataset as a starting point.
```
plink --bfile RefInd --extract Schlebusch_2012_6.bim --make-bed --out RefInd1_ext 
```
Make a list of CG and AT SNPs in your data:
```
sed 's/\t/ /g' RefInd1_ext.bim | grep " C G" >ATCGlist
sed 's/\t/ /g' RefInd1_ext.bim | grep " G C" >>ATCGlist
sed 's/\t/ /g' RefInd1_ext.bim | grep " A T" >>ATCGlist
sed 's/\t/ /g' RefInd1_ext.bim | grep " T A" >>ATCGlist
```
Exclude the CG and AT SNPs from both your RefInd and data:
```
plink  --bfile RefInd1_ext --exclude ATCGlist --make-bed --out RefInd1_ext2 
plink  --bfile Schlebusch_2012_6 --exclude ATCGlist --make-bed --out Schlebusch_2012_7
```
Merge with RefInd
```
plink --bfile RefInd1_ext2 --bmerge Schlebusch_2012_7.bed Schlebusch_2012_7.bim Schlebusch_2012_7.fam --make-bed --out MergeRef1  
```
*At this point (or some of the following ones, you might get an error message* - **error is generated because of the strand mismatches**. 

This will generate a file (MergeRef1.missnp) which contains the info on the SNPs where there are mismatches - you need to flip the strand of these SNPs in your data.

```
Schlebusch_2012_7 --flip MergeRef1-merge.missnp --make-bed --out  Schlebusch_2012_8  
```
Try merging again:
```
plink --bfile RefInd1_ext2 --bmerge Schlebusch_2012_8.bed Schlebusch_2012_8.bim Schlebusch_2012_8.fam --make-bed --out MergeRef2  
```
Now it works. No .nof file is generated which means none of your SNPs are monomorphic anymore

*Possible issue: You've done all the flipping and merging and you are still getting an error?*

```
Error: ~50 variants with 3+ alleles present. 
* If you believe this is due to strand inconsistency, try --flip with 
MergeRef2-merge.missnp. 
```
It could be that a few positions are causing this issue. In this case (assuming its just a few - check what the output says) you can exclude those problematic variants and remerge:

```
plink --bfile source1 --exclude merged.missnp --make-bed --out source1_tmp
plink --bfile source2 --exclude merged.missnp --make-bed --out source2_tmp
plink --bfile source1_tmp --bmerge source2_tmp --recode transpose --allow-no-sex --out Excluded_Remerged_Dataset
rm source1_tmp.*
rm source2_tmp.*
```
If unsure, at any time you can check Plink's documentation (https://www.cog-genomics.org/plink/1.9/data) on how to do this.

After that you can continue with merging the next two reference datasets (Patin & Gurdasani). Just as prevously, remove the A/T & C/G SNPs in each dataset & keep only the overlap between the three datasets (i.e. extract one from the other and then do the merging and flipping).

*Checkpoint* - You can name this dataset ```ALL_REFERENCE_DATASETS``` so that it's easier to follow.

### A.3 Merge the Unknown dataset to the all-three-reference-combined-dataset one
Finally, after you've merged all the reference datasets, we need to merge in the Unkowns as well. To do this, as previously, filter out the A/T & C/G SNPs and merge, flip and exclude (if you have to). The only difference here is that when adding the Unknowns, we just merge it to the rest (instead of keeping only the overlap, we keep everything at this point)  

*Checkpoint* - You can name this dataset ```ALL_COMBINED_DATASETS``` so that it's easier to follow.

### A.4 Remember to extract your fake (RefInd) from your data

```
plink --bfile ALL_COMBINED_DATASETS --remove RefInd1.fam --make-bed --out FINAL_DATASET  
```
This is the final files for the next exercise.

Now you have generated your input files for the next exercise which will deal with population structure analysis. You will look at the population structure of your unknown samples in comparison to the known reference populations from HapMap and HGDP.

*Again the official Plink documentation (https://www.cog-genomics.org/plink/1.9/data) has a good explanation on how to do all of the above.*

Make sure to remember what your two final datasets are called. The one where you merged all the reference datasets and the one where you did that but also added the Unknown. 

## Step 7 Remove all unnecessary iterations of the same data except .log files from the merging/filtering steps and of course the final dataset
As our storage on uppmax is quite limited, this is a good point to free up your space from all the data that was produced from each of the previous steps and just keep the final dataset.  

## Step 8 Reflect and review

How many SNPs are left for your analyses? 
What was the genotyping rate for the combined reference dataset and what was it for the unknowns? 
Finally what is the genotyping rate for the final merged dataset? 

This are the final files for the next exercise. Now you have generated your input files for the next part which will deal with population structure analysis. You will look at the population structure of your unknown samples in comparison to the known reference populations that came from HapMap and HGDP.

=============================================================================

# PART 2: Population structure inference 

**Word of advice - the following analyses might take quite some time to run. Therefore, plan accordingly and leave yourself enough room for error.**

Furthermore, make sure to doublecheck your script before running it. 90% of issues you will encounter when running scripts are going to be typos (*Sanchez,R., Smith,M., 2022*). Check if you have written the correct file names, the correct extensions, whether you have the right path to your file or whether you have typed in your correct project name. Since we have a limited CPU hours allocated for this course it is really important to be efficient!  

## Step 1 LD Pruning 

Before running PCA & ADMIXTURE it is advisable to prune the data to thin the marker set for linkage disequilibrium.
You can read up on how to prune for LD(https://dalexander.github.io/admixture/admixture-manual.pdf).

```plink --bfile FINAL_DATASET --indep-pairwise 10 10 0.1
   plink --bfile FINAL_DATASET --extract plink.prune.in --make-bed --out FINAL_DATASET_PRUNED

```
Check how much you're pruning out. Perhaps you can tweak the parameters.

## Step 2 Principal component Analysis with Eigensoft

The first population structure method we will look at is Principal Components Analysis with Eigensoft. 
There is a really good explanation on how to run smartPCA(https://gaworkshop.readthedocs.io/en/latest/contents/05_pca/pca.html). If you are feeling able, and still on time, it would be cool if you ran a smartPCA with only the *Reference* populations and one from the *Combined* dataset and see if you can spot a difference.
Look at the output PDF. How many of the PCs do you think contain useful information. What part of the variation is represented by each of the PCs. Can you see the percentage variation that each PC explains?

## Step 3 Projected PCA

There is a reason why we are doing a projected PCA instead of a "conventional" PCA. It would be cool if you thought about why we have decided to do this step. Maybe read up a bit on what it is and what it does.

For running a projected PCA, you first need to produce some files:

  - a tped file (your final merged dataset); 
  - a list of the modern populations in the dataset;
  - a list of the unknown populations in the dataset;

To produce the tped, you can use plink as you've done so far. 

And to produce the two population lists:

```
awk {'print $1'} FINAL_DATASET.tfam |sort|uniq > Sorted_list_of_the_populations.txt
```
Then, this sorted list of populations, you can subset to an ```unknown.txt``` and a ```modern.txt``` using any approach you like.
*Note Populations IDs in "modern reference pops" must exactly match the IDs in the tfam, while population IDs in "unkown pops to include" can be more general (e.g. HG_SE instead of HG_SE_SF12) - No error messages when this is violated but weird results may occur!

Hint: Check the ```pca_lsqpeoj.sh``` file if it is calling the following modules:

```
module load bioinfo-tools
module load plink/1.90b4.9
module load python/2.7.11
module load eigensoft/7.2.0
```
Also, to double-check the environment before getting too far into the job, adding the following statements after the module load statements will double-check the sources of the two errors that have been seen. The checks can be disabled by setting DOUBLECHECK_OK to something other than 0 at the top.
```
DOUBLECHECK_OK=0
if [[ $DOUBLECHECK_OK == 0 ]] ; then
test -e ped_add_pop.py && { (( ++DOUBLECHECK_OK )); } || { echo "ped_add_pop.py
not found in current directory."; }
which smartpca >/dev/null && { (( ++DOUBLECHECK_OK )); } || { echo "smartpca
not found."; }

[[ $DOUBLECHECK_OK == 2 ]] && { echo "passed doublechecks"; } || { echo "failed
doublechecks"; exit 1; }
fi
```

After you've sorted the abovementioned, make a ```tped``` version of the dataset and you can go ahead & submit a sbatch job:

```
sbatch -A uppmax2023-2-18 -M snowy -p core -n 2 -t 07:00:00 -J PCA pca_lsqproj.sh FINAL_DATASET_PRUNED_TPED MODERN.txt UNKNOWN.txt
```
This should take a few hours. You can check whether the job is still running with ``` sacct ``` or ```jobinfo```.

When it's done, assuming it has ended sucessfuly, copy the ```.evec``` and ```.eval``` files to your computer's R directory and plot the PCA.

There are multiple ways of doing this, but if you can't think of a better one here's a suggestion:

```
library(ggplot2)

pca<-read.table("FINAL_DATASET_PRUNED_TPED.popsubset.evec")
####### DIVIDE THE TABLE INTO MODERN & UNKNOWN - make sure to edit the right numbers! 

modern<-pca[c(1:1610),]
modern$V12<-factor(modern$V12)
ancient<-pca[c(1611:1626),]

###### PLOTTING THE MODERN SAMPLES
par(mfrow=c(1,1))
layout(matrix(c(1,1), ncol=1), heights=c(1,2))
par(mar=c(4,4,1,1))
with(pca,plot(pca$V2~pca$V3,pch=20,cex=0.5,col = rgb(0, 0, 0, 0.15),xlab="PCA2",ylab="PCA1"))

###### LABELING THE PLOT CREATED WITH THE MODERN SAMPLES 

mod_lab2<-aggregate(modern[,2:3],list(modern$V12),mean)
for(j in 1:100) 
{
  text(mod_lab2$V3[j],mod_lab2$V2[j],labels=mod_lab2$Group.1[j],cex=0.7)   
}

###### SPECIFYING THE UNKOWN SAMPLES - change the numbers accordingly based on the position in the final file! 

UNK1<-pca[c(1611),]
UNK2<-pca[c(1612),]
UNK3<-pca[c(1613),]
UNK4<-pca[c(1614),]
UNK5<-pca[c(1615),]
UNK6<-pca[c(1616),]
UNK7<-pca[c(1617),]
UNK8<-pca[c(1618),]
UNK9<-pca[c(1619),]
UNK10<-pca[c(1620),]
UNK11<-pca[c(1621),]
UNK12<-pca[c(1622),]
UNK13<-pca[c(1623),]
UNK14<-pca[c(1624),]
UNK15<-pca[c(1625),]
UNK16<-pca[c(1626),]


###### PLOTTING THE UNKNOWN SAMPLES 

UNK1$V12<-factor(UNK1$V12)
points(UNK1$V3,UNK1$V2,pch=17,cex=1.4,col="burlywood",bg="burlywood")
for(k in 1:35) 
  
{
  text(UNK1$V3[k],UNK1$V2[k],labels=UNK1$V12[k],pos=c(4),cex=0.9)   
}

####################################2
UNK2$V12<-factor(UNK2$V12)
points(UNK2$V3,UNK2$V2,pch=3,cex=1.4,col="yellow3",bg="yellow3")
for(k in 1:13) 
  
{
  text(UNK2$V3[k],UNK2$V2[k],labels=UNK2$V12[k],pos=c(3),cex=0.9)   
}

####################################3
UNK3$V12<-factor(UNK3$V12)
points(UNK3$V3,UNK3$V2,pch=4,cex=1.4,col="yellow2",bg="yellow2")
for(k in 1:13) 
  
{
  text(UNK3$V3[k],UNK3$V2[k],labels=UNK3$V12[k],pos=c(3),cex=0.9)   
}
####################################4
UNK4$V12<-factor(UNK4$V12)
points(UNK4$V3,UNK4$V2,pch=5,cex=1.4,col="palevioletred3",bg="palevioletred3")
for(k in 1:13) 
  
{
  text(UNK4$V3[k],UNK4$V2[k],labels=UNK4$V12[k],pos=c(1),cex=0.9)   
}
```
Most likely you can guess where this is going, so just use this as a reference and write the R script for the rest of your unknown individuals.
After running it, you should get your nice projected PCA. Look at the output PDF, where are the unknown samples being projected? What clues does the PCA give you?

## Step 4 ADMIXTURE 

ADMIXTURE is a similar tool to STRUCTURE but runs much quicker, especially on large datasets.
ADMIXTURE runs directly from .bed or .ped files and needs no extra parameters for file preparation. You do not specify burin and repeats, ADMIXTURE exits when it converged on a solution (Delta< minimum value)

First, you have to load the module:

```
module load bioinfo-tools
module load ADMIXTURE/1.3.0
```
A basic ADMIXTURE run looks like this:

```
admixture -s time FINAL_DATASET_PRUNED.bed 2
```

This command executes the program with a seed set from system clock time, it gives the input file (remember the extension) and the K value at which to run ADMIXTURE (2 in the previous command).

For ADMIXTURE you also need to run many iterations at each K value, thus a compute cluster and some scripting is useful.

Make a script from the code below to run Admixture for K = 2-6 with 3 iterations at each K value.

```
#!/bin/bash
for kval in {2..6}; do
   for itir in {1..3}; do
           (echo '#!/bin/bash -l'
           echo "
           mkdir ${kval}.${itir}
           cd ${kval}.${itir}
#working folder where admixture and the bed files are
module load bioinfo-tools
module load ADMIXTURE/1.3.0
admixture -j3 -s $RANDOM path/to/your/FINAL_DATASET_PRUNED.bed ${kval}
cd ../
rm -r ${kval}.${int}
#moves it to a different folder and renames it so you end up with multiple iterations
exit 0") |
               sbatch -M snowy -p core -n 3 -t 72:0:0 -A YOUR_UPPMAX_PROJECT_HERE -J ADMX.${kval}.${itir} -o ADMX.${kval}.${itir}.output -e ADMX${kval}.${itir}.output
       done
done
```

You will need to **replace the project name** and the **path** to your pruned .bed-file in the script then just run it:

```
bash NAME_OF_THE_SCRIPT.sh
```
After a successful Admixture run, you should be seeing a bunch of new folders being created. 
In each one of them you will find the output ```P``` and ```Q``` for each iteration. 

*If you decide to work on PONG locally, you need to edit each ```Q``` file by adding the right k & iteration number to each one. After doing this for all folders, you can move all the ```Q``` files to a directory of your liking.*

## Step 5 PONG 

When the admixture analysis is finally done we use PONG to combine the different iterations and visualize our results. 

To be able to run PONG we thus need to generate **three different files**.

The first being the ***filemap***. This is the only input that is strictly required to run PONG. It consists of three columns.
From the PONG manual: 


```
Column 1. The runID, a unique label for the Q matrix (e.g. the string “run5_K7”).

Column 2. The K value for the Q matrix. Each value of K between Kmin and Kmax must
be represented by at least one Q matrix in the filemap; if not, pong will abort.

Column 3. The path to the Q matrix, relative to the location of the filemap. 
```

In order to create what we need we can edit and run the loop below. Make sure to edit the right number of **k** and **i**.
Also check whether the prefix of the admixture output Q files is the same in both the loop and the Q files. As long as you are starting PONG in the same directory of all the Admixture output files you should be ok.

```
for i in {2..7};
do
    for j in {1..3};
    do
    echo -e "k${i}_r${j}\t${i}\t${i}.${j}/FINAL_DATASET_PRUNED.${i}.Q" >> unknown_FILEMAP.txt
    done
done

```
Example of what this file looks like, depending on how many K and how many iterations you run:

```
k2_r1	2	2.1/FINAL_DATASET_PRUNED.2.Q
k2_r2	2	2.2/FINAL_DATASET_PRUNED.2.Q
k2_r3	2	2.3/FINAL_DATASET_PRUNED.2.Q
k3_r1	3	3.1/FINAL_DATASET_PRUNED.3.Q
k3_r2	3	3.2/FINAL_DATASET_PRUNED.3.Q
k3_r3	3	3.3/FINAL_DATASET_PRUNED.3.Q
k4_r1	4	4.1/FINAL_DATASET_PRUNED.4.Q
k4_r2	4	4.2/FINAL_DATASET_PRUNED.4.Q
k4_r3	4	4.3/FINAL_DATASET_PRUNED.4.Q
k5_r1	5	5.1/FINAL_DATASET_PRUNED.5.Q
k5_r2	5	5.2/FINAL_DATASET_PRUNED.5.Q
k5_r3	5	5.3/FINAL_DATASET_PRUNED.5.Q
k6_r1	6	6.1/FINAL_DATASET_PRUNED.6.Q
k6_r2	6	6.2/FINAL_DATASET_PRUNED.6.Q
k6_r3	6	6.3/FINAL_DATASET_PRUNED.6.Q
k7_r1	7	7.1/FINAL_DATASET_PRUNED.7.Q
k7_r2	7	7.2/FINAL_DATASET_PRUNED.7.Q
k7_r3	7	7.3/FINAL_DATASET_PRUNED.7.Q
```

The next file we need to create is the ***ind2pop*** file. It is just a list of which population each individual belongs to.
We have this information in the `.fam` file so we can just cut out the field we need:

```
cut -f 1 -d " " FINAL_DATASET_PRUNED.fam > unknown_IND2POP.txt

``` 
The ind2pop file should look exactly like the first column of the fam.

The ***poporder*** file is a key between what your populations are called and what "Proper" name you want to show up in your final plot. 
Also as the name suggests, **gives the order in which you want the populations to be in**, so order them as you think it makes most sense.

*Note that the file needs to be **tab-delimited**, i.e separated by tabs (If you have spaces issues may occure). 
Also, make sure within the names in this file to not have any spaces. If you do, try replacing them with underscores (_).*
*Again, be careful when editing this file, typos, extra spaces etc might cause problems.*

Below is just an example of what a poporder file should look like:

```
aDNA_Rick	pickle
aDNA_Morty	m0rty
aDNA_Ragnar	rgnr
aDNA_martians	mrt99
Napoleon_Bonaparte	nb001
Netherlands_BA	Dutch
```
An easy way of doing this is by :
```
cut -f 1 -d " " FINAL_DATASET_PRUNED.fam | uniq > unknown_POPORDER_prep.txt
```
and then just use this simple script:
```
#!usr/bin/env bash

file=$(cat unknown_POPORDER_prep.txt)

for line in $file
do
    echo -e "$line\t$line" >> unknown_POPORDER.txt
done
```
This way you will get a file with two columns. You can edit the the names of the populations you are interested in by changing the names of the second column.

### Now we have all the files we need. Time to run PONG.

*It is by far best to run PONG on your local computer. All you need to do is install PONG and copy the ```.Q``` files to your computer together with the three files (filemap,ind2pop,poporder) and run it. If you don't know how to install PONG(https://github.com/ramachandran-lab/pong), it is available through the module system on Uppmax, although very often it is quite laggy. 

* Mac Users - you will need to install XQuartz (an open-source version of the X.Org X server, a component of the X Window System that runs on macOS)https://www.xquartz.org
* Windows Users - MobaXterm https://mobaxterm.mobatek.net
* Ubuntu Users - I'm sure you can figure it out yourselves

```
module load pong 
```
**Super important:**
Since we are several people who are going to run PONG at the same time we need to use a different port, otherwise, we will collide with each other. The default port for PONG is 4000. Any other free port will work, like 4001, 2, etc. Make sure you are using a **unique port** before proceeding. If multiple people are trying to run PONG on the same port you have problems, so talk to your classmates and find a unique port in the 4000s that people aren't using.

```
pong -m unknown_FILEMAP.txt -i unknown_IND2POP.txt -n unknown_POPORDER.txt -g --port YOUR_PORT_NUMBER_HERE
```

When PONG is done it will start hosting a webserver that displays the results at port 4000 by default:  http://localhost:4000. Pong needs to be running for you to look at the results, i.e. **if you close it it will not work.**
It will look like this:

```
-------------------------------------------------------------------
                            p o n g
      by A. Behr, K. Liu, T. Devlin, G. Liu-Fang, and S. Ramachandran
                       Version 1.5 (2021)
-------------------------------------------------------------------
-------------------------------------------------------------------

Parsing input and generating cluster network graph
Matching clusters within each K and finding representative runs
For K=2, there are 3 modes across 3 runs.
For K=3, there are 3 modes across 3 runs.
For K=4, there are 3 modes across 3 runs.
For K=5, there are 3 modes across 3 runs.
For K=6, there are 3 modes across 3 runs.
For K=7, there are 3 modes across 3 runs.
Matching clusters across K
Finding best alignment for all runs within and across K
match time: 0.64s
align time: 0.01s
total time: 1.21s
-----------------------------------------------------------
pong server is now running locally & listening on port 4005
Open your web browser and navigate to http://localhost:4005 to see the visualization
```

The web server is hosted on the **same login node as you were running pong**. In case you are unsure of which one that is you can use the command `hostname` to figure it out:

```
hostname
```

To view files interactively you need to have an X11 connection. So when you connect to rackham from a new tab do:

```
ssh -AY YOUR_USERNAME_HERE@rackham.uppmax.uu.se
```

Make sure that you connect to **the same Rackham (i.e 1,2,3 etc) as you got from hostname**.  
In a **new tab** (if you didn't put PONG in the background) type:

```
firefox http://localhost:YOUR_PORT_NUMBER
```


Once you are finished with looking at the PONG output you can click and save some output to file and then close the program by hitting `ctrl - c` in the terminal tab running it. 

### What do the results from Admixture & PCA tell you? 
What do you see? What information does this give you about your misterious unknown populations? 
Do the results of your population PCA correspond to the population structure results you got from the ADMIXTURE plots? Can you figure out who they are?

**After you have finished with your analyses, delete all unnecessary files from your folder.**

# Part 3: Review & reflect

Okay, hopefully this project wasn't too difficult, too boring or too long. Maybe it was even fun? And aside from that maybe you even learnt something. Regardless of how you felt about it, this project was constructed to provide you with a **real life** scenario of what most exploratory analysis our field of study look like, warts & all. More times than not you are working with something you don't know enough about (duh! - that's why it's called research) using techniques you're trying for the first time and failing on many attempts. Most of the manuals & how-toos you need to look for yourself, and references come in the shape of many scattered papers on the subject. And as you would do when you want to publish your findings in real life, the final step of the project is to write a 10-15 page report on your findings. Make sure to read up on the papers recommended above. They should provide you with enough of a backstory & will surely help you a lot in constructing the grand picture of things. You are free to also look into similar papers on the subject. Good luck & happy writing!
 
