# RNASeq DE Analysis

*Differential Expression analysis of RNAseq data in Loblolly Pine*

  * [Abstract](#abstract)
  * [Background of samples](#background-of-samples)
  * [Location of data](#location-of-data)
  * [Methods](#methods)
    + [Step 1 - Data Prep](#step-1---data-prep)
    + [Step 2 - Load Count Data](#step-2---load-count-data)
    + [Step 3 - Conduct pairwise DE tests among families](#step-3---conduct-pairwise-de-tests-among-families)

## Abstract

* The NCSU TIP has been established for over 60 years with the objective of increasing genetic value of family-level selections for production deployment. Given the size and resources required to conduct tree breeding, experiments and selections of families used for deployment occur across thousands of acres in the Southeastern United States. 

* ***This study examines the number and magnitude of differentially expressed transcripts among a set of 50+ families which span across a wide geographic range.**  Additionally, for those families which have phenotypes, we assess possible transcripts that may be differentially expressed between extreme ends of the phenotypic distribution.

## Background of samples

* A total of 144 biological replicates were grown from a group of families subset from the lower gulf elite population (LGEP) and 80 bioloigcal replicates were grown from the east-west diallel (EW) experiment.

* LGEP
   
   - 144 Biological Replicates x 3 technical replicates = 432 technical replicates

* EW

   - 80 Biological Replicates x 3 technical replicates = 240 technical replicates

* Note that some samples from the EW had to be dropped after sequencing due to contamination, so the total number of bioloigcal replicates analyzed do not add up to 80+144.

## Location of data

* All raw data are stored in the same locations as noted in the breeding value prediction project but are provided again here for consistancy.

Data Subject Type | Data File Type | Path | Notes
--- | --- | --- | ---
**raw read files**  | | `/media/disk6/ARF/RNASEQ/shared/rawreads/86kSalmon` | Raw files returned from GSL
| | *raw tar* | `./EWtarfiles or ./LGEPtarfiles`
| | *raw fasta* | `./EWfasta or ./LGEPfasta` 
**trimmed and filtered read files** | |`/media/disk6/ARF/RNASEQ/shared/trimmedfiltreads/86k` | Files post trim & adapater removal
|  |*EW* | `./EW/lane01 ... ./lane12` | 
|  |*LGEP* | `./LGEP/lane01 ... ./lane18` | 
**salmon count files** | |`/media/disk6/ARF/RNASEQ/shared/counts/86kSalmon` | Direcotries containing quant.sf files
|  |*EW bio reps* | `./bio_EW/Sample_<animal_id>/` | 
|  |*LGEP bio reps* | `./bio_LGEP/Sample_<animal_id>/` | 
**experimental data resources**  | | `/media/disk6/ARF/RNASEQ/shared/resources` | Experiment information
|  |*sequencing* | `./exptdesign/sequencing` | 
|  |*pedigree* | `./pedigree` | 
|  |*phenotypes* | `./phenos` | 


* All scripts and markdown files relating to this project are stored within this repository.

## Methods

### Step 1 - Data Prep

   Data prep includes everything from unpacking the original tar files recieved by the GSL up to estimating transcript abundance with Salmon. Additionally, this step includes identification of the indicies used within both batches and creates an experimental info matrix containing all meta data from both batches.
      
   See the [raw reads README](https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/rawreads/README.md) for step by step processing of files.  

### Step 2 - Load Count Data
      
   Once counts have been estimated, the next step involves reading in the aligned biological replicate counts using the tximport package.
      
   Additionally, the phenotype and other sample meta-data is constructed for normalization.
   
   To see this process for **biological reps**, navigate to: [load counts bio rep html file](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step2.loadcounts/load.counts.html) which contains the complete markdown and output.


### Step3 - Conduct pairwise DE tests among families

  An initial filter was applied so that only families which had **3 or more** bioloigcal replicates were kept for further analysis.  The 57 families which passed this threshold were then used with the R package `DESeq2` to test for differential epxression among all pairs of families.
  
 The result of this step is a list of 1596 comparisons.  Each comparison contains the results from DESeq2 between family X vs. family Y.
 
 [Generate pairwise DE tests: markdown file](https://github.com/arfesta/RNAseq-DE-analysis/blob/master/analyses/DESEQ2_bio.0.Rmd)
 
 
### Step4 - Conduct DE tests among State of origin

  Similar to conducting pairwise tests among all families, the same filter was applied so that only families which had **3 or more** bioloigcal replicates were kept for further analysis. However, an additional filter was applied here to remove all families which were not open-pollinated. The 41 families which passed this threshold were then used with the R package `DESeq2` to test for differential epxression among the following states `NC, SC, GA, FL, TX`
  
 The result of this step is a list of 840 comparisons.  Each comparison contains the results from DESeq2 between family X vs. family Y.
 
 [Generate State origin DE tests: markdown file](https://github.com/arfesta/RNAseq-DE-analysis/blob/master/analyses/DESEQ2_bio.0_states.Rmd)
 
 ### Step5 - Conduct DE contrast between high and low volume

  Families which had **3 or more** bioloigcal replicates and contained a volume breeding value were kept for this analysis.  The bottom and upper quartile of the volume phenotypic distribution were used with DESeq2 to generate a contrast between high and low volume breeding value families.
  
 The result of this step is a single comparison between 10 low volume families and 11 high volume families.
 
 [Generate High vs Low vol DE test: markdown file](https://github.com/arfesta/RNAseq-DE-analysis/blob/master/analyses/DESEQ2_high.vs.low.vol.Rmd)
