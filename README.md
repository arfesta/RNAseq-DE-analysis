# RNASeq DE Analysis

*Differential Expression analysis of RNAseq data in Loblolly Pine*

- [RNASeq DE Analysis](#RNASeq DE Analysis)
  * [Abstract](#abstract)
  * [Background of samples](#background-of-samples)
  * [Location of data](#location-of-data)
  * [Analyses](#analyses)
    + [Step 1 - Data Prep](#step-1---data-prep)
    + [Step 2 - Load Count Data](#step-2---load-count-data)

## Abstract

* The NCSU TIP has been established for over 60 years with the objective of increasing genetic value of family-level selections for production deployment. Given the size and resources required to conduct tree breeding, experiments and selections of families used for deployment occur across thousands of acres in the Southeastern United States. 

* ***This study examines the number and magnitude of differentially expressed transcripts among a set of 70 families which span across a wide geographic range.***  Additionally, for those families which have phenotypes, we assess possible transcripts that may be differentially expressed between extreme ends of the phenotypic distribution.

## Background of samples

* LGEP
   
   - 144 Biological Replicates x 3 technical replicates = 432 technical replicates

* EW

   - 80 Biological Replicates x 3 technical replicates = 240 technical replicates

## Location of data

Data Subject Type | Data File Type | Path | Notes
--- | --- | --- | ---
**raw read files**  | | `/media/disk6/ARF/RNASEQ/rawreads/86kSalmon` | Raw files returned from GSL
| | *raw tar* | `./EWtarfiles or ./LGEPtarfiles`
| | *raw fasta* | `./EWfasta or ./LGEPfasta` 
**trimmed and filtered read files** | |`/media/disk6/ARF/RNASEQ/trimmedfiltreads/86k` | Files post trim & adapater removal
|  |*EW* | `./EW/lane01 ... ./lane12` | 
|  |*LGEP* | `./LGEP/lane01 ... ./lane18` | 
**salmon count files** | |`/media/disk6/ARF/RNASEQ/counts/86kSalmon` | Direcotries containing quant.sf files
|  |*EW tech reps* | `./EW/lane01 ... ./lane12` | 
|  |*LGEP tech reps* | `./LGEP/lane01 ... ./lane18` |
|  |*EW bio reps* | `./bio_EW/Sample_<animal_id>/` | 
|  |*LGEP bio reps* | `./bio_LGEP/Sample_<animal_id>/` | 
**experimental data resources**  | | `/media/disk6/ARF/RNASEQ/Breeding-Value-Prediction/disk6directory/resources` | Experiment information
|  |*sequencing* | `./exptdesign/sequencing` | `./EWtarfiles or ./LGEPtarfiles`
|  |*pedigree* | `./pedigree` | `./EWfasta or ./LGEPfasta`
|  |*phenotypes* | `./phenos` | `./EWfasta or ./LGEPfasta`


## Analyses

### Step 1 - Data Prep

   Data prep includes everything from unpacking the original tar files recieved by the GSL up to estimating transcript               
      abundance with Salmon. Additionally, this step includes identification of the indicies used within both batches and creates an experimental info matrix containing all meta data from both batches.
      
   See the [raw reads README](https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/rawreads/README.md) for step by step processing of files.  

### Step 2 - Load Count Data
      
   Once counts have been estimated, the next step involves reading in the aligned biological replicate counts using the tximport package.
      
   Additionally, the phenotype and other sample meta-data is constructed for normalization.
   
   To see this process for **biological reps**, navigate to: [load counts bio rep html file](http://htmlpreview.github.com/?https://github.com/arfesta/Breeding-Value-Prediction/blob/master/disk6directory/analyses/step2.loadcounts/load.counts.html) which contains the complete markdown and output.

