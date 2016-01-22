# cosopt
Bioinformatics / detection of periodic patterns in gene expression data


## Release Info
Realsed in April 6, 2015.
Updated in November 18, 2015. Validated on Mac OS 10.8.5 and R 3.2.0. 

##Introduction
COSOPT is a statistical method to distinguish a period and phase for a transcript in time-series gene expression analysis.  The source code is provided by C language and R. 

IMPORTANT: This COSOPT has a small difference from the original COSOPT described in Straume et al. 2004. In the original paper, SD within replicates of a probe was used to generate surrogates with white noise,  however we fixed SD to 0.1 in our analysis. If you want to change SD values yourself in the R package, please modify <sigma> variable as follows.


## R version 
### Installation 

You can find a compressed package such as "cosopt_0.3.tar.gz" in code_R/ directory. Just use install.packages command as below: 

```R:
> install.packages("cosopt_0.3.tar.gz",repos = NULL, type = "source")
```

### Running 

We recommend to run COSOPT against a small number of samples on R.

```R:
> timepoints <- c(0,8,16,24,32,40,48)     # time points of gene expression data
> data <- c(-1,0,1,0,-1,0,1)  # signal intensity of each probe in gene expression data
> sigma <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1)  # sigma (standard deviation of each probe in the original source code) 
> cosopt(data,sigma,timepoints,plotting=TRUE) # run COSOPT with plotting image
gene=1 timecourse=7
plotting=1
minimum_freq=0.031250
```

## C version 
The original source code is written in C. You can find it from  code_C/.

### Running COSOPT by C with MPI environment
Here is MPI_COSOPT which is optimized for parallel computation, MPI. Running MPI_COSOPT requires knowledge of parallel computing. You can also run MPI_COSOPT on SunGridEngine(SGE/UGE). 

## Supplementary data
All supplementary data are available from [the github wiki page](https://github.com/mhiromi/cosopt/wiki/Supplementary-data) or our paper. 


## Citation 
Please cite the following paper if you used this program.

* Hiromi Matsumae, Ryosuke Ishiwata, Toshifumi Minamoto, Norio Ishida, Soichi Ogishima, and Hiroshi Tanaka. "Detection of periodic patterns in microarray data reveals novel oscillating transcripts of biological rhythms in <I>Ciona intestinalis</I>" 2015. Artificial Life and Robotics, Springer. doi:10.​1007/​s10015-015-0237-6
Link to the paper: [http://link.springer.com/article/10.1007/s10015-015-0237-6](http://link.springer.com/article/10.1007/s10015-015-0237-6)

## Contact 
Please report a bug and ask a question to the following contacts. 

* Hiromi Matsumae: hiromi [dot] matsumae [atmark] gmail [dot] com
* Ryosuke Ishiwata: ishiwata [atmark] phys [dot] cs [dot] is [dot] nagoya-u [dot] ac [dot] jp

