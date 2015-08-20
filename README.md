# cosopt
Bioinformatics / detection of periodic patterns in gene expression data


## Release Info
Realsed in April 6, 2015.
Updated in Augst 20, 2015. Validated on Mac OS 10.8.5 and R 3.2.0. 

##Introduction
COSOPT is a statistical method to distinguish a period and phase for a transcript in time-series gene expression analysis.  The source code is provided by C language and R. 

IMPORTANT: This COSOPT has a small difference from the original COSOPT described in Straume et al. 2004. In the original paper, SD within replicates of a probe was used to generate surrogates with white noise,  however we fixed SD to 0.1 in our analysis. If you want to change SD values yourself in the R package, please modify <sigma> variable as follows.


## Source Code

### original source code written in C
The code is in code_R/.

### Running COSOPT on R
It is recommended to run COSOPT against a small number of samples on R.

timepoints <- c(0,8,16,24,32,40,48)     # time points of gene expression data

data <- c(-1,0,1,0,-1,0,1)  # signal intensity of each probe in gene expression data

sigma <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1)  # sigma (standard deviation of each probe in the original source code) 

cosopt(data,sigma,timepoints,plotting=TRUE) # run COSOPT

### Running in C with MPI environment
Here is MPI_COSOPT which is optimized for parallel computation, MPI. Running MPI_COSOPT requires knowledge of parallel computing. You can also run MPI_COSOPT on SunGridEngine(SGE/UGE). 

## Supplementary data
All supplementary data are available from [the github wiki page](https://github.com/mhiromi/cosopt/wiki/Supplementary-data).


## Citation 
Please cite the following paper if you used this program.

* Hiromi Matsumae, Ryosuke Ishiwata, Toshifumi Minamoto, Norio Ishida, Soichi Ogishima, and Hiroshi Tanaka. <I> submitted</I>.


## Contact 
Please report a bug and ask a question to the following contacts. 

* Hiromi Matsumae: hiromi [dot] matsumae [atmark] gmail [dot] com
* Ryosuke Ishiwata: ishiwata [atmark] phys [dot] cs [dot] is [dot] nagoya-u [dot] ac [dot] jp

