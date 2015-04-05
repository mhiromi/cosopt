# cosopt
Bioinformatics / detection of periodic patterns in gene expression data


# Release Info
Realsed in March 8, 2015.

#Introduction
COSOPT is a statistical method to distinguish a period and phase for a transcript in time-series gene expression analysis.  The source code is provided by C language and R. 

IMPORTANT: This COSOPT has a small difference from the original COSOPT described in Straume et al. 2004. In the original paper, SD within replicates of a probe was used to generate surrogates with white noise,  however we fixed SD to 0.1 in our analysis. If you want to change SD values yourself in the R package, please modify <sigma> variable as follows.

# An example in R 

timepoints <- c(0,8,16,24,32,40,48)     # time points of gene expression data

data <- c(-1,0,1,0,-1,0,1)  # signal intensity of each probe in gene expression data

sigma <- c(0.1,0.1,0.1,0.1,0.1,0.1,0.1)  # sigma (standard deviation of each probe in the original source code) 

cosopt(data,sigma,timepoints,plotting=TRUE) # run COSOPT

# Citation 
Please cite the following paper if you used this program.

* Hiromi Matsumae, Ryosuke Ishiwata, Toshifumi Minamoto, Norio Ishida, Soichi Ogishima,  and Hiroshi Tanaka. <I> submitted </I>.


# Contact 
Please report a bug and ask a question to the following contacts. 

* Hiromi Matsumae: hiromi [dot] matsumae [atmark] gmail [dot] com
* Ryosuke Ishiwata: ishiwata [atmark] phys [dot] cs [dot] is [dot] nagoya-u [dot] ac [dot] jp

