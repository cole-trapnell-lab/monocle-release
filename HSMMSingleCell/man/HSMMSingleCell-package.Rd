\name{HSMMSingleCell-package}
\alias{HSMMSingleCell-package}
\alias{HSMMSingleCell}
\docType{package}
\title{
A single-cell RNA-Seq study of differentiating human skeletal muscle myoblasts  
}
\description{
Skeletal myoblasts undergo a well-characterized sequence of morphological and transcriptional changes during differentiation. In this experiment, primary human skeletal muscle myoblasts (HSMM) were expanded under high mitogen conditions (GM) and then differentiated by switching to low-mitogen media (DM).  RNA-Seq libraries were sequenced from each of several hundred cells taken over a time-course of serum-induced differentiation. Between 49 and 77 cells were captured at each of four time points (0, 24, 48, 72 hours) following serum switch using the Fluidigm C1 microfluidic system. RNA from each cell was isolated and used to construct mRNA-Seq libraries, which were then sequenced to a depth of approximately 4 million reads per library, resulting in a complete gene expression profile for each cell. 

For single-cell mRNA sequencing, dissociated cells were captured and processed with the C1 Single-Cell Auto Prep System (Fluidigm) following manufacturer's protocol 100-5950.  
Starting with a suspension of cells at a concentration of approximately 250 cells per microliter, up to 96 single cells are captured in each C1 microfluidic device. 
In this study, we used one C1 capture chip at 0, 24, 48, and 72 hours after switching to differentiation medium, for a total of four independent captures. 
After imaging with a microscope to identify which sites have captured a single cell, processing of the cells occurs within the C1 instrument to perform the steps of cell lysis, cDNA synthesis with reverse transcriptase, and PCR amplification of each cDNA library.  
The cDNA synthesis and PCR use reagents from the SMARTer Ultra Low RNA Kit for Illumina Sequencing (Clontech 634936).  
The SMARTer chemistry utilizes a strand-switching mechanism so that both the 1st and 2nd strands of cDNA are synthesized in a single reaction.  Following harvest from the C1 microfluidic device, each cDNA library is subjected to tagmentation (simultaneous fragmentation and tagging with sequencing adapters) using the Nextera XT DNA Sample Preparation Kit (Illumina FC-131-1096) as described in Fluidigm protocol 100-5950.  
PCR amplification of the tagmented cDNA uses Index Primers from the Nextera XT DNA Sample Preparation Index Kit (Illumina FC-131-1002).  
Following PCR, the cDNA libraries from individual cells are pooled and purified using AMPure XP beads (Agencourt Bioscience Corp A63880) as described in Fluidigm protocol 100-5950.
Libraries that contained fewer than 1 million reads or for which less than 80\% of fragments mapped to non-mitochondrial protein coding genes were excluded.
}
\details{
\tabular{ll}{
Package: \tab monocle\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2013-11-19\cr
License: \tab Artistic-2.0\cr
}

}
\author{
Cole Trapnell
Maintainer: Cole Trapnell <cole@broadinstitute.org>
}
\references{
Trapnell, Cacchiarelli, et al
}
\keyword{ package }
