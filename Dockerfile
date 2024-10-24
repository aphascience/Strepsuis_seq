FROM ubuntu:24.04

####################### METADATA ##########################
LABEL Maintainer="Richard Ellis <richard.ellis@apha.gov.uk>"
LABEL base.image=ubuntu:24.04
LABEL software="Strepsuis_seq"
LABEL about.documentation="https://github.com/APHA-CSU/Strepsuis_seq"

################## INSTALL DEPENDANCIES ###################

## intsall trimming tool (bbduk?)...
## install bowtie
## install samtools 
## install srst2 (need to update to py3...)