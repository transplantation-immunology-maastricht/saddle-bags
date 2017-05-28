# Bhast
Ben's HLA Allele Submission Tool

A tool for generating an EMBL-formatted submission of a standard novel HLA allele. 

## Download the executable
Download an executable for Windows from the Release page.
[Download Bhast for Windows Here](https://github.com/transplantation-immunology/EMBL-HLA-Submission/releases)

## Run using Python
Alternatively, you can run this program using Python 2.7. This works for Mac and Linux users (or Windows). There are prerequesites, you can install them inside an Anaconda environment.
```
python AlleleSubmissionEMBL.py
```

## To configure Anaconda
Anaconda uses separate environments to run your programs in.  
Install Anaconda for python 2.7.  
https://www.continuum.io/downloads  
To set up the environment in anaconda:  

Linux/Mac:  
```
conda create --name AlleleSubEnvironment biopython six  
source activate AlleleSubEnvironment  
pip install pyinstaller packaging  
source deactivate  
```  
Windows:  
```  
conda create --name AlleleSubEnvironment biopython six pywin32  
call activate AlleleSubEnvironment && pip install pyinstaller packaging && call deactivate  
```

## Run using a bash or .bat script using anaconda
You can execute the following scripts to run the tool inside of Anaconda:  
Linux/Mac:  
```
bash Run_allele_submission.sh  
```
Windows:  
```
Run_allele_submission.bat
```

## Input Data
Input data sequence must consist of  
1) 5' UTR  
2) Any number of Exons, with an introns distributed between them.  
3) 3' UTR

Introns, Exons, and UTRs are to be distinguished using capital or lowercase nuclotide letters.  
Exons are capital, while introns and UTRs are lowercase.

Like this:  
fiveprimeutrEXONONEintrononeEXONTWOintrontwoEXONTHREEthreeprimeutr  
agctagctagctAGCTAGCtagctagctAGCTAGCtagctagctAGCTAGCTAgctagctagctag

All spaces, tabs, and newlines are not interpreted by this program,  They are removed.
If you wish, you may use those characters for convenient visualization of splice sites, like this:

agctagctagct  
AGCTAGC  
tagctagct  
AGCTAGC  
tagctagct  
AGCTAGCTA  
gctagctagctag  

## Annotating your exons
Identifying the exons in your HLA sequence is a nontrivial challenge. We have provided a list of common sequences surrounding exon boundaries. See the [Release Page](https://github.com/transplantation-immunology/EMBL-HLA-Submission/releases) for a .pdf reference.  These sequences can probably be found within your HLA consensus sequence at exon boundary sites. It may also help to use the [IMGT/HLA sequence alignment tool](http://www.ebi.ac.uk/ipd/imgt/hla/align.html) for more information on common exon patterns.

## Output Data
The resulting report is in the form of an EMBL HLA Novel Allele submission flatfile.  You can submit this to EMBL as a new HLA allele



## For more information on EMBL's ENA format:  
http://www.ebi.ac.uk/ena/submit/sequence-submission  
http://www.ebi.ac.uk/ena/submit/entry-upload-templates  
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt  
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html  
http://www.ebi.ac.uk/ena/software/flat-file-validator  
