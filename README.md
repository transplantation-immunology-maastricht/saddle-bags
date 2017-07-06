# Saddlebags
A tool for generating allele submissions for the EMBL and IMGT nucleotide databases. 

## Download the executable
Download an executable for Windows from the Release page.
[Download Saddlebags for Windows Here](https://github.com/transplantation-immunology/saddle-bags/releases)

## Run using Python
Alternatively, you can run this program using Python 2.7. This works for Mac and Linux users (or Windows). There are prerequesites, I recommend you install them inside an Anaconda environment. See Run_allele_submission.sh and Run_allele_submission.bat for an example of this in Linux and Windows environments, respectively.

```
python AlleleSubmissionMain.py
```

## To configure Anaconda
Anaconda uses separate environments to run your programs in.  
Install Anaconda for python 2.7.  
https://www.continuum.io/downloads  
To set up the environment in anaconda:  

Linux/Mac:  
```
conda create --name AlleleSubEnvironment biopython six pycurl 
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
Identifying the exons in your HLA sequence is a nontrivial challenge. We have provided a list of common sequences surrounding exon boundaries. See the [Release Page](https://github.com/transplantation-immunology/saddle-bags/releases) for a .pdf reference.  These sequences can probably be found within your HLA consensus sequence at exon boundary sites. It may also help to use the [IMGT/HLA sequence alignment tool](http://www.ebi.ac.uk/ipd/imgt/hla/align.html) for more information on common exon patterns.

## EMBL Metadata Input format
Sample ID: Specified by the submitting laboratory, you may use a value that is informative to you.
Gene:
Class I or II:
Allele Local Name:

Saddlebags will submit to EMBL Test environment by default, you must specify that the software target the Live / Production environment.

EMBL sequence submissions must be associated with a Study/Project. You may specify the accession number of an existing EMBL study (Get this accession number from [EMBL Webin](https://www.ebi.ac.uk/ena/submit/sra/#home) ), or Saddlebags can create a new project to your specifications.


## For more information on EMBL's ENA format:  
http://www.ebi.ac.uk/ena/submit/sequence-submission  
http://www.ebi.ac.uk/ena/submit/entry-upload-templates  
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt  
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html  
http://www.ebi.ac.uk/ena/software/flat-file-validator  

## For more information on IMGT metadata:

TODO: Put a description of the IMGT metadata form. There is lots of information that goes in here, and much of it is confusing.  How should Primers and Sequencing methodology be provided? What are the options for Ethnic Origin or Sex or Cosanguineous?
