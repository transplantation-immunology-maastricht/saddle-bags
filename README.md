#Bhast
Ben's HLA Allele Submission Tool
A tool for generating an EMBL-formatted submission of a standard novel HLA allele. 

##Download the executable
When the compiled executable is available, i'll put a link for the file download here.

##Run using Python
Run this program using standard python call:  
```
python AlleleSubmissionEMBL.py
```

##To configure Anaconda
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

##Run using a bash or .bat script using anaconda
You can execute the following scripts to run the tool inside of Anaconda:  
Linux/Mac:  
```
bash Run_allele_submission.sh  
```
Windows:  
```
Run_allele_submission.bat
```

##Input Data
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

##Output Data
The resulting report is in the form of an EMBL HLA Novel Allele submission flatfile.  You can submit this to EMBL as a new HLA allele.

##For more information on EMBL's ENA format:  
http://www.ebi.ac.uk/ena/submit/sequence-submission  
http://www.ebi.ac.uk/ena/submit/entry-upload-templates  
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt  
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html  
http://www.ebi.ac.uk/ena/software/flat-file-validator  
