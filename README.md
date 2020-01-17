# Saddlebags
A tool for generating allele submissions for the EMBL nucleotide database. 

## Download the executable
Download an executable for Windows from the Release page.
[Download Saddlebags for Windows Here](https://github.com/transplantation-immunology-maastricht/saddle-bags/releases)

This executable requires java runtime in order to run. The EMBL ENA submission is performed using java.

## Run using Python
Alternatively, you can run this program using Python 3.6. The required packages can be installed within a python virtual environment. This works for Mac and Linux users, or Windows. See Run_allele_submission.bat and Run_allele_submission.sh for an example of launching Saddlebags in Windows, or Linux/Mac, respectively..

## To configure Virtual Environment.
The general installation instructions are similar for Windows and Linux/Mac environments.

Install the Git commandline. In Windows, there is an option to include the git bash console, I recommend doing this, as it is very useful.

Install Python 3.6, 64 or 32 bit.
[Download Python](https://www.python.org/downloads/release/python-363/)

The python installer includes pip, which is necessary for the next steps.

Install Virtual Environment. This requires a couple steps, because there is a windows wrapper to get virtual environment working with powershell. 

```
pip install virtualenv
pip install virtualenvwrapper-win
```
[More info](http://timmyreilly.azurewebsites.net/python-pip-virtualenv-installation-on-windows/)

Create the virtual environment. I typically use an environment called "minionvenv":
```
mkvirtualenv minionvenv
```
The mkvirtualenv command should automatically activate the new environment. If not, activate the environment.
```
C:\Users\ben\Envs\minionvenv\Scripts\activate
```

Is your environment activated? Then use pip to install the packages that saddlebags needs:
```
pip install biopython six pywin32 pyinstaller packaging pycurl
```

I found that Installing pycurl was a bit difficult inside of virtualenv.

In Ubuntu:

https://stackoverflow.com/questions/37669428/error-in-installation-pycurl-7-19-0
sudo apt-get install libgnutls-dev

In linux, i needed Needed the "dev" verson of python 3.6
sudo apt-get install python3.6-dev

In windows:

I downloaded the pip wheel file manually, and installed it. Wheel file can be found here:
https://www.lfd.uci.edu/~gohlke/pythonlibs/
and I installed manually using pip:

```
cd C:\Curl
pip install pycurl-(...).whl
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

## Logging
Logging is disabled by default, but it is possible to turn on the verbose logging feature by changing a configuration value.

Find the file "Saddlebags.Config.xml". When you exit Saddlebags, it creates this configuration file to store commonly used values. It should be placed in your home directory, in a subfolder called "saddlebags_temp".

For example, you might find it in one of these locations in Windows or Ubuntu, respetively:
```
C:\Users\ben\saddlebags_temp
/home/ben/saddlebags_temp
```
Early verions of saddlebags created configs in slightly different locations. If there is any confusion, it is safe to delete Saddlebags.Config.xml and allow Saddlebags to recreate it.

Edit Saddlebags.Config.xml, using Notepad or your favorite text editor. Turn logging on by changing the 0 to a 1, like this

from
```
<logging>0</logging>
```
to
```
<logging>1</logging>
```
Save Saddlebags.Config.xml, and restart Saddlebags. It will now create a log file in the same directory, like this:
```
C:\Users\ben\saddlebags_temp\Saddlebags.Log.txt
/home/ben/saddlebags_temp/Saddlebags.Log.txt
```
The text of this log file may help with identifying problems.

## Annotating your exons
This service uses the [NMDP BeTheMatch Allele Calling Tool](https://github.com/nmdp-bioinformatics/service-act) to automatically annotate the genomic features of a full length HLA sequence. 

If the service is unavailable, or you are annotating nonstandard sequences, it may be necessary to identify your genomic features manually. We have provided a list of common sequences surrounding exon boundaries to assist with this. See the [Release Page](https://github.com/transplantation-immunology/saddle-bags/releases) for a .pdf reference.  These sequences can probably be found within your HLA consensus sequence at exon boundary sites. It may also help to use the [IMGT/HLA sequence alignment tool](http://www.ebi.ac.uk/ipd/imgt/hla/align.html) for more information on common exon patterns.

## Sequence Format for Manual Annotation

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


## EMBL Metadata Input format
Sample ID: Specified by the submitting laboratory, you may use a value that is informative to you, such as a sample or experiment number.
Gene: Specify the gene name.  For example, 'HLA-A'
Class I or II: This tool requires a full-length Class I or II HLA sequence.
Allele Local Name: Specify a identifier for this sequence. You may use the name of the next-closest HLA allele.

Saddlebags will submit to EMBL Test environment by default, you must specify that the software target the Live / Production environment.

EMBL sequence submissions must be associated with a Study/Project. You may specify the accession number of an existing EMBL study (Get this accession number from [EMBL Webin](https://www.ebi.ac.uk/ena/submit/sra/#home) ), or Saddlebags can create a new project to your specifications. This sequence metadata is collected by Saddlebags and included in the EMBL submission.


## For more information on EMBL's ENA format:  
http://www.ebi.ac.uk/ena/submit/sequence-submission  
http://www.ebi.ac.uk/ena/submit/entry-upload-templates  
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt  
ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html  
http://www.ebi.ac.uk/ena/software/flat-file-validator  

## Contributing:  
More info on contributing to this project on the [MUMC Transplantation Immunology Wiki](https://github.com/transplantation-immunology-maastricht/General-information/wiki/Contributing-to-Github)


