# This file is part of saddle-bags.
#
# saddle-bags is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# saddle-bags is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Publsic License
# along with saddle-bags. If not, see <http://www.gnu.org/licenses/>.

import sys
from sys import exc_info

from datetime import datetime

try:
    from sys import _MEIPASS
except Exception:
    print ('No MEIPASS Directory. This is not running from a compiled EXE file.')

from os import makedirs
from os.path import join, expanduser, isfile, abspath, isdir, split

from tkinter import messagebox

from xml.etree import ElementTree as ET
from xml.dom import minidom as MD

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from io import StringIO
from json import loads

from saddlebags.HlaGene import HlaGene, GeneLocus

def getClosestAllele(alleleCallWithGFE):
    logEvent ('I am determining the closest known allele here. This is the ACT/GFE that was passed in:\n\n')
    #print (str(alleleCallWithGFE))
    
    # TODO: Obviously I need to be smart about this, not return a hard-coded value.
    return ('HLA-B*07:96:01')

def getAlleleDescription(alleleCallWithGFE):
    # TODO: When I fetch a GFE, it will include the next closest allele. 
    # Should I write an allele description for it?  I can generate a description like:
    # "The closest allele found is A*01:01 but with a polymorphism blah blah."
    # "IMGT/HLA requires a description of the next closest allele, Should I use this one?"
    
    logEvent ('I am generating an Allele Description here. This is the ACT/GFE that was passed in:\n\n')
    #print (str(alleleCallWithGFE))
        
    return ('The next closest allele is A*01, blah blah, with a polymorphism\n'
    + 'at an important locus.')
    
    

def assignIcon(tkRootWindow):
    logEvent ('Assigning Icon for the GUI.')

    # Find window location inside executable
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        #base_path = _MEIPASS
        iconFileLocation = resourcePath('images\\horse_image_icon.ico')
        #print ('I am assigning this icon:' + imageFileLocation)
        tkRootWindow.wm_iconbitmap(iconFileLocation)
    except Exception:
        #base_path = os.path.abspath(".")
        logEvent('Could not assign icon based on path inside executable.')
        
        logEvent (exc_info())

    # Linux
    # I have given up on setting an icon in linux. I can't seem to load up any file format.
  
 
def resourcePath(relativePath):
    # Where will I find my resources? This should work in, or outside, a compiled EXE
    if hasattr(sys, '_MEIPASS'):
        return join(_MEIPASS, relativePath)
    return join(abspath('.'), relativePath)

# This is a short wrapper method to use biopython's translation method. 
# Most of this code is just checking for things that went wrong

def translateSequence(inputSequence):

    proteinSequence = ''
    
    try:
        # Do nothing if the input sequence is blank.
        if( len(inputSequence) > 0 ):
            
            coding_dna = Seq(inputSequence, generic_dna)        
            proteinSequence = str(coding_dna.translate())   
            logEvent ('Exon Sequence before translation:' + coding_dna)     
            logEvent ('Translated Protein:' + proteinSequence)
            
            # Perform Sanity Checks.
            # Stop codon *should* be at the end of the protein.  
            # Here we seek out the first instance of a stop codon, 
            # and remove the peptides afterwards.
            # because that's what happens in real life.
            stopCodonLocation = proteinSequence.find('*')
            
            # If no stop codon was found
            if (stopCodonLocation == -1):
                assignConfigurationValue('is_pseudo_gene','1')
                logEvent ('No Stop Codon found. This is a "pseudo-gene".')
                # If multiple of three (correct codon length)
                if(len(coding_dna) % 3 == 0):
                    messagebox.showinfo('No Stop Codon Found', 
                        'The translated protein does not contain a stop codon.\n' + 
                        'This is indicated by a /pseudo flag in the sequence submission.'
                         )
                    
                # Wrong Codon Length
                else:
                    messagebox.showinfo('No Stop Codon Found', 
                        'The translated protein does not contain a stop codon.\n' + 
                        'The coding nucleotide sequence length (' + str(len(coding_dna))  + ') is not a multiple of 3.\n' + 
                        'This is indicated by a /pseudo flag in the sequence submission.')

            # If Stop Codon is in the end of the protein (This is expected and correct)
            elif (stopCodonLocation == len(proteinSequence) - 1):                
                assignConfigurationValue('is_pseudo_gene','0')
                
                # If multiple of three (correct codon length)
                if(len(coding_dna) % 3 == 0):
                    # Everything is fine in this case.  Trim off the stop codon
                    logEvent ('The stop codon is in the correct position. This is not a "pseudo-gene".')
                    proteinSequence = proteinSequence[0:stopCodonLocation]
                    pass 
                # Wrong Codon Length
                else:
                    logEvent ('The stop codon is in the correct position, but there are extra nucleotides. This is not a "pseudo-gene".')
                    messagebox.showinfo('Extra Nucleotides After the Stop Codon', 
                        'The stop codon is at the correct position in the protein, but ' + 
                        'The coding nucleotide sequence length (' + str(len(coding_dna))  + ') is not a multiple of 3.\n\n' +
                        'Please double check your sequence.')
                    proteinSequence = proteinSequence[0:stopCodonLocation]
                                        
            # Else Stop Codon is premature (before the end of the protein) 
            else:
                logEvent ('A premature stop codon was found. This is a "pseudo-gene".')
                assignConfigurationValue('is_pseudo_gene','1')
                
                # If multiple of three (correct codon length)
                if(len(coding_dna) % 3 == 0):
                    messagebox.showinfo('Premature Stop Codon Detected',
                        'Premature stop codon found:\nProtein Position (' + 
                        str(stopCodonLocation + 1) + '/' +
                        str(len(proteinSequence)) + ')\n\n' + 
                        'This is indicated by a /pseudo flag in the sequence submission.\n' +
                        'Double check your protein sequence,\n' + 
                        'this might indicate a missense mutation.\n\n' + 
                        'Translated Protein:\n' + proteinSequence + 
                        '\n\nProtein in EMBL Submission:\n' + proteinSequence[0:stopCodonLocation] + 
                        '\n'
                        )
                    proteinSequence = proteinSequence[0:stopCodonLocation]


                # Wrong Codon Length
                else:
                    messagebox.showinfo('Premature Stop Codon Detected',
                        'Premature stop codon found:\nProtein Position (' + 
                        str(stopCodonLocation + 1) + '/' +
                        str(len(proteinSequence)) + ')\n\n' + 
                        'This is indicated by a /pseudo flag in the sequence submission.\n' +
                        'Nucleotide count is not a multiple of 3,\n' +
                        'Double check your protein sequence,\n' + 
                        'this might indicate a missense mutation.\n\n' + 
                        'Translated Protein:\n' + proteinSequence + 
                        '\n\nProtein in EMBL Submission:\n' + proteinSequence[0:stopCodonLocation] + 
                        '\n'
                        )
                    proteinSequence = proteinSequence[0:stopCodonLocation]
        else:
            logEvent('Translating a nucleotide sequence of length 0.  That was easy.')
            pass

        return proteinSequence
    
    except Exception:
        logEvent ('Problem when translating protein:')
        logEvent (str(exc_info()))
        messagebox.showinfo('Protein Translation Error', 
            'I could not translate your protein:\n' +  str(exc_info()))
        
        raise
    
def collectAndValidateRoughSequence(guiSequenceInputObject):
    try:
        roughNucleotideSequence = collectRoughSequence(guiSequenceInputObject)
        annotatedSequence = None

        #Is this sequence in Fasta format?        
        try:
            #print('Checking if sequence is fasta format.')    
            fileHandleObject = StringIO(roughNucleotideSequence)
            fastaSeqList = list(SeqIO.parse(fileHandleObject, 'fasta'))
            #print ('The length of the fasta seq list is:' + str(len(fastaSeqList)))
            if(len(fastaSeqList) == 1):
                annotatedSequence = cleanSequence(str(fastaSeqList[0].seq))
                #print ('found exactly 1 fasta sequence:' + annotatedSequence)
                logEvent ('The input sequence is in .fasta format.')
            else:
                logEvent('This sequence is not in .fasta format.') 
        except Exception:
            logEvent('This sequence is not in .fasta format: ' + str(exc_info()))            
            
        #Is this sequence in Fastq format?        
        try:
            #print('Checking if sequence is fastq format.')            
            fileHandleObject = StringIO(roughNucleotideSequence)
            fastqSeqList = list(SeqIO.parse(fileHandleObject, 'fastq'))
            #print ('The length of the fasta seq list is:' + str(len(fastaSeqList)))
            if(len(fastqSeqList) == 1):
                annotatedSequence = cleanSequence(str(fastqSeqList[0].seq))
                #print ('found exactly 1 fasta sequence:' + annotatedSequence)
                logEvent ('The input sequence is in .fastq format.')
            else:
                logEvent('This sequence is not in .fastq format.') 
        except Exception:
            logEvent('This sequence is not in .fastq format: ' + str(exc_info()))            

    
            
        # TODO: If this file is xml what should we do?  Just give up i suppose.
        # We want to accept HML.  But there are too many xml formats.
        # Yeah I dunno about HML, we will not implement that right now.
        
        #else Is XML?
        #    Warn user that XML isn't supported
        #    Return the entire xml text as the sequence
        
        
        # If we haven't found an annotated sequence yet, this is not fasta or fastq.
        if(annotatedSequence is None):
            annotatedSequence = cleanSequence(roughNucleotideSequence)
        #    Rough text = all of the gui text.
            
    
        #Are we using any nonstandard / ambiguous nucleotides?
        for nucleotideCharacter in annotatedSequence:
            if(nucleotideCharacter not in ('A','G','C','T','a','g','c','t')):
                messagebox.showerror('Nonstandard Nucleotide' 
                    , 'I found a non-standard\n'
                    + 'character in your nucleotide\n'
                    + 'sequence: ' 
                    + str(nucleotideCharacter) + '\n'
                    + 'You should use standard nucleotide\n'
                    + 'characters in your submission.\n'
                    + 'I will attempt to continue.')
                break
    

    
        # TODO: Fix this last-ditch effort.
        
        #annotatedSequence = roughNucleotideSequence
        return annotatedSequence
            
    except Exception:
        #except Exception, e:
        messagebox.showerror('Error Reading Input Sequence.'
            , str(exc_info()))
        raise
                

def collectRoughSequence(guiSequenceInputObject):
    # This method gets the text from a gui object and returns it.
    # There is no validation performed here.
    try:
        roughNucleotideSequence = guiSequenceInputObject.get('1.0', 'end')
        return roughNucleotideSequence
            
    except Exception:
        #except Exception, e:
        messagebox.showerror('Error Reading Input Sequence.'
            , str(exc_info()))
        raise
                


# 
def isSequenceAlreadyAnnotated(inputSequenceText):
    # The easy case.
    if ('a' in inputSequenceText and
        'g' in inputSequenceText and
        'c' in inputSequenceText and
        't' in inputSequenceText and
        'A' in inputSequenceText and
        'G' in inputSequenceText and
        'C' in inputSequenceText and
        'T' in inputSequenceText
        ):
        return True
    
    # TODO: This isn't perfect. It must have all 8 nucleotides to return true.
    # Circle back on this one later.
    return False

def parseExons(roughFeatureSequence, alleleCallWithGFEJson):
    
    # TODO: Parse the JSON in the alleleCallWithGFEJson.  
    # There should be some information about the exons in here.
    # ex = uppercase
    # utr = lowercase
    
    try:

        fivePrimeSequence = ''
        threePrimeSequence = ''
        exonDictionary = {}
        intronDictionary = {}
        
        parsedJson = loads(alleleCallWithGFEJson)

        if (len(parsedJson.keys()) > 1):
                        
            if 'status' in parsedJson.keys():
                queryStatus = parsedJson['status']
                if(str(queryStatus) == '200'):
                    logEvent ('Response from Allele call was status 200. That is quite normal.')
                elif(str(queryStatus) == '500'):
                    logEvent ('500 status found. Unknown server error.')
                    logEvent ('JSON Results:' + str(parsedJson))
                    #messagebox.showerror('Error Annotating Sequence.',
                    #    'The ACT service returned a 500 status, there was an unknown server problem. I have no annotation results.')
                    raise Exception ('The ACT service returned a 500 status, there was an unknown server problem. I have no annotation results, you can try to annotate the sequence manually.')
                    
            else:
                logEvent ('JSON results have no "status" element. This is not a problem.')
               
            # Is sequence novel?
            if 'typing_status' in parsedJson.keys():
                logEvent ('typing_status element was found.')
            else:
                logEvent ('JSON Results:' + str(parsedJson))
                raise Exception ('No typing_status element was found in Json results. Cannot continue.')
            
            typingStatusDictionary = parsedJson['typing_status']
            
            # Loop through the recognized Features
            if 'features' in parsedJson.keys():
                # We found features.
                featureList = parsedJson['features']
                logEvent ('I found this many Known Features:' + str(len(featureList)))
                
                for featureDictionary in featureList:
                    
                    term=str(featureDictionary['term'])
                    rank=str(featureDictionary['rank'])
                    sequence= str(featureDictionary['sequence'])
                    
                    #print ('Known Feature' 
                    #    + ':' + term
                    #    + ':' + rank
                    #    + ':' + sequence)
                    
                    if(term == 'five_prime_UTR'):
                        fivePrimeSequence = sequence.lower()
                    elif(term == 'three_prime_UTR'):
                        threePrimeSequence = sequence.lower()
                    elif(term == 'exon'):    
                        exonDictionary[rank] = sequence.upper()
                    elif(term == 'intron'):    
                        intronDictionary[rank] = sequence.lower()
                    else:
                        raise Exception('Unknown Feature Term, expected exon or intron:' + term)
                
                #print ('fivePrimeSequence:\n' + str(fivePrimeSequence))
                #print ('threePrimeSequence:\n' + str(threePrimeSequence))
                #print ('exonCount:\n' + str(len(exonDictionary.keys())))
                #print ('intronCount:\n' + str(len(intronDictionary.keys())))
            else:
                raise Exception ('Unable to identify any HLA exon features, unable to annotate sequence.')
                # no features found
                #return roughFeatureSequence
                
            # Loop through the Novel Features
            if 'novel_features' in typingStatusDictionary.keys():
                logEvent ('Novel Features were found.')
                # We found features.
                featureList = typingStatusDictionary['novel_features']
                #print 'This many Novel Features:' + str(len(featureList))
                
                for featureDictionary in featureList:

                    term=str(featureDictionary['term'])
                    rank=str(featureDictionary['rank'])
                    sequence= str(featureDictionary['sequence'])
                    
                    #print ('Novel Feature' 
                    #    + ':' + term
                    #    + ':' + rank
                    #    + ':' + sequence)
                    
                    if(term == 'five_prime_UTR'):
                        fivePrimeSequence = sequence.lower()
                        logEvent ('Novel 5Prime Sequence:' + fivePrimeSequence)
                    elif(term == 'three_prime_UTR'):
                        threePrimeSequence = sequence.lower()
                        logEvent ('Novel 3Prime Sequence:' + threePrimeSequence)
                    elif(term == 'exon'):    
                        exonDictionary[rank] = sequence.upper()
                        logEvent ('Novel Exon ' + str(rank) + ' Sequence:' + exonDictionary[rank])
                    elif(term == 'intron'):    
                        intronDictionary[rank] = sequence.lower()
                        logEvent ('Novel Intron ' + str(rank) + ' Sequence:' + intronDictionary[rank])
                    else:
                        raise Exception('Unknown Feature Term, expected exon or intron:' + term)
                
            else:
                logEvent ('No novel features were found. Presumably this is a known HLA allele. No problem.')
                #raise Exception ('Unable to identify any HLA exon features, unable to annotate sequence.')
            
            if (len(fivePrimeSequence) < 1):
                logEvent ('I cannot find a five prime UTR.')
                logEvent ('Rough Sequence:\n' + cleanSequence(roughFeatureSequence).upper())
                #print ('Annotated Sequence:\n' + cleanSequence(annotatedSequence).upper())
                raise Exception('GFE service did not find a 5\' UTR sequence. You will need to annotate the genomic features manually.')
            # What if the reported 5' UTR is less than what is returned by GFE?
            elif cleanSequence(fivePrimeSequence).upper() in cleanSequence(roughFeatureSequence).upper():
                # This means that we provided a longer sequence than what is available in the GFE service.
                #messagebox.showinfo('Short 5\' Sequence', 
                #    'The 5\' sequence from the GFE service is shorter than your provided sequence.\n'
                #    + 'I will use your sequence instead.'
                #     )
                beginIndex = cleanSequence(roughFeatureSequence).upper().find(cleanSequence(fivePrimeSequence).upper())
                endIndex = beginIndex + len(fivePrimeSequence)
                logEvent ('GFE sequence exists in rough sequence, at index: (' + str(beginIndex) + ':' + str(endIndex) + ')')
                logEvent ('previous fivePrime Sequence=\n' + fivePrimeSequence)
                fivePrimeSequence = cleanSequence(roughFeatureSequence)[0:endIndex].lower()
                logEvent ('new fivePrime Sequence=\n' + fivePrimeSequence)
            #print ('FOUND THIS ANNOTATED SEQUENCE:\n' + str(annotatedSequence))

            
            logEvent ('Annotating 5\' UTR:' + str(fivePrimeSequence))
            annotatedSequence = fivePrimeSequence + '\n'

            # arbitrarily choose 50.
            # TODO: this loop range is arbitrary. 
            # Maybe indexString should just loop through the exon and intron dictionarys
            for i in range(1,50):
                indexString = str(i)
                if indexString in exonDictionary.keys():
                    logEvent ('Annotating exon#' + indexString + ':' + str(exonDictionary[indexString]))
                    annotatedSequence += (str(exonDictionary[indexString]) + '\n')
                    
                if indexString in intronDictionary.keys():
                    logEvent ('Annotating intron#' + indexString + ':' + str(intronDictionary[indexString]))
                    annotatedSequence += (str(intronDictionary[indexString]) + '\n')
                
            logEvent ('Annotating 3\' UTR:' + str(threePrimeSequence))
            
            if (len(threePrimeSequence) < 1):
                #print ('Rough Sequence:\n' + cleanSequence(roughFeatureSequence).upper())
                #print ('Annotated Sequence:\n' + cleanSequence(annotatedSequence).upper())

                #raise Exception('GFE service did not find a 3\' UTR sequence. You will need to annotate the genomic features manually.')
                logEvent('There is no three prime sequence.')
                
                # if sequence so far is in the rough sequence
                if cleanSequence(annotatedSequence).upper() in cleanSequence(roughFeatureSequence).upper():
                    # use rest of the sequence as the UTR.
                    beginIndex = cleanSequence(roughFeatureSequence).upper().find(cleanSequence(annotatedSequence).upper()) + len(cleanSequence(annotatedSequence))
                    #endIndex = len(roughFeatureSequence)
                    threePrimeSequence = cleanSequence(roughFeatureSequence)[beginIndex:].lower()
                    logEvent('Using the rest of the sequence as the 3\' UTR:\n' + threePrimeSequence)
                    annotatedSequence += threePrimeSequence
                
                
            else:
            
                annotatedSequence += threePrimeSequence
            
            # TODO: I need to have better checks here.  What is missing?
            # Do the annotated sequence and rough sequence match?
            if(cleanSequence(annotatedSequence).upper() == cleanSequence(roughFeatureSequence).upper()):
            
                return annotatedSequence
            
            else:
                logEvent ('Rough Sequence:\n' + cleanSequence(roughFeatureSequence).upper())
                logEvent ('Annotated Sequence:\n' + cleanSequence(annotatedSequence).upper())

                raise Exception('Annotated sequence and rough sequence do not match.')
        
    
            
        else:
            raise Exception ('No keys found in the JSON Dictionary, unable to annotate sequence.')
            # no keys in JSON dictionary.
            #return roughFeatureSequence

    
        raise Exception ('Reached end of parsing without returning a value.')
    
    
    except Exception:
        logEvent ('Exception when parsing exons:')
        logEvent (str((exc_info())))
        messagebox.showinfo('Exon Parsing Error', 
            'I had trouble annotating your sequence:\n' 
            +  str(str(exc_info()) 
            + '. You will have to annotate manually.')) 
        return roughFeatureSequence
        
        #raise
        

def cleanSequence(inputSequenceText):
    # Trim out any spaces, tabs, newlines. 
    cleanedSequence = inputSequenceText.replace(' ','').replace('\n','').replace('\t','').replace('\r','')
    return cleanedSequence

# The input file should be a string of nucleotides, with capital letters to identify exons and introns.
# Annotations are expected and read in this format:
# fiveprimeutrEXONONEintrononeEXONTWOintrontwoEXONTHREEthreeprimeutr
# agctagctagctAGCTAGCtagctagctAGCTAGCtagctagctAGCTAGCTAgctagctagctag
# All spaces, line feeds, and tabs are removed and ignored.  
def identifyGenomicFeatures(inputSequenceText):

    # TODO: I should accept a Fasta Input. 
    # Disregard the header line completely. Is there still sequence?
    resultGeneLoci = HlaGene()
    
    cleanedGene = cleanSequence(inputSequenceText)
    
    # Capitalize, so I can store a copy of the full unannotated sequence.
    unannotatedGene = cleanedGene.upper()
    resultGeneLoci.fullSequence = unannotatedGene
    logEvent('Total Sequence Length = ' + str(len(unannotatedGene)))

    # Loop through the cleaned and annotated input sequence, 
    # capitals and lowercase letters to determine exon start and end
    if(len(cleanedGene) > 0):
        
        # Is the first feature an exon or an intron?
        # If we begin in an Exon
        if( cleanedGene[0] in ('A','G','C','T')):                
            insideAnExon = True
        # If we begin in an Intron/UTR
        elif( cleanedGene[0] in ('a','g','c','t')):  
            insideAnExon = False
        else:
            # Nonstandard nucleotide? I should start panicking.
            #raise Exception('Nonstandard Nucleotide, not sure how to handle it')
            logEvent('Nonstandard Nucleotide at the beginning of the sequence, not sure how to handle it')
            insideAnExon = False
        
        
        locusBeginPosition = 0
        for x in range(0, len(cleanedGene)):
            currentChar = cleanedGene[x]
            
            # Is this a standard nucleotide character?
            if(currentChar.upper() in ('A','G','C','T')):

                if(currentChar.isupper()):
                    if(insideAnExon):
                        #We're STILL in an exon.  In this case, I should just do nothing and continue.  
                        pass
                    else:
                        #In this case, we're just starting an EXON.
                        #Store the last Intron in the list.
                        currentIntron = GeneLocus()
                        currentIntron.sequence = cleanedGene[locusBeginPosition:x].upper()
                        currentIntron.exon = False
                        resultGeneLoci.loci.append(currentIntron)                    
                        insideAnExon=True
                        locusBeginPosition = x
                        pass
                        
                else:
                    if not (insideAnExon):
                        #We're STILL in an intron.  Continue.
                        pass
                    else:
                        #Starting a new Intron.
                        # Store an Exon in the list.
                        currentExon = GeneLocus()
                        currentExon.sequence = cleanedGene[locusBeginPosition:x].upper()
                        currentExon.exon = True
                        resultGeneLoci.loci.append(currentExon)     
                        insideAnExon = False
                        locusBeginPosition=x
                        pass
            else:
                logEvent('Nonstandard nucleotide detected at position ' + str(x) + ' : ' + currentChar 
                    + '.  If this is a wildcard character, you might be ok.')

        # We've reached the end of the loop and we still need to store the last feature.
        # Should be a 3' UTR, but I can't be sure, people like to put in weird sequences.
        currentIntron = GeneLocus()
        currentIntron.sequence = cleanedGene[locusBeginPosition:len(cleanedGene)].upper()
        currentIntron.exon = insideAnExon
        resultGeneLoci.loci.append(currentIntron)    

        # Annotate the loci (name them) and print the results of the read file.
        resultGeneLoci.annotateLoci()
        resultGeneLoci.printGeneSummary()
    
    # If the sequence is empty
    else:
        logEvent('Empty sequence, I don\'t have anything to do.')
        
    return resultGeneLoci    
    #self.sequenceAnnotation = resultGeneLoci

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        logEvent('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput

# I'm storing global variables in a dictionary for now. 
def initializeGlobalVariables():    
    global globalVariables 
    
    if not ("globalVariables" in globals()):
        globalVariables={}
        
def assignConfigurationValue(configurationKey, configurationValue):
    # Overwrite config value without question.
    initializeGlobalVariables()
    globalVariables[configurationKey] = configurationValue
    
def assignIfNotExists(configurationKey, configurationValue):
    # Use this assigner if we want to declare important, new configuration values.
    # Using this method, we will not overwrite custom values
    # But we will provide critical new config values.
    initializeGlobalVariables()
    if configurationKey not in globalVariables.keys():
        assignConfigurationValue(configurationKey, configurationValue)
    
    
def getConfigurationValue(configurationKey):
    if configurationKey in globalVariables.keys():
        return globalVariables[configurationKey]
    else:
        logEvent ('Configuration Key Not Found:' + configurationKey)
        #raise KeyError('Key Not Found:' + configurationKey)
        return None


def logEvent(eventText):
    fullLogMessage = str(datetime.now()) + ' : ' + str(eventText)
    print(fullLogMessage)
        
    if 'logging' in globalVariables.keys():
        
        if (getConfigurationValue('logging') == '1'
            or getConfigurationValue('logging').upper() == 'TRUE'
            or getConfigurationValue('logging').upper() == 'T'
            ):
            
            logFileName = join(join(
                expanduser("~"),'saddlebags_temp'),'Saddlebags.Log.txt')
            logFile = open(logFileName, 'a')
        
            logFile.write(fullLogMessage + '\n')
            logFile.close()   

def assignConfigName():
    # Join together the working directory, a subfolder called "saddlebags_temp", and the config name.
    assignConfigurationValue('config_file_location',join(join(
        expanduser("~"),'saddlebags_temp'),'Saddlebags.Config.xml'))
    
def writeConfigurationFile():
    assignConfigName()
    logEvent ('Writing a config file to:\n' + globalVariables['config_file_location'])
    
    root = ET.Element("config")
    
    # TODO: Potential problem: I want to store primer descriptions.
    # This is a list of primers. Can i save or load a list of strings in here?
    for key in globalVariables.keys():
        # Some config values I don't want to store.

        if(key not in [
            'embl_password'
            ,'imgt_password'
            , 'sequence'
            , 'source_hla'
            ]):
            ET.SubElement(root, key).text = globalVariables[key]

    xmlText = ET.tostring(root, encoding='utf8', method='xml')
    prettyXmlText = MD.parseString(xmlText).toprettyxml()
    
    xmlOutput = createOutputFile(globalVariables['config_file_location'])
    xmlOutput.write(prettyXmlText)
    xmlOutput.close()

def loadConfigurationFile():
    assignConfigName()
    
    if not isfile(globalVariables['config_file_location']):
        logEvent ('The config file does not exist yet. I will not load it:\n' + globalVariables['config_file_location'])
        
    else:
        logEvent ('The config file already exists, I will load it:\n' + globalVariables['config_file_location'])
        
        tree = ET.parse(globalVariables['config_file_location'])
        root = tree.getroot()
        
        for child in root:
            assignConfigurationValue(child.tag, child.text)
            
    # Here is where I assign the common/critical configuration values
    # test_submission indicates if we should use the "test" values.
    # I think I'll use this value for both EMBL and IMGT submissions, if it applies.
    assignIfNotExists('test_submission', '1')
    
    # Logging is turned off by default. Users can change this to 1 to turn on a logfile.
    assignIfNotExists('logging','0')
    
    # I'm storing FTP without the ftp:// identifier, because it is not necessary.
    # The test and prod ftp sites have the same address. This is intentional, embl doesn't have a test ftp
    assignIfNotExists('embl_ftp_upload_site_test', 'webin.ebi.ac.uk')
    assignIfNotExists('embl_ftp_upload_site_prod', 'webin.ebi.ac.uk')
    assignIfNotExists('embl_rest_address_test', 'https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/')
    assignIfNotExists('embl_rest_address_prod', 'https://www.ebi.ac.uk/ena/submit/drop-box/submit/')
    assignIfNotExists('nmdp_act_rest_address', 'http://act.b12x.org/act' )
    
            
    
       