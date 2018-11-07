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
from json import loads, dumps # what dumb method names.

from saddlebags.AlleleSubmission import HlaGene, GeneFeature, SubmissionBatch, AlleleSubmission

#from gfe_client import TypeSeqAPI, typeseq_get
# TODO: I don't like to use global imports like this, can i change it to only import a function or something, Yes.
import gfe_client


# TODO I dont think I'm using this locus, it is not necessary anymore.
def fetchSequenceAlleleCallWithGFE(rawSequence, locus):

    cleanedSequence = cleanSequence(rawSequence.upper())

    # TODO: get the act address into a configuration file. I want to be able to use a local service. Which i already have?
    config = gfe_client.Configuration()
    config.host = 'http://act.b12x.org'
    api = gfe_client.ApiClient(configuration=config)

    ann_api = gfe_client.TypeSeqApi(api_client=api)

    responseText = ann_api.typeseq_get(sequence=cleanedSequence, imgthla_version="3.31.0", locus='HLA-A')

    annotation = responseText.to_dict()
    jsonResponse = dumps(annotation, indent=4)

    return jsonResponse



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
    # Circle back on this one later, should I store a variable somewhere if the sequence has been annotated?
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

        # The json should be a String, but it is returned from the NMDP ACT API as a "Typing" object. I convert it to a String to make everyone happy.
        jsonString = str(alleleCallWithGFEJson)
        print('THIS IS THE JSON STRING\n:' + jsonString)
        parsedJson = loads(jsonString)

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
            #if 'typing_status' in parsedJson.keys():
            #    logEvent ('typing_status element was found.')
            #else:
            #    logEvent ('JSON Results:' + str(parsedJson))
            #    raise Exception ('No typing_status element was found in Json results. Cannot continue.')
            # TODO We no longer have a typing_status element in the returned JSON. I will find the typing status here...

            # Is sequence novel?
            if 'status' in parsedJson.keys():
                logEvent ('status element was found.')
            else:
                logEvent ('JSON Results:' + str(parsedJson))
                raise Exception ('No status element was found in Json results. Cannot continue.')

            typingStatusNode = parsedJson['status']

            print('I will try to loop through the keys of the typing status dictionary.')
            print('the typing status dictionary looks like this:')
            print(typingStatusNode)

            
            
            if(typingStatusNode == 'documented'):
                logEvent ('This is a known/documented allele.')
            elif(typingStatusNode == 'novel'):
                logEvent ('This is a novel allele.')                
            elif(typingStatusNode == 'novel_combination'):
                logEvent ('This is a novel combination of gene features.')
            else:
                print ('I do not understand the status of this allele. Expected either "documented" or "novel" or "novel_combination":')
                print(typingStatusNode)
                raise Exception('Unknown Typing status, expected documented or novel:' + typingStatusNode)

            # TODO: This format is different. Re-work this logic until it is better.

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

            # TODO: Once we know the genomic features, I can add this next-closest allele description.

            alleleDescription = ''
            if 'hla' in parsedJson.keys():
                closestAllele = parsedJson['hla']
                alleleDescription += (closestAllele)

                # if seqdiff is available, parse that list and specify the allele differences.
                # This has the most info but its missing from old versions of ACT.
                # else if novel_features are available, at least I can tell what feature the modification is in.
                # else, assume the sequence is known, not novel.
                if 'seqdiff' in parsedJson.keys():
                    seqDiffList = parsedJson['seqdiff']

                    # seqDiffList can be None, if there are no known sequence differences.
                    if(seqDiffList is not None):
                        for seqDiffDictionary in seqDiffList:

                            alleleDescription += ('\nFT                  Position '
                                + str(seqDiffDictionary['location']) + ' in '
                                + str(seqDiffDictionary['term']) + ' : '
                                + str(seqDiffDictionary['ref']) + '->' + str(seqDiffDictionary['inseq'])
                                )
                    else:
                        logEvent('No unknown features. This allele is documented. Cool.')



                # Loop Novel Features
                elif 'novel_features' in parsedJson.keys():
                    novelFeatureList = parsedJson['novel_features']

                    # In this case I only have info about the feature it is in. Darn.
                    # TODO Make this a bit more descriptive
                    for featureDictionary in novelFeatureList:
                        alleleDescription += ('\nNovel '
                            + str(featureDictionary['term']) + ' ' + str(featureDictionary['rank']))


                else:
                    # Add the Sequence Differences.
                    alleleDescription += ('No novel features identified. This sequence is not novel?')

            else:
                alleleDescption += 'Could not determine closest HLA allele, please provide a detailed description of the novel sequence.'

            assignConfigurationValue('closest_allele_written_description', alleleDescription)

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

            
            #logEvent ('Annotating 5\' UTR:' + str(fivePrimeSequence))
            annotatedSequence = fivePrimeSequence + '\n'

            # arbitrarily choose 50.
            # TODO: this loop range is arbitrary, i'm using a for loop so i can access the index. Someone could calculate a max index.
            # Maybe indexString should just loop through the exon and intron dictionarys
            for i in range(1,50):
                indexString = str(i)
                if indexString in exonDictionary.keys():
                    #logEvent ('Annotating exon#' + indexString + ':' + str(exonDictionary[indexString]))
                    annotatedSequence += (str(exonDictionary[indexString]) + '\n')
                    
                if indexString in intronDictionary.keys():
                    #logEvent ('Annotating intron#' + indexString + ':' + str(intronDictionary[indexString]))
                    annotatedSequence += (str(intronDictionary[indexString]) + '\n')
                
            #logEvent ('Annotating 3\' UTR:' + str(threePrimeSequence))
            
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
                        currentIntron = GeneFeature()
                        currentIntron.sequence = cleanedGene[locusBeginPosition:x].upper()
                        currentIntron.exon = False
                        resultGeneLoci.features.append(currentIntron)
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
                        currentExon = GeneFeature()
                        currentExon.sequence = cleanedGene[locusBeginPosition:x].upper()
                        currentExon.exon = True
                        resultGeneLoci.features.append(currentExon)
                        insideAnExon = False
                        locusBeginPosition=x
                        pass
            else:
                logEvent('Nonstandard nucleotide detected at position ' + str(x) + ' : ' + currentChar 
                    + '.  If this is a wildcard character, you might be ok.')

        # We've reached the end of the loop and we still need to store the last feature.
        # Should be a 3' UTR, but I can't be sure, people like to put in weird sequences.
        currentIntron = GeneFeature()
        currentIntron.sequence = cleanedGene[locusBeginPosition:len(cleanedGene)].upper()
        currentIntron.exon = insideAnExon
        resultGeneLoci.features.append(currentIntron)

        # Annotate the features (name them) and print the results of the read file.
        resultGeneLoci.annotateFeatures()
        #resultGeneLoci.printGeneSummary()
    
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

# Clear my configuration, fearlessly and without hesitation.
def clearGlobalVariables():
    global globalVariables
    globalVariables = {}

# I'm storing global variables in a dictionary for now. 
def initializeGlobalVariables():    
    global globalVariables 
    
    if not ("globalVariables" in globals()):
        globalVariables={}
        
def assignConfigurationValue(configurationKey, configurationValue):
    # Overwrite config value without question.
    initializeGlobalVariables()

    # Lists or Strings are handled well by the configuration serializer.
    if(type(configurationValue) is list or type(configurationValue) is str ):
        globalVariables[configurationKey] = serializeConfigValue(configurationValue)
    else:
        globalVariables[configurationKey] = configurationValue

    globalVariables[configurationKey] = configurationValue
    #logEvent ('Just stored configuration key ' + configurationKey + ' which is ' + str(configurationValue) + ' of type ' + str(type(configurationValue)))

def assignIfNotExists(configurationKey, configurationValue):
    # Use this assigner if we want to declare important, new configuration values.
    # Using this method, we will not overwrite custom values
    # But we will provide critical new config values.
    initializeGlobalVariables()
    if configurationKey not in globalVariables.keys():
        assignConfigurationValue(configurationKey, configurationValue)

def getConfigurationValue(configurationKey):
    if configurationKey in globalVariables.keys():

        configurationValue = globalVariables[configurationKey]

        #logEvent ('Retrieving configuration key ' + configurationKey + ' which is ' + str(configurationValue) + ' of type ' + str(type(configurationValue)))

        if (type(configurationValue) is str):
            return deserializeConfigValue(configurationValue)
        else:
            return configurationValue
    else:
        logEvent ('Configuration Key Not Found:' + configurationKey)
        #raise KeyError('Key Not Found:' + configurationKey)
        return None

def logEvent(eventText):
    initializeGlobalVariables()
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

# Encode semicolons within a config value.
# Store lists as a string separated by a semicolon.
def serializeConfigValue(configListObject):
    #print('serializing a config value')
    if (type(configListObject) == str):
        return configListObject.replace(';','@@@')
    elif(type(configListObject) == list):
        serializedString = ''
        # TODO: what if the list elements are not strings? Just don't do that. Lol.
        # TODO: shouldn't there be built in methods for serializing? probably.
        # Encode semicolons in the list elements.
        for configString in configListObject:
            serializedString = serializedString + str(configString).replace(';','@@@') + ';'
        serializedString = serializedString[:-1]
        return serializedString
    else:
        raise(TypeException('Unknown configuration type, can not serialize:' + str(type(configListObject))))


# Split strings containing semicolon into a list.
# Decode @@@ back into a semicolon.
def deserializeConfigValue(serializedConfigString):
    #print('deserializing a config value')
    if (';' in serializedConfigString):
        configList = serializedConfigString.split(';')
        # Maybe a for loop is not necessary but I want to make sure i can change the string in the list.
        for i in range (0,len(configList)):
            configList[i] = configList[i].replace('@@@',';')
        return configList
    else:
        # TODO What if it's not actually a string? Could be an arbitrary object...
        return serializedConfigString.replace('@@@',';')

def writeConfigurationFile():
    assignConfigName()
    logEvent ('Writing a config file to:\n' + getConfigurationValue('config_file_location'))
    
    root = ET.Element("config")

    # Root node stores "normal" configuration keys, stuff related to software
    # and not necessarily HLA or sequence submission.
    configElement = ET.SubElement(root, 'config_file_location').text = getConfigurationValue('config_file_location')


    # The purpose of this loop is to make sure I don't miss any keys but I can just do that myself.
    for key in globalVariables.keys():

        # "normal" configuration keys, stuff related to software and not necessarily HLA or sequence submission.
        # Some config values I don't want to store. I can add more to this list if i want.
        # Don't store passwords.
        # I handle the submission batch manually.
        if(key not in [
            'embl_password'
            , 'imgt_password'
            , 'submission_batch'
            ]):

            # getConfigurationValue will handle serializing and encoding semicolons.
            ET.SubElement(root, key).text = getConfigurationValue(key)

    # Add keys for "each" batch of submissions.
    # TODO: i may add functionality for multiple batches later. Put this in a loop.
    submissionBatch = getConfigurationValue('submission_batch')
    
    # If the config is not already initiated, this can be None. Make a new one.
    if (submissionBatch is None):
        submissionBatch = SubmissionBatch()
    
    
    # Create a node object, most of this stuff can be parameters on the node.
    submissionBatchElement = ET.SubElement(root, 'submission_batch')
    submissionBatchElement.set('imgtsubmitterid', serializeConfigValue(submissionBatch.imgtSubmitterId))
    submissionBatchElement.set('imgtsubmittername', serializeConfigValue(submissionBatch.imgtSubmitterName))
    submissionBatchElement.set('imgtaltcontact', serializeConfigValue(submissionBatch.imgtAltContact))
    submissionBatchElement.set('imgtsubmitteremail', serializeConfigValue(submissionBatch.imgtSubmitterEmail))
    submissionBatchElement.set('laboforigin', serializeConfigValue(submissionBatch.labOfOrigin))
    submissionBatchElement.set('labcontact', serializeConfigValue(submissionBatch.labContact))

    # Keys for each submission.
    for hlaAllele in submissionBatch.submissionBatch:
        print ('I found an HLA allele.')
        submissionElement = ET.SubElement(submissionBatchElement, 'submission')

        # Most of this stuff is attrbutes. Store the Sequence as the text of this element.
        submissionElement.text = hlaAllele.submittedGene.fullSequence
        submissionElement.set('genelocus',                       serializeConfigValue(hlaAllele.submittedGene.geneLocus))
        submissionElement.set('localallelename'                , serializeConfigValue(hlaAllele.localAlleleName))
        submissionElement.set('closestallelewrittendescription', serializeConfigValue(hlaAllele.closestAlleleWrittenDescription))
        submissionElement.set('imgtsubmissionidentifier'       , serializeConfigValue(hlaAllele.imgtSubmissionIdentifier))
        submissionElement.set('imgtsubmissionversion', serializeConfigValue(hlaAllele.imgtSubmissionVersion))
        submissionElement.set('cellid', serializeConfigValue(hlaAllele.cellId))
        submissionElement.set('ethnicorigin', serializeConfigValue(hlaAllele.ethnicOrigin))
        submissionElement.set('sex', serializeConfigValue(hlaAllele.sex))
        submissionElement.set('consanguineous', serializeConfigValue(hlaAllele.consanguineous))
        submissionElement.set('homozygous', serializeConfigValue(hlaAllele.homozygous))
        submissionElement.set('typedalleles', serializeConfigValue(hlaAllele.typedAlleles))
        submissionElement.set('materialavailability', serializeConfigValue(hlaAllele.materialAvailability))
        submissionElement.set('cellbank', serializeConfigValue(hlaAllele.cellBank))
        submissionElement.set('primarysequencingmethodology', serializeConfigValue(hlaAllele.primarySequencingMethodology))
        submissionElement.set('secondarysequencingmethodology', serializeConfigValue(hlaAllele.secondarySequencingMethodology))
        submissionElement.set('primertype', serializeConfigValue(hlaAllele.primerType))
        submissionElement.set('primers', serializeConfigValue(hlaAllele.primers))
        submissionElement.set('sequencedinisolation', serializeConfigValue(hlaAllele.sequencedInIsolation))
        submissionElement.set('numofreactions', serializeConfigValue(hlaAllele.noOfReactions))
        submissionElement.set('methodcomments', serializeConfigValue(hlaAllele.methodComments))
        submissionElement.set('citations', serializeConfigValue(hlaAllele.citations))

    xmlText = ET.tostring(root, encoding='utf8', method='xml')
    prettyXmlText = MD.parseString(xmlText).toprettyxml()
    
    xmlOutput = createOutputFile(globalVariables['config_file_location'])
    xmlOutput.write(prettyXmlText)
    xmlOutput.close()

def loadConfigurationFile():
    # TODO: should I clear my configuration first? I have a method to purge my globals. I don't know right now, but probably not.
    assignConfigName()

    if not isfile(globalVariables['config_file_location']):
        logEvent ('The config file does not exist yet. I will not load it:\n' + globalVariables['config_file_location'])
    else:
        logEvent ('The config file already exists, I will load it:\n' + globalVariables['config_file_location'])

        tree = ET.parse(globalVariables['config_file_location'])
        root = tree.getroot()
        
        for child in root:
            #print ('The child tag is:' + child.tag)

            # If the child node is a submission batch
            # TODO: Think about how I can handle multiple submission batches.
            if(child.tag == 'submission_batch'):
                submissionBatch = SubmissionBatch()

                # Assign some information about this batch of submissions.
                submissionBatch.imgtSubmitterId = deserializeConfigValue(child.attrib['imgtsubmitterid'])
                submissionBatch.imgtSubmitterName = deserializeConfigValue(child.attrib['imgtsubmittername'])
                submissionBatch.imgtAltContact = deserializeConfigValue(child.attrib['imgtaltcontact'])
                submissionBatch.imgtSubmitterEmail = deserializeConfigValue(child.attrib['imgtsubmitteremail'])
                submissionBatch.labOfOrigin = deserializeConfigValue(child.attrib['laboforigin'])
                submissionBatch.labContact = deserializeConfigValue(child.attrib['labcontact'])

                # Loop the children, they are submission objects. Load up their information.
                for submissionChild in child:
                    print ('The submission child tag is:' + submissionChild.tag)
                    #print ('This submission has the text:' + submissionChild.text)
                    # Add a few submissions to this batch.
                    # Submission # 1
                    submission = AlleleSubmission()
                    submission.submittedGene.fullSequence = submissionChild.text
                    submission.submittedGene.geneLocus = deserializeConfigValue(submissionChild.attrib['genelocus'])
                    submission.localAlleleName = deserializeConfigValue(submissionChild.attrib['localallelename'])
                    submission.closestAlleleWrittenDescription = deserializeConfigValue(submissionChild.attrib['closestallelewrittendescription'])
                    submission.imgtSubmissionIdentifier = deserializeConfigValue(submissionChild.attrib['imgtsubmissionidentifier'])
                    submission.imgtSubmissionVersion = deserializeConfigValue(submissionChild.attrib['imgtsubmissionversion'])
                    submission.cellId = deserializeConfigValue(submissionChild.attrib['cellid'])
                    submission.ethnicOrigin = deserializeConfigValue(submissionChild.attrib['ethnicorigin'])
                    submission.sex = deserializeConfigValue(submissionChild.attrib['sex'])
                    submission.consanguineous = deserializeConfigValue(submissionChild.attrib['consanguineous'])
                    submission.homozygous = deserializeConfigValue(submissionChild.attrib['homozygous'])
                    submission.typedAlleles = deserializeConfigValue(submissionChild.attrib['typedalleles'])
                    submission.materialAvailability = deserializeConfigValue(submissionChild.attrib['materialavailability'])
                    submission.cellBank = deserializeConfigValue(submissionChild.attrib['cellbank'])
                    submission.primarySequencingMethodology = deserializeConfigValue(submissionChild.attrib['primarysequencingmethodology'])
                    submission.secondarySequencingMethodology = deserializeConfigValue(submissionChild.attrib['secondarysequencingmethodology'])
                    submission.primerType = deserializeConfigValue(submissionChild.attrib['primertype'])
                    submission.primers = deserializeConfigValue(submissionChild.attrib['primers'])
                    submission.sequencedInIsolation = deserializeConfigValue(submissionChild.attrib['sequencedinisolation'])
                    submission.numOfReactions = deserializeConfigValue(submissionChild.attrib['numofreactions'])
                    submission.methodComments = deserializeConfigValue(submissionChild.attrib['methodcomments'])
                    submission.citations = deserializeConfigValue(submissionChild.attrib['citations'])
                    submissionBatch.submissionBatch.append(submission)

                # Store my submission batch in the global variables.
                assignConfigurationValue('submission_batch', submissionBatch)

            else:
                # Any arbitrary configuration value, just store it.
                assignConfigurationValue(child.tag, child.text)


        # Here is where I assign the common/critical configuration values
        # test_submission indicates if we should use the "test" values.
        # I think I'll use this value for both EMBL and IMGT submissions, if it applies.
        assignIfNotExists('test_submission', '1')

        # Logging is turned off by default. Users can change this to 1 to turn on a logfile.
        assignIfNotExists('logging','0')

        # I'm storing FTP without the ftp:// identifier, because it is not necessary.
        # The test and prod ftp sites have the same address. This is intentional, embl doesn't have a test ftp
        # TODO : I still need this stuff? Probably. I think the act service does not need the method name anymore though.
        assignIfNotExists('embl_ftp_upload_site_test', 'webin.ebi.ac.uk')
        assignIfNotExists('embl_ftp_upload_site_prod', 'webin.ebi.ac.uk')
        assignIfNotExists('embl_rest_address_test', 'https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/')
        assignIfNotExists('embl_rest_address_prod', 'https://www.ebi.ac.uk/ena/submit/drop-box/submit/')
        #assignIfNotExists('nmdp_act_rest_address', 'http://act.b12x.org/type_align')
        assignIfNotExists('nmdp_act_rest_address', 'http://act.b12x.org')
        # TODO: i should use the nmdp configuration value when I call the GFE/ACT services. It is currently hardcoded AFAIK
