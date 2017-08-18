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
# You should have received a copy of the GNU Lesser General Public License
# along with saddle-bags. If not, see <http://www.gnu.org/licenses/>.

import xml.etree.ElementTree as ET
import xml.dom.minidom

from os.path import isdir, split
from os import makedirs

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# TODO: Maybe I shouldn't have GUI methods in here.
# But it's hard to refactor those out.
# I can use the HlaSequenceException class.  And catch the exception further upstream
# That's the strategy.
import tkMessageBox

import sys
from os.path import join, isfile, expanduser

from HlaGene import HlaGene, GeneLocus

# This is a short wrapper method to use biopython's translation method. 
# Most of this code is just checking for things that went wrong
def translateSequence(inputSequence):

    proteinSequence = ''
    
    try:
        # Do nothing if the input sequence is blank.
        if( len(inputSequence) > 0 ):
            
            coding_dna = Seq(inputSequence, generic_dna)        
            proteinSequence = str(coding_dna.translate())   
            print ('Exon Sequence before translation:' + coding_dna)     
            print ('Translated Protein:' + proteinSequence)
            
            # Perform Sanity Checks.
            # Stop codon *should* be at the end of the protein.  
            # Here we seek out the first instance of a stop codon, 
            # and remove the peptides afterwards.
            # because that's what happens in real life.
            stopCodonLocation = proteinSequence.find('*')
            
            # If no stop codon was found
            if (stopCodonLocation == -1):
                assignConfigurationValue('is_pseudo_gene','1')
                # If multiple of three (correct codon length)
                if(len(coding_dna) % 3 == 0):
                    tkMessageBox.showinfo('No Stop Codon Found', 
                        'The translated protein does not contain a stop codon.\n' + 
                        'This is indicated by a /pseudo flag in the sequence submission.'
                         )
                    
                # Wrong Codon Length
                else:
                    tkMessageBox.showinfo('No Stop Codon Found', 
                        'The translated protein does not contain a stop codon.\n' + 
                        'The coding nucleotide sequence length (' + str(len(coding_dna))  + ') is not a multiple of 3.\n' + 
                        'This is indicated by a /pseudo flag in the sequence submission.')

            # If Stop Codon is in the end of the protein (This is expected and correct)
            elif (stopCodonLocation == len(proteinSequence) - 1):
                assignConfigurationValue('is_pseudo_gene','0')
                
                # If multiple of three (correct codon length)
                if(len(coding_dna) % 3 == 0):
                    # Everything is fine in this case.  Trim off the stop codon
                    proteinSequence = proteinSequence[0:stopCodonLocation]
                    pass 
                # Wrong Codon Length
                else:
                    tkMessageBox.showinfo('Extra Nucleotides After the Stop Codon', 
                        'The stop codon is at the correct position in the protein, but ' + 
                        'The coding nucleotide sequence length (' + str(len(coding_dna))  + ') is not a multiple of 3.\n\n' +
                        'Please double check your sequence.')
                    proteinSequence = proteinSequence[0:stopCodonLocation]
                                        
            # Else Stop Codon is premature (before the end of the protein) 
            else:
                assignConfigurationValue('is_pseudo_gene','1')
                
                # If multiple of three (correct codon length)
                if(len(coding_dna) % 3 == 0):
                    tkMessageBox.showinfo('Premature Stop Codon Detected',
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
                    tkMessageBox.showinfo('Premature Stop Codon Detected',
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
            print('Translating a nucleotide sequence of length 0.  That was easy.')
            pass

        return proteinSequence
    
    except Exception:
        print 'Problem when translating protein:'
        print sys.exc_info()[1]
        tkMessageBox.showinfo('Protein Translation Error', 
            'I could not translate your protein:\n' +  str(sys.exc_info()[1]))
        
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

def parseExons(roughFeatureSequence, alleleCallWithGFE):
    
    # TODO: Parse the JSON in the alleleCallWithGFE.  
    # There should be some information about the exons in here.
    # ex = uppercase
    # utr = lowercase
    
    annotatedSequence = roughFeatureSequence   
    annotatedSequence = 'aaaCCCgggTTTaaacgttga'
    return annotatedSequence
    

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
    print('Total Sequence Length = ' + str(len(unannotatedGene)))

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
            print('Nonstandard Nucleotide at the beginning of the sequence, not sure how to handle it')
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
                print('Nonstandard nucleotide detected at position ' + str(x) + ' : ' + currentChar 
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
        print('Empty sequence, I don\'t have anything to do.')
        
    return resultGeneLoci    
    #self.sequenceAnnotation = resultGeneLoci

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput

# I'm storing global variables in a dictionary for now. 
def initializeGlobalVariables():    
    global globalVariables 
    
    if not ("globalVariables" in globals()):
        globalVariables={}
        
def assignConfigurationValue(configurationKey, configurationValue):
    initializeGlobalVariables()
    globalVariables[configurationKey] = configurationValue
    
def getConfigurationValue(configurationKey):
    if configurationKey in globalVariables.keys():
        return globalVariables[configurationKey]
    else:
        print ('Configuration Key Not Found:' + configurationKey)
        #raise KeyError('Key Not Found:' + configurationKey)
        return None

def assignConfigName():
    assignConfigurationValue('config_file_location',join(expanduser("~"),'Saddlebags.Config.xml'))
    
def writeConfigurationFile():
    assignConfigName()
    print ('Writing a config file to:\n' + globalVariables['config_file_location'])
    
    root = ET.Element("config")
    
    for key in globalVariables.keys():
        # Some config values I don't want to store.
        # Add to this: Sequence Text, EMBL Submission Text, IMGT Submission Text
        if(key not in [
            'embl_password'
            ,'imgt_password'
            , 'sequence'
            ]):
            ET.SubElement(root, key).text = globalVariables[key]

    xmlText = ET.tostring(root, encoding='utf8', method='xml')
    prettyXmlText = xml.dom.minidom.parseString(xmlText).toprettyxml()
    
    xmlOutput = createOutputFile(globalVariables['config_file_location'])
    xmlOutput.write(prettyXmlText)
    xmlOutput.close()

def loadConfigurationFile():
    assignConfigName()
    
    if not isfile(globalVariables['config_file_location']):
        print ('The config file does not exist yet. I will not load it:\n' + globalVariables['config_file_location'])
        
        # Here is where I assign the common configuration values
        # test_submission indicates if we should use the "test" values.
        # I think I'll use this value for both EMBL and IMGT submissions, if it applies.
        assignConfigurationValue('test_submission', '1')
        
        # I'm storing FTP without the ftp:// identifier, because it is not necessary.
        # The test and prod ftp sites have the same address. This is intentional, embl doesn't have a test ftp
        assignConfigurationValue('embl_ftp_upload_site_test', 'webin.ebi.ac.uk')
        assignConfigurationValue('embl_ftp_upload_site_prod', 'webin.ebi.ac.uk')
        assignConfigurationValue('embl_rest_address_test', 'https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/')
        assignConfigurationValue('embl_rest_address_prod', 'https://www.ebi.ac.uk/ena/submit/drop-box/submit/')
        assignConfigurationValue('nmdp_act_rest_address', 'http://act.b12x.org/act' )

    else:
        print ('The config file already exists, I will load it:\n' + globalVariables['config_file_location'])
        
        tree = ET.parse(globalVariables['config_file_location'])
        root = tree.getroot()
        
        for child in root:
            assignConfigurationValue(child.tag, child.text)
            
       