# This file is part of EMBL-HLA-Submission.
#
# EMBL-HLA-Submission is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EMBL-HLA-Submission is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with EMBL-HLA-Submission. If not, see <http://www.gnu.org/licenses/>.

# Version 1.0

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys
import tkMessageBox

import math

from HLAGene import *

# The AlleleGenerator class contains logic to generate an EMBL HLA allele submission 
# In ENA format.  
class AlleleGenerator():   
    
    def __init__(self):
 
        self.inputFileName = ''
        self.outputFileName = ''
        self.sequenceAnnotation = HLAGene()
        self.inputCellNummer = 0
        self.inputGene = ''
        self.inputAllele = '' 

    # This is a short wrapper method to use biopython's translation method. 
    # Most of this code is just checking for things that went wrong
    def translateSequence(self,inputSequence):

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
                    # If multiple of three (correct codon length)
                    if(len(coding_dna) % 3 == 0):
                        tkMessageBox.showinfo('No Stop Codon Found', 
                            'The translated protein does not contain a stop codon.' )
                        
                    # Wrong Codon Length
                    else:
                        tkMessageBox.showinfo('No Stop Codon Found', 
                            'The translated protein does not contain a stop codon.\n' + 
                            'The coding nucleotide sequence length (' + str(len(coding_dna))  + ') is not a multiple of 3.')

                # If Stop Codon is in the end of the protein (This is expected and correct)
                elif (stopCodonLocation == len(proteinSequence) - 1):
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
                    # If multiple of three (correct codon length)
                    if(len(coding_dna) % 3 == 0):
                        tkMessageBox.showinfo('Premature Stop Codon Detected',
                            'Premature stop codon found:\nProtein Position (' + 
                            str(stopCodonLocation + 1) + '/' +
                            str(len(proteinSequence)) + ')\n\n' +
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

    # The input file should be a string of nucleotides, with capital letters to identify exons and introns.
    # Annotations are expected and read in this format:
    # fiveprimeutrEXONONEintrononeEXONTWOintrontwoEXONTHREEthreeprimeutr
    # agctagctagctAGCTAGCtagctagctAGCTAGCtagctagctAGCTAGCTAgctagctagctag
    # All spaces, line feeds, and tabs are removed and ignored.  
    def processInputSequence(self, inputSequenceText):

        resultGeneLoci = HLAGene()
        
        # Trim out any spaces, tabs, newlines.  Uppercase.
        cleanedGene = inputSequenceText.replace(' ','').replace('\n','').replace('\t','').replace('\r','')

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
            
        self.sequenceAnnotation = resultGeneLoci

    

    def printHeader(self):
        headerText = ''
        
        # Print header
        headerText += 'ID   XXX; XXX; linear; genomic DNA; XXX; XXX; ' + str(self.sequenceAnnotation.totalLength()) + ' BP.\n'
        headerText += 'XX\n'
        # A valid document should have an AC (Accession Number) and DE (Description) field.
        # I don't have an AC number available, so it's blank.
        headerText += 'AC   \n'
        headerText += 'XX\n'
        headerText += 'DE   Human Leukocyte Antigen\n'
        headerText += 'XX\n'

        # Print key
        headerText += ('FH   Key             Location/Qualifiers\n')
        headerText += ('FH\n')
        
        # Print source
        # It's from a human.
        headerText += ('FT   source          1..' + str(self.sequenceAnnotation.totalLength()) + '\n')
        headerText += ('FT                   /organism="Homo sapiens"\n')
        headerText += ('FT                   /db_xref="taxon:9606"\n')
        headerText += ('FT                   /mol_type="genomic DNA"\n')
        headerText += ('FT                   /chromosome="6"\n')
        headerText += ('FT                   /isolate="' + str(self.inputCellNummer) + '"\n')    
        
        return headerText
    
    def printMRNA(self):
        mRNAText = ''
        # Print mRNA
        mRNAText += ('FT   mRNA            join(')
        
        # Iterate through the indices of the UTRs and exons.
        # The 3' and 5' UTR are included in the mRNA
        for x in range(0,len(self.sequenceAnnotation.loci)):
            geneLocus = self.sequenceAnnotation.loci[x]
            # If it is an exon or UTR
            if (geneLocus.exon or 'UT' in geneLocus.name):
                mRNAText += str(geneLocus.beginIndex) + '..' + str(geneLocus.endIndex) + ','

        # Trim off the last comma and add a parenthese
        mRNAText = mRNAText[0:len(mRNAText)-1] + ')\n'

        mRNAText += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
        mRNAText += ('FT                   /allele="' + str(self.inputAllele) + '"\n')  
        mRNAText += ('FT                   /product=\"MHC class I antigen\"\n')  
        
        return mRNAText
    
    
    def printCDS(self):
        cdsText = ''
        
        # Print CDS
        # CDS is the coding sequence.  It should include the exons, but not the UTRs/Introns
        # The range 1:featureCount-1 will exclude the UTRs.
        cdsText += ('FT   CDS             join(') 
        for x in range(0,len(self.sequenceAnnotation.loci)):
            geneLocus = self.sequenceAnnotation.loci[x]
            if (geneLocus.exon):
                cdsText += str(geneLocus.beginIndex) + '..' + str(geneLocus.endIndex)
                if not x==len(self.sequenceAnnotation.loci)-2:
                    cdsText += ','
                else:
                    cdsText += ')\n'

        cdsText += ('FT                   /transl_table=1\n')
        cdsText += ('FT                   /codon_start=1\n')
        cdsText += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
        cdsText += ('FT                   /allele="' + str(self.inputAllele) + '"\n')  
        cdsText += ('FT                   /product=\"MHC class I antigen\"\n')  
        cdsText += ('FT                   /translation=\"')

        # Some simple formatting for the peptide sequence, making it human and computer readable.  
        # 80 peptides per line.  Except the first line, which is 66.
        # 66 is 80-14, where 14 is the length of { /translation=" }
        peptideSequence = self.translateSequence(self.sequenceAnnotation.getExonSequence())
        if(len(peptideSequence) < 66):
            cdsText += (peptideSequence) + '\"\n'
        else:
            cdsText += peptideSequence[0:66] + '\n'
            i=66
            while (i < len(peptideSequence)):
                cdsText += 'FT                   ' + peptideSequence[i:i+80] + '\n'   
                i += 80
                
        return cdsText
    
    def printFeatures(self):
        featureText = ''
        
        exonIndex = 1
        intronIndex = 1
        
        geneHas3UTR = False
        geneHas5UTR = False
            
        for x in range(0,len(self.sequenceAnnotation.loci)):
            currentFeature = self.sequenceAnnotation.loci[x]

            # 3' UTR
            if(currentFeature.name == '3UT'):
                featureText += ('FT   3\'UTR           ' + str(currentFeature.beginIndex) + '..' + str(currentFeature.endIndex) + '\n')
                featureText += ('FT                   /note=\"3\'UTR\"\n')
                featureText += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
                featureText += ('FT                   /allele="' + str(self.inputAllele) + '"\n')
                geneHas3UTR = True  
                
            # 5' UTR
            elif(currentFeature.name == '5UT'):
                featureText += ('FT   5\'UTR           ' + str(currentFeature.beginIndex) + '..' + str(currentFeature.endIndex) + '\n')
                featureText += ('FT                   /note=\"5\'UTR\"\n')
                featureText += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
                featureText += ('FT                   /allele="' + str(self.inputAllele) + '"\n')
                geneHas5UTR = True   
            
            # Exon
            elif(currentFeature.exon):
                featureText += ('FT   exon            ' + str(currentFeature.beginIndex) 
                    + '..' + str(currentFeature.endIndex) + '\n')
                featureText += ('FT                   /number=' + str(exonIndex) + '\n') 
                featureText += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
                featureText += ('FT                   /allele="' + str(self.inputAllele) + '"\n')  
                exonIndex += 1
            
            # Intron
            else:
                featureText += ('FT   intron          ' + str(currentFeature.beginIndex) 
                    + '..' + str(currentFeature.endIndex) + '\n')   
                featureText += ('FT                   /number=' + str(intronIndex) + '\n') 
                featureText += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
                featureText += ('FT                   /allele="' + str(self.inputAllele) + '"\n') 
                intronIndex += 1

       
        featureText += ('XX\n')
        
        # Do a quick sanity check.  If we are missing either UTR I should warn the user.
        # But move on with your life, this is not worth getting upset over.
        if (not geneHas3UTR and not geneHas5UTR):
            tkMessageBox.showinfo('Missing UTRs', 
                'This sequence has no 5\' or 3\' UTR.\n\n' + 
                'Use lowercase nucleotides at the\n' + 
                'beginning and end of your DNA\n' +
                'sequence to specify the 5\' and 3\' UTRs.' )
        elif (not geneHas5UTR):
            tkMessageBox.showinfo('Missing 5\' UTR', 
                'This sequence has no 5\' UTR.\n\n' + 
                'Use lowercase nucleotides at the\n' + 
                'beginning and end of your DNA\n' +
                'sequence to specify the 5\' and 3\' UTRs.' )            
        elif (not geneHas3UTR):
            tkMessageBox.showinfo('Missing 3\' UTR', 
                'This sequence has no 3\' UTR.\n\n' + 
                'Use lowercase nucleotides at the\n' + 
                'beginning and end of your DNA\n' +
                'sequence to specify the 5\' and 3\' UTRs.' )    
        else:
            print('The UTRs look fine.')
            pass
        
        return featureText
    
    def printSequence(self):
        sequenceText = ''
 
        completeSequence = self.sequenceAnnotation.getCompleteSequence().upper()
                
        cCount = completeSequence.count('C')
        gCount = completeSequence.count('G')
        tCount = completeSequence.count('T')
        aCount = completeSequence.count('A')
        otherCount = self.sequenceAnnotation.totalLength() - (cCount + gCount + tCount + aCount)

        sequenceText += ('SQ   Sequence ' + str(self.sequenceAnnotation.totalLength()) + ' BP; ' 
            + str(aCount) + ' A; ' + str(cCount) + ' C; ' 
            + str(gCount) + ' G; ' + str(tCount) + ' T; ' 
            + str(otherCount) + ' other;\n')

        # Here's some logic to print the sequence information in groups of 10.
        # This format is specified in the User manual specified by EMBL.
        currentSeqIndex = 0

        while (currentSeqIndex < self.sequenceAnnotation.totalLength()):
            # The character code for a sequence region is two blank spaces,
            # followed by three blank spaces, for a total of 5 blanks.
            sequenceText += '     '
            sequenceRow = self.sequenceAnnotation.getCompleteSequence()[currentSeqIndex : currentSeqIndex + 60]

            # A sequenceChunk is 10 nucleotides in this context.
            # Format specifies up to six "chunks" per line.
            for i in range(0,6):
                sequenceChunk = sequenceRow[i*10 : (i+1)*10]
                sequenceText += sequenceChunk + ' '

            # If line is complete (=60 bp), we can print the nucleotide index and move on to the next row.
            if(len(sequenceRow) == 60):
                sequenceText += str(currentSeqIndex + 60) + '\n'
            # but if line is not complete (this is more likely, and more complicated.)
            else:
                # Fill with spaces to align the nucleotide indices at the end of the sequence.
                numberSpaces = 60-len(sequenceRow)
                for n in range (0, numberSpaces):
                    sequenceText += ' '
                sequenceText += (str(len(sequenceRow) + currentSeqIndex) + '\n')
            
            # The next row of the sequence
            currentSeqIndex += 60     
            
        return sequenceText
            
        
    # Create the text submission based on the ENA format.
    def buildENASubmission(self):    
        
        # ENA format is the preferred submission type for EMBL.  More information:
        # http://www.ebi.ac.uk/ena/submit/sequence-submission
        # http://www.ebi.ac.uk/ena/submit/entry-upload-templates
        # ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
        # ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html
        # http://www.ebi.ac.uk/ena/software/flat-file-validator

        documentBuffer = ''

        totalLength = self.sequenceAnnotation.totalLength()
        print('total calculated length = ' + str(totalLength))
        
        if(totalLength > 0):

            # These are the main sections of the ENA submission.
            documentBuffer += self.printHeader()
            documentBuffer += self.printMRNA()
            documentBuffer += self.printCDS()
            documentBuffer += self.printFeatures()
            documentBuffer += self.printSequence()
    
            # Print entry terminator.  The last line of an ENA entry.
            documentBuffer += ('//\n')
            
        else: 
            tkMessageBox.showinfo('No HLA Sequence Found', 
                'The HLA sequence is empty.\nPlease fill in an annotated HLA sequence\nbefore generating the submission.' )
            
            pass
        

        return documentBuffer

    # Simple method to write the results to a file on your computer.
    def outputENASubmissionToFile(self, outputText): 

        outputFileObject = open(self.outputFileName, 'w')  
        outputFileObject.write(outputText)
        outputFileObject.close()

