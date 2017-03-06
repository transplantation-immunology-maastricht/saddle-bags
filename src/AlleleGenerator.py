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
import tkMessageBox

import math

from HLAGene import *

# The AlleleGenerator class contains logic to generate an EMBL HLA allele submission 
# In ENA format.  
class AlleleGenerator():   
 
    inputFileName = ''
    outputFileName = ''
    sequenceAnnotation = HLAGene()
    inputCellNummer = 0#12345
    inputGene = ''#HLA-C'
    inputAllele = ''#C0316ext'  

    # This is a short wrapper method to use biopython's translation method. 
    def translateSequence(self,inputSequence):

        coding_dna = Seq(inputSequence, generic_dna)
        
        peptideSequence = str(coding_dna.translate())
        print ('Translated Protein:' + peptideSequence)
        
        #Stop codon *should* be at the end of the protein.  
        stopCodonLocation = peptideSequence.find('*')
        
        if (stopCodonLocation == -1):
            if(len(coding_dna) % 3 == 0):
                tkMessageBox.showinfo('No Stop Codon Found', 
                    'The translated protein does not contain a stop codon.' )
            else:
                tkMessageBox.showinfo('No Stop Codon Found', 
                    'The translated protein does not contain a stop codon.\n' + 
                    'It looks like a frame shift,\n' +
                    'The coding nucleotide sequence has length: ' + str(len(coding_dna)) + '.')
        else:
            if (stopCodonLocation == len(peptideSequence) - 1):
                #Stop codon is the last character in the peptide sequence.  That's just fine, but trim off the stop codon.
                peptideSequence = peptideSequence[0:stopCodonLocation]
                pass
            else:
                tkMessageBox.showinfo('Premature Stop Codon Detected',
                    'Premature stop codon found:\nPeptide Position (' + 
                    str(stopCodonLocation + 1) + '/' +
                    str(len(peptideSequence)) + ')\n\n' +
                    'Double check your peptide sequence,\n' + 
                    'Some aminos from the 3\' / C-Terminus\nwere spliced out.\n\n' + 
                    'Before : ' + peptideSequence + 
                    '\nAfter  : ' + peptideSequence[0:stopCodonLocation] + 
                    '\n'
                    )
                peptideSequence = peptideSequence[0:stopCodonLocation]

        return peptideSequence

    # 
    def readInputSequence(self):

        print('Reading file: ' + self.inputFileName)

        fileObject = open(self.inputFileName, 'r')        
        fullFile = fileObject.read()

        self.processInputSequence(fullFile)        

    # The input file should be a string of nucleotides, with capital letters to identify exons and introns.
    # Annotations are expected and read in this format:
    # fiveprimeutrEXONONEintrononeEXONTWOintrontwoEXONTHREEthreeprimeutr
    # agctagctagctAGCTAGCtagctagctAGCTAGCtagctagctAGCTAGCTAgctagctagctag
    # All spaces, line feeds, and tabs are removed and ignored.  
    def processInputSequence(self, inputSequenceText):

        resultGeneLoci = HLAGene()
        # Why do I need to initialize loci array?  I would have thought calling Gene() would do it.  
        resultGeneLoci.loci = []

        # Trim out any spaces, tabs, newlines.  Uppercase.
        cleanedGene = inputSequenceText.replace(' ','').replace('\n','').replace('\t','').replace('\r','')

        # Trim out the annotation marks, and capitalize, so I have a copy of the full sequence.
        unannotatedGene = cleanedGene.upper()
        resultGeneLoci.fullSequence = unannotatedGene
        print('Total Sequence Length = ' + str(len(unannotatedGene)))

        # Loop through the cleaned and annotated input sequence, 
        # to search for exon annotation characters ( '[' and ']' )
        # I assume that the first and last loci are the 5' and 3' UTR, 
        # I assume that Exons and Introns will alternate beyond that.
        # It no longer uses ( '[' and ']' ) to specify exons.  I check for 
        # capitals and lowercase letters to determine exon start and end
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
                        # Store an Exon in the list.
                        currentExon = GeneLocus()
                        currentExon.sequence = cleanedGene[locusBeginPosition:x].upper()
                        currentExon.exon = True
                        resultGeneLoci.loci.append(currentExon)     
                        insideAnExon = False
                        locusBeginPosition=x

                        #Starting a new Intron.
                        pass
            else:
                print('Nonstandard nucleotide detected at position ' + str(x) + ' : ' + currentChar 
                    + '.  If this is a wildcard character, you might be ok.')

        # Store the last(3') UTR as an intron.
        currentIntron = GeneLocus()
        currentIntron.sequence = cleanedGene[locusBeginPosition:len(cleanedGene)].upper()
        currentIntron.exon = False
        resultGeneLoci.loci.append(currentIntron)    

        # Annotate the loci (name them) and print the results of the read file.
        resultGeneLoci.annotateLoci()
        resultGeneLoci.printGeneSummary()

        self.sequenceAnnotation = resultGeneLoci

    # Create the text submission based on the ENA format.
    def buildENASubmission(self):    
        
        # ENA format is the preferred submission type for EMBL.  More information:
        # http://www.ebi.ac.uk/ena/submit/sequence-submission
        # http://www.ebi.ac.uk/ena/submit/entry-upload-templates
        # ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt
        # ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html
        # http://www.ebi.ac.uk/ena/software/flat-file-validator

        documentBuffer = ''

        # These variables are for test data, they should be filled in by GUI.
        #self.inputCellNummer = 23445
        #self.inputGene = 'HLA-C'
        #self.inputAllele = 'C0316ext'        

        completeSequence = self.sequenceAnnotation.getCompleteSequence()
        exonSequence = self.sequenceAnnotation.getExonSequence()
        totalLength = self.sequenceAnnotation.totalLength()
        featureCount = len(self.sequenceAnnotation.loci)
        print('total calculated length = ' + str(totalLength))

        # Print header
        documentBuffer += ('ID   XXX; XXX; linear; genomic DNA; XXX; XXX; ' + str(totalLength) + ' BP.\n')
        documentBuffer += ('XX\n')

        # A valid document should have an AC (Accession Number) and DE (Description) field.
        # I don't have an AC number available, so it's blank.
        documentBuffer += ('AC   \n')
        documentBuffer += ('XX\n')

        documentBuffer += ('DE   Human Leukocyte Antigen\n')
        documentBuffer += ('XX\n')

        # Print key
        documentBuffer += ('FH   Key             Location/Qualifiers\n')
        documentBuffer += ('FH\n')

        # Print source
        # It's from a human.
        documentBuffer += ('FT   source          1..' + str(totalLength) + '\n')
        documentBuffer += ('FT                   /organism="Homo sapiens"\n')
        documentBuffer += ('FT                   /db_xref="taxon:9606"\n')
        documentBuffer += ('FT                   /mol_type="genomic DNA"\n')
        documentBuffer += ('FT                   /chromosome="6"\n')
        documentBuffer += ('FT                   /isolate="' + str(self.inputCellNummer) + '"\n')    

        # Print mRNA
        documentBuffer += ('FT   mRNA            join(')
        # Iterate through the indices of the UTRs and exons.
        # The 3' and 5' UTR are included in the mRNA
        # But not in the CDS (coding sequence), since they're untranslated.
        documentBuffer += (str(self.sequenceAnnotation.loci[0].beginIndex) 
            + '..' + str(self.sequenceAnnotation.loci[0].endIndex) + ',')

        for x in range(1,featureCount-1):
            geneLocus = self.sequenceAnnotation.loci[x]
            if (geneLocus.exon):
                documentBuffer += str(geneLocus.beginIndex) + '..' + str(geneLocus.endIndex) + ','

        documentBuffer += (str(self.sequenceAnnotation.loci[featureCount-1].beginIndex) 
            + '..' + str(self.sequenceAnnotation.loci[featureCount-1].endIndex) + ')\n')

        documentBuffer += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
        documentBuffer += ('FT                   /allele="' + str(self.inputAllele) + '"\n')  
        documentBuffer += ('FT                   /product=\"MHC class I antigen\"\n')  
                
        # Print CDS
        # CDS is the coding sequence.  It should include the exons, but not the UTRs/Introns
        # The range 1:featureCount-1 will exclude the UTRs.
        documentBuffer += ('FT   CDS             join(') 
        for x in range(1,featureCount-1):
            geneLocus = self.sequenceAnnotation.loci[x]
            if (geneLocus.exon):
                documentBuffer += str(geneLocus.beginIndex) + '..' + str(geneLocus.endIndex)
                if not x==featureCount-2:
                    documentBuffer += ','
                else:
                    documentBuffer += ')\n'

        documentBuffer += ('FT                   /transl_table=1\n')
        documentBuffer += ('FT                   /codon_start=1\n')
        documentBuffer += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
        documentBuffer += ('FT                   /allele="' + str(self.inputAllele) + '"\n')  
        documentBuffer += ('FT                   /product=\"MHC class I antigen\"\n')  
        documentBuffer += ('FT                   /translation=\"')

        # Some simple formatting for the peptide sequence, making it human and computer readable.  
        # 80 peptides per line.  Except the first line, which is 66.
        # 66 is 80-14, where 14 is the length of { /translation=" }
        peptideSequence = self.translateSequence(exonSequence)
        if(len(peptideSequence) < 66):
            documentBuffer += (peptideSequence) + '\"\n'
        else:
            documentBuffer += peptideSequence[0:66] + '\n'
            i=66
            while (i < len(peptideSequence)):
                documentBuffer += 'FT                   ' + peptideSequence[i:i+80] + '\n'   
                i += 80
 
        # Print 5'UTR        
        utr = self.sequenceAnnotation.loci[0]
        documentBuffer += ('FT   5\'UTR           ' + str(utr.beginIndex) + '..' + str(utr.endIndex) + '\n')
        documentBuffer += ('FT                   /note=\"5\'UTR\"\n')
        documentBuffer += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
        documentBuffer += ('FT                   /allele="' + str(self.inputAllele) + '"\n')  

        # Print alternating Ex/Int/Ex
        for x in range(1,featureCount-1):
            currentFeature = self.sequenceAnnotation.loci[x]

            if(currentFeature.exon):
                documentBuffer += ('FT   exon            ' + str(currentFeature.beginIndex) 
                    + '..' + str(currentFeature.endIndex) + '\n')
            else:
                documentBuffer += ('FT   intron          ' + str(currentFeature.beginIndex) 
                    + '..' + str(currentFeature.endIndex) + '\n')

            geneNumber = int(math.ceil(x / 2.0))
            documentBuffer += ('FT                   /number=' + str(geneNumber) + '\n') 
            documentBuffer += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
            documentBuffer += ('FT                   /allele="' + str(self.inputAllele) + '"\n')  
                

        # Print 3'UTR
        utr = self.sequenceAnnotation.loci[len(self.sequenceAnnotation.loci)-1]
        documentBuffer += ('FT   3\'UTR           ' + str(utr.beginIndex) + '..' + str(utr.endIndex) + '\n')
        documentBuffer += ('FT                   /note=\"3\'UTR\"\n')
        documentBuffer += ('FT                   /gene="' + str(self.inputGene) + '"\n') 
        documentBuffer += ('FT                   /allele="' + str(self.inputAllele) + '"\n')  
        documentBuffer += ('XX\n')

        # Print sequence
        # There's a sweet biopython method which can count the nucleotides.
        # Bio.Seq.count('A')
        # I didn't use it. 
        cCount = 0
        gCount = 0
        tCount = 0
        aCount = 0
        otherCount = 0
        for nucleotide in completeSequence:
            if nucleotide == 'C':
                cCount+=1
            elif nucleotide == 'G':
                gCount+=1
            elif nucleotide == 'T':
                tCount+=1
            elif nucleotide == 'A':
                aCount+=1
            else:
                otherCount+=1

        documentBuffer += ('SQ   Sequence ' + str(totalLength) + ' BP; ' 
            + str(aCount) + ' A; ' + str(cCount) + ' C; ' 
            + str(gCount) + ' G; ' + str(tCount) + ' T; ' 
            + str(otherCount) + ' other;\n')

        # Here's some logic to print the sequence information in groups of 10.
        # This format is specified in the User manual specified by EMBL.
        rowCount = 0
        columnCount = 0
        currentSeqIndex = 0

        while (currentSeqIndex < totalLength):
            # The character code for a sequence region is two blank spaces,
            # followed by three blank spaces, for a total of 5 blanks.
            documentBuffer += '     '
            sequenceRow = completeSequence[currentSeqIndex : currentSeqIndex + 60]

            # A sequenceChunk is 10 nucleotides in this context.
            # Format specifies up to six "chunks" per line.
            for i in range(0,6):
                sequenceChunk = sequenceRow[i*10 : (i+1)*10]
                documentBuffer += sequenceChunk + ' '

            # If line is complete (=60 bp), we can print the nucleotide index and move on to the next row.
            if(len(sequenceRow) == 60):
                documentBuffer += str(currentSeqIndex + 60) + '\n'
            # but if line is not complete (this is more likely, and more complicated.)
            else:
                # Fill with spaces to align the nucleotide indices at the end of the sequence.
                numberSpaces = 60-len(sequenceRow)
                for n in range (0, numberSpaces):
                    documentBuffer += ' '
                documentBuffer += (str(len(sequenceRow) + currentSeqIndex) + '\n')
            
            # The next row of the sequence
            currentSeqIndex += 60     
   
        # Print entry terminator.  The last line of an ENA entry.
        documentBuffer += ('//\n')

        return documentBuffer

    # Simple method to write the results to a file on your computer.
    def outputENASubmissionToFile(self, outputText): 

        outputFileObject = open(self.outputFileName, 'w')  
        outputFileObject.write(outputText)
        outputFileObject.close()

