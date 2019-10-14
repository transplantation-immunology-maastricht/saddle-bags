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

#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

from saddlebags.AlleleSubmission import HlaGene, AlleleSubmission, SubmissionBatch
from saddlebags.AlleleSubCommon import translateSequence
from saddlebags.HlaSequenceException import HlaSequenceException

import logging

# The AlleleGenerator class contains logic to generate an EMBL HLA allele submission 
# In ENA format.  
class EnaSubGenerator():

    def __init__(self):
        #self.sequenceAnnotation = HlaGene()
        # Instead of storing an HlaGene, im now storing a ALleleSubmission object.
        # A ALleleSubmission object has it's own HLA gene, so this is better.
        # I think I need to also store the submission batch because it needs that info too. This is a bit redundant, but oh well.
        self.submission = AlleleSubmission()
        self.submissionBatch = SubmissionBatch()


    def printHeader(self):      
        #print('The EMBL Print Header Method.')
        headerText = ''
        
        # Print header
        totalLength = self.submission.submittedGene.totalLength()
        headerText += 'ID   XXX; XXX; linear; genomic DNA; XXX; XXX; ' + str(totalLength) + ' BP.\n'
        headerText += 'XX\n'
        # A valid document should have an AC (Accession Number) and DE (Description) field.
        # I don't have an AC number available, so it's blank.
        headerText += 'AC   \n'
        headerText += 'XX\n'
        #headerText += 'DE   Human Leukocyte Antigen\n'
        #Requested change to the DE line.  It should look like:
        #Homo sapiens HLA-B gene for MHC class I antigen, allele "/allele name"
        headerText += ('DE   Homo sapiens ' + str(self.submission.submittedGene.geneLocus)
            + ' gene for MHC class ' + str(self.submission.submittedGene.hlaClass)
            + ' antigen, allele "' + str(self.submission.localAlleleName) + '"\n')
        headerText += 'XX\n'

        # Print key
        headerText += ('FH   Key             Location/Qualifiers\n')
        headerText += ('FH\n')
        
        # Print source
        # It's from a human.
        headerText += ('FT   source          1..' + str(self.submission.submittedGene.totalLength()) + '\n')
        headerText += ('FT                   /organism="Homo sapiens"\n')
        headerText += ('FT                   /db_xref="taxon:9606"\n')
        headerText += ('FT                   /mol_type="genomic DNA"\n')
        headerText += ('FT                   /chromosome="6"\n')
        headerText += ('FT                   /isolate="' + str(self.submission.cellId) + '"\n')
        
        return headerText
    
    def printMRNA(self):
        mRNAText = ''
        # Print mRNA
        mRNAText += ('FT   mRNA            join(')
        
        # Iterate through the indices of the UTRs and exons.
        # The 3' and 5' UTR are included in the mRNA
        for x in range(0, len(self.submission.submittedGene.features)):
            geneLocus = self.submission.submittedGene.features[x]
            # If it is an exon or UTR
            if (geneLocus.exon or 'UT' in geneLocus.name):
                mRNAText += str(geneLocus.beginIndex) + '..' + str(geneLocus.endIndex) + ','

        # Trim off the last comma and add a parenthese
        mRNAText = mRNAText[0:len(mRNAText)-1] + ')\n'

        mRNAText += ('FT                   /gene="' + str(self.submission.submittedGene.geneLocus) + '"\n')
        mRNAText += ('FT                   /allele="' + str(self.submission.localAlleleName) + '"\n')
        # TODO: Is this sufficient to allow I, II, 1, and 2?
        # TODO: What if it's class III?
        mRNAText += ('FT                   /product=\"MHC class ' + str(('I' if ('1'==str(self.submission.submittedGene.hlaClass)) else 'II')) + ' antigen\"\n')
        
        return mRNAText
    
    
    def printCDS(self):
        # I need to perform the translation first, so I know if this is a "pseudogene" or not
        peptideSequence = translateSequence(self.submission)
        
        cdsText = ''
        
        # Print CDS
        # CDS is the coding sequence.  It should include the exons, but not the UTRs/Introns
        # The range 1:featureCount-1 will exclude the UTRs.
        cdsText += ('FT   CDS             join(') 
        for x in range(0, len(self.submission.submittedGene.features)):
            geneLocus = self.submission.submittedGene.features[x]
            if (geneLocus.exon):
                cdsText += str(geneLocus.beginIndex) + '..' + str(geneLocus.endIndex)
                if not x == len(self.submission.submittedGene.features) - 2:
                    cdsText += ','
                else:
                    cdsText += ')\n'

        cdsText += ('FT                   /transl_table=1\n')
        cdsText += ('FT                   /codon_start=1\n')
        
        # If this sequence has premature stop codon, add the "/pseudo" flag.
        # This indicates the gene is a /pseudo gene, not a complete protein.
        if(self.submission.isPseudoGene):
            logging.info("putting pseudo in the submission")
            cdsText += ('FT                   /pseudo\n')
        else:
            logging.info("not putting pseudo in the submission")
            pass
        
        
        cdsText += ('FT                   /gene="' + str(self.submission.submittedGene.geneLocus) + '"\n')
        cdsText += ('FT                   /allele="' + str(self.submission.localAlleleName) + '"\n')
        cdsText += ('FT                   /product=\"MHC class ' + str(('I' if ('1'==str(self.submission.submittedGene.hlaClass)) else 'II')) + ' antigen\"\n')
        cdsText += ('FT                   /translation=\"')

        # Some simple formatting for the peptide sequence, making it human and computer readable.  
        # 80 peptides per line.  Except the first line, which is 66.
        # 66 is 80-14, where 14 is the length of { /translation=" }
        
        # The translation is commented out here. I had to move it to the top of this method.
        #peptideSequence = self.translateSequence(self.submission.submittedGene.getExonSequence())
        if(len(peptideSequence) < 66):
            cdsText += (peptideSequence) + '\"\n'
        else:
            cdsText += peptideSequence[0:66] + '\n'
            i=66
            while (i < len(peptideSequence)):
                cdsText += 'FT                   ' + peptideSequence[i:i+80]   
                i += 80
                
                # If we're not yet at the end of the sequence, go to the next line
                if(i < len(peptideSequence)):
                    cdsText += '\n'
                # We're at the end. close the quote and new line.
                else:
                    cdsText += '\"\n'
                
        return cdsText
    
    def printFeatures(self):
        featureText = ''
        
        exonIndex = 1
        intronIndex = 1
        
        geneHas3UTR = False
        geneHas5UTR = False
            
        for x in range(0, len(self.submission.submittedGene.features)):
            currentFeature = self.submission.submittedGene.features[x]

            # 3' UTR
            if(currentFeature.name == '3UT'):
                featureText += ('FT   3\'UTR           ' + str(currentFeature.beginIndex) + '..' + str(currentFeature.endIndex) + '\n')
                featureText += ('FT                   /note=\"3\'UTR\"\n')
                featureText += ('FT                   /gene="' + str(self.submission.submittedGene.geneLocus) + '"\n')
                featureText += ('FT                   /allele="' + str(self.submission.localAlleleName) + '"\n')
                geneHas3UTR = True  
                
            # 5' UTR
            elif(currentFeature.name == '5UT'):
                featureText += ('FT   5\'UTR           ' + str(currentFeature.beginIndex) + '..' + str(currentFeature.endIndex) + '\n')
                featureText += ('FT                   /note=\"5\'UTR\"\n')
                featureText += ('FT                   /gene="' + str(self.submission.submittedGene.geneLocus) + '"\n')
                featureText += ('FT                   /allele="' + str(self.submission.localAlleleName) + '"\n')
                geneHas5UTR = True   
            
            # Exon
            elif(currentFeature.exon):
                featureText += ('FT   exon            ' + str(currentFeature.beginIndex) 
                    + '..' + str(currentFeature.endIndex) + '\n')
                featureText += ('FT                   /number=' + str(exonIndex) + '\n') 
                featureText += ('FT                   /gene="' + str(self.submission.submittedGene.geneLocus) + '"\n')
                featureText += ('FT                   /allele="' + str(self.submission.localAlleleName) + '"\n')
                exonIndex += 1
            
            # Intron
            else:
                featureText += ('FT   intron          ' + str(currentFeature.beginIndex) 
                    + '..' + str(currentFeature.endIndex) + '\n')   
                featureText += ('FT                   /number=' + str(intronIndex) + '\n') 
                featureText += ('FT                   /gene="' + str(self.submission.submittedGene.geneLocus) + '"\n')
                featureText += ('FT                   /allele="' + str(self.submission.localAlleleName) + '"\n')
                intronIndex += 1

       
        featureText += ('XX\n')
        
        # Do a quick sanity check.  If we are missing either UTR I should warn the user.
        # But move on with your life, this is not worth getting upset over.
        if (not geneHas3UTR and not geneHas5UTR):
            raise HlaSequenceException( 
               'This sequence has no 5\' or 3\' UTR.\n\n' + 
                'Use lowercase nucleotides at the\n' + 
                'beginning and end of your DNA\n' +
                'sequence to specify the 5\' and 3\' UTRs.' )
        elif (not geneHas5UTR):
            raise HlaSequenceException(  
                'This sequence has no 5\' UTR.\n\n' + 
                'Use lowercase nucleotides at the\n' + 
                'beginning and end of your DNA\n' +
                'sequence to specify the 5\' and 3\' UTRs.' )           
        elif (not geneHas3UTR):
            raise HlaSequenceException(   
                'This sequence has no 3\' UTR.\n\n' + 
                'Use lowercase nucleotides at the\n' + 
                'beginning and end of your DNA\n' +
                'sequence to specify the 5\' and 3\' UTRs.' ) 
        else:
            logging.info('The UTRs look fine.')
            pass
        
        return featureText
    
    def printSequence(self):
        sequenceText = ''
 
        completeSequence = self.submission.submittedGene.getCompleteSequence().upper()
                
        cCount = completeSequence.count('C')
        gCount = completeSequence.count('G')
        tCount = completeSequence.count('T')
        aCount = completeSequence.count('A')
        otherCount = self.submission.submittedGene.totalLength() - (cCount + gCount + tCount + aCount)

        sequenceText += ('SQ   Sequence ' + str(self.submission.submittedGene.totalLength()) + ' BP; '
            + str(aCount) + ' A; ' + str(cCount) + ' C; ' 
            + str(gCount) + ' G; ' + str(tCount) + ' T; ' 
            + str(otherCount) + ' other;\n')

        # Here's some logic to print the sequence information in groups of 10.
        # This format is specified in the User manual specified by EMBL.
        currentSeqIndex = 0

        while (currentSeqIndex < self.submission.submittedGene.totalLength()):
            # The character code for a sequence region is two blank spaces,
            # followed by three blank spaces, for a total of 5 blanks.
            sequenceText += '     '
            sequenceRow = self.submission.submittedGene.getCompleteSequence()[currentSeqIndex : currentSeqIndex + 60]

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

        try:
            documentBuffer = ''
    
            totalLength = self.submission.submittedGene.totalLength()
            logging.info('total calculated length = ' + str(totalLength))
            
            if(totalLength > 0 and self.validateInputs()):
    
                # These are the main sections of the ENA submission.
                documentBuffer += self.printHeader()
                documentBuffer += self.printMRNA()
                documentBuffer += self.printCDS()
                documentBuffer += self.printFeatures()
                documentBuffer += self.printSequence()
        
                # Print entry terminator.  The last line of an ENA entry.
                documentBuffer += ('//\n')
                
            else: 
                #tkMessageBox.showinfo('No HLA Sequence Found', 
                #    'The HLA sequence is empty.\nPlease fill in an annotated HLA sequence\nbefore generating the submission.' )
                # TODO: Make a specific exception type for this.

                if(not totalLength > 0):
                    raise Exception('The HLA sequence is empty.\nPlease fill in an annotated HLA sequence\nbefore generating the submission.' )
                else:
                    raise Exception ('Inputs do not look valid. Not sure what is missing.')
                return None
    
    
            return documentBuffer
        except HlaSequenceException:
            # Throw a fit.
            raise
    
    # Return True if our input values are all present and accomodated for.
    # If something is missing, then throw a fit and give up.
    # TODO: I should probably not raise these exceptions actually.
    # Instead, I should have the GUI Automatically open the choose options screen
    
    # TODO: Maybe I should delete this method, and add error handling to the generate methods.
    def validateInputs(self):
        #raise Exception ('Validate Inputs Method is being used, after all.')
        
        if (self.submission.cellId is None or len(self.submission.cellId) < 1):
            logging.error('Invalid Sequence ID:' + str(self.submission.cellId))
            return False
        
        elif (self.submission.submittedGene is None):
            logging.error('Invalid Submitted Gene:' + str(self.submission.submittedGene))
            return False
        
        elif (self.submission.submittedGene.geneLocus is None or len(self.submission.submittedGene.geneLocus) < 1):
            logging.error('Invalid Input Gene:' + str(self.submission.submittedGene.geneLocus))
            return False
        
        elif (self.submission.localAlleleName is None or len(self.submission.localAlleleName) < 1):
            logging.error('Invalid Input Allele:' + str(self.submission.localAlleleName))
            return False

        elif (self.submission.submittedGene.hlaClass is None or len(self.submission.submittedGene.hlaClass) < 1):
            logging.error('Invalid Input Class:' + str(self.submission.submittedGene.hlaClass))
            return False
        
        else:
            return True

