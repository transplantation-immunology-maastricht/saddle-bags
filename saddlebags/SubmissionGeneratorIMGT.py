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

import sys

import datetime
import tkMessageBox

from AlleleSubCommon import *
#import math

from HLAGene import *

# The AlleleGenerator class contains logic to generate an IMGT HLA allele submission 
# In ENA format.  
class SubmissionGeneratorIMGT():   
    
    def __init__(self):
        self.sequenceAnnotation = HLAGene()
    
    # Create the text submission based on the IMGT format.
    def buildIMGTSubmission(self):    
        
        documentBuffer = ''

        totalLength = self.sequenceAnnotation.totalLength()
        print('total calculated length = ' + str(totalLength))
        
        if(totalLength > 0):
            
            print ('im gonna add the header in here:')

            # These are the main sections of the ENA submission.
            documentBuffer += self.printHeader()
            documentBuffer += self.printSubmitter()
            documentBuffer += self.printSource()
            documentBuffer += self.printMethods()
            documentBuffer += self.printFeatures()
            documentBuffer += self.printSequence()

            # Print entry terminator.  The last line of an ENA entry.
            documentBuffer += ('//\n')
            
        else: 
            tkMessageBox.showinfo('No HLA Sequence Found', 
                'The HLA sequence is empty.\nPlease fill in an annotated HLA sequence\nbefore generating the submission.' )
            
            pass
        

        return documentBuffer


    def printHeader(self):
        
        headerText = ''
        
        # TODO: Get these values from IMGT, they shouldn't be hardcoded.  
        # Maybe it should be an unknown identifier with 
        imgtIdentifier = 'HWS10012345'
        imgtIdentifierWithVersion = 'HWS10012345.1'        
        currentSubmissionDate = '{:%d/%m/%Y}'.format(datetime.datetime.now())

        headerText += 'ID   ' + str(imgtIdentifier) + '; Sequence Submission; Confidential; ' + str(self.sequenceAnnotation.totalLength()) + ' BP.\n'
        headerText += 'XX\n'
        headerText += 'AC   ' + str(imgtIdentifier) + ';\n'
        headerText += 'XX\n'
        headerText += 'SV   ' + str(imgtIdentifierWithVersion) + '\n'
        headerText += 'XX\n'
        headerText += 'DT   ' + str(currentSubmissionDate) + ' (Submitted)\n'
        headerText += 'DT   ' + str(getConfigurationValue('embl_release_date')) + ' (Release)\n'
        headerText += 'XX\n'

        # TODO: I'm using the local allele name that is assigned by the user.
        # Maybe this allele name should be based on the closest allele.
        # Do I want the allele name, or should I generate a new one based on the closest allele?
                
        headerText += 'DE   ' + str(getConfigurationValue('allele_name')) + '\n'
        headerText += 'XX\n'
        headerText += 'KW   HLA WEB SUBMISSION;\n'
        headerText += 'XX\n'
        
        # The new allele description is split into multiple lines. I should add a new 'CC' line for each part of the description.
        rawDescription = str(getConfigurationValue('closest_allele_written_description'))
        rawDescriptionLineTokens = rawDescription.split('\n')
        for lineToken in rawDescriptionLineTokens:
            headerText += 'CC   ' + lineToken + '\n'
            
        #headerText += 'CC   A*03:01:01:01new is identical to A*03:01:01:01 except for position 382 is a A\n'
        #headerText += 'CC   in the new allele. This result in an amino change from W to stopcodon.\n'
        
        headerText += 'XX\n'
        headerText += 'OS   Homo sapiens (human);\n'
        headerText += 'OC   Eukaryota; Metazoa; Chordata; Vertebrata; Mammalia; Eutheria; Primates;\n'
        headerText += 'OC   Catarrhini; Hominidae; Homo.\n'
        headerText += 'XX\n'
        # TODO: Our submission says GENBANK, but we're using EMBL Numbers.  Also what does that [1] mean?
        headerText += 'DR   GENBANK; ' + str(getConfigurationValue('embl_sequence_accession')) + '.\n'
        headerText += 'XX\n'
        headerText += 'RN   [1]\n'
        # TODO: This submission is Unpublished.  What if it is published?
        # Ask James what a published study looks like. I need to include study name etc.
        headerText += 'RC   Unpublished.\n'
        headerText += 'XX\n'
        headerText += 'FH   Key            Location/Qualifier\n'
        headerText += 'FH\n'
        
        return headerText

    def printSubmitter(self):
        submitterText = ''
        
        # TODO: I don't know any of this data.  Should it be int he form?
        # Maybe I just need the submitter ID, and i can or can not get the rest?
        # I should be able to calculate the indices, at least.

        submitterText += 'FT   submittor      1..' + str(self.sequenceAnnotation.totalLength()) + '\n'
        submitterText += 'FT                  /ID="**IMGT_SUBMITTER_EMAIL_ID**"\n'
        submitterText += 'FT                  /name="**IMGT_SUBMITTER_NAME**"\n'
        submitterText += 'FT                  /alt_contact=""\n'
        submitterText += 'FT                  /email="**IMGT_SUBMITTER_EMAIL_ADDRESS**"\n'
        
        return submitterText

    def printSource(self):
        sourceText = ''
        
        # TODO: Submitting Laboratory Information. Can this be fetched from IMGT?

        sourceText += 'FT   source         1..' + str(self.sequenceAnnotation.totalLength()) + '\n'
        sourceText += 'FT                  /cell_id="' + str(getConfigurationValue('sample_id')) + '"\n'
        sourceText += 'FT                  /ethnic_origin="' + str(getConfigurationValue('ethnic_origin')) + '"\n'
        sourceText += 'FT                  /sex="' + str(getConfigurationValue('sex')) + '"\n'
        sourceText += 'FT                  /consanguineous="' + str(getConfigurationValue('consanguineous')) + '"\n'
        sourceText += 'FT                  /homozygous="Yes"\n'
        sourceText += 'FT                  /lab_of_origin="**IMGT_SUBMITTING_LAB_NAME**"\n'
        sourceText += 'FT                  /lab_contact="**IMGT_SUBMITTER_NAME**"\n'
        
        # TODO: No Material Available.  What if Material is available?
        # I think I need to add this to the form still.
        # Same story with "cell_bank"
        
        sourceText += 'FT                  /material_available="No Material Available"\n'
        sourceText += 'FT                  /cell_bank="Not Available"\n'
        
        # TODO: James suggested that I only allow valid fully-sequenced alleles.
        # Should I validate this, or should I leave that work to IMGT?
        
        sourceText += 'FT                  /HLA-A*="02:01,03new"\n'
        sourceText += 'FT                  /HLA-C*="07,-"\n'
        sourceText += 'FT                  /HLA-B*="07,-"\n'
        sourceText += 'FT                  /HLA-DRB1*="15:01,-"\n'
        
        return sourceText
    
    def printMethods(self):
        methodsText = ''
        
        # TODO: Get primer info from the form.  Make sure this all is correct        
        
        methodsText += 'FT   method         1..' + str(self.sequenceAnnotation.totalLength()) + '\n'
        
        # TODO: What are the options for sequencing methodology?
        # I can provide an open-text field.        
        
        methodsText += 'FT                  /primary_sequencing="Direct sequencing of PCR product from DNA (SBT)"\n'
        methodsText += 'FT                  /secondary_sequencing="Direct sequencing of PCR product from DNA (SBT)"\n'
        methodsText += 'FT                  /type_of_primer="Both allele and locus specific"\n'
        methodsText += 'FT                  /sequenced_in_isolation="Yes"\n'
        
        # TODO Add these primers dynamically
        # A primer has these pieces of information
        # "ID" ("primer_1") "Sequence" "Feature" "locus/indices" 
        # locus seems to be genomic index, from the beginning of the sequence.
        # I suppose this has to be locations in the reference sequence?
        # I should store a dictionary of primers in the configuration. 
        # Errr, nodes underneath the Primer nodes.
        # They put a "tab" character between some of this data.  Why? Because Tabs, sigh.
                
        methodsText += 'FT                  /primer_1="97022    GAGCCCCGCTTCAACGCC    E2    257-274"\n'
        methodsText += 'FT                  /primer_2="09148    CCAGGCGTGGCTCTCAGA    5UT    -265--248"\n'
        methodsText += 'FT                  /primer_3="09152    AACCTACGTAGGGTCCTTCA    5UT    -161--142"\n'
        methodsText += 'FT                  /primer_4="09154    AGTGTCGTCGCGGTCGCT    5UT    -72--55"\n'
        methodsText += 'FT                  /primer_5="09167    CAGACSCCGAGGATGGCC    5UT    -12-6"\n'
        methodsText += 'FT                  /primer_6="09162    AACACCCAACACACATTAGGT    I7    2745-2765"\n'
        methodsText += 'FT                  /primer_7="09168    GGGAGCACAGGTCAGCGTGGGAAG    3UT    3075-3098"\n'
        methodsText += 'FT                  /primer_8="98008    GTTTAGGCCAAAAATYCCCCC    I2    635-655"\n'
        methodsText += 'FT                  /no_of_reactions="3"\n'
        methodsText += 'FT                  /sequencing_direction="Both"\n'
        
        
        # TODO: There's something up with these primers.
        # Why are they in the comments? Did we run out of space?
        
        methodsText += 'FT                  /method_comments="98021 GTCCAGGCTGGTGTCTGG I3 1432-1449\n' 
        methodsText += 'FT                  01026seq GGGGAGAAGCAASGGGC I1 108-124  02100seq\n'
        methodsText += 'FT                  CCGCACGCACCCACCG 5UT -44--29  03026 GAGGTTCCTCTAGGACCTTAA I5\n'
        methodsText += 'FT                  2439-2459  03052 TAAGGAGGGAGAYGGGGGT I4 1847-1865  03055\n'
        methodsText += 'FT                  CTGCYGTGAKGTGGAGGAG E5 2035-2053  14256 GAATCCTCCTGGGTTTCCAG\n'
        methodsText += 'FT                  I3 1115-1134  97094seq TGTCGTCCACGTAGC E2 279-293  98070\n'
        methodsText += 'FT                  GGCCTAAACTGAAAATGAAACC I2 622-643  00029 GGTCCCAATTGTCTCCCCTC\n'
        methodsText += 'FT                  I3 1055-1074  02038seq GGCCAGCAATGATGC E5 1981-1995  03017\n'
        methodsText += 'FT                  CCTTTGCAGAAACAAAGTCAGGGT 3UT 2970-2993  03050\n'
        methodsText += 'FT                  TTAAGGTCCTAGAGGAACCTC I5 2439-2459  14019 CCAGACACCAGCCTGGAC\n'
        methodsText += 'FT                  I3 1432-1449      Exons and introns are defined as in regular\n'
        methodsText += 'FT                  HLA genes although in this allele a stopcodon is present in\n'
        methodsText += 'FT                  Exon 2."\n'
        
        # This is the "closest allele, right?"
        
        methodsText += 'FT                  /alignment="' + str(getConfigurationValue('closest_known_allele')) + '"\n'
        
        return methodsText

    def printFeatures(self):
        featureText = ''
        
        # TODO: I might double check with James Robinson about the backslashes before "number". 
        # Seems inconsistent.
        
        
        #featureText += 'FT   CDS            join(248..320,457..720,962..1237,1816..2091,2194..2310,\n'
        #featureText += 'FT                  2753..2785,2928..2975,3145..3149)\n'
        # Coding sequence is just the exons.  Print out each exon.
        # Ignoring line-breaks for now, this might create a really wide line. Ok?
        featureText += ('FT   CDS            join(') 
        for x in range(0,len(self.sequenceAnnotation.loci)):
            geneLocus = self.sequenceAnnotation.loci[x]
            if (geneLocus.exon):
                featureText += str(geneLocus.beginIndex) + '..' + str(geneLocus.endIndex)
                if not x==len(self.sequenceAnnotation.loci)-2:
                    featureText += ','
                else:
                    featureText += ')\n'
                    
                    
        exonIndex = 1
        intronIndex = 1
        
        geneHas3UTR = False
        geneHas5UTR = False
            
        for x in range(0,len(self.sequenceAnnotation.loci)):
            currentFeature = self.sequenceAnnotation.loci[x]

            # 3' UTR
            if(currentFeature.name == '3UT'):
                featureText += ('FT   3\' UTR         ' + str(currentFeature.beginIndex) + '..' + str(currentFeature.endIndex) + '\n')
                geneHas3UTR = True  
                
            # 5' UTR
            elif(currentFeature.name == '5UT'):
                featureText += ('FT   5\' UTR         ' + str(currentFeature.beginIndex) + '..' + str(currentFeature.endIndex) + '\n')
                geneHas5UTR = True   
            
            # Exon
            elif(currentFeature.exon):
                featureText += ('FT   Exon           ' + str(currentFeature.beginIndex) 
                    + '..' + str(currentFeature.endIndex) + '\n')
                featureText += ('FT                  \\number="' + str(exonIndex) + '"\n') 
                exonIndex += 1
            
            # Intron
            else:
                featureText += ('FT   Intron         ' + str(currentFeature.beginIndex) 
                    + '..' + str(currentFeature.endIndex) + '\n')   
                featureText += ('FT                  \\number="' + str(intronIndex) + '"\n')
                intronIndex += 1

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


    # Return True if our input values are all present and accomodated for.
    # If something is missing, then throw a fit and give up.
    # TODO: I should probably not raise these exceptions actually.
    # Instead, I should have the GUI Automatically open the choose options screen
    def validateInputs(self):
        
        # TODO: I'm using the self. values.  These should mostly be configuration values, load them from there instead.
        
        # TODO: This method is not being used?  Right, I should just delete this method
        # Instead of this method, maybe I should consider adding more robust error handling to the sequence generator.
        
        raise Exception ('Validate Inputs Method is being used, after all.')
        
        if (self.inputSampleID is None or len(self.inputSampleID) < 1):
            raise Exception ('Invalid Sequence ID:' + str(self.inputSampleID))
            return False
        
        elif (self.sequenceAnnotation is None):
            raise Exception ('Invalid Sequence Annotation:' + str(self.sequenceAnnotation))
            return False
        
        elif (getConfigurationValue('gene') is None or len(getConfigurationValue('gene')) < 1):
            raise Exception ('Invalid Input Gene:' + str(getConfigurationValue('gene')))
            return False
        
        elif (self.inputAllele is None or len(self.inputAllele) < 1):
            raise Exception ('Invalid Input Allele:' + str(self.inputAllele))
            return False
        
        elif (self.inputClass is None or len(self.inputClass) < 1):
            raise Exception ('Invalid Input Class:' + str(self.inputClass))
            return False
        
        else:
            return True


