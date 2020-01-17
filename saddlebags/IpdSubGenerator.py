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

from datetime import datetime

from saddlebags.HlaSequence import HlaSequence
from saddlebags.AlleleSubmission import  AlleleSubmission, SubmissionBatch

#from saddlebags.AcademicCitation import AcademicCitation
# TODO: I removed AcademicCitation because I'm pretty sure we don't actually need that in the submission. James and Dom agree this isn't necessary.


import logging

# The IpdSubGenerator class contains logic to generate an IPD-IMGT/HLA allele submission flatfile
class IpdSubGenerator():
    
    def __init__(self):
        # Instead of storing an HlaGene, im now storing a ALleleSubmission object.
        # A ALleleSubmission object has it's own HLA gene, so this is better.
        # I think I need to also store the submission batch because it needs that info too. This is a bit redundant, but oh well.
        #self.sequenceAnnotation = HlaGene()
        self.submission = AlleleSubmission()
        self.submissionBatch = SubmissionBatch()
    
    # Create the text submission based on the IPD format.
    def buildIpdSubmission(self):
        
        documentBuffer = ''

        totalLength = self.submission.submittedAllele.totalLength()
        logging.info('total calculated length = ' + str(totalLength))
        
        if(totalLength > 0):

            # These are the main sections of the IPD-IMGT/HLA submission.
            documentBuffer += self.printHeader()
            documentBuffer += self.printCitations()
            documentBuffer += self.printKeyHeader()
            documentBuffer += self.printSubmitter()
            documentBuffer += self.printSource()
            documentBuffer += self.printMethods()
            documentBuffer += self.printFeatures()
            documentBuffer += self.printSequence()

            # Print entry terminator.  The last line of an ENA entry.
            documentBuffer += ('//\n')
            
        else:
            # TODO: Remove the messagebox stuff from this class. Put it in the GUI.

            #messagebox.showinfo('No HLA Sequence Found',
            #    'The HLA sequence is empty.\nPlease fill in an annotated HLA sequence\nbefore generating the submission.' )
            raise Exception('The HLA sequence is empty. Please fill in an annotated HLA sequence before generating the submission.')
            
            pass
        

        return documentBuffer


    def printHeader(self):
        
        headerText = ''
        
        # TODO: just use the ENA identifier, or something else.
        ipdIdentifier = self.submission.ipdSubmissionIdentifier
        ipdIdentifierWithVersion = str(self.submission.ipdSubmissionIdentifier + '.' + self.submission.ipdSubmissionVersion)
        #currentSubmissionDate = '{:%d/%m/%Y}'.format(datetime.now())   # I don't need this?


        # TODO: I'm removing the IMGT/HLA identifier from the input file, because I think I don't need it.
        # Check my notes, i think I was going to just provide an ENA identifier. IMGT/HLA can give us an "HS" identifier, I shouldn't assign those.
        # It's assigned by IMGT/HLA?
        headerText += 'ID   ' + str(ipdIdentifier) + '; Sequence Submission; Confidential; ' + str(self.submission.submittedAllele.totalLength()) + ' BP.\n'
        headerText += 'XX\n'
        headerText += 'AC   ' + str(ipdIdentifier) + ';\n'
        headerText += 'XX\n'
        headerText += 'SV   ' + str(ipdIdentifierWithVersion) + '\n'
        headerText += 'XX\n'

        # TODO: The DT fields refer to versions of the IMGT/HLA  database.
        # The manual.md on github describes it better, but i think these are DB versions assigned by IMGT/HLA.
        # Check our previous submissions and check with James. But I think we can skip this section...Not sure.
        #headerText += 'DT   ' + str(currentSubmissionDate) + ' (Submitted)\n'
        #headerText += 'DT   ' + str(self.submission.ena_release_date')) + ' (Release)\n'
        #headerText += 'XX\n'

        # TODO: I'm using the local allele name that is assigned by the user.
        # Maybe this allele name should be based on the closest allele.
        # Do I want the allele name, or should I generate a new one based on the closest allele?
        headerText += 'DE   ' + str(self.submission.localAlleleName) + '\n'
        headerText += 'XX\n'
        headerText += 'KW   HLA WEB SUBMISSION;\n'
        headerText += 'XX\n'
        
        # The new allele description is split into multiple lines. I should add a new 'CC' line for each part of the description.
        rawDescription = str(self.submission.closestAlleleWrittenDescription)
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
        # TODO: Our submission says GENBANK, but we're using ENA Numbers. This works fine but maybe I should write ENA instead?
        headerText += 'DR   GENBANK; ' + str(self.submission.enaAccessionIdentifier) + '.\n'

        return headerText

    def printCitations(self):
        citationText = ''
        citationText += 'XX\n'

        citationTextList = self.submission.citations
        #
        #str(self.submission.ipd_submission_identifier'))

        citationText += 'RN   [1]\n'
        citationText += 'RC   Unpublished.\n'
        citationText += 'XX\n'

        # TODO: Fix the citation input. The citation input is a free form input String
        # Actually, don't. Based on discussions with james and Dom I can probably assume most submissions have no citations.
        # They are able to add any relevent citations later.
        # Pretty sure I don't need to parse the @@@ characters. Those are in the configuration file.
        # I planned to just split by newline character, but maybe James has a suggestion.

        """
        if (citationTextList is None or len(citationTextList) < 1):
            citationText += 'RN   [1]\n'
            citationText += 'RC   Unpublished.\n'
            citationText += 'XX\n'
        else:

            #citationTextList = citationText.split(';')
            #citationTextEnumeration = enumerate(citationText)

            # Each citation is 5 tokens, separated by semicolons.
            # Is this a multiple of 5?
            if(len(citationTextList) % 5 == 0):

                # A weird for loop to parse the info.
                citationCount = int(len(citationTextList) / 5)

                for x in range(0, citationCount):
                    citationObject = AcademicCitation()

                    # I don't like how i'm indexing here. Moduluses are confusing.
                    # Also, super-cheap way of escaping my semicolons here
                    citationObject.referencePosition = citationTextList[0 + 5 * x].replace('@@@', ';')
                    citationObject.referenceCrossReference = citationTextList[1 + 5 * x].replace('@@@', ';')
                    citationObject.referenceAuthors = citationTextList[2 + 5 * x].replace('@@@', ';')
                    citationObject.referenceTitle = citationTextList[3 + 5 * x].replace('@@@', ';')
                    citationObject.referenceLocation = citationTextList[4 + 5 * x].replace('@@@', ';')

                    citationText += citationObject.printCitation(x + 1)


            else:
                raise Exception('Invalid Citation. The Citation information is malformed in the configuration file. Please double check your citation information.')

        """
        return citationText

    def printKeyHeader(self):
        keyHeaderText = ''
        keyHeaderText += 'FH   Key            Location/Qualifier\n'
        keyHeaderText += 'FH\n'

        return keyHeaderText


    def printSubmitter(self):
        submitterText = ''
        
        # TODO: I don't know any of this data.  Should it be int he form?
        # Maybe I just need the submitter ID, and i can or can not get the rest?
        # I should be able to calculate the indices, at least.

        submitterText += 'FT   submittor      1..' + str(self.submission.submittedAllele.totalLength()) + '\n'
        submitterText += 'FT                  /ID="' + str(self.submissionBatch.ipdSubmitterId) + '"\n'
        submitterText += 'FT                  /name="' + str(self.submissionBatch.ipdSubmitterName) + '"\n'
        submitterText += 'FT                  /alt_contact="' + str(self.submissionBatch.ipdAltContact) + '"\n'
        submitterText += 'FT                  /email="' + str(self.submissionBatch.ipdSubmitterEmail) + '"\n'
        
        return submitterText

    def printSource(self):
        sourceText = ''
        
        # TODO: Submitting Laboratory Information. Can this be fetched from imgt, or use some sort of lookup

        sourceText += 'FT   source         1..' + str(self.submission.submittedAllele.totalLength()) + '\n'
        sourceText += 'FT                  /cell_id="' + str(self.submission.cellId) + '"\n'
        sourceText += 'FT                  /ethnic_origin="' + str(self.submission.ethnicOrigin) + '"\n'
        sourceText += 'FT                  /sex="' + str(self.submission.sex) + '"\n'
        sourceText += 'FT                  /consanguineous="' + str(self.submission.consanguineous) + '"\n'
        sourceText += 'FT                  /homozygous="' + str(self.submission.homozygous) + '"\n'
        sourceText += 'FT                  /lab_of_origin="' + str(self.submissionBatch.labOfOrigin) + '"\n'
        sourceText += 'FT                  /lab_contact="' + str(self.submissionBatch.labContact) + '"\n'
        
        # TODO: No Material Available.  What if Material is available?
        # I think I need to add this to the form still.
        # Same story with "cell_bank"
        
        sourceText += 'FT                  /material_available="' + str(self.submission.materialAvailability) + '"\n'
        sourceText += 'FT                  /cell_bank="' + str(self.submission.cellBank) + '"\n'
        
        # TODO: James suggested that I only allow valid fully-sequenced alleles.
        # Should I validate this, or should I leave that work to ipd?
        # Yeah he's better at validating that than I am, I can't maintain a list of alleles, James can.
        # But this should be a dropdown list. I can store these in a list in the config file but why

        # The HLA alleles in the ipd submission file are stored in sections labeled as HLA-A and HLA-B.
        # That's dumb, since we could just include that in a single string. They want me to make other alleles.
        # I think I'll just store the proper allele name in saddlebags and split based on the * character.

        typedHlaAlleles = self.submission.typedAlleles
        if typedHlaAlleles is not None:
            # I see a problem with this. If the hlaAlleleList only has one element, i think I'll try to split up that string.
            # That seems bad and will crash this.

            for loci in sorted(typedHlaAlleles.keys()):
            #for index, hlaAllele in enumerate(typedHlaAlleles):

                #gene, nomenclatureFields = hlaAllele.split('*')
                #methodsText += 'FT                  /primer_' + str(index + 1) + '="' + hlaAllele + '"\n'
                sourceText += 'FT                  /' + str(loci) + '*="' + typedHlaAlleles[loci] + '"\n'

        else:
            # TODO: What do if there are no other HLA alleles?
            pass

        #sourceText += 'FT                  /HLA-A*="01"\n'
        #sourceText += 'FT                  /HLA-C*="07"\n'
        #sourceText += 'FT                  /HLA-B*="07"\n'
        #sourceText += 'FT                  /HLA-DRB1*="15:01,-"\n'
        
        return sourceText
    
    def printMethods(self):
        methodsText = ''
        
        # TODO: Get primer info from the form.  Make sure this all is correct        
        
        methodsText += 'FT   method         1..' + str(self.submission.submittedAllele.totalLength()) + '\n'
        
        # TODO: What are the options for sequencing methodology?
        # I can provide an open-text field.        
        
        methodsText += 'FT                  /primary_sequencing="' + str(self.submission.primarySequencingMethodology) + '"\n'
        methodsText += 'FT                  /secondary_sequencing="' + str(self.submission.secondarySequencingMethodology) + '"\n'
        methodsText += 'FT                  /type_of_primer="' + str(self.submission.primerType) + '"\n'
        methodsText += 'FT                  /sequenced_in_isolation="' + str(self.submission.sequencedInIsolation) + '"\n'
        
        # TODO Add these primers dynamically
        # A primer has these pieces of information
        # "ID" ("primer_1") "Sequence" "Feature" "locus/indices" 
        # locus seems to be genomic index, from the beginning of the sequence.
        # I suppose this has to be locations in the reference sequence?
        # I should store a dictionary of primers in the configuration. 
        # Errr, nodes underneath the Primer nodes.
        # They put a "tab" character between some of this data.  Why? Because Tabs, sigh.
        # Loop: For each line in Primer text. I'll need to split it by newline character. That is fine.
        # methodsText += 'FT                  /primer_1=""\n'
        # methodsText += 'FT                  /primer_2=""\n'

        primerList = self.submission.primers.split(';')


        #print('I will add the primers to the ipd submission:' + str(primerList))

        if primerList is not None:
            for index, primer in enumerate(primerList):
                methodsText += 'FT                  /primer_' + str(index + 1) + '="' + primer + '"\n'
        else:
            # TODO: What do I do if we have not provided any primers? Im not sure right now.
            pass


        methodsText += 'FT                  /no_of_reactions="' + str(self.submission.numOfReactions) + '"\n'
        methodsText += 'FT                  /sequencing_direction="' + str(self.submission.sequencingDirection) + '"\n'
        
        # TODO: There's something up with these primers.
        # Why are they in the comments? Did we run out of space?
        
        # Don't put primers in comments. Make this a loop. Multiline.
        methodsText += 'FT                  /method_comments="' + str(self.submission.methodComments) + '"\n'

        # TODO the alignment seems to be rejected by the IMGT/HLA Validator. For now I will not include this, because this information is included in the header.
        # I can set this when the GFE/ACT was performed.
        #methodsText += 'FT                  /alignment="' + str(self.submission.closest_allele_written_description')) + '"\n'
        
        return methodsText

    def printFeatures(self):
        featureText = ''
        
        # TODO: I might double check with James about the backslashes before "number".
        # Seems inconsistent.
        
        
        #featureText += 'FT   CDS            join(248..320,457..720,962..1237,1816..2091,2194..2310,\n'
        #featureText += 'FT                  2753..2785,2928..2975,3145..3149)\n'
        # Coding sequence is just the exons.  Print out each exon.
        # Ignoring line-breaks for now, this might create a really wide line. Ok?
        featureText += ('FT   CDS            join(') 
        for x in range(0, len(self.submission.submittedAllele.features)):
            geneLocus = self.submission.submittedAllele.features[x]
            if (geneLocus.exon):
                featureText += str(geneLocus.beginIndex) + '..' + str(geneLocus.endIndex)
                if not x == len(self.submission.submittedAllele.features) - 2:
                    featureText += ','
                else:
                    featureText += ')\n'
                    
                    
        exonIndex = 1
        intronIndex = 1
        
        geneHas3UTR = False
        geneHas5UTR = False
            
        for x in range(0, len(self.submission.submittedAllele.features)):
            currentFeature = self.submission.submittedAllele.features[x]

            # test, temp
            #print('this feature is named:' + currentFeature.name + '\n')
            #print('Total of ' + str(len(self.submission.submittedAllele.loci)) + ' features')

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
        # But move on with your life, it is not necessary to become upset.
        if (not geneHas3UTR and not geneHas5UTR):
            messagebox.showinfo('Missing UTRs', 
                'This sequence has no 5\' or 3\' UTR.\n\n' + 
                'Use lowercase nucleotides at the\n' + 
                'beginning and end of your DNA\n' +
                'sequence to specify the 5\' and 3\' UTRs.' )
        elif (not geneHas5UTR):
            messagebox.showinfo('Missing 5\' UTR', 
                'This sequence has no 5\' UTR.\n\n' + 
                'Use lowercase nucleotides at the\n' + 
                'beginning and end of your DNA\n' +
                'sequence to specify the 5\' and 3\' UTRs.' )            
        elif (not geneHas3UTR):
            messagebox.showinfo('Missing 3\' UTR', 
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
        
        completeSequence = self.submission.submittedAllele.getAnnotatedSequence(includeLineBreaks=False).upper()
                
        cCount = completeSequence.count('C')
        gCount = completeSequence.count('G')
        tCount = completeSequence.count('T')
        aCount = completeSequence.count('A')
        otherCount = self.submission.submittedAllele.totalLength() - (cCount + gCount + tCount + aCount)

        sequenceText += ('SQ   Sequence ' + str(self.submission.submittedAllele.totalLength()) + ' BP; '
                         + str(aCount) + ' A; ' + str(cCount) + ' C; '
                         + str(gCount) + ' G; ' + str(tCount) + ' T; '
                         + str(otherCount) + ' other;\n')

        # Here's some logic to print the sequence information in groups of 10.
        # This format is specified in the User manual specified by ENA.
        currentSeqIndex = 0

        while (currentSeqIndex < self.submission.submittedAllele.totalLength()):
            # The character code for a sequence region is two blank spaces,
            # followed by three blank spaces, for a total of 5 blanks.
            sequenceText += '     '
            sequenceRow = self.submission.submittedAllele.getAnnotatedSequence(includeLineBreaks=False).upper()[currentSeqIndex: currentSeqIndex + 60]

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


# TODO: I suppose this method should be in an IPD SubGenerator file.
def createIPDZipFile(zipFileName):
    logging.debug('Saving Zip File:' + str(zipFileName))

    # create a temp working directory in the current folder.
    zipDirectory = getSaddlebagsDirectory()
    workingDirectory = join(zipDirectory, 'submission_temp')
    #makedirs(workingDirectory)


    # Loop through my submission batch.
    submissionBatch = getConfigurationValue('submission_batch')
    if (submissionBatch == None):
        logging.warning ('There is no submission batch, I cannot create a .zip file.')

    submissionFileList = []

    submissioncount =0
    for submissionObject in submissionBatch.submissionBatch:

        print ('Generating Submission #' + str(submissioncount))
        submissioncount += 1

        # Create a submission for each entry in the batch.
        allGen = IpdSubGenerator()
        allGen.submission = submissionObject
        allGen.submissionBatch = submissionBatch
        ipdSubmission = allGen.buildIpdSubmission()

        submissionLocalFileName = str(submissionObject.localAlleleName) + '_submission.txt'
        submissionFileName = join(workingDirectory, submissionLocalFileName)
        submissionFileList.append(submissionLocalFileName)

        submissionFileObject = createOutputFile(submissionFileName)
        submissionFileObject.write(ipdSubmission)
        submissionFileObject.close()

        print ('I just saved this file: ' + submissionFileName)

    # create a zip file from the list of files.
    zipFileName = join(zipDirectory,zipFileName)

    with ZipFile(zipFileName,'w') as zip:
        # writing each file one by one
        for fileName in submissionFileList:
            fullPath = join(workingDirectory,fileName)
            zip.write(fullPath,fileName)

    # Delete the files.
    # TODO: Report results and submission identifiers. I should create a report, maybe save our submitted sequence .zip.
    # Give a report of what was submitted. It should be clear to the user what was submitted.
    # At the very least, just tell the user where the .zip file is. I can check what was submitted myself.
    try:
        successfulDeletion = True
        for fileName in submissionFileList:
            fullPath = join(workingDirectory, fileName)
            if isfile(fullPath):
                remove(fullPath)
            else:
                successfulDeletion = False
                raise Exception('ERROR when deleting file, it does not seem to exist:' + fullPath)

        # delete the working directory
        if(successfulDeletion):
            if isdir(workingDirectory):
                rmdir(workingDirectory)
            else:
                raise Exception('ERROR when deleting workign directory, it does not seem to exist:' + workingDirectory)

    except Exception as error:
        logging.error('ERROR when removing working directory and submission files:' + str(error))
        raise()



# TODO: This method is not being used?  Right, I should just delete this method, check if i need the logic anywhere.
"""
    # Return True if our input values are all present and accomodated for.
    # If something is missing, then throw a fit and give up.
    # TODO: I should probably not raise these exceptions actually.
    # Instead, I should have the GUI Automatically open the choose options screen
    def validateInputs(self):
        
        # TODO: I'm using the self. values.  These should mostly be configuration values, load them from there instead.
        
        
        
   
        # Instead of this method, maybe I should consider adding more robust error handling to the sequence generator.
        
        raise Exception ('Validate Inputs Method is being used, after all.')
        
        if (self.inputSampleID is None or len(self.inputSampleID) < 1):
            raise Exception ('Invalid Sequence ID:' + str(self.inputSampleID))
            return False
        
        elif (self.submission.submittedAllele is None):
            raise Exception ('Invalid Sequence Annotation:' + str(self.submission.submittedAllele))
            return False
        
        elif (self.submission.submittedAllele.geneLocus is None or len(self.submission.submittedAllele.geneLocus) < 1):
            raise Exception ('Invalid Input Gene:' + str(self.submission.submittedAllele.geneLocus))
            return False
        
        #elif (self.inputAllele is None or len(self.inputAllele) < 1):
        #   raise Exception ('Invalid Input Allele:' + str(self.inputAllele))
        #   return False
        
        #elif (self.inputClass is None or len(self.inputClass) < 1):
            raise Exception ('Invalid Input Class:' + str(self.inputClass))
            return False
        
        else:
            return True
            
"""


