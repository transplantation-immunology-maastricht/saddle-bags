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

from sys import exc_info
from io import StringIO, BytesIO
from urllib.parse import urlencode
from pycurl import Curl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from json import loads

from saddlebags.AlleleSubCommon import showInfoBox

import logging

def translateSequence(submission):
    # This is a short wrapper method to use biopython's translation method.
    # Most of this code is just checking for things that went wrong
    # TODO: I'll need to fix the imports for these biopython methods. Don't worry, it'll break when the time comes.
    # TODO: This method should be a class method of AlleleSubmission. Move it there.
    inputSequence = submission.submittedAllele.getExonSequence()
    proteinSequence = ''
    alleleLocalName = submission.localAlleleName

    try:
        # Do nothing if the input sequence is blank.
        if (len(inputSequence) > 0):

            coding_dna = Seq(inputSequence, generic_dna)
            proteinSequence = str(coding_dna.translate())
            logging.debug('Translating allele:' + alleleLocalName)
            logging.debug('Exon Sequence before translation:' + coding_dna)
            logging.debug('Translated Protein:' + proteinSequence)

            # Perform Sanity Checks.
            # Stop codon *should* be at the end of the protein.
            # Here we seek out the first instance of a stop codon,
            # and remove the peptides afterwards.
            # because that's what happens in real life.
            stopCodonLocation = proteinSequence.find('*')

            # If no stop codon was found
            if (stopCodonLocation == -1):
                submission.isPseudoGene = True
                # assignConfigurationValue('is_pseudo_gene','1')
                logging.info('No Stop Codon found. This is a "pseudo-gene".')
                # If multiple of three (correct codon length)
                if (len(coding_dna) % 3 == 0):
                    showInfoBox('No Stop Codon Found',
                                        'The translated protein does not contain a stop codon.\n' +
                                        'This is indicated by a /pseudo flag in the sequence submission.'
                                        )

                # Wrong Codon Length
                else:
                    showInfoBox('No Stop Codon Found',
                                        'The translated protein does not contain a stop codon.\n' +
                                        'The coding nucleotide sequence length (' + str(
                                            len(coding_dna)) + ') is not a multiple of 3.\n' +
                                        'This is indicated by a /pseudo flag in the sequence submission.')

            # If Stop Codon is in the end of the protein (This is expected and correct)
            elif (stopCodonLocation == len(proteinSequence) - 1):
                submission.isPseudoGene = False
                # assignConfigurationValue('is_pseudo_gene','0')

                # If multiple of three (correct codon length)
                if (len(coding_dna) % 3 == 0):
                    # Everything is fine in this case.  Trim off the stop codon
                    logging.info('The stop codon is in the correct position. This is not a "pseudo-gene".')
                    proteinSequence = proteinSequence[0:stopCodonLocation]
                    pass
                    # Wrong Codon Length
                else:
                    logging.info(
                        'The stop codon is in the correct position, but there are extra nucleotides. This is not a "pseudo-gene".')
                    showInfoBox('Extra Nucleotides After the Stop Codon',
                                        'The stop codon is at the correct position in the protein, but ' +
                                        'The coding nucleotide sequence length (' + str(
                                            len(coding_dna)) + ') is not a multiple of 3.\n\n' +
                                        'Please double check your sequence.')
                    proteinSequence = proteinSequence[0:stopCodonLocation]

            # Else Stop Codon is premature (before the end of the protein)
            else:
                logging.info('A premature stop codon was found. This is a "pseudo-gene".')
                submission.isPseudoGene = True
                # assignConfigurationValue('is_pseudo_gene','1')

                # If multiple of three (correct codon length)
                if (len(coding_dna) % 3 == 0):
                    showInfoBox('Premature Stop Codon Detected',
                                        'Premature stop codon found:\nProtein Position (' +
                                        str(stopCodonLocation + 1) + '/' +
                                        str(len(proteinSequence)) + ')\n\n' +
                                        'This is indicated by a /pseudo flag in the sequence submission.\n' +
                                        'Double check your protein sequence,\n' +
                                        'this might indicate a missense mutation.\n\n' +
                                        'Translated Protein:\n' + proteinSequence +
                                        '\n\nProtein in ENA Submission:\n' + proteinSequence[0:stopCodonLocation] +
                                        '\n'
                                        )
                    proteinSequence = proteinSequence[0:stopCodonLocation]


                # Wrong Codon Length
                else:
                    showInfoBox('Premature Stop Codon Detected',
                                        'Premature stop codon found:\nProtein Position (' +
                                        str(stopCodonLocation + 1) + '/' +
                                        str(len(proteinSequence)) + ')\n\n' +
                                        'This is indicated by a /pseudo flag in the sequence submission.\n' +
                                        'Nucleotide count is not a multiple of 3,\n' +
                                        'Double check your protein sequence,\n' +
                                        'this might indicate a missense mutation.\n\n' +
                                        'Translated Protein:\n' + proteinSequence +
                                        '\n\nProtein in ENA Submission:\n' + proteinSequence[0:stopCodonLocation] +
                                        '\n'
                                        )
                    proteinSequence = proteinSequence[0:stopCodonLocation]
        else:
            logging.warning('Translating a nucleotide sequence of length 0. Done. That was easy.')
            pass

        return proteinSequence

    except Exception:
        logging.error('Problem when translating protein:')
        logging.error(str(exc_info()))
        showInfoBox('Protein Translation Error','I could not translate your protein:\n' + str(exc_info()))

        raise

def collectAndValidateRoughSequence(roughNucleotideSequence):
    # This method can clean up a sequence input.
    # Should work for fasta and fastq inputs. XML in the future???
    try:
        cleanedSequence = None

        # Is this sequence in Fasta format?
        try:
            logging.debug('Checking if sequence is fasta format.')

            fileHandleObject = StringIO(roughNucleotideSequence)
            fastaSeqList = list(SeqIO.parse(fileHandleObject, 'fasta'))
            logging.debug('The length of the fasta seq list is:' + str(len(fastaSeqList)))
            if (len(fastaSeqList) == 1):
                cleanedSequence = cleanSequence(str(fastaSeqList[0].seq))
                logging.debug('The input sequence is in .fasta format.')
            else:
                logging.debug('This sequence is not in .fasta format.')
        except Exception:
            logging.error('Exception when determining file type: ' + str(exc_info()))

        # Is this sequence in Fastq format?
        try:
            logging.debug('Checking if sequence is fastq format.')
            fileHandleObject = StringIO(roughNucleotideSequence)
            fastqSeqList = list(SeqIO.parse(fileHandleObject, 'fastq'))
            logging.debug('The length of the fasta seq list is:' + str(len(fastaSeqList)))
            if (len(fastqSeqList) == 1):
                cleanedSequence = cleanSequence(str(fastqSeqList[0].seq))
                logging.debug('The input sequence is in .fastq format.')
            else:
                logging.debug('This sequence is not in .fastq format.')
        except Exception:
            logging.error('Exception when determining file type: ' + str(exc_info()))

        # TODO: If this file is xml what should we do?  Just give up i suppose.
        # TODO: I could warn the user that XML isn't supported yet...
        # We want to accept HML.  But there are too many xml formats.
        # Yeah I dunno about HML, we will not implement that right now.

        # If we haven't found an annotated sequence yet, this is not fasta or fastq.
        if (cleanedSequence is None):
            cleanedSequence = cleanSequence(roughNucleotideSequence)

        # Are we using any nonstandard / ambiguous nucleotides?
        for nucleotideCharacter in cleanedSequence:
            if (nucleotideCharacter not in ('A', 'G', 'C', 'T', 'a', 'g', 'c', 't')):
                showInfoBox('Nonstandard Nucleotide'
                                     , 'I found a non-standard\n'
                                     + 'character in your nucleotide\n'
                                     + 'sequence: '
                                     + str(nucleotideCharacter) + '\n'
                                     + 'You should use standard nucleotide\n'
                                     + 'characters in your submission.\n'
                                     + 'I will attempt to continue.')
                break

        return cleanedSequence

    except Exception:
        # except Exception, e:
        showInfoBox('Error Reading Input Sequence.'
                             , str(exc_info()))
        raise

def cleanSequence(inputSequenceText):
    # Trim out any spaces, tabs, newlines.
    if (inputSequenceText is None):
        logging.warning('Attempting to clean an input sequence that is None! I will return None.')
        return None
    else:
        logging.debug('Cleaning Input Sequence.')

    cleanedSequence = inputSequenceText.replace(' ', '').replace('\n', '').replace('\t', '').replace('\r', '')
    return cleanedSequence

class GeneFeature():
    
    def __init__(self):
        self.name = None
        self.sequence = None
        self.exon = False
        self.beginIndex = 0
        self.endIndex = 0

    def length(self):
        return 1 + self.endIndex - self.beginIndex     


class HlaSequence():
    # The HlaSequence class represents an entire HLA alleles, consisting of a series of loci.

    def __init__(self):
        self.rawSequence = None
        self.features = []
        self.geneLocus = None
        self.hlaClass = None

    def totalLength(self):
        #logging.info('Calculating the total length. It is:' + str(len(self.getCompleteSequence())))
        #logging.info('I have this many features: ' + str(len(self.features)))
        fullSeq = self.getAnnotatedSequence(includeLineBreaks=False)
        return 0 if fullSeq is None else len(fullSeq)

    def getAnnotatedSequence(self, includeLineBreaks=True):
        # Combine the UTRs, Exons, and Introns into a contiguous sequence.
        # If there is no annotation, return the raw sequence.
        # TODO: This assumes the features are in the correct order. They probably are but what if they aren't?
        if(len(self.features) < 1):
            logging.warning('There is no stored annotation, so I will return the raw sequence as an annotated sequence.')
            return self.rawSequence
        else:
            #logging.debug('getAnnotatedSequence: I identified ' + str(len(self.features)) + ' features.')
            sequence=''
            for i in range(0, len(self.features)):
                currentFeature = self.features[i]
                if(currentFeature.exon):
                    sequence += currentFeature.sequence.upper()
                else:
                    sequence += currentFeature.sequence.lower()
                if(includeLineBreaks):
                    sequence += '\n'
            return sequence

    def getExonSequence(self):
        # Combine the Exons into a contiguous sequence
        sequence=''
        for i in range(0, len(self.features)):
            if(self.features[i].exon):
                sequence += self.features[i].sequence
        return sequence
 
    def nameAnnotatedFeatures(self):
        # This method names the UTRs, Exons, and Introns, and records their indices.
        # A HLA gene is always expected to have the pattern
        # # 5UT -> EX1 -> IN1 -> EX2 -> IN2 -> ... -> EXN -> 3UT
        # TODO This may not be a safe assumption when we are dealing with Class II submissions.

        # An index for the nucleotide position
        lociBeginIndex = 1
        
        # Start with the names 'EX1' and 'IN1'        
        exonIndex = 1
        intronIndex = 1
        
        lociCount = len(self.features)
        
        if(lociCount == 0):
            logging.warning('No loci to annotate. This sequence does not have any gene features.')
            # TODO: do I need to raise an exception or just ignore the problem?
            
        else:
            for lociIndex in range(0, len(self.features)):

                if (self.features[lociIndex].exon):
                    self.features[lociIndex].name = 'EX' + str(exonIndex)
                    exonIndex += 1
                else :
                    # Is it a UTR or an Intron?
                    if (lociIndex == 0):
                        self.features[lociIndex].name = '5UT'
                    elif (lociIndex == lociCount - 1):
                        self.features[lociIndex].name = '3UT'
                    else:
                        self.features[lociIndex].name = 'I' + str(intronIndex)
                        intronIndex += 1

                # Determine start and end indices of these exons.
                # Attempting to make index that looks like:
                #5UT:	1-65
                #EX1:	66-137
                # I1:	138-267
                self.features[lociIndex].beginIndex = lociBeginIndex
                lociBeginIndex += len(self.features[lociIndex].sequence)
                self.features[lociIndex].endIndex = lociBeginIndex - 1

    def printGeneSummary(self):
        # Print a summary of the inputted sequence to console.
        logging.info('\nPrinting Gene Summary')
        for x in range(0, len(self.features)):
            currentLocus = self.features[x]
            logging.info(currentLocus.name + ":\t"
            + str(currentLocus.beginIndex) + '-' + str(currentLocus.endIndex)
            + '\n' + currentLocus.sequence
            )
        logging.info('')
 
    #TODO: Might Need to detect if a sequence is already annotated. This is probably not the way to do it.
    # def isSequenceAlreadyAnnotated(self):
    #     # The easy case.
    #     if ('a' in inputSequenceText and
    #             'g' in inputSequenceText and
    #             'c' in inputSequenceText and
    #             't' in inputSequenceText and
    #             'A' in inputSequenceText and
    #             'G' in inputSequenceText and
    #             'C' in inputSequenceText and
    #             'T' in inputSequenceText
    #     ):
    #         return True
    #
    #     # TODO: This isn't perfect. It must have all 8 nucleotides to return true.
    #     # Circle back on this one later, should I store a variable somewhere if the sequence has been annotated?
    #     return False

    def annotateSequenceUsingService(self, rawRequestURL=None):
        # Just a wrapper method, call this when I want to annotate using the service.
        sequenceAnnotation = self.fetchAnnotationJson(rawRequestURL=rawRequestURL)
        self.identifyFeaturesFromJson(sequenceAnnotation)

    def fetchAnnotationJson(self, rawRequestURL=None):
        try:
            postData = {'sequence': self.rawSequence}

            # Using configuration here causes circular dependency. So I'll just pass it in.
            if(rawRequestURL is None):
                logging.error('You must pass a rawRequestURL to fetchAnnotationJson.')
                return
            else:
                requestURL = rawRequestURL + '?' + urlencode(postData)

            resultsIoObject = BytesIO()

            curlObject = Curl()
            curlObject.setopt(curlObject.URL, requestURL)
            curlObject.setopt(curlObject.WRITEDATA, resultsIoObject)

            curlObject.perform()
            curlObject.close()

            getBody = resultsIoObject.getvalue().decode('utf8')

            logging.debug('JSON Request Body:\n' + getBody)

            # TODO:
            # Detect error <head><title>414 Request-URI Too Large</title></head>
            # For larger DRB alleles the webserver fails.
            # Detect error if the result is not json.
            # Maybe this error detection happens in parseExons. But i maybe need to detect server errors here.
            # Simple case is an empty string.
            if(getBody is None or len(getBody)<1):
                logging.error('The JSON results were an empty string. Is there a problem with the ACT server?:' + str(requestURL))
                showInfoBox('Problem Accessing Annotation Service','The JSON results were an empty string. Is there a problem with the ACT server?')
                return None

            # If it's an html error we can respond nicely.
            if(getBody[0:5]=='<html>'):
                # TODO: this might not work if i get some other kind of html.
                errorCode = getBody[getBody.find('<title>'):getBody.find('</title>')]
                logging.error('The annotation JSON results are html, this probably indicates an issue with the annotation webserver:\n' + str(requestURL))
                showInfoBox('Problem Accessing Annotation Service', 'The annotation results are HTML, not JSON, probably an issue with the ACT webserver:\n' + str(errorCode))
                return None

            return getBody

        except Exception:
            logging.error('Exception when performing CURL:\n')
            logging.error(str(exc_info()))
            logging.error('URL:' + str(requestURL))

            raise

    def identifyFeaturesFromJson(self, sequenceAnnotationJson):
        # This method parses the Json text from the ACT service, and identifies the genomic features.
        # It performs some sanity checks and then sets the according features in this HlaGene object.

        # The json should be a String, but it is returned from the NMDP ACT API as a "Typing" object. I convert it to a String to make everyone happy.
        jsonString = str(sequenceAnnotationJson)
        if(jsonString is not None and len(jsonString) > 1):

            try:
                self.features = []
                fivePrimeSequence = ''
                threePrimeSequence = ''



                #logging.debug('THIS IS THE JSON STRING\n:' + jsonString)
                parsedJson = loads(jsonString)

                if (len(parsedJson.keys()) > 1):
                    # Loop through the recognized Features
                    if 'features' in parsedJson.keys():
                        # We found features.
                        featureList = parsedJson['features']
                        logging.info('I found this many Known Features:' + str(len(featureList)))

                        for featureDictionary in featureList:

                            term = str(featureDictionary['term'])
                            rank = str(featureDictionary['rank'])
                            sequence = str(featureDictionary['sequence'])

                            logging.debug('Known Feature'
                                          + ':' + term
                                          + ':' + rank
                                          + ':' + sequence)

                            currentFeature = GeneFeature()

                            if (term == 'five_prime_UTR'):
                                currentFeature.sequence = sequence.lower()
                                currentFeature.exon = False
                                fivePrimeSequence = currentFeature.sequence
                            elif (term == 'three_prime_UTR'):
                                currentFeature.sequence = sequence.lower()
                                currentFeature.exon = False
                                threePrimeSequence = currentFeature.sequence
                            elif (term == 'exon'):
                                currentFeature.sequence = sequence.upper()
                                currentFeature.exon = True
                            elif (term == 'intron'):
                                currentFeature.sequence = sequence.lower()
                                currentFeature.exon = False
                            else:
                                raise Exception('Unknown Feature Term, expected exon or intron:' + term)

                            self.features.append(currentFeature)

                    else:
                        raise Exception('Unable to identify any HLA exon features, unable to annotate sequence.')

                    if (len(fivePrimeSequence) < 1):
                        logging.warning('I cannot find a five prime UTR.')
                        logging.info('Rough Sequence:\n' + cleanSequence(self.rawSequence).upper())
                        logging.info('Annotated Sequence:\n' + cleanSequence(annotatedSequence).upper())
                        raise Exception('GFE service did not find a 5\' UTR sequence. You will need to annotate the genomic features manually.')

                    elif cleanSequence(fivePrimeSequence).upper() in cleanSequence(self.rawSequence).upper():
                        # What if the reported 5' UTR is less than what is returned by GFE?
                        # TODO: I don't know if this code is working. Hard to debug.
                        beginIndex = cleanSequence(self.rawSequence).upper().find(cleanSequence(fivePrimeSequence).upper())
                        endIndex = beginIndex + len(fivePrimeSequence)
                        logging.info('GFE sequence exists in rough sequence, at index: (' + str(beginIndex) + ':' + str(endIndex) + ')')
                        logging.info('previous fivePrime Sequence=\n' + fivePrimeSequence)
                        fivePrimeSequence = cleanSequence(self.rawSequence)[0:endIndex].lower()
                        logging.info('new fivePrime Sequence=\n' + fivePrimeSequence)

                    self.nameAnnotatedFeatures()

                    # Final check: Do the annotated sequence and rough sequence match?
                    if (cleanSequence(self.getAnnotatedSequence(includeLineBreaks=False)).upper() == cleanSequence(self.rawSequence).upper()):
                        logging.info('Successful annotation.')
                        pass

                    else:
                        logging.error('Rough Sequence:\n' + cleanSequence(self.rawSequence).upper())
                        logging.error('Annotated Sequence:\n' + cleanSequence(self.getAnnotatedSequence(includeLineBreaks=True)).upper())
                        raise Exception('Annotated sequence and rough sequence do not match. Something went wrong in Annotation.')

                else:
                    raise Exception('No keys found in the JSON Dictionary, unable to annotate sequence.')

            except Exception:
                logging.error(str((exc_info())))
                showInfoBox('Exon Parsing Error', 'I had trouble annotating your sequence:\n' + str(exc_info()) + '. You will have to annotate manually.')
        else:
            logging.error('JSON Parse is empty.')

    def identifyFeaturesFromFormattedSequence(self):
        # The input file should be a string of nucleotides, with capital letters to identify exons and introns.
        # Annotations are expected and read in this format:
        # fiveprimeutrEXONONEintrononeEXONTWOintrontwoEXONTHREEthreeprimeutr
        # agctagctagctAGCTAGCtagctagctAGCTAGCtagctagctAGCTAGCTAgctagctagctag
        # All spaces, line feeds, and tabs are removed and ignored.
        # Return an HlaGene object

        inputSequenceText = self.rawSequence
        self.features = []

        if (inputSequenceText is None):
            logging.warning('Attempting to Identify Genomic Features on an input sequence that is None.')
            inputSequenceText = ''

        logging.debug('Identifying Genomic Features.')
        # TODO: I should accept a Fasta Input. Think i did that already. think that's fine.
        # TODO: I probably need to change the call to cleanSequence to collectAndValidateRoughSequence
        # Disregard the header line completely. Is there still sequence?

        cleanedInputText = cleanSequence(inputSequenceText)

        # Capitalize, so I can store a copy of the full unannotated sequence.
        unannotatedGene = cleanedInputText.upper()
        self.rawSequence = unannotatedGene
        logging.info('Total Unannotated Sequence Length = ' + str(len(unannotatedGene)))

        # Loop through the cleaned and annotated input sequence,
        # capitals and lowercase letters to determine exon start and end
        if (len(cleanedInputText) > 0):

            # Is the first feature an exon or an intron?
            # If we begin in an Exon
            if (cleanedInputText[0] in ('A', 'G', 'C', 'T')):
                insideAnExon = True
            # If we begin in an Intron/UTR
            elif (cleanedInputText[0] in ('a', 'g', 'c', 't')):
                insideAnExon = False
            else:
                # Nonstandard nucleotide? I should start panicking.
                # raise Exception('Nonstandard Nucleotide, not sure how to handle it')
                logging.error('Nonstandard Nucleotide at the beginning of the sequence, not sure how to handle it', 'ERROR')
                insideAnExon = False

            locusBeginPosition = 0
            for x in range(0, len(cleanedInputText)):
                currentChar = cleanedInputText[x]

                # Is this a standard nucleotide character?
                if (currentChar.upper() in ('A', 'G', 'C', 'T')):

                    if (currentChar.isupper()):
                        if (insideAnExon):
                            # We're STILL in an exon.  In this case, I should just do nothing and continue.
                            pass
                        else:
                            # In this case, we're just starting an EXON.
                            # Store the last Intron in the list.
                            currentIntron = GeneFeature()
                            currentIntron.sequence = cleanedInputText[locusBeginPosition:x].upper()
                            currentIntron.exon = False
                            self.features.append(currentIntron)
                            insideAnExon = True
                            locusBeginPosition = x
                            pass

                    else:
                        if not (insideAnExon):
                            # We're STILL in an intron.  Continue.
                            pass
                        else:
                            # Starting a new Intron.
                            # Store an Exon in the list.
                            currentExon = GeneFeature()
                            currentExon.sequence = cleanedInputText[locusBeginPosition:x].upper()
                            currentExon.exon = True
                            self.features.append(currentExon)
                            insideAnExon = False
                            locusBeginPosition = x
                            pass
                else:
                    logging.warning('Nonstandard nucleotide detected at position ' + str(x) + ' : ' + currentChar
                                    + '.  If this is a wildcard character, you might be ok.')

            # We've reached the end of the loop and we still need to store the last feature.
            # Should be a 3' UTR, but I can't be sure, people like to put in weird sequences.
            currentFeature = GeneFeature()
            currentFeature.sequence = cleanedInputText[locusBeginPosition:len(cleanedInputText)].upper()
            currentFeature.exon = insideAnExon
            self.features.append(currentFeature)

            # Annotate the features (name them) and print the results of the read file.
            self.nameAnnotatedFeatures()
            # resultGeneLoci.printGeneSummary()

        # If the sequence is empty
        else:
            logging.warning('Empty sequence during gene annotation, I don\'t have anything to do.')

        # self.sequenceAnnotation = resultGeneLoci

