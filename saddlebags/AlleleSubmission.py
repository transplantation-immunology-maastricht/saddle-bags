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

# The GeneLocus class specifies a locus on a Gene, Either an Exon, intron, or UTR.
# An HlaGene object describes an HLA sequence.
# An AlleleSubmission object contains much of the submission-specific metadata necessary to store the objects.
# An AlleleSubmission has an HlaGene, which has multiple GeneLocus.
# Multiple submissions are stored in a SubmissionBatch object.

class SubmissionBatch():
    def __init__(self):
        # It's dumb to reuse a variable name, should think of a better name here.
        self.submissionBatch = []
        self.enaUserName = None
        self.enaPassword = None
        self.ipdSubmitterId = None
        self.ipdSubmitterName = None
        self.ipdAltContact = None
        self.ipdSubmitterEmail = None
        self.labOfOrigin = None
        self.labContact = None
        self.studyAccession = None
        self.chooseStudy = "2" # 2 = new study. 1 = existing study, use the studyaccession number. Study=Project
        self.studyId = None
        self.studyShortTitle = None
        self.studyAbstract = None


class AlleleSubmission():

    def __init__(self):
        self.submittedGene=HlaGene()
        self.localAlleleName = None
        self.closestAlleleWrittenDescription = None
        self.ipdSubmissionIdentifier = None
        self.ipdSubmissionVersion = None
        self.enaAccessionIdentifier = None
        # TODO: i think this column is intended for use identifying cell line names, if we are submitting the HLA types of cell lines.
        self.cellId = None
        self.ethnicOrigin = None
        self.sex = None
        self.consanguineous = None
        self.homozygous = None
        # Necessary = A,B, DRB1. The rest are extra, and they help James trust the submitted sequence.
        # I store the typed alleles as a dictionary. Key is the Locus (HLA-A) and the value is a String with the alleles, separated by a comma (02:01,03:01:14)
        self.typedAlleles = {}
        self.materialAvailability = None
        self.cellBank = None
        self.primarySequencingMethodology = None
        self.secondarySequencingMethodology = None
        self.primerType = None
        self.primers = None
        self.sequencedInIsolation = None
        self.sequencingDirection = None
        self.numOfReactions = None
        self.methodComments = None
        self.citations = None
        self.enaSubmissionText = None
        self.ipdSubmissionText = None
        self.isPseudoGene = False # It's considered a pseudogene if length of the coding sequence is not a multiple of 3.

class GeneFeature():
    
    def __init__(self):
        self.name = None
        self.sequence = None
        self.exon = False
        self.beginIndex = 0
        self.endIndex = 0

    def length(self):
        return 1 + self.endIndex - self.beginIndex     

# The Gene class represents an entire HLA Gene, consisting of a series of loci.
class HlaGene():

    def __init__(self):
        self.fullSequence = None
        self.features = []
        self.geneLocus = None
        # TODO: I don't think I'm generating the class anywhere. I'm putting a ? mark here so I can spot mistakes. Just make it Blank or None or something
        self.hlaClass = None

    def totalLength(self):
        #print('Calculating the total length. It is:' + str(len(self.getCompleteSequence())))
        #print('I have this many features: ' + str(len(self.features)))
        return len(self.getCompleteSequence())

    # Combine the UTRs, Exons, and Introns into a contiguous sequence.
    def getCompleteSequence(self):
        sequence=''
        for i in range(0, len(self.features)):
            currentFeature = self.features[i]
            if(currentFeature.exon):
                sequence += currentFeature.sequence.upper()
            else:
                sequence += currentFeature.sequence.lower()
        return sequence

    # Combine the Exons into a contiguous sequence
    def getExonSequence(self):

        sequence=''
        for i in range(0, len(self.features)):
            if(self.features[i].exon):
                sequence += self.features[i].sequence
        return sequence
 
    # This method names the UTRs, Exons, and Introns, and records their indices.
    # A HLA gene is always expected to have the pattern
    # # 5UT -> EX1 -> IN1 -> EX2 -> IN2 -> ... -> EXN -> 3UT
    # TODO This may not be a safe assumption when we are dealing with Class II submissions.
    def annotateFeatures(self):
        
        #print('Annotating Gene Now')

        # An index for the nucleotide position
        lociBeginIndex = 1
        
        # Start with the names 'EX1' and 'IN1'        
        exonIndex = 1
        intronIndex = 1
        
        lociCount = len(self.features)
        
        if(lociCount == 0):
            print('No loci to annotate.  We are done here.')
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
            
            
    # Print a summary of the inputted sequence to console.
    # Todo: Log this info instead of printing it. I can't because of a weird circular dependency. Look into that.
    def printGeneSummary(self):
        #print('\nPrinting Gene Summary')
        for x in range(0, len(self.features)):
            currentLocus = self.features[x]
            print(currentLocus.name + ":\t"
            + str(currentLocus.beginIndex) + '-' + str(currentLocus.endIndex)
            + '\n' + currentLocus.sequence
            )
        print('')
 


