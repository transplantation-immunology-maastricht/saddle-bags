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

# TODO: I can't use the logEvent here because it's making a circular dependancy. Think about that, it would be nice to log this information.
#from saddlebags.AlleleSubCommon import logEvent

class SubmissionBatch():
    def __init__(self):
        # It's dumb to reuse a variable name, should think of a better name here.
        self.submissionBatch = []
        self.imgtSubmitterId = ''
        self.imgtSubmitterName = ''
        self.imgtAltContact = ''
        self.imgtSubmitterEmail = ''
        self.labOfOrigin = ''
        self.labContact = ''


class AlleleSubmission():

    def __init__(self):
        self.submittedGene=HlaGene()
        self.localAlleleName = ''
        self.closestAlleleWrittenDescription = ''
        self.imgtSubmissionIdentifier = ''
        self.imgtSubmissionVersion = ''
        # TODO: i think this column is intended for use identifying cell line names, if we are submitting the HLA types of cell lines.
        self.cellId = ''
        self.ethnicOrigin = ''
        self.sex = ''
        self.consanguineous = ''
        self.homozygous = ''
        # Necessary = A,B, DRB1. The rest are extra, and they help James trust the submitted sequence.
        self.typedAlleles = ''
        self.materialAvailability = ''
        self.cellBank = ''
        self.primarySequencingMethodology = ''
        self.secondarySequencingMethodology = ''
        self.primerType = ''
        self.primers = ''
        self.sequencedInIsolation = ''
        self.noOfReactions = ''
        self.methodComments = ''
        self.citations = ''

class GeneFeature():
    
    def __init__(self):
        self.name = ''
        self.sequence = ''
        self.exon = False
        self.beginIndex = 0
        self.endIndex = 0

    def length(self):
        return 1 + self.endIndex - self.beginIndex     

# The Gene class represents an entire HLA Gene, consisting of a series of loci.
class HlaGene():

    def __init__(self):
        self.fullSequence = ''
        self.features = []
        self.geneLocus = ''

    def totalLength(self):
        return len(self.getCompleteSequence())

    # Combine the UTRs, Exons, and Introns into a contiguous sequence.
    def getCompleteSequence(self):
        sequence=''
        for i in range(0, len(self.features)):
            sequence += self.features[i].sequence
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
        
        print('Annotating Gene Now')

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
    # Todo: Log this info instead of printing it. I can't because of a circular dependency. Look into that.
    def printGeneSummary(self):
        print('\nPrinting Gene Summary')
        for x in range(0, len(self.features)):
            currentLocus = self.features[x]
            print(currentLocus.name + ":\t"
            + str(currentLocus.beginIndex) + '-' + str(currentLocus.endIndex)
            + '\n' + currentLocus.sequence
            )
        print('')
 
