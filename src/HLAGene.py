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

# The GeneLocus class specifies a locus on a Gene, 
# Either an Exon, intron, or UTR.
class GeneLocus():
    
    def __init__(self):
        self.name = ''
        self.sequence = ''
        self.exon = False
        self.beginIndex = 0
        self.endIndex = 0

    def length(self):
        return 1 + self.endIndex - self.beginIndex     

# The Gene class represents an entire HLA Gene, consisting of a series of loci.
class HLAGene():

    def __init__(self):
        self.fullSequence = ''
        self.loci = []

    def totalLength(self):
        return len(self.getCompleteSequence())

    # Combine the UTRs, Exons, and Introns into a contiguous sequence.
    def getCompleteSequence(self):
        sequence=''
        for i in range(0, len(self.loci)):
            sequence += self.loci[i].sequence
        return sequence

    # Combine the Exons into a contiguous sequence
    def getExonSequence(self):

        sequence=''
        for i in range(0, len(self.loci)):
            if(self.loci[i].exon):
                sequence += self.loci[i].sequence
        return sequence
 
    # This method names the UTRs, Exons, and Introns, and records their indices.
    # A HLA gene is always expected to have the pattern
    # # 5UT -> EX1 -> IN1 -> EX2 -> IN2 -> ... -> EXN -> 3UT
    # But I can't assume this so I will check.
    def annotateLoci(self):
        
        print('Annotating Gene Now')

        # An index for the nucleotide position
        lociBeginIndex = 1
        
        # Start with the names 'EX1' and 'IN1'        
        exonIndex = 1
        intronIndex = 1
        
        lociCount = len(self.loci)
        
        if(lociCount == 0):
            print('No loci to annotate.  We are done here.')
            # TODO: do I need to raise an exception or just ignore the problem?
            
        else:
            for lociIndex in range(0, len(self.loci)):

                if (self.loci[lociIndex].exon):
                    self.loci[lociIndex].name = 'EX' + str(exonIndex)
                    exonIndex += 1
                else :
                    # Is it a UTR or an Intron?
                    if (lociIndex == 0):
                        self.loci[lociIndex].name = '5UT'
                    elif (lociIndex == lociCount - 1):
                        self.loci[lociIndex].name = '3UT'
                    else:
                        self.loci[lociIndex].name = 'I' + str(intronIndex)
                        intronIndex += 1

                # Determine start and end indices of these exons.
                # Attempting to make index that looks like:
                #5UT:	1-65
                #EX1:	66-137
                # I1:	138-267
                self.loci[lociIndex].beginIndex = lociBeginIndex
                lociBeginIndex += len(self.loci[lociIndex].sequence)
                self.loci[lociIndex].endIndex = lociBeginIndex - 1
            
            
    # Print a summary of the inputted sequence to console.  
    def printGeneSummary(self):
        print('\nPrinting Gene Summary')
        for x in range(0, len(self.loci)):
            currentLocus = self.loci[x]
            print(currentLocus.name + ":\t"
            + str(currentLocus.beginIndex) + '-' + str(currentLocus.endIndex)
            + '\n' + currentLocus.sequence
            )
        print('')
 
