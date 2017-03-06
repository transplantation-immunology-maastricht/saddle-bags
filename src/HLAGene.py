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

# The GeneLocus class specifies a locus on a Gene, 
# Either an Exon, intron, or UTR.
class GeneLocus(): 
  
    name = ''
    sequence = ''
    exon = False
    beginIndex = 0
    endIndex = 0

    def length(self):
        return 1 + self.endIndex - self.beginIndex     

# The Gene class represents an entire HLA Gene, consisting of a series of loci.
class HLAGene():

    fullSequence = ''
    loci = []

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
        for i in range(1, len(self.loci)-1):
            if(self.loci[i].exon):
                sequence += self.loci[i].sequence
        return sequence
 
    # This method names the UTRs, Exons, and Introns, and records their indices.
    # A HLA gene is always expected to have the pattern
    # # 5UT -> EX1 -> IN1 -> EX2 -> IN2 -> ... -> EXN -> 3UT
    def annotateLoci(self):
        
        print('Annotating Gene Now')

        lociBeginIndex = 1
        if(len(self.loci) > 2):
            for x in range(0, len(self.loci)):

                # Determine the name of this loci.  
                # 5UT -> EX1 -> IN1 -> EX2 -> IN2 -> ... -> EXN -> 3UT
                if(x==0):
                    self.loci[x].name = '5UT'
                elif(x==len(self.loci)-1):
                    self.loci[x].name = '3UT'
                elif(x%2 == 1):
                    self.loci[x].name = 'EX' + str(x/2 + 1)
                else:
                    self.loci[x].name = 'I' + str(x/2)

                # Determine start and end indices of these exons.
                # Attempting to make index that looks like:
                #5UT:	1-65
                #EX1:	66-137
                # I1:	138-267
                self.loci[x].beginIndex = lociBeginIndex
                lociBeginIndex += len(self.loci[x].sequence)
                self.loci[x].endIndex = lociBeginIndex - 1
                

        else:
            print('I expected at least three loci in order to annotate them.  Please double check your input file.')

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
 
