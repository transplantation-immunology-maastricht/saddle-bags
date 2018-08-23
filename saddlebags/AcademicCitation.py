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

# A literature citation has 5 pieces of information, all is text for flexibility.

class AcademicCitation():
    
    def __init__(self):
        self.referencePosition = ''
        self.referenceCrossReference = ''
        self.referenceAuthors = ''
        self.referenceTitle = ''
        self.referenceLocation = ''

    # format based on suggestions in the IMGT HLA Manual
    # https://github.com/ANHIG/IMGTHLA/blob/Latest/Manual.md
    def printCitation(self, citationNumber):
        citationText = ''

        citationText += 'RN   [' + str(citationNumber) + ']\n'
        citationText += 'RP   ' + self.referencePosition + '\n'

        if (self.referenceCrossReference is not None and len(self.referenceCrossReference) > 0):
            citationText += 'RX   ' + self.referenceCrossReference + '\n'

        citationText += 'RA   ' + self.referenceAuthors + '\n'
        citationText += 'RT   ' + self.referenceTitle + '\n'
        citationText += 'RL   ' + self.referenceLocation + '\n'

        citationText += 'XX\n'

        return citationText

