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

from saddlebags.HlaSequence import HlaSequence

import logging

class SubmissionBatch():
    def __init__(self, includeInitialSubmission):

        if(includeInitialSubmission):
            # Starting with a single empty submission in the batch.
            self.submissionBatch = [AlleleSubmission()]
        else:
            # Starting with an empty batch
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
        self.submittedAllele=HlaSequence()
        self.localAlleleName = None
        self.closestAlleleWrittenDescription = None
        self.ipdSubmissionIdentifier = None
        self.ipdSubmissionVersion = None
        self.enaAccessionIdentifier = None
        # TODO: i think this column is intended for use identifying cell line names, if we are submitting the HLA types of cell lines. I'm just using it as a sample ID or cellnum. Lets see where that breaks.
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
        self.isPseudoGene = False # A null allele uses pseudogene if length of the coding sequence is not a multiple of 3.

