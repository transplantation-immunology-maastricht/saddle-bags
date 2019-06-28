# The logging module should be loaded first, because subsequent imports will supercede my own logging settings.
# I don't like the logging module.
import logging


#from hf_cypher.cypherQuery import CypherQuery
from nose.tools import assert_equal, assert_true

from saddlebags.AlleleSubCommon import assignConfigurationValue, writeConfigurationFile,\
    identifyGenomicFeatures, parseExons, fetchSequenceAlleleCallWithGFE, clearGlobalVariables,\
    initializeGlobalVariables, getConfigurationValue, loadConfigurationFile, assignIfNotExists,\
    initializeLog,cleanSequence,loadFromCSV, createIPDZipFile, parseTypedAlleleInput, showYesNoBox

from saddlebags.EnaSubGenerator import EnaSubGenerator
from saddlebags.IpdSubGenerator import IpdSubGenerator
from saddlebags.AlleleSubmission import SubmissionBatch, AlleleSubmission
from saddlebags.IpdGoogleDriveUpload import uploadZipToIpdHla
#from saddlebags.EnaSub import performFullEnaSubmission
from saddlebags.EnaSub import performBatchEnaSubmission

from os.path import join, expanduser

from json import dumps

initializeLog()

#json_file = os.path.join(os.path.dirname(__file__), '')
#file_data = open(json_file).read()
#expected = json.loads(file_data)

# To run this, activate the environment
# source /home/ben/minionvenv/bin/activate
# Browse to the project folder and run nosetests.
# cd /home/ben/Github/saddlebags
# nosetests
# to see the console output (probably useful, nosetests hides all output by default.)
# nosetests --nocapture

#These tests seem to work, commenting to isolate if there is a problem in later tests or if we just overloaded the server

def testLoadConfigAndBatchEnaSubmission():
    # TODO: Comment this test, it's just for development.
    # This is a simple test to load submission, or to perform direct submission to ENA.
    # This is not really a standalone test. I am assuming you have a decent configuration file with proper configuration values.
    # I'm also assuming you have a .csv file full of alleles to submit.

    logging.info('Loading Config file...')
    # Clear previous configuration, just to be sure.
    clearGlobalVariables()
    assert (getConfigurationValue('submission_batch') is None)
    loadConfigurationFile()
    assert (getConfigurationValue('submission_batch') is not None)
    submissionBatch = getConfigurationValue('submission_batch')
    assert submissionBatch.submissionBatch is not None
    logging.info('Clearing submissions from submission batch...')
    logging.info('Before clearing the batch, it had a length of:' + str(len(submissionBatch.submissionBatch)))
    submissionBatch.submissionBatch = []
    assert (len(submissionBatch.submissionBatch) == 0)

    logging.info('Loading sequences from .csv file...')
    csvFileLocation = '/home/ben/BensDRASequences/DRA_InputCSV.csv'
    logging.info('Loading sequence submission values from this file:' + str(csvFileLocation))
    loadFromCSV(csvFileLocation)
    submissionBatch = getConfigurationValue('submission_batch')
    assert_true(len(submissionBatch.submissionBatch) >= 1) # 4 is a bit arbitrary. I have...21 sequences in my batch.

    logging.info('Writing Config File...')
    writeConfigurationFile()

    logging.info('Starting batch submission at ENA')

    for submission in submissionBatch.submissionBatch:
        allGen = EnaSubGenerator()
        allGen.submission = submission
        allGen.submissionBatch = submissionBatch
        submission.enaSubmissionText = allGen.buildENASubmission()
        assert_true(submission.enaSubmissionText is not None)
        assert_true(len(submission.enaSubmissionText) > 3)
        # logging.info('Submission Generated:\n' + str(submission.enaSubmissionText))

    assert (getConfigurationValue('submission_batch') is not None)
    submissionBatch = getConfigurationValue('submission_batch')
    performBatchEnaSubmission(submissionBatch)

    # Go ahead and save everything.
    writeConfigurationFile()


# def testWriteConfigFile():
#     # This will generate a full configuration and write it to file.
#     # It will overwrite your previous configuration, sorry.
#     assignIfNotExists('test_submission', '1')
#     assignIfNotExists('logging', 'DEBUG')
#     assignIfNotExists('proxy', '')
#     # TODO: Is it true that we don't upload any files to FTP anymore? using the commandline submission, I don't do any FTP submission anymore. I can probably delete these ftp options.
#     assignIfNotExists('ena_ftp_upload_site_test', 'webin.ebi.ac.uk')
#     assignIfNotExists('ena_ftp_upload_site_prod', 'webin.ebi.ac.uk')
#     assignIfNotExists('ena_rest_address_test', 'https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/')
#     assignIfNotExists('ena_rest_address_prod', 'https://www.ebi.ac.uk/ena/submit/drop-box/submit/')
#     assignIfNotExists('nmdp_act_rest_address', 'http://act.b12x.org')
#
#     submissionBatch = SubmissionBatch()
#
#     # Assign some information about this batch of submissions.
#     submissionBatch.ipdSubmitterId = 'IPD_Submitter_ID'
#     submissionBatch.ipdSubmitterName = 'Ben Bioinformaticist'
#     submissionBatch.ipdAltContact = 'IPD_Alt_Contact'
#     submissionBatch.ipdSubmitterEmail = 'Ben@BioinformaticsResearch.com'
#     submissionBatch.labOfOrigin = 'Maastricht University Medical Center'
#     submissionBatch.labContact = 'Marspel JR Spilanus'
#
#     # Add a few submissions to this batch.
#     # Submission # 1
#     submission1 = AlleleSubmission()
#
#     annotatedSequence = cleanSequence("""cagaagcagaggggtcagggcgaagtcccagggccccaggcgtggctctcagggtctcaggccccgaaggcggtgtatggattggggagtcccagccttggggattccccaactccgcagtttcttttctccctctcccaacctatgtagggtccttcttcctggatactcacgacgcggacccagttctcactcccattgggtgtcgggtttccagagaagccaatcagtgtcgtcgcggtcgcggttctaaagtccgcacgcacccaccgggactcagattctccccagacgccgagg
# ATGGCCGTCATGGCGCCCCGAACCCTCGTCCTGCTACTCTCGGGGGCTCTGGCCCTGACCCAGACCTGGGCGG
# gtgagtgcggggtcgggagggaaacggcctctgtggggagaagcaacgggcccgcctggcgggggcgcaggacccgggaagccgcgccgggaggagggtcgggcgggtctcagccactcctcgtccccag
# GCTCTCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGTCCGGAGTATTGGGACGGGGAGACACGGAAAGTGAAGGCCCACTCACAGACTCACCGAGTGGACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCG
# gtgagtgaccccggcccggggcgcaggtcacgacctctcatcccccacggacgggccaggtcgcccacagtctccgggtccgagatccgccccgaagccgcgggaccccgagacccttgccccgggagaggcccaggcgcctttacccggtttcattttcagtttaggccaaaaatccccccaggttggtcggggcggggcggggctcgggggaccgggctgaccgcggggtccgggccag
# GTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGG
# gtaccaggggccacggggcgcctccctgatcgcctgtagatctcccgggctggcctcccacaaggaggggagacaattgggaccaacactagaatatcgccctccctctggtcctgagggagaggaatcctcctgggtttccagatcctgtaccagagagtgactctgaggttccgccctgctctctgacacaattaagggataaaatctctgaaggaatgacgggaagacgatccctcgaatactgatgagtggttccctttgacacacacaggcagcagccttgggcccgtgacttttcctctcaggccttgttctctgcttcacactcaatgtgtgtgggggtctgagtccagcacttctgagtccttcagcctccactcaggtcaggaccagaagtcgctgttccctcttcagggactagaattttccacggaataggagattatcccaggtgcctgtgtccaggctggtgtctgggttctgtgctcccttccccatcccaggtgtcctgtccattctcaagatagccacatgtgtgctggaggagtgtcccatgacagatgcaaaatgcctgaatgatctgactcttcctgacag
# ACGCCCCCAAAACGCATATGACTCACCACGCTGTCTCTGACCATGAAGCCACCCTGAGGTGCTGGGCCCTGAGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGACAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTTTGCCCAAGCCCCTCACCCTGAGATGGG
# gtaaggagggagacgggggtgtcatgtcttttagggaaagcaggagcctctctgacctttagcagggtcagggcccctcaccttcccctcttttcccag
# AGCCGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCTTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCTGTGATGTGGAGGAGGAAGAGCTCAG
# gtggggaaggggtgaagggtgggtctgagatttcttgtctcactgagggttccaagacccaggtagaagtgtgccctgcctcgttactgggaagcaccacccacaattatgggcctacccagcctgggccctgtgtgccagcacttactcttttgtaaagcacctgttaaaatgaaggacagatttatcaccttgattacagcggtgatgggacctgatcccagcagtcacaagtcacaggggaaggtccctgaggaccttcaggagggcggttggtccaggacccacacctgctttcttcatgtttcctgatcccgccctgggtctgcagtcacacatttctggaaacttctctgaggtccaagacttggaggttcctctaggaccttaaggccctgactcctttctggtatctcacaggacattttcttcccacag
# ATAGAAAAGGAGGGAGCTACTCTCAGGCTGCAA
# gtaagtatgaaggaggctgatgcctgaggtccttgggatattgtgtttgggagcccatgggggagctcacccaccccacaattcctcctctagccacatcttctgtgggatctgaccaggttctgtttttgttctaccccag
# GCAGTGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGTAAAG
# gtgagagcctggagggcctgatgtgtgttgggtgttgggcggaacagtggacacagctgtgctatggggtttctttccattggatgtattgagcatgcgatgggctgtttaaagtgtgacccctcactgtgacagatacgaatttgttcatgaatatttttttctatag
# TGTGA
# gacagctgccttgtgtgggactgagaggcaagagttgttcctgcccttccctttgtgacttgaagaaccctgactttgtttctgcaaaggcacctgcatgtgtctgtgttcgtgtaggcataatgtgaggaggtggggagaccaccccacccccatgtccaccatgaccctcttcccacgctgacctgtgctccctccccaatcatctttcctgttccagagaggtggggctgaggtgtctccatctctgtctcaacttcatggtgcactgagctgtaacttcttccttccctattaaaa
# """)
#
#     submission1.submittedGene = identifyGenomicFeatures(annotatedSequence)
#
#     #submission1.submittedGene.annotateFeatures() # TODO: I think this won't work, because, the sequence is not annotated (capital and lowercase).
#     submission1.submittedGene.geneLocus = 'HLA-DRB1'
#     submission1.submittedGene.hlaClass = 'II'
#     submission1.localAlleleName = 'A02_new44'
#     submission1.closestAlleleWrittenDescription = 'DRB1:11 with a SNP in exon 2'
#     submission1.ipdSubmissionIdentifier = 'HWS100234567'
#     submission1.ipdSubmissionVersion = '1'
#     submission1.enaAccessionIdentifier = 'ENAID42423'
#     submission1.cellId = 'LMS_ID_14444'
#     submission1.ethnicOrigin = 'Caucasoid - Dutch/Irish/English, Europe'
#     submission1.sex = 'M'
#     submission1.consanguineous = 'Unknown'
#     submission1.homozygous = 'Unknown'
#     submission1.typedAlleles = parseTypedAlleleInput('HLA-A*01:01:01:01;HLA-B*40:72:01;HLA-DRB1*07:75;HLA-A*03:73;HLA-B*44:87;HLA-DRB1*14:90')
#     submission1.materialAvailability = 'No Material Available'
#     submission1.cellBank = 'Not Available'
#     submission1.primarySequencingMethodology = 'Direct sequencing of PCR product from DNA (SBT)'
#     submission1.secondarySequencingMethodology = 'MinION 2D Amplicon Sequencing'
#     submission1.primerType = 'Both allele and locus specific'
#     submission1.primers = '03PID03     CCCAAAGGGTTTCCCGGGAAATTT 3UT 3015-3042;04PID04     AAAGGGTTTCCCGGGAAATTTCCC 5UT 4015-4042'
#     submission1.sequencedInIsolation = 'Yes'
#     submission1.sequencingDirection = 'Both'
#     submission1.numOfReactions = '3'
#     submission1.methodComments = 'Sanger SBT with confirmatory MinION Long Range sequencing.'
#     submission1.citations = '1-3503;PUBMED; 9349617.;Laforet M, Froelich N, Parissiadis A, Pfeiffer B, Schell A, Faller B, Woehl-Jaegle ML, Cazenave JP, Tongio MM;A nucleotide insertion in exon 4 is responsible for the absence of expression of an HLA-A*01 allele;Tissue Antigens 50:347 - 50(1997).;1-3503;PUBMED; 9349617.;Laforet M, Froelich N, Parissiadis A, Pfeiffer B, Schell A, Faller B, Woehl-Jaegle ML, Cazenave JP, Tongio MM;A nucleotide insertion in exon 4 is responsible for the absence of expression of an HLA-A*01 allele;Tissue Antigens 50:347 - 50(1997).'
#     submissionBatch.submissionBatch.append(submission1)
#
#     # Submission # 2
#     submission2 = AlleleSubmission()
#     annotatedSequence = cleanSequence("""cagaagcagaggggtcagggcgaagtcccagggccccaggcgtggctctcagggtctcaggccccgaaggcggtgtatggattggggagtcccagccttggggattccccaactccgcagtttcttttctccctctcccaacctatgtagggtccttcttcctggatactcacgacgcggacccagttctcactcccattgggtgtcgggtttccagagaagccaatcagtgtcgtcgcggtcgcggttctaaagtccgcacgcacccaccgggactcagattctccccagacgccgagg
# ATGGCCGTCATGGCGCCCCGAACCCTCGTCCTGCTACTCTCGGGGGCTCTGGCCCTGACCCAGACCTGGGCGG
# gtgagtgcggggtcgggagggaaacggcctctgtggggagaagcaacgggcccgcctggcgggggcgcaggacccgggaagccgcgccgggaggagggtcgggcgggtctcagccactcctcgtccccag
# GCTCTCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGTCCGGAGTATTGGGACGGGGAGACACGGAAAGTGAAGGCCCACTCACAGACTCACCGAGTGGACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCG
# gtgagtgaccccggcccggggcgcaggtcacgacctctcatcccccacggacgggccaggtcgcccacagtctccgggtccgagatccgccccgaagccgcgggaccccgagacccttgccccgggagaggcccaggcgcctttacccggtttcattttcagtttaggccaaaaatccccccaggttggtcggggcggggcggggctcgggggaccgggctgaccgcggggtccgggccag
# GTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGG
# gtaccaggggccacggggcgcctccctgatcgcctgtagatctcccgggctggcctcccacaaggaggggagacaattgggaccaacactagaatatcgccctccctctggtcctgagggagaggaatcctcctgggtttccagatcctgtaccagagagtgactctgaggttccgccctgctctctgacacaattaagggataaaatctctgaaggaatgacgggaagacgatccctcgaatactgatgagtggttccctttgacacacacaggcagcagccttgggcccgtgacttttcctctcaggccttgttctctgcttcacactcaatgtgtgtgggggtctgagtccagcacttctgagtccttcagcctccactcaggtcaggaccagaagtcgctgttccctcttcagggactagaattttccacggaataggagattatcccaggtgcctgtgtccaggctggtgtctgggttctgtgctcccttccccatcccaggtgtcctgtccattctcaagatagccacatgtgtgctggaggagtgtcccatgacagatgcaaaatgcctgaatgatctgactcttcctgacag
# ACGCCCCCAAAACGCATATGACTCACCACGCTGTCTCTGACCATGAAGCCACCCTGAGGTGCTGGGCCCTGAGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGACAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTTTGCCCAAGCCCCTCACCCTGAGATGGG
# gtaaggagggagacgggggtgtcatgtcttttagggaaagcaggagcctctctgacctttagcagggtcagggcccctcaccttcccctcttttcccag
# AGCCGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCTTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCTGTGATGTGGAGGAGGAAGAGCTCAG
# gtggggaaggggtgaagggtgggtctgagatttcttgtctcactgagggttccaagacccaggtagaagtgtgccctgcctcgttactgggaagcaccacccacaattatgggcctacccagcctgggccctgtgtgccagcacttactcttttgtaaagcacctgttaaaatgaaggacagatttatcaccttgattacagcggtgatgggacctgatcccagcagtcacaagtcacaggggaaggtccctgaggaccttcaggagggcggttggtccaggacccacacctgctttcttcatgtttcctgatcccgccctgggtctgcagtcacacatttctggaaacttctctgaggtccaagacttggaggttcctctaggaccttaaggccctgactcctttctggtatctcacaggacattttcttcccacag
# ATAGAAAAGGAGGGAGCTACTCTCAGGCTGCAA
# gtaagtatgaaggaggctgatgcctgaggtccttgggatattgtgtttgggagcccatgggggagctcacccaccccacaattcctcctctagccacatcttctgtgggatctgaccaggttctgtttttgttctaccccag
# GCAGTGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGTAAAG
# gtgagagcctggagggcctgatgtgtgttgggtgttgggcggaacagtggacacagctgtgctatggggtttctttccattggatgtattgagcatgcgatgggctgtttaaagtgtgacccctcactgtgacagatacgaatttgttcatgaatatttttttctatag
# TGTGA
# gacagctgccttgtgtgggactgagaggcaagagttgttcctgcccttccctttgtgacttgaagaaccctgactttgtttctgcaaaggcacctgcatgtgtctgtgttcgtgtaggcataatgtgaggaggtggggagaccaccccacccccatgtccaccatgaccctcttcccacgctgacctgtgctccctccccaatcatctttcctgttccagagaggtggggctgaggtgtctccatctctgtctcaacttcatggtgcactgagctgtaacttcttccttccctattaaaa
# """)
#
#     submission2.submittedGene = identifyGenomicFeatures(annotatedSequence)
#     submission2.submittedGene.geneLocus = 'HLA-A'
#     submission2.submittedGene.hlaClass = 'I'
#
#     print('The second submission has this many features:' + str(len(submission2.submittedGene.features)))
#     assert_true(len(submission2.submittedGene.features) == 17)
#
#
#     #submission2.submittedGene.annotateFeatures()
#     submission2.localAlleleName = 'A_02_new36'
#     submission2.closestAlleleWrittenDescription = 'This is an A02 allele and it has interesting snps like 604G->T'
#     submission2.ipdSubmissionIdentifier = 'HWS104234667'
#     submission2.ipdSubmissionVersion = '1'
#     submission2.enaAccessionIdentifier = 'ENAID12323'
#     submission2.cellId = '12345_LMS_ID'
#     submission2.ethnicOrigin = 'Caucasoid - Dutch/Irish/English, Europe'
#     submission2.sex = 'M'
#     submission2.consanguineous = 'Unknown'
#     submission2.homozygous = 'Unknown'
#     # Necessary = A,B, DRB1. The rest are extra, and they help James trust the submitted sequence.
#     submission2.typedAlleles = parseTypedAlleleInput('HLA-A*01:01:01:01;HLA-B*40:72:01;HLA-DRB1*07:75;HLA-A*03:73;HLA-B*44:87;HLA-DRB1*14:90')
#     submission2.materialAvailability = 'No Material Available'
#     submission2.cellBank = 'Not Available'
#     submission2.primarySequencingMethodology = 'Direct sequencing of PCR product from DNA (SBT)'
#     submission2.secondarySequencingMethodology = 'MinION 1D2 Amplicon Sequencing'
#     submission2.primerType = 'Both allele and locus specific'
#     submission2.primers = '03PID03     CCCAAAGGGTTTCCCGGGAAATTT 3UT 3015-3042;04PID04     AAAGGGTTTCCCGGGAAATTTCCC 5UT 4015-4042'
#     submission2.sequencedInIsolation = 'Yes'
#     submission2.sequencingDirection = 'Both'
#     submission2.numOfReactions = '3'
#     submission2.methodComments = 'Sanger SBT with confirmatory MinION Long Range sequencing.'
#     submission2.citations = '1-3503;PUBMED; 9349617.;Laforet M, Froelich N, Parissiadis A, Pfeiffer B, Schell A, Faller B, Woehl-Jaegle ML, Cazenave JP, Tongio MM;A nucleotide insertion in exon 4 is responsible for the absence of expression of an HLA-A*01 allele;Tissue Antigens 50:347 - 50(1997).;1-3503;PUBMED; 9349617.;Laforet M, Froelich N, Parissiadis A, Pfeiffer B, Schell A, Faller B, Woehl-Jaegle ML, Cazenave JP, Tongio MM;A nucleotide insertion in exon 4 is responsible for the absence of expression of an HLA-A*01 allele;Tissue Antigens 50:347 - 50(1997).'
#
#     submissionBatch.submissionBatch.append(submission2)
#
#     # Store my submission batch in the global variables.
#     assignConfigurationValue('submission_batch', submissionBatch)
#
#     print('just assigned my submission batch, it looks like this:' + str(submissionBatch) + ' and has a submitter id of: ' + str(submissionBatch.ipdSubmitterId))
#
#
#     writeConfigurationFile()
#
#
# def testReadConfigFile():
#     # Read config file and validate that all values we expect are here.
#
#     # Clear previous configuration, just to be sure.
#     # This assertion depends on the previous test. I wrote a bunch of configuration values and
#     # then stored them, there should be something in submission_batch right now.
#     assert (getConfigurationValue('submission_batch') is not None)
#
#     submissionBatch = getConfigurationValue('submission_batch')
#     #print('Before clearing variables, my submission has has a submitter id of: ' + str(submissionBatch.ipdSubmitterId))
#
#     clearGlobalVariables()
#
#
#     assert (getConfigurationValue('submission_batch') is None)
#
#     #assert (1==0)
#     # Read configuration file
#     loadConfigurationFile()
#     assert (getConfigurationValue('submission_batch') is not None)
#
#     submissionBatch = getConfigurationValue('submission_batch')
#     #print('After clearing variables, my submission has has a submitter id of: ' + str(submissionBatch.ipdSubmitterId))
#
#     # Assert that the values exist and are what we expect.
#     assert (len(getConfigurationValue('submission_batch').submissionBatch) == 2)





# # TODO: I think this one should generally be permanently commented. Uncomment for development, I don't need to submit, willy nilly.
# def testEMBLENASubmission():
#     assert (getConfigurationValue('submission_batch') is not None)
#
#     submissionBatch = getConfigurationValue('submission_batch')
#
#     performBatchEnaSubmission(submissionBatch)




# # TODO: the Google upload / zip file was working but i disabled it while i work on the ENA stuff.
#
#
# def testIPDSubmissionFlatfileWithoutAnnotating():
#     #Can I still create an submission file, even though annotation is not working?
#     print ('Test: Creating IPD Flatfile')
#
#     initializeLog()
#
#     ipdGenerator = IpdSubGenerator()
#
#
#
#     ipdGenerator.submission = getConfigurationValue('submission_batch').submissionBatch[1]
#     ipdGenerator.submissionBatch = getConfigurationValue('submission_batch')
#
#     #ipdGenerator.sequenceAnnotation = identifyGenomicFeatures(annotatedSequence)
#     print ('I gave the ipdGenerator this submission: ' + str(ipdGenerator.submission))
#     print ('the submission has this local allele name:' + str(ipdGenerator.submission.localAlleleName))
#
#     ipdSubmission = ipdGenerator.buildIpdSubmission()
#
#     print('IPD SUBMISSION:\n' + ipdSubmission)
#
#     assert_true(len(ipdSubmission) > 20)
#     assert_true(ipdSubmission is not None)
#
#     # Maybe I don't want to safe my configuration with all these test values but maybe I don't care :D
#     writeConfigurationFile()
#
#     #Print out the submission to a file, for easier validation.
#     outputFileObject = open('Test_IPD_Submission_HLAA02010101_NoAnnotation', 'w')
#     submissionText = ipdSubmission
#     outputFileObject.write(submissionText)
#
#
#
# def testInputSequenceFromCSV():
#     csvFileLocation = 'testsequences/TestInputCSV.csv'
#     print ('Loading sequence submission values from this file:' + str(csvFileLocation))
#     loadFromCSV(csvFileLocation)
#
#     # we should have at least 4 submissions in the batch now.
#
#     submissionBatch = getConfigurationValue('submission_batch')
#     assert_true(len(submissionBatch.submissionBatch) >= 4)
#
#
#
#
# def testCreateZipFile():
#     print('Creating a zip file of our submissions.')
#     # There should at least be 4 submissions in the batch right now.
#
#     zipFileLocation = 'TestSubmissionZipFile.zip'
#     createIPDZipFile(zipFileLocation)
#
#
#
# def testGoogleUpload():
#     # TODO: go through the whole pipeline here. Implement the methods to check credentials.
#     # TODO: Implement upload test helloworld, wait 5 sec, and check if it exists.
#     # TODO: Implement send email to self, wait 5 sec, check if it exists.
#     print('Uploading a zip file to the google drive:')
#     zipFileName = 'TestSubmissionZipFile.zip'
#
#     uploadZipToIpdHla(zipFileName)
#
#






# #TODO: Annotation is not working. I can try to revert to the SeqAnn code, or fix the gfe service.
# If the gfe service needs help, talk to Martin.
# def testCreateIPDSubmissionFlatfileA():
#     print ('Test: Creating IPD Flatfile')
#
#     initializeLog()
#
#     mySubmissionBatch = getConfigurationValue('submission_batch').submissionBatch
#
#     roughFeatureSequence = mySubmissionBatch[1].submittedGene.fullSequence
#
#     # roughFeatureSequence = 'GAAAGACCTGAAAGATCACGGTGCCTTCATTTCAACTGTGAGACATGATGTAATTTTCACAAATCTACAACAGTAAGATATAGTGCAACAGGACCAGATTAAGGTCTCCTGGTTTGCAACCATGTCCCCTCCATCTCCTTTACTCCTGAACACACTCACTCCTGCAAACAGTTCTCTTGTCAAGTGGGAAATGAATGCTCTTACAAGGCTCAAACTTGTGAACACATCACTGACCAGCACAGAGCTAAAATAATTGGGGCTAAAAATACCGCCCCAATTAAAGTGTTTTACATGCAACTGGTTCAAACCTTTCAAGTACTAAAAACAATCCTGTAAAGAAGGAAATTCTGTTTCAGAAGAGGACCTTCATACAGCATCTCTGACCAGCAACTGATGATGCTATTGAACTCAGATGCTGATTCGTTCTCCAACACTAGATTACCCAATCCAGGAGCAAGGAAATCAGTAACTTACTCCCTATAACTTGGAATGTGGGTGGAGGGGTTCATAGTTCTCCCTGAGTGAGACTTGCCTGCTGCTCTGGCCCCTGGTCCTGTCCTGTTCTCCAGCATGGTGTGTCTGAGGCTCCCTGGAGGCTCCTGCATGGCAGTTCTGACAGTGACACTGATGGTGCTGAGCTCCCCACTGGCTTTGGCTGGGGACACCAGACGTAAGTGCACATTGTGGGTGCTGAGCTACTATGGGGTGGGGAAAATAGGGAGTTTTGTTAACATTGTGCCCAGGCCATGTACCTTAAGAAATTGTGACGTATTTTTCAGAGATTGCCCATCTTTATCATATGGATCCCAAATAATTTCCCCACCACAAAAGGAGCTTGGCTACTTGCCCACTCCATGAGACTTGTGTAAGGGGCCTCCATACAGGTCATTTCTACTCAAATCTCCACCAATAAAACCTTTGCATCACATGTCCTCAGGGTCTTTAGAGGATTTAGAAATAAGGATGCTAAAATAAATTCCCCATACAGCACTTGCCTTTATTATGTTGACTTATGTCAGACAAAGGAGGTTTTTTCTGAAAATTTTGTGGAAGTCAAGGGAATTTAAAGGGTCTCTCCTAAACGATCCTGGGTTATGTCACCCACAGGACCTTTGGTGTTGGCCCCTCTTCCTCATATGTGAGGATGGACCCAGTGGCCTCCCCATTATCTACTTTCTTTTCTTTCTGAACTCCAATGTTTATAAAGCCTGTACCCCTGTAGTGTATGTAGGTTGTCTGACAGAAGTTATACTTAGTGCTCTTTCTTTCTTGTGGGGAAAAATCCCTGGAACTGAAGCTGAGATCTTTAGTACTTGGAGTCACCTTACAGATACAGAGCATTTATGAGGTATTCTTTGGTGCCTAAAGAACTTAAGGCATCCTCTGAAAACCTGGCCCAGGTTAGTGTTTATTATGAATCTCTTTTAACCTTTCTATACTTGTTTCTCCTACATCTCCTAAGTGCTCCAACTAGACATGACAGAAGAGATTTAACTAACGTAGTATGAGTTATATAAAATTCTATTTTTGTAAGTCAAAAATAATCAAATATCAAAAATTTAATAATGTTCAAACTATATACTCTGTGTGGGGTTACCGAGACAATGTGGACATTGTTCACATCTCATAGGGCTGAAAGTCAATGGGCAAGTCCTGGGAACTCATTGTCTTACTGGGGTCTTGTCCTAAATTTCATAGGTTCACCCATCATGCCCTCAACTTTCCTTAATTAGCCATGTCTGCTTACCTCTTCCTCCAGTTTCCCTCTTATTTTTCCCCAGCTATGTTGTTATCATTTCCAGAAATCTCTAAAGCTTGCACAGATCCTTAGCACTATGAGATCCATTGAAAGAGATAGTTTTTTTCTTTTTGAGATAGGGCCTGGCTCTGTCACCCAGGCTGTAGCTCAGTGGTGCGATCGAGGCTCACTGCAACCTCTGCCTCCCACGCTCAAGTGATCCTCCCTCCTCAGGCTCCAGAGTAGCTGGGAATACAGGCAGGCAACCACGCCCAGCTAATTTTTGTAATTTTGGTAGAAATGAGATTTTGCCATGTTGCCCAGGCTGGTCTTAAACTGCTGGACTCAAGCAATCCTCCTGCCTTGGCCTCCCAACATGCTAGGATTATAGATGTGAGCCACTGTGCCCAGGCAAAAAGAGATGACTCTCAATAAAAAAAAGTCCTTTTTCTTAAATCACTGTTTCTTTATCTGTGAATTCTTCTTCCAACTAGAAGGAGGAGAAAGAAGTTTGCCTGTATTTCTCACCAGGAGGAGAAGGGGTCTAGTGTGACATCAGAATGAAAGAGTGCTGGAGTTTGAGCCCCTTCTTGCTTTCCAAGATCCCCACAGTTATCACTTCCCATACCCTGGTTTATTCATGTAAACCACACTTATTTTTCTTAGCAGCTACTGTGTACTCGGCTCCATTCTAGGTTCAGATCATTCTATTTGATTAAGACAGAGAGGGTCCCGACTCTCATGGAAGTTACACAACAATAGAGGAGACAGACACTAACCCAATAAGCATTTAACAAAGAATAAAATATTAGAGAGTCATAGTGCACTGAAGAAAAGACATCAGGTTTGTGAAGAAGAGAGACATGGATTCACCTACTTTAGTTCATATGTTTAGGGAGCTCTACCTGAGAAAGTGACATTCAGCTGAGACAACAAAATAAGTAGACAGTCATGAAGATCTAATGGACGAAAGCTCCAGAGAGACCAAATGGGGGGAAAGCCCTGGTGTGGGAAATTATGTGGAGAGAGAGAAAGACGGCTAGAGGGGCTGATGTATAGAAAGTAAGGAAATGGAGAGGCAGAAGATGAGGTAGGACACAGAGAGAAAGTCAGGAGCCTCATCATTATAGGCTCTGATGTCCACGGTAAAAAATTTGAATTTTATTTTATTTATTTATTTATTTATTTATTTATTTATTTTATTTATTTATTTATTTTTTTAAGATGGAGTCTCGCTTTGTTGCCCAGGCTGGAGTGCAGTGGCGTGATCTCAGCTCACTGCAAGCTCCACCTCCTGGGTTCATGCCATTCTCCTGCCTCAGCCTCCCTAGTAGCTGGGACTACAGGCACCTGCCACCACGCCTGGCTAATTTTTAGTATTTTTAGTAGTGAGGGCATTTTGTCATGTTAACCAGGGTGGTCTCCATCTCCTGACCTCGTGATCCACCTGCCTCAGCCTCCCAAAGTGCTGAGACTACAGGTGTGAGCCACCACGCCTGGCCTATTTTTTTTTTTTTGAGACGGAGTTTCGTTCTTGTTGCCCAGGCTGGAGTGCAATGGTGCAATCTCGGCTCACTGCAACCTTCGCCTCCCTGGTTCAAGTGATTCTCCTGCCTCAGCTTCCCAAGTATCTGAGATTACAGGCACCCACCACCATACCTGGCTAATTTTTATTTTTTTGTATTTTTAGTAGACATGGGGTATCACCATGTTGACCATGCTGGTCTCGAACTCCTGACCTCAGATAATCTGTCTGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCATGAGCCACCGCGGCCGGACTGAATTTTATTTAAATAGATACGAGAAGCTACTGTATGGTTACAAGGAGAGTCAATTTATATTCAATTTATTTTTATTTTATTTTATTTTTTTGGTGATGGGGTCTTGCTCTGTTGCCCAGGCTAGATTGTAGTGGCACAATCTCGGCTCACTGCAACCTCTGCCTTCTGGGTTCAGGCGATTCTCCTGACTCAGCCTCCAGAGTAGCTGGGACCACAGGTACATGCCACCACACCTGGCTAAGTTTTTGTATTTTTTAGTAGAGACAGGATTTCACCATGTTAGCCAGGATGGTCTCGATCTCCTGATCTCGTGATCCACCCACCTTGGCCTCCCAAAGTGCAGGGATTACAGGCCTGAGCCACTGCGCCCGGCCTATATTCAATTTTTAAAACTAATTCTAGCTACTCTGTGGGGATTGGAATGTTGGGGTTCACAAGTGGTCAGGAAGACTATTTAGGAGCACAGCAGGGAATTCTCCAGCGAAAACAGGCTTGTGGCTTCATGGAGTGCATTAGTGATAAAGACGGTGAAAAAGATAAAGTGGACAGACTTGGCATGTATTTTTCCTTAGCTTGTTAATGAATTACTGTAAAGGGGGTAGAACAATCAAGCTTATTCCTAAGGATTTTGTTTTGACAAATAAGTGGGTGGTAGTGTTGTTTATTGAGATAGGAAAAACTATGGGAGGAAATTATTTGAAGTGGGTGGTTGGAAATAAAAGTTTTGTTTAAATTTGAGATGATTTATTGACATTTATGTGGAGCAATCAGAAGGTCAATGGCATTTAAGAGACTCATGGTGAGGCTAGGGCTTCAAGTATTTATGTTGGCGGCATCAATACGTGTAGTGTGTTAAATTCCAGGGAGTGGAAGAGGATACATAGGGAGATGGATTGTGTGGAGAAAAAAGAACAGGGCACAGGCCAGCAAAGGGGGCTGAGAAAGAGCCCAGGGATGTTGGAGAAAAATCAAGAGAACATGATGCATGTAAGTCAAGGAAAATAGATTTTTTTCAAGGAGAAGGGAGAGGCCAATTGTGGTGAGTACCACTAAGCGGAGGGGGAAGTGAGAACGTGACAGAGAAGCAAGTGCTGGGTTTGGTGGAGTTGATATTTGCAGTCAGTGGAGTATCCAGGGAGGAAACTGGATTGGACAATTTGAAGAGCGAGTAGAAGTGAGGATGAGGTTAAGGTTGACTGTTTTGAGTAGAGAGGTTCAGGGAAGGACTGCACTCTGGGTTCAGGGAGCCAGCTGGATCAAAAGGAAAAGGCTAAAGAGGCTGAAGAGAAGCAGGAGGACCTGTGAACCAGAGATGCTCAGTCATTATTAGCGAGGAAATACTAGAAAGCCCCTGTGTGCAGTGATGACTACTCATGCAGAAGGTCACACAGCCAATATTTAACACAGCCAGTATTTCACACAGCTAATATTTATTAGTGACATAGAATATACCAGTTATTACTCTAGGTCATGAGAATGGAGTGATAAATAAAATGAATCCGGTCGCCATCAGTATATGCCATGTAACATTTTGCAGTGACTGTGTACCAGGCCTGTGAATTTCAGTATGCAATTTCAATAATGATCCTGCTGTATCTGTGGTGTTTAAAAACATATACATCTCTGGAATCTAAAATTGAGAGGTTATAAGTAAAACCCAGTATTACAAATTGAGTGCTGGAAATCAGATTGCAGTTTAAATCTGAGCATATAGAAAGTCCCTTTCTTCTATGTCAGCAGATGCCTTTTGTGTGAGGTTTAGCTGGACTGCATTATTAGACATAAACCAGTGTTTCTGCCCTATGTTTTCAGAATGACAATTCTTTATGAAACTCATAGAAGAACAGAAGACAACTGCAAAATCATGATGAAGATAGTAATTGCTTTAGAATTAAGGAATACAAAAAATAATGTGAGCTGTAGTTATAGGGATCATAAAAGTTTAAATGGGAATGTATTTGAGTATGTGATCAGTGCTAAGAAGAGTCATCATTTAATTTTACACTTAACAGTAATCTCGTGAGGATTACGCTATTATTAAATGCATTTGATAGATTACAAAAAGGCTTATGGTTGGTAAAAATTGACCCAAGTAGAAGAGATCATGTTTTTATTCAGGTTTTCTGATTCTAGAGTTTGAGAGTTTGTCCATCATTAGTGAGTAGTGACTATATTGTGTCTGAATTATTGACAGAATTTCTGATATTCATATGTACCAGGTTGTTTCTTAGAGTGGGGATAGAGATGCAAGGGCTGCTAGTTCCGATGTATTGGGGAAACTTTCATTCATTTTGCATTTATCATTTTAAAGTTCTGTATGTCTATAATGGTCATGTGTTGAAGAACACAAGGAAGTATTAAATCACTCCTTCTTCTAAGGTTTGACTAGCAAGTTGGGCTAGAGTTACCAAATAAAATACATGTTCCTAGTGAAATCTGAATTTCAGATACAAAACCATAATTTATTGAAAATCCAAATTTAACTGGGCATCCTCTGGTTTTATTTGCCACATCTGTCAACCCTAAGTGCGACACATGGACATGGATTACAGTGCTAACCATGCAAGCCACGGTGACAGCAACTTCACACATGTTTATTTTTAACTTTCTCTGTAAGAAAGTGCTTAGATAATTTAGGGATAAAAAGATAGACATTGCTTGATCCAGGGTGCACACCTCTCTGCCACCATTTCTAAAGGCAAAGGGAGATTTCTGCAGGTCTTGCTCACAGTCTGGGGAGCTGCTCATTTTTGTAAAGTGTCTGTATGAGAATGTCATTTTCTTGGTTTCCTCCTTTCCGAGGGGACTTGACTACAAAACCAAGAGTTCTGCCTCTGGCCAAGGCTGGTAATTTGATGCCTGCTAGTATTGTTGGGAGTGGGAGACTGAAAGAAATGAGTTAGTTGGGGCATTTAACGGGAATAAAATAGCTGTGGTTGTGACTCATTACTACAGATAATTAGTGGACCAGTGGCAGAGAAATTAAGAAAAAAGATGATGTGAATGATAAATGATATGATTAGTGACTGCTTGGTAAGGCAAGGAAATCATTAAATCTTGGTTCTCATCAAGTTCATTTTCTGGAAAGATAGCACTGTATTGGGAGCAGAATTCTACAAAACCTTCTTTTTACATAGGACCAAGATTTTCAACAAATATTTTTCAATGCAATTCTCAGCTGCTCCATAACTAATAGTAGCTTGTTCAACACAGATTTTTTCAGATGATTCACACCTGTGGTACTTACCCAGGGATGGTTCACCACCCCTCCCTTCCCTCTCATCATCGTTGGGGAACGGTGACAATGTTTGGAACAATTTTTGGTTGTCACAAACAGGGGTTTCTTCTGATATTCAATGAGTAGAAGCCAGGGACACTGCTAGAGAACCCACAATGTTCAGAACAGCCTCTGCCATCAACAAGGAATTATCTGGTCCAAAATGTCAATAGTGCTGAGGCTAAGAGCACTGGTTCACACTGTGCTCTTTCTGAAAATTCTAGACTCACATCTGTTATACACTCACCACACAGTTTAGTCTTTTATTTTTGCTTGTTTCATTATAAACAATTAGACAGTTGTATAAATTCAACCACTTTCTTGTTGAATCCATTTAGTCAATGCAAGCTCAATATTTTCATATTTATTTTTTGCCTTATGCAATATTTTTCAACATTTTCATGAGTTGTCGGTCATCACTATCTCTATTAACTTTCAACAACTTGCCCTTGTAAGTCACAAATAGTGCTGCTGCTGAAATTATTTCTCACTAACATGCCTCAGATTTCTGTAGTGATTCTACATTTGATATTATTCACAATGTAAAATGCTTCTATTTATTCATTTCGCTTTTACCCAAGGATTATTTTTAAGTTATTTTTGTCATTTTCACACTTCAAACATAAAGACAAAAACATCAAAAATATAGTGTTTTACATATGTGCATATTTTCACACATATATGTATGTATATTTATATGTATTGAAACTACAGAAGCACATGTCACCAATAAGAGCTCTGAGACACCTTCGACCACTTACCCTTATCAGATGAGTTGTGGAAACAAGTTTTTTTTAACTGAATTTCTGAGCTTTGTGGATTTAGAAATGCAAAGGAAGGTTTGTGGACATTCACAGGGATCATGATTTTATTCTCCTTAAAACTCTTCTGTACTTTCCAATTGTCCTTAGTATAAATCCAAAATCCTAACACCACCCAAGAGGCTTTTCAATACCTGGCTCCTGTGATTTCTCCTGTCTAATCTTTTACCCTCCTTCCCCTCAGCCTCTCTGCTTTAGTGAACTTTCTCCTAGTTTTTTGAAGAAGTTCATCAATTCAAGCTTTTGTACATGGGATTTCCTAAACCTGAAATGTGCCTCCCGTTTTGTCCAAACAGACACGGGCTCCACTCTGCCCCCTGGCTCACACCTGCTTAACCTGTCAAGTCACATCTGAACCGTCACTATTAAGAGGGTCCTTCTCTGGCACCCTAATGTAATTGAGATCATCCTATTATTCTCTGTTCTAGAACTCCACACTTCCGACATTTCTCATTCCTGTCTAAGCTCTTGTGTGTTTGGTGTTGGGCCATCACTTTCACTGCTCTTTAAGCTCCCCCAGCGGAGTGGAGAGGTCTGTTTTCCCTCGTTTGGATTCCTAGAGGCAGCGCAGACCAGGCACAAGGTCAGCACTAAGGAAGGGTTCACAGGATGAACGCGGTGGGTGCTGTTTAAGGAACCGGTAAACGTGTGGGATGAGAGAAGGAGCAGAGTGTCTTTGGGGTGGAGGCTCCCAGGAGGAGGCGGCGCGGGCTGCGGTGCGGGGCGGATCCTCCTCCAGCTCCTGCTTGGAGGTCTCCAGAACAGGCTGGAGGTAGGGAGGGGGGTCCCAAAAGCCTGGGGATCAGACGTGGTTTTCCCGCCTGGTCCCCCAGGCCCCCTTTCGCCTCAGGAAGACAGAGGAGGAGCCCCTGGGCTGCAGGTGGTGGGCGTTGCGGCGGGGGCCGGTTAAGGTTCCCAGTGCCCGCACCCGGCCCAGGGAGCCCCGGATGGCGGCGTCACTGTCAGTGTCTTCTCAGGAGGCCGCCTGTGTGACTGGATCGTTCGTGTCCCCACAGCACGTTTCTTGGAGTACTCTACGTCTGAGTGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCTATAACCAAGAGGAGTACGTGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGATGAGGAGTACTGGAACAGCCAGAAGGACATCCTGGAAGACGAGCGGGCCGCGGTGGACACCTACTGCAGACACAACTACGGGGTTGTGGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCGCGGCGCGGGGCGGGGCCTGAGTCCCTGTGAGCTGGGAATCTGAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAGAGAGAGACAGAGAGACAGAGAGAGAGAGAGCGCCATCTGTGAGCATTTAGAATCCTCTCTATCCTGAGCAAGGAGTTCTGAGGGCACAGGTGTGTGTGTAGAGTGTGGATTTGTCTGTGTCTGTGAGGCTGTTGTGGGAGGGGAGGCAGGAGGGGGCTGCTTCTTATTCTTGGAGGACTCTGTGGGGAGGTGACAAGGGAGGTGGGTGCGGGCGGCTGGAGAGAGAGGTGACCTTGATTGTCTCGGGTCCTTAGAGATGCAGGGAAGGGAAATGTAAGGGGTGTGTGGTTGGGGTGAAGGTTTAGGGGAGGAGAGCTGAGGGGTAAGGAAGGTTTGGGATAATGTGAGGAGGCCAGTTCCAGACTGTCCCTGGCACACACCCTTCATGTAATCTCTGAAATAAAAGTGTGTGCTGTTTGTTTGTAAAAGCATTAGATTAATTTCTAGGGGAATTGAGGAGACCTCTGAGGCATCTCTGAAGCTTCTTTAGGTCTAAATTTCTTGCTAGTTTTTTGTTTTTTATTGTGTATATTTTTACATAGTAGAAATGACTGTGAAACTAACTTTTTGAATTAAAGTTTTAACACAGTTACTATTTTATTATAATGCTAATAGTTTTCTAGTAGTTACATATTATTCTTTTATATATAATAGTTGTGACACAACTTACCTCACTTTCCCCTTTGTTGACCTTTATTATGACATTCACCAAAATTTGAAAATGTATGTTTCTGGTTAATTTTTAATTTATATTTTTTTCATTTATAATTCTTTTGAATTATTTTGACCTATTTATTGGCCAGTTTTAATAACTGCTGTAAGAATTCCCTATTGTATTTGGTAGGGAATGGACAATGATCTACTGCCTAATATCTCGAGGGCTTAGTATTTTTCTCAGTGACTTTGTGGGTTCTTTGTACTGTGAGATTATTAACACTTTATTGATATTTGATTCAGCATTTGCTCCAGTTTGTGGTTTGTATGTTGATTTTGAAAATTCTTTTCCATGTTAAGAATTTGAACATTTTTATATAATAAAATATGTTGCAAAATTTTTATTAATGATTTACAATCCATCTTAAATCTGCCATTTTGTGGTATTGTTGTCTCCAGGTTTCTCCTTACTTCTAAAAAAAATTGCATTTATTGAGAGTCTGCTAGTGTTAGGGATTTTCCTGGGCATAAGCACCCCAAGTGACGAGTCCCAGACACTGCCTTAATCCAAATGTGATTCTGGAAAGAAAAATCATTTTACAATGATAGGCCTAATAATAATTAAGCTTGTGTTGCATGGGAGATGCATTGATCAGCTAAATGTAAATATAAGAACTTTCAAAACTAAAATGACGTTCCTTAATCCTTCTCTCTGCTTTATGACTCATGCTTTTCTGGGAAAGTAAAAATTTGGAGAATCATTTCTGTCTGTCCCACCTTCCCAGGGGCAGAACCATTTCTGTGGTGTTCTAAGGTGTGAGTGCATGGCGGTAGTATTCCTAAAAATTCATATTCGGTTTCGTCATGTACCCAACTCTGTCCCGTTATCTATCAACATTGTTTTAAATCATATATTTCTGTCAAGGTGTACAAGGATGATAAATAGGTGCCAAGTGGAGCACCCAAGTGTGATGAGCCCCCTCACAGTGGAATGGAGTGTGAAGCTTTATGACCTCATAAATTGAAGGTTATCTTCAGTCATTGTTTTATATATTTTACATGCATTAATCCTCATATAATCCCAAGAGGTAAATTAGTATAATTATCCTTCATTATAGGTGACAAAGTTGAGACACAGAAGAATCGAACTCTTAAGGCAGACCTTGGATTTGAACCAGGCAACCTGGCTCAGATATCAGTTTTAATTACTACACTCTGTACTTTCAAAGATTTGTAAACACTTTGACAATGCATGACAATTTCAAGCTATGAAGAAACAAACACAATTTTTCACAATATCTCTCAAATCTAATAGGTCCTCACTATCAAGATTAAGTTCCAGGCTGATGACACTGTAAGGCCACATGGCCAGCTGTGCTGGAGGCCTGGTCAAGGTCAGAGCCTGGGTTTGCAGAGAAGCAGACAAACAGCCAAACAAGGAGACTTACTCTGTCTTCATGACTCATTCCCTCTACCTTTTTTCTCCTAGTCCATCCTAAGGTGACTGTGTATCCTTCAAAGACCCAGCCCCTGCAGCACCACAACCTCCTGGTCTGTTCTGTGAGTGGTTTCTATCCAGGCAGCATTGAAGTCAGGTGGTTCCGGAATGGCCAGGAAGAGAAGACTGGGGTGGTGTCCACAGGCCTGATCCACAATGGAGACTGGACCTTCCAGACCCTGGTGATGCTGGAAACAGTTCCTCGGAGTGGAGAGGTTTACACCTGCCAAGTGGAGCACCCAAGCGTGACAAGCCCTCTCACAGTGGAATGGAGTGAGCAGCTTTCTGACTTCATAAATTTCTCACCCACCAAGAAGGGGACTGTGCTCATCCCTGAGTGTCAGGTTTCTCCTCTCCGACATCCTATTTTCATTTGCTCCATGTTCTCATCTCCATCAGCACAGGTCACTGGGGGTAGCCCTGTAGGTGTTTCTAGAAACACCTGTACCTCCTGGAGAAGCAGTCTCGCCTGCCAGGCAGGAGAGGCTGTCCCTCTTTTGAACCTCCCCATGATGTCACAGGTCAGGGTCACCCACCCTCCCCGGGCTCCAGGCACTGCCTCTGGGTCTGAGACTGAGTTTCTGGTGCTGTTGATCTGAGTTATTTGTTGTGATCTGGGAAGAGGAGAAGTGTAGGGGCCTTCCTGACATGAGGGGAGTCCAATCTCAGCTCTGCCTTTTATTAGCTCTGTCACTCTAGACAAACTACTTAGCCTCATTGAGTCTCAGGCTTTCTGTGGATCAGATGTTGAACTCTTGCCTTACATCAAGGCTGTAATATTTGAATGAGTTTGATGTCTGAACCTTGTAACTGTTCAGTGTGATTTGAAATCCTTTTTTTCTCCAGAAATGGCTAGTTATTTTAGTTCTTGTGGGGCAGACTTCTTCCCCATTTTCAAAGCTCTGAATCTTAGAGTCTCAATTAAAGAGGTTCAATTTGGAATAAACACTAAACCTGGCTTCCTCTCTCAGGAGCACGGTCTGAATCTGCACAGAGCAAGATGCTGAGTGGAGTCGGGGGCTTTGTGCTGGGCCTGCTCTTCCTTGGGGCCGGGCTGTTCATCTACTTCAGGAATCAGAAAGGTGAGGAGCCTTTGGTAGCTGGCTCTCTCCATAGGCTTTTCTGGAGGAGGAACTATGGCTTTGCTGAGGTTAGTTCTCAGTATATGAGTGGCCCTGAATAAAGCCTTTCTTTCCCCAAACGGCTCTAATGTCCTGCTAATCCAGAAATCATCAGTGCATGGTTACTATGTGAAAGCATAATAGCTTGTGGCCTGCAGAGACAAGAGGAAGGTTAACAAGTAGGGGTCCTTTGGTTTGAGATCTTGGAGCAGATTAAGGAAGAGCCACTAAGACTAATGGAATTACACTGGATCCTGTGACAGACACTTCACCCTTCATGGGTCACATGGTCTGTTTCTGCTCCTCTCTGCCCTGGCTGGTGTGGGTTGTAGTGACAGAGAACTCTCCGGTGGGAGATCTGGGGCTGGGACATTGTGTTGGAAGACAGATTTGCTTCCATAAATTTTAAGTGTATATATTTTCCTCTTTTTCCCAGGACACTCTGGACTTCAGCCAAGAGGTAATACCTTTTAATCCTCTTTTAGAAACAGATACGGTTTCCCTAGTGAGAGGTGAAGCCAGCTGGACTTCTGGGTCGGGTAGGGACTTGCAGAACTTTCCTGTCTTAGGAGAGGTTTCTAAATGCACCAATCAGTGCTCTGTAAAAACACACCAATTGGCACTCTGTGGCTAGATAGATGTTTGTAAAATGGACTAATCAGCACTCTGTAAAATGGAGCAATCCACACTCTGTAAAATGGACCAATCAATGCTCTTTAAAATGGACCAATCAGCAGGACATGGGCGGGGACAAATAAGGGAATACAAGCTGGCCACCCCAGCCAGCAGCAGCAACCCGCTCAGGTCGCCTTCCATGCTGTGGAAGCTTTGTTCTTTTGCTCTTCACAATAAATCTTGCTGTTGCTCACTCTTCGGGTCTGTGCCACCTTTAAGAGCTGTAACACTCACTGTGAAGATTCGCGGCTTCATTCTTGAAGTCAGCGAAACCACGAACCCACCGGAAGGAACAAACTCTGGACACACTAGAATTGATGGTAGAGGTGATAAGGCATGAGACAGAAATAATAGGAAAGACTTTGGATCCAAATTTCTGATCAGGCAATTTACACCAAAACTCCTCCTCTCCACTTAGAAAAGGCCTGTGCTCTGTGGGACTATTGGCTCTGGGAGACTCAGGAACTTGTTTTTCTTCTTCCTGCAGTGCTCTCATCTGAGTCCCTGAAAGAGAGGAAAAGAAACTGTTAGTAGAGTCAGGTTGAAAACAACACTCTCCTCTGTCTTTTGCAGGATTCCTGAGCTGAAGTGCAGATGACACATTCAAAGAAGAACTTTCTGCCCCAGCTTTGCAGGATGAAAAGCTTTCCCTCCTGGCTGT'
#     hlaLocus = mySubmissionBatch[1].submittedGene.geneLocus
#     print('I stored the gene name successfully:' + str(hlaLocus))
#     # roughFeatureSequence = 'aag\nCGTCGT\nccg\nGGCTGA\naat'
#     alleleCallWithGFE = fetchSequenceAlleleCallWithGFE(roughFeatureSequence, hlaLocus)
#     assert_true(len(alleleCallWithGFE) > 3)
#
#     annotatedSequence = parseExons(roughFeatureSequence, alleleCallWithGFE)
#     assert_true(len(annotatedSequence) > 3)
#
#     ipdGenerator = IpdSubGenerator()
#
#     ipdGenerator.sequenceAnnotation = identifyGenomicFeatures(annotatedSequence)
#
#     ipdSubmission = ipdGenerator.buildIPDSubmission()
#
#     print('IPD SUBMISSION:\n' + ipdSubmission)
#
#     assert_true(len(ipdSubmission) > 3)
#     assert_true(ipdSubmission is not None)
#
#     # Maybe I don't want to safe my configuration with all these test values but maybe I don't care :D
#     writeConfigurationFile()
#
#     #Print out the submission to a file, for easier validation.
#     outputFileObject = open('Test_IPD_Submission_HLA_DRB1.txt', 'w')
#     submissionText = ipdSubmission
#     outputFileObject.write(submissionText)




# Last time i checked, the DRB1 submission was not working, because the request is too big.


# def testCreateIPDSubmissionFlatfileDRB1():
#     print ('Test: Creating IPD Flatfile')
#
#     mySubmissionBatch = getConfigurationValue('submission_batch').submissionBatch
#
#     roughFeatureSequence = mySubmissionBatch[0].submittedGene.fullSequence
#     # roughFeatureSequence = 'GAAAGACCTGAAAGATCACGGTGCCTTCATTTCAACTGTGAGACATGATGTAATTTTCACAAATCTACAACAGTAAGATATAGTGCAACAGGACCAGATTAAGGTCTCCTGGTTTGCAACCATGTCCCCTCCATCTCCTTTACTCCTGAACACACTCACTCCTGCAAACAGTTCTCTTGTCAAGTGGGAAATGAATGCTCTTACAAGGCTCAAACTTGTGAACACATCACTGACCAGCACAGAGCTAAAATAATTGGGGCTAAAAATACCGCCCCAATTAAAGTGTTTTACATGCAACTGGTTCAAACCTTTCAAGTACTAAAAACAATCCTGTAAAGAAGGAAATTCTGTTTCAGAAGAGGACCTTCATACAGCATCTCTGACCAGCAACTGATGATGCTATTGAACTCAGATGCTGATTCGTTCTCCAACACTAGATTACCCAATCCAGGAGCAAGGAAATCAGTAACTTACTCCCTATAACTTGGAATGTGGGTGGAGGGGTTCATAGTTCTCCCTGAGTGAGACTTGCCTGCTGCTCTGGCCCCTGGTCCTGTCCTGTTCTCCAGCATGGTGTGTCTGAGGCTCCCTGGAGGCTCCTGCATGGCAGTTCTGACAGTGACACTGATGGTGCTGAGCTCCCCACTGGCTTTGGCTGGGGACACCAGACGTAAGTGCACATTGTGGGTGCTGAGCTACTATGGGGTGGGGAAAATAGGGAGTTTTGTTAACATTGTGCCCAGGCCATGTACCTTAAGAAATTGTGACGTATTTTTCAGAGATTGCCCATCTTTATCATATGGATCCCAAATAATTTCCCCACCACAAAAGGAGCTTGGCTACTTGCCCACTCCATGAGACTTGTGTAAGGGGCCTCCATACAGGTCATTTCTACTCAAATCTCCACCAATAAAACCTTTGCATCACATGTCCTCAGGGTCTTTAGAGGATTTAGAAATAAGGATGCTAAAATAAATTCCCCATACAGCACTTGCCTTTATTATGTTGACTTATGTCAGACAAAGGAGGTTTTTTCTGAAAATTTTGTGGAAGTCAAGGGAATTTAAAGGGTCTCTCCTAAACGATCCTGGGTTATGTCACCCACAGGACCTTTGGTGTTGGCCCCTCTTCCTCATATGTGAGGATGGACCCAGTGGCCTCCCCATTATCTACTTTCTTTTCTTTCTGAACTCCAATGTTTATAAAGCCTGTACCCCTGTAGTGTATGTAGGTTGTCTGACAGAAGTTATACTTAGTGCTCTTTCTTTCTTGTGGGGAAAAATCCCTGGAACTGAAGCTGAGATCTTTAGTACTTGGAGTCACCTTACAGATACAGAGCATTTATGAGGTATTCTTTGGTGCCTAAAGAACTTAAGGCATCCTCTGAAAACCTGGCCCAGGTTAGTGTTTATTATGAATCTCTTTTAACCTTTCTATACTTGTTTCTCCTACATCTCCTAAGTGCTCCAACTAGACATGACAGAAGAGATTTAACTAACGTAGTATGAGTTATATAAAATTCTATTTTTGTAAGTCAAAAATAATCAAATATCAAAAATTTAATAATGTTCAAACTATATACTCTGTGTGGGGTTACCGAGACAATGTGGACATTGTTCACATCTCATAGGGCTGAAAGTCAATGGGCAAGTCCTGGGAACTCATTGTCTTACTGGGGTCTTGTCCTAAATTTCATAGGTTCACCCATCATGCCCTCAACTTTCCTTAATTAGCCATGTCTGCTTACCTCTTCCTCCAGTTTCCCTCTTATTTTTCCCCAGCTATGTTGTTATCATTTCCAGAAATCTCTAAAGCTTGCACAGATCCTTAGCACTATGAGATCCATTGAAAGAGATAGTTTTTTTCTTTTTGAGATAGGGCCTGGCTCTGTCACCCAGGCTGTAGCTCAGTGGTGCGATCGAGGCTCACTGCAACCTCTGCCTCCCACGCTCAAGTGATCCTCCCTCCTCAGGCTCCAGAGTAGCTGGGAATACAGGCAGGCAACCACGCCCAGCTAATTTTTGTAATTTTGGTAGAAATGAGATTTTGCCATGTTGCCCAGGCTGGTCTTAAACTGCTGGACTCAAGCAATCCTCCTGCCTTGGCCTCCCAACATGCTAGGATTATAGATGTGAGCCACTGTGCCCAGGCAAAAAGAGATGACTCTCAATAAAAAAAAGTCCTTTTTCTTAAATCACTGTTTCTTTATCTGTGAATTCTTCTTCCAACTAGAAGGAGGAGAAAGAAGTTTGCCTGTATTTCTCACCAGGAGGAGAAGGGGTCTAGTGTGACATCAGAATGAAAGAGTGCTGGAGTTTGAGCCCCTTCTTGCTTTCCAAGATCCCCACAGTTATCACTTCCCATACCCTGGTTTATTCATGTAAACCACACTTATTTTTCTTAGCAGCTACTGTGTACTCGGCTCCATTCTAGGTTCAGATCATTCTATTTGATTAAGACAGAGAGGGTCCCGACTCTCATGGAAGTTACACAACAATAGAGGAGACAGACACTAACCCAATAAGCATTTAACAAAGAATAAAATATTAGAGAGTCATAGTGCACTGAAGAAAAGACATCAGGTTTGTGAAGAAGAGAGACATGGATTCACCTACTTTAGTTCATATGTTTAGGGAGCTCTACCTGAGAAAGTGACATTCAGCTGAGACAACAAAATAAGTAGACAGTCATGAAGATCTAATGGACGAAAGCTCCAGAGAGACCAAATGGGGGGAAAGCCCTGGTGTGGGAAATTATGTGGAGAGAGAGAAAGACGGCTAGAGGGGCTGATGTATAGAAAGTAAGGAAATGGAGAGGCAGAAGATGAGGTAGGACACAGAGAGAAAGTCAGGAGCCTCATCATTATAGGCTCTGATGTCCACGGTAAAAAATTTGAATTTTATTTTATTTATTTATTTATTTATTTATTTATTTATTTTATTTATTTATTTATTTTTTTAAGATGGAGTCTCGCTTTGTTGCCCAGGCTGGAGTGCAGTGGCGTGATCTCAGCTCACTGCAAGCTCCACCTCCTGGGTTCATGCCATTCTCCTGCCTCAGCCTCCCTAGTAGCTGGGACTACAGGCACCTGCCACCACGCCTGGCTAATTTTTAGTATTTTTAGTAGTGAGGGCATTTTGTCATGTTAACCAGGGTGGTCTCCATCTCCTGACCTCGTGATCCACCTGCCTCAGCCTCCCAAAGTGCTGAGACTACAGGTGTGAGCCACCACGCCTGGCCTATTTTTTTTTTTTTGAGACGGAGTTTCGTTCTTGTTGCCCAGGCTGGAGTGCAATGGTGCAATCTCGGCTCACTGCAACCTTCGCCTCCCTGGTTCAAGTGATTCTCCTGCCTCAGCTTCCCAAGTATCTGAGATTACAGGCACCCACCACCATACCTGGCTAATTTTTATTTTTTTGTATTTTTAGTAGACATGGGGTATCACCATGTTGACCATGCTGGTCTCGAACTCCTGACCTCAGATAATCTGTCTGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCATGAGCCACCGCGGCCGGACTGAATTTTATTTAAATAGATACGAGAAGCTACTGTATGGTTACAAGGAGAGTCAATTTATATTCAATTTATTTTTATTTTATTTTATTTTTTTGGTGATGGGGTCTTGCTCTGTTGCCCAGGCTAGATTGTAGTGGCACAATCTCGGCTCACTGCAACCTCTGCCTTCTGGGTTCAGGCGATTCTCCTGACTCAGCCTCCAGAGTAGCTGGGACCACAGGTACATGCCACCACACCTGGCTAAGTTTTTGTATTTTTTAGTAGAGACAGGATTTCACCATGTTAGCCAGGATGGTCTCGATCTCCTGATCTCGTGATCCACCCACCTTGGCCTCCCAAAGTGCAGGGATTACAGGCCTGAGCCACTGCGCCCGGCCTATATTCAATTTTTAAAACTAATTCTAGCTACTCTGTGGGGATTGGAATGTTGGGGTTCACAAGTGGTCAGGAAGACTATTTAGGAGCACAGCAGGGAATTCTCCAGCGAAAACAGGCTTGTGGCTTCATGGAGTGCATTAGTGATAAAGACGGTGAAAAAGATAAAGTGGACAGACTTGGCATGTATTTTTCCTTAGCTTGTTAATGAATTACTGTAAAGGGGGTAGAACAATCAAGCTTATTCCTAAGGATTTTGTTTTGACAAATAAGTGGGTGGTAGTGTTGTTTATTGAGATAGGAAAAACTATGGGAGGAAATTATTTGAAGTGGGTGGTTGGAAATAAAAGTTTTGTTTAAATTTGAGATGATTTATTGACATTTATGTGGAGCAATCAGAAGGTCAATGGCATTTAAGAGACTCATGGTGAGGCTAGGGCTTCAAGTATTTATGTTGGCGGCATCAATACGTGTAGTGTGTTAAATTCCAGGGAGTGGAAGAGGATACATAGGGAGATGGATTGTGTGGAGAAAAAAGAACAGGGCACAGGCCAGCAAAGGGGGCTGAGAAAGAGCCCAGGGATGTTGGAGAAAAATCAAGAGAACATGATGCATGTAAGTCAAGGAAAATAGATTTTTTTCAAGGAGAAGGGAGAGGCCAATTGTGGTGAGTACCACTAAGCGGAGGGGGAAGTGAGAACGTGACAGAGAAGCAAGTGCTGGGTTTGGTGGAGTTGATATTTGCAGTCAGTGGAGTATCCAGGGAGGAAACTGGATTGGACAATTTGAAGAGCGAGTAGAAGTGAGGATGAGGTTAAGGTTGACTGTTTTGAGTAGAGAGGTTCAGGGAAGGACTGCACTCTGGGTTCAGGGAGCCAGCTGGATCAAAAGGAAAAGGCTAAAGAGGCTGAAGAGAAGCAGGAGGACCTGTGAACCAGAGATGCTCAGTCATTATTAGCGAGGAAATACTAGAAAGCCCCTGTGTGCAGTGATGACTACTCATGCAGAAGGTCACACAGCCAATATTTAACACAGCCAGTATTTCACACAGCTAATATTTATTAGTGACATAGAATATACCAGTTATTACTCTAGGTCATGAGAATGGAGTGATAAATAAAATGAATCCGGTCGCCATCAGTATATGCCATGTAACATTTTGCAGTGACTGTGTACCAGGCCTGTGAATTTCAGTATGCAATTTCAATAATGATCCTGCTGTATCTGTGGTGTTTAAAAACATATACATCTCTGGAATCTAAAATTGAGAGGTTATAAGTAAAACCCAGTATTACAAATTGAGTGCTGGAAATCAGATTGCAGTTTAAATCTGAGCATATAGAAAGTCCCTTTCTTCTATGTCAGCAGATGCCTTTTGTGTGAGGTTTAGCTGGACTGCATTATTAGACATAAACCAGTGTTTCTGCCCTATGTTTTCAGAATGACAATTCTTTATGAAACTCATAGAAGAACAGAAGACAACTGCAAAATCATGATGAAGATAGTAATTGCTTTAGAATTAAGGAATACAAAAAATAATGTGAGCTGTAGTTATAGGGATCATAAAAGTTTAAATGGGAATGTATTTGAGTATGTGATCAGTGCTAAGAAGAGTCATCATTTAATTTTACACTTAACAGTAATCTCGTGAGGATTACGCTATTATTAAATGCATTTGATAGATTACAAAAAGGCTTATGGTTGGTAAAAATTGACCCAAGTAGAAGAGATCATGTTTTTATTCAGGTTTTCTGATTCTAGAGTTTGAGAGTTTGTCCATCATTAGTGAGTAGTGACTATATTGTGTCTGAATTATTGACAGAATTTCTGATATTCATATGTACCAGGTTGTTTCTTAGAGTGGGGATAGAGATGCAAGGGCTGCTAGTTCCGATGTATTGGGGAAACTTTCATTCATTTTGCATTTATCATTTTAAAGTTCTGTATGTCTATAATGGTCATGTGTTGAAGAACACAAGGAAGTATTAAATCACTCCTTCTTCTAAGGTTTGACTAGCAAGTTGGGCTAGAGTTACCAAATAAAATACATGTTCCTAGTGAAATCTGAATTTCAGATACAAAACCATAATTTATTGAAAATCCAAATTTAACTGGGCATCCTCTGGTTTTATTTGCCACATCTGTCAACCCTAAGTGCGACACATGGACATGGATTACAGTGCTAACCATGCAAGCCACGGTGACAGCAACTTCACACATGTTTATTTTTAACTTTCTCTGTAAGAAAGTGCTTAGATAATTTAGGGATAAAAAGATAGACATTGCTTGATCCAGGGTGCACACCTCTCTGCCACCATTTCTAAAGGCAAAGGGAGATTTCTGCAGGTCTTGCTCACAGTCTGGGGAGCTGCTCATTTTTGTAAAGTGTCTGTATGAGAATGTCATTTTCTTGGTTTCCTCCTTTCCGAGGGGACTTGACTACAAAACCAAGAGTTCTGCCTCTGGCCAAGGCTGGTAATTTGATGCCTGCTAGTATTGTTGGGAGTGGGAGACTGAAAGAAATGAGTTAGTTGGGGCATTTAACGGGAATAAAATAGCTGTGGTTGTGACTCATTACTACAGATAATTAGTGGACCAGTGGCAGAGAAATTAAGAAAAAAGATGATGTGAATGATAAATGATATGATTAGTGACTGCTTGGTAAGGCAAGGAAATCATTAAATCTTGGTTCTCATCAAGTTCATTTTCTGGAAAGATAGCACTGTATTGGGAGCAGAATTCTACAAAACCTTCTTTTTACATAGGACCAAGATTTTCAACAAATATTTTTCAATGCAATTCTCAGCTGCTCCATAACTAATAGTAGCTTGTTCAACACAGATTTTTTCAGATGATTCACACCTGTGGTACTTACCCAGGGATGGTTCACCACCCCTCCCTTCCCTCTCATCATCGTTGGGGAACGGTGACAATGTTTGGAACAATTTTTGGTTGTCACAAACAGGGGTTTCTTCTGATATTCAATGAGTAGAAGCCAGGGACACTGCTAGAGAACCCACAATGTTCAGAACAGCCTCTGCCATCAACAAGGAATTATCTGGTCCAAAATGTCAATAGTGCTGAGGCTAAGAGCACTGGTTCACACTGTGCTCTTTCTGAAAATTCTAGACTCACATCTGTTATACACTCACCACACAGTTTAGTCTTTTATTTTTGCTTGTTTCATTATAAACAATTAGACAGTTGTATAAATTCAACCACTTTCTTGTTGAATCCATTTAGTCAATGCAAGCTCAATATTTTCATATTTATTTTTTGCCTTATGCAATATTTTTCAACATTTTCATGAGTTGTCGGTCATCACTATCTCTATTAACTTTCAACAACTTGCCCTTGTAAGTCACAAATAGTGCTGCTGCTGAAATTATTTCTCACTAACATGCCTCAGATTTCTGTAGTGATTCTACATTTGATATTATTCACAATGTAAAATGCTTCTATTTATTCATTTCGCTTTTACCCAAGGATTATTTTTAAGTTATTTTTGTCATTTTCACACTTCAAACATAAAGACAAAAACATCAAAAATATAGTGTTTTACATATGTGCATATTTTCACACATATATGTATGTATATTTATATGTATTGAAACTACAGAAGCACATGTCACCAATAAGAGCTCTGAGACACCTTCGACCACTTACCCTTATCAGATGAGTTGTGGAAACAAGTTTTTTTTAACTGAATTTCTGAGCTTTGTGGATTTAGAAATGCAAAGGAAGGTTTGTGGACATTCACAGGGATCATGATTTTATTCTCCTTAAAACTCTTCTGTACTTTCCAATTGTCCTTAGTATAAATCCAAAATCCTAACACCACCCAAGAGGCTTTTCAATACCTGGCTCCTGTGATTTCTCCTGTCTAATCTTTTACCCTCCTTCCCCTCAGCCTCTCTGCTTTAGTGAACTTTCTCCTAGTTTTTTGAAGAAGTTCATCAATTCAAGCTTTTGTACATGGGATTTCCTAAACCTGAAATGTGCCTCCCGTTTTGTCCAAACAGACACGGGCTCCACTCTGCCCCCTGGCTCACACCTGCTTAACCTGTCAAGTCACATCTGAACCGTCACTATTAAGAGGGTCCTTCTCTGGCACCCTAATGTAATTGAGATCATCCTATTATTCTCTGTTCTAGAACTCCACACTTCCGACATTTCTCATTCCTGTCTAAGCTCTTGTGTGTTTGGTGTTGGGCCATCACTTTCACTGCTCTTTAAGCTCCCCCAGCGGAGTGGAGAGGTCTGTTTTCCCTCGTTTGGATTCCTAGAGGCAGCGCAGACCAGGCACAAGGTCAGCACTAAGGAAGGGTTCACAGGATGAACGCGGTGGGTGCTGTTTAAGGAACCGGTAAACGTGTGGGATGAGAGAAGGAGCAGAGTGTCTTTGGGGTGGAGGCTCCCAGGAGGAGGCGGCGCGGGCTGCGGTGCGGGGCGGATCCTCCTCCAGCTCCTGCTTGGAGGTCTCCAGAACAGGCTGGAGGTAGGGAGGGGGGTCCCAAAAGCCTGGGGATCAGACGTGGTTTTCCCGCCTGGTCCCCCAGGCCCCCTTTCGCCTCAGGAAGACAGAGGAGGAGCCCCTGGGCTGCAGGTGGTGGGCGTTGCGGCGGGGGCCGGTTAAGGTTCCCAGTGCCCGCACCCGGCCCAGGGAGCCCCGGATGGCGGCGTCACTGTCAGTGTCTTCTCAGGAGGCCGCCTGTGTGACTGGATCGTTCGTGTCCCCACAGCACGTTTCTTGGAGTACTCTACGTCTGAGTGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCTATAACCAAGAGGAGTACGTGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGATGAGGAGTACTGGAACAGCCAGAAGGACATCCTGGAAGACGAGCGGGCCGCGGTGGACACCTACTGCAGACACAACTACGGGGTTGTGGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCGCGGCGCGGGGCGGGGCCTGAGTCCCTGTGAGCTGGGAATCTGAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAGAGAGAGACAGAGAGACAGAGAGAGAGAGAGCGCCATCTGTGAGCATTTAGAATCCTCTCTATCCTGAGCAAGGAGTTCTGAGGGCACAGGTGTGTGTGTAGAGTGTGGATTTGTCTGTGTCTGTGAGGCTGTTGTGGGAGGGGAGGCAGGAGGGGGCTGCTTCTTATTCTTGGAGGACTCTGTGGGGAGGTGACAAGGGAGGTGGGTGCGGGCGGCTGGAGAGAGAGGTGACCTTGATTGTCTCGGGTCCTTAGAGATGCAGGGAAGGGAAATGTAAGGGGTGTGTGGTTGGGGTGAAGGTTTAGGGGAGGAGAGCTGAGGGGTAAGGAAGGTTTGGGATAATGTGAGGAGGCCAGTTCCAGACTGTCCCTGGCACACACCCTTCATGTAATCTCTGAAATAAAAGTGTGTGCTGTTTGTTTGTAAAAGCATTAGATTAATTTCTAGGGGAATTGAGGAGACCTCTGAGGCATCTCTGAAGCTTCTTTAGGTCTAAATTTCTTGCTAGTTTTTTGTTTTTTATTGTGTATATTTTTACATAGTAGAAATGACTGTGAAACTAACTTTTTGAATTAAAGTTTTAACACAGTTACTATTTTATTATAATGCTAATAGTTTTCTAGTAGTTACATATTATTCTTTTATATATAATAGTTGTGACACAACTTACCTCACTTTCCCCTTTGTTGACCTTTATTATGACATTCACCAAAATTTGAAAATGTATGTTTCTGGTTAATTTTTAATTTATATTTTTTTCATTTATAATTCTTTTGAATTATTTTGACCTATTTATTGGCCAGTTTTAATAACTGCTGTAAGAATTCCCTATTGTATTTGGTAGGGAATGGACAATGATCTACTGCCTAATATCTCGAGGGCTTAGTATTTTTCTCAGTGACTTTGTGGGTTCTTTGTACTGTGAGATTATTAACACTTTATTGATATTTGATTCAGCATTTGCTCCAGTTTGTGGTTTGTATGTTGATTTTGAAAATTCTTTTCCATGTTAAGAATTTGAACATTTTTATATAATAAAATATGTTGCAAAATTTTTATTAATGATTTACAATCCATCTTAAATCTGCCATTTTGTGGTATTGTTGTCTCCAGGTTTCTCCTTACTTCTAAAAAAAATTGCATTTATTGAGAGTCTGCTAGTGTTAGGGATTTTCCTGGGCATAAGCACCCCAAGTGACGAGTCCCAGACACTGCCTTAATCCAAATGTGATTCTGGAAAGAAAAATCATTTTACAATGATAGGCCTAATAATAATTAAGCTTGTGTTGCATGGGAGATGCATTGATCAGCTAAATGTAAATATAAGAACTTTCAAAACTAAAATGACGTTCCTTAATCCTTCTCTCTGCTTTATGACTCATGCTTTTCTGGGAAAGTAAAAATTTGGAGAATCATTTCTGTCTGTCCCACCTTCCCAGGGGCAGAACCATTTCTGTGGTGTTCTAAGGTGTGAGTGCATGGCGGTAGTATTCCTAAAAATTCATATTCGGTTTCGTCATGTACCCAACTCTGTCCCGTTATCTATCAACATTGTTTTAAATCATATATTTCTGTCAAGGTGTACAAGGATGATAAATAGGTGCCAAGTGGAGCACCCAAGTGTGATGAGCCCCCTCACAGTGGAATGGAGTGTGAAGCTTTATGACCTCATAAATTGAAGGTTATCTTCAGTCATTGTTTTATATATTTTACATGCATTAATCCTCATATAATCCCAAGAGGTAAATTAGTATAATTATCCTTCATTATAGGTGACAAAGTTGAGACACAGAAGAATCGAACTCTTAAGGCAGACCTTGGATTTGAACCAGGCAACCTGGCTCAGATATCAGTTTTAATTACTACACTCTGTACTTTCAAAGATTTGTAAACACTTTGACAATGCATGACAATTTCAAGCTATGAAGAAACAAACACAATTTTTCACAATATCTCTCAAATCTAATAGGTCCTCACTATCAAGATTAAGTTCCAGGCTGATGACACTGTAAGGCCACATGGCCAGCTGTGCTGGAGGCCTGGTCAAGGTCAGAGCCTGGGTTTGCAGAGAAGCAGACAAACAGCCAAACAAGGAGACTTACTCTGTCTTCATGACTCATTCCCTCTACCTTTTTTCTCCTAGTCCATCCTAAGGTGACTGTGTATCCTTCAAAGACCCAGCCCCTGCAGCACCACAACCTCCTGGTCTGTTCTGTGAGTGGTTTCTATCCAGGCAGCATTGAAGTCAGGTGGTTCCGGAATGGCCAGGAAGAGAAGACTGGGGTGGTGTCCACAGGCCTGATCCACAATGGAGACTGGACCTTCCAGACCCTGGTGATGCTGGAAACAGTTCCTCGGAGTGGAGAGGTTTACACCTGCCAAGTGGAGCACCCAAGCGTGACAAGCCCTCTCACAGTGGAATGGAGTGAGCAGCTTTCTGACTTCATAAATTTCTCACCCACCAAGAAGGGGACTGTGCTCATCCCTGAGTGTCAGGTTTCTCCTCTCCGACATCCTATTTTCATTTGCTCCATGTTCTCATCTCCATCAGCACAGGTCACTGGGGGTAGCCCTGTAGGTGTTTCTAGAAACACCTGTACCTCCTGGAGAAGCAGTCTCGCCTGCCAGGCAGGAGAGGCTGTCCCTCTTTTGAACCTCCCCATGATGTCACAGGTCAGGGTCACCCACCCTCCCCGGGCTCCAGGCACTGCCTCTGGGTCTGAGACTGAGTTTCTGGTGCTGTTGATCTGAGTTATTTGTTGTGATCTGGGAAGAGGAGAAGTGTAGGGGCCTTCCTGACATGAGGGGAGTCCAATCTCAGCTCTGCCTTTTATTAGCTCTGTCACTCTAGACAAACTACTTAGCCTCATTGAGTCTCAGGCTTTCTGTGGATCAGATGTTGAACTCTTGCCTTACATCAAGGCTGTAATATTTGAATGAGTTTGATGTCTGAACCTTGTAACTGTTCAGTGTGATTTGAAATCCTTTTTTTCTCCAGAAATGGCTAGTTATTTTAGTTCTTGTGGGGCAGACTTCTTCCCCATTTTCAAAGCTCTGAATCTTAGAGTCTCAATTAAAGAGGTTCAATTTGGAATAAACACTAAACCTGGCTTCCTCTCTCAGGAGCACGGTCTGAATCTGCACAGAGCAAGATGCTGAGTGGAGTCGGGGGCTTTGTGCTGGGCCTGCTCTTCCTTGGGGCCGGGCTGTTCATCTACTTCAGGAATCAGAAAGGTGAGGAGCCTTTGGTAGCTGGCTCTCTCCATAGGCTTTTCTGGAGGAGGAACTATGGCTTTGCTGAGGTTAGTTCTCAGTATATGAGTGGCCCTGAATAAAGCCTTTCTTTCCCCAAACGGCTCTAATGTCCTGCTAATCCAGAAATCATCAGTGCATGGTTACTATGTGAAAGCATAATAGCTTGTGGCCTGCAGAGACAAGAGGAAGGTTAACAAGTAGGGGTCCTTTGGTTTGAGATCTTGGAGCAGATTAAGGAAGAGCCACTAAGACTAATGGAATTACACTGGATCCTGTGACAGACACTTCACCCTTCATGGGTCACATGGTCTGTTTCTGCTCCTCTCTGCCCTGGCTGGTGTGGGTTGTAGTGACAGAGAACTCTCCGGTGGGAGATCTGGGGCTGGGACATTGTGTTGGAAGACAGATTTGCTTCCATAAATTTTAAGTGTATATATTTTCCTCTTTTTCCCAGGACACTCTGGACTTCAGCCAAGAGGTAATACCTTTTAATCCTCTTTTAGAAACAGATACGGTTTCCCTAGTGAGAGGTGAAGCCAGCTGGACTTCTGGGTCGGGTAGGGACTTGCAGAACTTTCCTGTCTTAGGAGAGGTTTCTAAATGCACCAATCAGTGCTCTGTAAAAACACACCAATTGGCACTCTGTGGCTAGATAGATGTTTGTAAAATGGACTAATCAGCACTCTGTAAAATGGAGCAATCCACACTCTGTAAAATGGACCAATCAATGCTCTTTAAAATGGACCAATCAGCAGGACATGGGCGGGGACAAATAAGGGAATACAAGCTGGCCACCCCAGCCAGCAGCAGCAACCCGCTCAGGTCGCCTTCCATGCTGTGGAAGCTTTGTTCTTTTGCTCTTCACAATAAATCTTGCTGTTGCTCACTCTTCGGGTCTGTGCCACCTTTAAGAGCTGTAACACTCACTGTGAAGATTCGCGGCTTCATTCTTGAAGTCAGCGAAACCACGAACCCACCGGAAGGAACAAACTCTGGACACACTAGAATTGATGGTAGAGGTGATAAGGCATGAGACAGAAATAATAGGAAAGACTTTGGATCCAAATTTCTGATCAGGCAATTTACACCAAAACTCCTCCTCTCCACTTAGAAAAGGCCTGTGCTCTGTGGGACTATTGGCTCTGGGAGACTCAGGAACTTGTTTTTCTTCTTCCTGCAGTGCTCTCATCTGAGTCCCTGAAAGAGAGGAAAAGAAACTGTTAGTAGAGTCAGGTTGAAAACAACACTCTCCTCTGTCTTTTGCAGGATTCCTGAGCTGAAGTGCAGATGACACATTCAAAGAAGAACTTTCTGCCCCAGCTTTGCAGGATGAAAAGCTTTCCCTCCTGGCTGT'
#     hlaLocus = mySubmissionBatch[0].submittedGene.geneLocus
#     print('I stored the gene name successfully:' + str(hlaLocus))
#     # roughFeatureSequence = 'aag\nCGTCGT\nccg\nGGCTGA\naat'
#     alleleCallWithGFE = fetchSequenceAlleleCallWithGFE(roughFeatureSequence, hlaLocus)
#     assert_true(len(alleleCallWithGFE) > 3)
#
#     # print 'ALLELE CALL RESULTS:'
#     # print alleleCallWithGFE
#
#     annotatedSequence = parseExons(roughFeatureSequence, alleleCallWithGFE)
#     #print ('ANNOTATED SEQUENCE:')
#     #print (annotatedSequence)
#     assert_true(len(annotatedSequence) > 3)
#
#     ipdGenerator = IpdSubGenerator()
#
#     ipdGenerator.sequenceAnnotation = identifyGenomicFeatures(annotatedSequence)
#
#
#     ipdSubmission = ipdGenerator.buildIPDSubmission()
#
#     print('IPD SUBMISSION:\n' + ipdSubmission)
#
#     assert_true(len(ipdSubmission) > 3)
#     assert_true(ipdSubmission is not None)
#
#     # Maybe I don't want to safe my configuration with all these test values but maybe I don't care :D
#     writeConfigurationFile()
#
#     #Print out the submission to a file, for easier validation.
#     outputFileObject = open('Test_IPD_Submission_HLA_DRB1.txt', 'w')
#     submissionText = ipdSubmission
#     outputFileObject.write(submissionText)





# def testCreateIPDSubmissionFlatfileC():
#     logEvent ('Test: Creating IPD Flatfile','INFO')
#
#     initializeGlobalVariables()
#
#     submission1 = AlleleSubmission()
#     submission1.submittedGene.fullSequence = 'AATCAGGACGAAGTCCCAGGTCCCGGGCGGGGCTCTCAGGGTCTCAGGCTCCAAGGGCCGTGTCTGCACTGGGGAGGCGCCGCGTTGAGGATTCTCCACTCCCCTGAGTTTCACTTCTTCTCCCAACCTGCGTCGGGTCCTTCTTCCTGAATACTCATGACGCGTCCCCAATTCCCACTCCCATTGGGTGTCGGATTCTAGAGAAGCCAATCAGCGTCTCCGCAGTCCCGGTTCTGAAGTCCCCAGTCACCCACCCGGACTCAGATTCTCCCCAGACGCCGAGATGCGGGTCATGGCGCCCCGAACCCTCATCCTGCTGCTCTCGGGAGCCCTGGCCCTGACCGAGACCTGGGCCGGTGAGTGCGGGGTTGGGAGGGAATCGGCCTCTGCGGAGAGGAGCGAGGGGCCCGCCCGGCGAGGGCGCAGGACCCGGGGAGCCGCGCAGGGAGGAGGGTCGGGCGGGTCTCAGCCCCTCCTCGCCCCCAGGCTCCCACTCCATGAGGTATTTCTACACCGCTGTGTCCCGGCCCGGCCGCGGGGAGCCCCACTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGTCCGAGAGGGGAGCCGCGGGCGCCGTGGGTGGAGCAGGAGGGGCCGGAGTATTGGGACCGGGAGACACAGAAGTACAAGCGCCAGGCACAGACTGACCGAGTGAGCCTGCGGAACCTGCGCGGCTACTACAACCAGAGCGAGGCCAGTGAGTGACCCCGGCCCGGGGCGCAGGTCACGACCCCTCCTCATCCCCCACGGACGGCCCGGGTCGCCCCAAGTCTCCCGGTCTGAGATCCACCCCGAGGCTGCGGAACCCGCCCAGACCCTCGACCGGAGAGAGCCCCAGTCACCTTTACCCGGTTTCATTTTCAGTTTAGGCCAAAATCCCCGCGGGTTGGTCGGGGCGGGGCGGGGCTCGGGGGACGGGGCTGACCGCGGGGGCGGGGCCAGGGTCTCACATCATCCAGAGGATGTATGGCTGCGACGTGGGGCCCGACGGGCGCCTCCTCCGCGGGTATGACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGATCTGCGCTCCTGGACCGCCGCGGACACGGCGGCTCAGATCACCCAGCGTAAGTGGGAGGCGGCCCGTGAGGCGGAGCAGCTGAGAGCCTACCTGGAGGGCCTGTGCGTGGAGTGGCTCCGCAGATACCTGAAGAATGGGAAGGAGACGCTGCAGCGCGCGGGTACCAGGGGCAGTGGGGAGCCTTCCCCATCTCCTGTAGATCTCCCGGGATGGCCTCCCACGAGGAGGGGAGGAAAATGGGATCAGCGCTAGAATATCGCCCTCCCTTGAATGGAGAATGGGATGAGTTTTCCTGAGTTTCCTCTGAGGGCCCCCTCTGCTCTCTAGGACAATTAAGGGATGAAGTCCTTGAGGAAATGGAGGGGAAGACAGTCCCTAGAATACTGATCAGGGGTCCCCTTTGACCACTTTGACCACTGCAGCAGCTGTGGTCAGGCTGCTGACCTTTCTCTCAGGCCTTGTTCTCTGCCTCACGCTCAATGTGTTTGAAGGTTTGATTCCAGCTTTTCTGAGTCCTTCGGCCTCCACTCAGGTCAGGACCAGAAGTCGCTGTTCCTCCCTCAGAGACTAGAACTTTCCAATGAATAGGAGATTATCCCAGGTGCCTGTGTCCAGGCTGGCGTCTGGGTTCTGTGCCCCCTTCCCCACCCCAGGTGTCCTGTCCGTTCTCAGGATGGTCACATGGGCGCTGTTGGAGTGTCGCAAGAGAGATACAAAGTGTCTGAATTTTCTGACTCTTCCCGTCAGAACACCCAAAGACACACGTGACCCACCATCCCGTCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGTGGGATGGGGAGGACCAAACTCAGGACACTGAGCTTGTGGAGACCAGGCCAGCAGGAGATGGAACCTTCCAGAAGTGGGCAGCTGTGGTGGTGCCTTCTGGAGAAGAGCAGAGATACACGTGCCATGTGCAGCACGAGGGGCTGCCGGAGCCCCTCACCCTGAGATGGGGTAAGGAGGGGGATGAGGGGTGATGTGTCTTCTCAGGGAAAGCAGAAGTCCTGGAGCCCTTCAGCCAGGTCAGGGCTGAGGCTTGGGGGTCAGGGCCCCTCACCTTCCCCTCCTTTCCCAGAGCCGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCGTTGCTGGCCTGGCTGTCCTGGCTGTCCTAGCTGTCCTAGGAGCTGTGGTGGCTGTTGTGATGTGTAGGAGGAAGAGCTCAGGTAGGGAAGGGGTGAGGAGTGGGGTCTGGGTTTTCTTGTTCCACTGGGAGTTTCAAGCCCCAGGTAGAAGTGTGCCCCACCTCGTTACTGGAAGCACCATCCACACATGGGCCATCCCAGCCTGGGACCCTGTGTGCCAGCACTTACTCTGTTGTGAAGCACATGACAATGAAGGACAGATGTATCACCTTGATGATTATGGTGTTGGGGTCCTTGATTCCAGCATTCATGAGTCAGGGGAAGGTCCCTGCTAAGGACAGACCTTAGGAGGGCAGTTGCTCCAGAACCCACAGCTGCTTTCCCCGTGTTTCCTGATCCTGCCCTGGGTCTGCAGTCATAGTTCTGGAAACTTCTCTTGGGTCCAAGACTAGGAGGTTCCCCTAAGATCGCATGGCCCTGACTCCTCCCTGTCCCCTCACAGGGCATTTTCTTCCCACAGGTGGAAAAGGAGGGAGCTGCTCTCAGGCTGCGTGTAAGTGATGGCGGTGGGCGTGTGGAGGAGCTCACCCACCCCATAATTCCTCTTGTCCCACATCTCCTGCGGGCTCTGACCAGGTCTTTTTTTTTGTTCTACCCCAGCCAGCAACAGTGCCCAGGGCTCTGATGAGTCTCTCATCGCTTGTAAAGGTGAGATTCTGGGGAGCTGAAGTGGTCGGGGGTGGGGCAGAGGGAAAAGGCCTAGGTAATGGGGATCCTTTGATTGGGACGTTTCGAATGTGTGGTGAGCTGTTCAGAGTGTGATCACTTACCATGACTGACCTGAATTTGTTCATGACTATTGTGTTCTGTAGCCTGAGACAGCTGCCTGTGTGGGACTGAGATGCAGGATTTCTTCACACCTCTCCTTTGTGACTTCAAGAGCCTCTGGCATCTCTTTCTGCAAAGGCATCTGAATGTGTCTGCGTTCCTGTTAGCATAATGTGAGGAGGTGGACAGACAGCCCACCCCCGTGTCCACCGTGACCCCT'
#     submission1.submittedGene.geneLocus = 'HLA-C'
#     submission1.localAlleleName = 'HLA-C_NEW_1'
#     submission1.closestAlleleWrittenDescription = 'Some C allele'
#     submission1.ipdSubmissionIdentifier = 'HWS177477467'
#     submission1.ipdSubmissionVersion = '1'
#     submission1.cellId = 'Cosang_Cell_Line_1'
#     submission1.ethnicOrigin = 'Caucasoid - Dutch/Irish/English, Europe'
#     submission1.sex = 'M'
#     submission1.consanguineous = 'Unknown'
#     submission1.homozygous = 'Unknown'
#     submission1.typedAlleles = 'HLA-A*01:01:01:01;HLA-B*40:72:01;HLA-DRB1*07:75;HLA-A*03:73;HLA-B*44:87;HLA-DRB1*14:90'
#     submission1.materialAvailability = 'No Material Available'
#     submission1.cellBank = 'Not Available'
#     submission1.primarySequencingMethodology = 'Direct sequencing of PCR product from DNA (SBT)'
#     submission1.secondarySequencingMethodology = 'MinION 2D Amplicon Sequencing'
#     submission1.primerType = 'Both allele and locus specific'
#     submission1.primers = '03PID03     CCCAAAGGGTTTCCCGGGAAATTT 3UT 3015-3042;04PID04     AAAGGGTTTCCCGGGAAATTTCCC 5UT 4015-4042'
#     submission1.sequencedInIsolation = 'Yes'
#     submission1.numOfReactions = '3'
#     submission1.methodComments = 'Sanger SBT with confirmatory MinION Long Range sequencing.'
#     submission1.citations = '1-3503;PUBMED; 9349617.;Laforet M, Froelich N, Parissiadis A, Pfeiffer B, Schell A, Faller B, Woehl-Jaegle ML, Cazenave JP, Tongio MM;A nucleotide insertion in exon 4 is responsible for the absence of expression of an HLA-A*01 allele;Tissue Antigens 50:347 - 50(1997).;1-3503;PUBMED; 9349617.;Laforet M, Froelich N, Parissiadis A, Pfeiffer B, Schell A, Faller B, Woehl-Jaegle ML, Cazenave JP, Tongio MM;A nucleotide insertion in exon 4 is responsible for the absence of expression of an HLA-A*01 allele;Tissue Antigens 50:347 - 50(1997).'
#
#     submissionBatch = getConfigurationValue('submission_batch')
#     submissionBatch.submissionBatch.append(submission1)
#
#     roughFeatureSequence = submission1.submittedGene.fullSequence
#     # roughFeatureSequence = 'GAAAGACCTGAAAGATCACGGTGCCTTCATTTCAACTGTGAGACATGATGTAATTTTCACAAATCTACAACAGTAAGATATAGTGCAACAGGACCAGATTAAGGTCTCCTGGTTTGCAACCATGTCCCCTCCATCTCCTTTACTCCTGAACACACTCACTCCTGCAAACAGTTCTCTTGTCAAGTGGGAAATGAATGCTCTTACAAGGCTCAAACTTGTGAACACATCACTGACCAGCACAGAGCTAAAATAATTGGGGCTAAAAATACCGCCCCAATTAAAGTGTTTTACATGCAACTGGTTCAAACCTTTCAAGTACTAAAAACAATCCTGTAAAGAAGGAAATTCTGTTTCAGAAGAGGACCTTCATACAGCATCTCTGACCAGCAACTGATGATGCTATTGAACTCAGATGCTGATTCGTTCTCCAACACTAGATTACCCAATCCAGGAGCAAGGAAATCAGTAACTTACTCCCTATAACTTGGAATGTGGGTGGAGGGGTTCATAGTTCTCCCTGAGTGAGACTTGCCTGCTGCTCTGGCCCCTGGTCCTGTCCTGTTCTCCAGCATGGTGTGTCTGAGGCTCCCTGGAGGCTCCTGCATGGCAGTTCTGACAGTGACACTGATGGTGCTGAGCTCCCCACTGGCTTTGGCTGGGGACACCAGACGTAAGTGCACATTGTGGGTGCTGAGCTACTATGGGGTGGGGAAAATAGGGAGTTTTGTTAACATTGTGCCCAGGCCATGTACCTTAAGAAATTGTGACGTATTTTTCAGAGATTGCCCATCTTTATCATATGGATCCCAAATAATTTCCCCACCACAAAAGGAGCTTGGCTACTTGCCCACTCCATGAGACTTGTGTAAGGGGCCTCCATACAGGTCATTTCTACTCAAATCTCCACCAATAAAACCTTTGCATCACATGTCCTCAGGGTCTTTAGAGGATTTAGAAATAAGGATGCTAAAATAAATTCCCCATACAGCACTTGCCTTTATTATGTTGACTTATGTCAGACAAAGGAGGTTTTTTCTGAAAATTTTGTGGAAGTCAAGGGAATTTAAAGGGTCTCTCCTAAACGATCCTGGGTTATGTCACCCACAGGACCTTTGGTGTTGGCCCCTCTTCCTCATATGTGAGGATGGACCCAGTGGCCTCCCCATTATCTACTTTCTTTTCTTTCTGAACTCCAATGTTTATAAAGCCTGTACCCCTGTAGTGTATGTAGGTTGTCTGACAGAAGTTATACTTAGTGCTCTTTCTTTCTTGTGGGGAAAAATCCCTGGAACTGAAGCTGAGATCTTTAGTACTTGGAGTCACCTTACAGATACAGAGCATTTATGAGGTATTCTTTGGTGCCTAAAGAACTTAAGGCATCCTCTGAAAACCTGGCCCAGGTTAGTGTTTATTATGAATCTCTTTTAACCTTTCTATACTTGTTTCTCCTACATCTCCTAAGTGCTCCAACTAGACATGACAGAAGAGATTTAACTAACGTAGTATGAGTTATATAAAATTCTATTTTTGTAAGTCAAAAATAATCAAATATCAAAAATTTAATAATGTTCAAACTATATACTCTGTGTGGGGTTACCGAGACAATGTGGACATTGTTCACATCTCATAGGGCTGAAAGTCAATGGGCAAGTCCTGGGAACTCATTGTCTTACTGGGGTCTTGTCCTAAATTTCATAGGTTCACCCATCATGCCCTCAACTTTCCTTAATTAGCCATGTCTGCTTACCTCTTCCTCCAGTTTCCCTCTTATTTTTCCCCAGCTATGTTGTTATCATTTCCAGAAATCTCTAAAGCTTGCACAGATCCTTAGCACTATGAGATCCATTGAAAGAGATAGTTTTTTTCTTTTTGAGATAGGGCCTGGCTCTGTCACCCAGGCTGTAGCTCAGTGGTGCGATCGAGGCTCACTGCAACCTCTGCCTCCCACGCTCAAGTGATCCTCCCTCCTCAGGCTCCAGAGTAGCTGGGAATACAGGCAGGCAACCACGCCCAGCTAATTTTTGTAATTTTGGTAGAAATGAGATTTTGCCATGTTGCCCAGGCTGGTCTTAAACTGCTGGACTCAAGCAATCCTCCTGCCTTGGCCTCCCAACATGCTAGGATTATAGATGTGAGCCACTGTGCCCAGGCAAAAAGAGATGACTCTCAATAAAAAAAAGTCCTTTTTCTTAAATCACTGTTTCTTTATCTGTGAATTCTTCTTCCAACTAGAAGGAGGAGAAAGAAGTTTGCCTGTATTTCTCACCAGGAGGAGAAGGGGTCTAGTGTGACATCAGAATGAAAGAGTGCTGGAGTTTGAGCCCCTTCTTGCTTTCCAAGATCCCCACAGTTATCACTTCCCATACCCTGGTTTATTCATGTAAACCACACTTATTTTTCTTAGCAGCTACTGTGTACTCGGCTCCATTCTAGGTTCAGATCATTCTATTTGATTAAGACAGAGAGGGTCCCGACTCTCATGGAAGTTACACAACAATAGAGGAGACAGACACTAACCCAATAAGCATTTAACAAAGAATAAAATATTAGAGAGTCATAGTGCACTGAAGAAAAGACATCAGGTTTGTGAAGAAGAGAGACATGGATTCACCTACTTTAGTTCATATGTTTAGGGAGCTCTACCTGAGAAAGTGACATTCAGCTGAGACAACAAAATAAGTAGACAGTCATGAAGATCTAATGGACGAAAGCTCCAGAGAGACCAAATGGGGGGAAAGCCCTGGTGTGGGAAATTATGTGGAGAGAGAGAAAGACGGCTAGAGGGGCTGATGTATAGAAAGTAAGGAAATGGAGAGGCAGAAGATGAGGTAGGACACAGAGAGAAAGTCAGGAGCCTCATCATTATAGGCTCTGATGTCCACGGTAAAAAATTTGAATTTTATTTTATTTATTTATTTATTTATTTATTTATTTATTTTATTTATTTATTTATTTTTTTAAGATGGAGTCTCGCTTTGTTGCCCAGGCTGGAGTGCAGTGGCGTGATCTCAGCTCACTGCAAGCTCCACCTCCTGGGTTCATGCCATTCTCCTGCCTCAGCCTCCCTAGTAGCTGGGACTACAGGCACCTGCCACCACGCCTGGCTAATTTTTAGTATTTTTAGTAGTGAGGGCATTTTGTCATGTTAACCAGGGTGGTCTCCATCTCCTGACCTCGTGATCCACCTGCCTCAGCCTCCCAAAGTGCTGAGACTACAGGTGTGAGCCACCACGCCTGGCCTATTTTTTTTTTTTTGAGACGGAGTTTCGTTCTTGTTGCCCAGGCTGGAGTGCAATGGTGCAATCTCGGCTCACTGCAACCTTCGCCTCCCTGGTTCAAGTGATTCTCCTGCCTCAGCTTCCCAAGTATCTGAGATTACAGGCACCCACCACCATACCTGGCTAATTTTTATTTTTTTGTATTTTTAGTAGACATGGGGTATCACCATGTTGACCATGCTGGTCTCGAACTCCTGACCTCAGATAATCTGTCTGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCATGAGCCACCGCGGCCGGACTGAATTTTATTTAAATAGATACGAGAAGCTACTGTATGGTTACAAGGAGAGTCAATTTATATTCAATTTATTTTTATTTTATTTTATTTTTTTGGTGATGGGGTCTTGCTCTGTTGCCCAGGCTAGATTGTAGTGGCACAATCTCGGCTCACTGCAACCTCTGCCTTCTGGGTTCAGGCGATTCTCCTGACTCAGCCTCCAGAGTAGCTGGGACCACAGGTACATGCCACCACACCTGGCTAAGTTTTTGTATTTTTTAGTAGAGACAGGATTTCACCATGTTAGCCAGGATGGTCTCGATCTCCTGATCTCGTGATCCACCCACCTTGGCCTCCCAAAGTGCAGGGATTACAGGCCTGAGCCACTGCGCCCGGCCTATATTCAATTTTTAAAACTAATTCTAGCTACTCTGTGGGGATTGGAATGTTGGGGTTCACAAGTGGTCAGGAAGACTATTTAGGAGCACAGCAGGGAATTCTCCAGCGAAAACAGGCTTGTGGCTTCATGGAGTGCATTAGTGATAAAGACGGTGAAAAAGATAAAGTGGACAGACTTGGCATGTATTTTTCCTTAGCTTGTTAATGAATTACTGTAAAGGGGGTAGAACAATCAAGCTTATTCCTAAGGATTTTGTTTTGACAAATAAGTGGGTGGTAGTGTTGTTTATTGAGATAGGAAAAACTATGGGAGGAAATTATTTGAAGTGGGTGGTTGGAAATAAAAGTTTTGTTTAAATTTGAGATGATTTATTGACATTTATGTGGAGCAATCAGAAGGTCAATGGCATTTAAGAGACTCATGGTGAGGCTAGGGCTTCAAGTATTTATGTTGGCGGCATCAATACGTGTAGTGTGTTAAATTCCAGGGAGTGGAAGAGGATACATAGGGAGATGGATTGTGTGGAGAAAAAAGAACAGGGCACAGGCCAGCAAAGGGGGCTGAGAAAGAGCCCAGGGATGTTGGAGAAAAATCAAGAGAACATGATGCATGTAAGTCAAGGAAAATAGATTTTTTTCAAGGAGAAGGGAGAGGCCAATTGTGGTGAGTACCACTAAGCGGAGGGGGAAGTGAGAACGTGACAGAGAAGCAAGTGCTGGGTTTGGTGGAGTTGATATTTGCAGTCAGTGGAGTATCCAGGGAGGAAACTGGATTGGACAATTTGAAGAGCGAGTAGAAGTGAGGATGAGGTTAAGGTTGACTGTTTTGAGTAGAGAGGTTCAGGGAAGGACTGCACTCTGGGTTCAGGGAGCCAGCTGGATCAAAAGGAAAAGGCTAAAGAGGCTGAAGAGAAGCAGGAGGACCTGTGAACCAGAGATGCTCAGTCATTATTAGCGAGGAAATACTAGAAAGCCCCTGTGTGCAGTGATGACTACTCATGCAGAAGGTCACACAGCCAATATTTAACACAGCCAGTATTTCACACAGCTAATATTTATTAGTGACATAGAATATACCAGTTATTACTCTAGGTCATGAGAATGGAGTGATAAATAAAATGAATCCGGTCGCCATCAGTATATGCCATGTAACATTTTGCAGTGACTGTGTACCAGGCCTGTGAATTTCAGTATGCAATTTCAATAATGATCCTGCTGTATCTGTGGTGTTTAAAAACATATACATCTCTGGAATCTAAAATTGAGAGGTTATAAGTAAAACCCAGTATTACAAATTGAGTGCTGGAAATCAGATTGCAGTTTAAATCTGAGCATATAGAAAGTCCCTTTCTTCTATGTCAGCAGATGCCTTTTGTGTGAGGTTTAGCTGGACTGCATTATTAGACATAAACCAGTGTTTCTGCCCTATGTTTTCAGAATGACAATTCTTTATGAAACTCATAGAAGAACAGAAGACAACTGCAAAATCATGATGAAGATAGTAATTGCTTTAGAATTAAGGAATACAAAAAATAATGTGAGCTGTAGTTATAGGGATCATAAAAGTTTAAATGGGAATGTATTTGAGTATGTGATCAGTGCTAAGAAGAGTCATCATTTAATTTTACACTTAACAGTAATCTCGTGAGGATTACGCTATTATTAAATGCATTTGATAGATTACAAAAAGGCTTATGGTTGGTAAAAATTGACCCAAGTAGAAGAGATCATGTTTTTATTCAGGTTTTCTGATTCTAGAGTTTGAGAGTTTGTCCATCATTAGTGAGTAGTGACTATATTGTGTCTGAATTATTGACAGAATTTCTGATATTCATATGTACCAGGTTGTTTCTTAGAGTGGGGATAGAGATGCAAGGGCTGCTAGTTCCGATGTATTGGGGAAACTTTCATTCATTTTGCATTTATCATTTTAAAGTTCTGTATGTCTATAATGGTCATGTGTTGAAGAACACAAGGAAGTATTAAATCACTCCTTCTTCTAAGGTTTGACTAGCAAGTTGGGCTAGAGTTACCAAATAAAATACATGTTCCTAGTGAAATCTGAATTTCAGATACAAAACCATAATTTATTGAAAATCCAAATTTAACTGGGCATCCTCTGGTTTTATTTGCCACATCTGTCAACCCTAAGTGCGACACATGGACATGGATTACAGTGCTAACCATGCAAGCCACGGTGACAGCAACTTCACACATGTTTATTTTTAACTTTCTCTGTAAGAAAGTGCTTAGATAATTTAGGGATAAAAAGATAGACATTGCTTGATCCAGGGTGCACACCTCTCTGCCACCATTTCTAAAGGCAAAGGGAGATTTCTGCAGGTCTTGCTCACAGTCTGGGGAGCTGCTCATTTTTGTAAAGTGTCTGTATGAGAATGTCATTTTCTTGGTTTCCTCCTTTCCGAGGGGACTTGACTACAAAACCAAGAGTTCTGCCTCTGGCCAAGGCTGGTAATTTGATGCCTGCTAGTATTGTTGGGAGTGGGAGACTGAAAGAAATGAGTTAGTTGGGGCATTTAACGGGAATAAAATAGCTGTGGTTGTGACTCATTACTACAGATAATTAGTGGACCAGTGGCAGAGAAATTAAGAAAAAAGATGATGTGAATGATAAATGATATGATTAGTGACTGCTTGGTAAGGCAAGGAAATCATTAAATCTTGGTTCTCATCAAGTTCATTTTCTGGAAAGATAGCACTGTATTGGGAGCAGAATTCTACAAAACCTTCTTTTTACATAGGACCAAGATTTTCAACAAATATTTTTCAATGCAATTCTCAGCTGCTCCATAACTAATAGTAGCTTGTTCAACACAGATTTTTTCAGATGATTCACACCTGTGGTACTTACCCAGGGATGGTTCACCACCCCTCCCTTCCCTCTCATCATCGTTGGGGAACGGTGACAATGTTTGGAACAATTTTTGGTTGTCACAAACAGGGGTTTCTTCTGATATTCAATGAGTAGAAGCCAGGGACACTGCTAGAGAACCCACAATGTTCAGAACAGCCTCTGCCATCAACAAGGAATTATCTGGTCCAAAATGTCAATAGTGCTGAGGCTAAGAGCACTGGTTCACACTGTGCTCTTTCTGAAAATTCTAGACTCACATCTGTTATACACTCACCACACAGTTTAGTCTTTTATTTTTGCTTGTTTCATTATAAACAATTAGACAGTTGTATAAATTCAACCACTTTCTTGTTGAATCCATTTAGTCAATGCAAGCTCAATATTTTCATATTTATTTTTTGCCTTATGCAATATTTTTCAACATTTTCATGAGTTGTCGGTCATCACTATCTCTATTAACTTTCAACAACTTGCCCTTGTAAGTCACAAATAGTGCTGCTGCTGAAATTATTTCTCACTAACATGCCTCAGATTTCTGTAGTGATTCTACATTTGATATTATTCACAATGTAAAATGCTTCTATTTATTCATTTCGCTTTTACCCAAGGATTATTTTTAAGTTATTTTTGTCATTTTCACACTTCAAACATAAAGACAAAAACATCAAAAATATAGTGTTTTACATATGTGCATATTTTCACACATATATGTATGTATATTTATATGTATTGAAACTACAGAAGCACATGTCACCAATAAGAGCTCTGAGACACCTTCGACCACTTACCCTTATCAGATGAGTTGTGGAAACAAGTTTTTTTTAACTGAATTTCTGAGCTTTGTGGATTTAGAAATGCAAAGGAAGGTTTGTGGACATTCACAGGGATCATGATTTTATTCTCCTTAAAACTCTTCTGTACTTTCCAATTGTCCTTAGTATAAATCCAAAATCCTAACACCACCCAAGAGGCTTTTCAATACCTGGCTCCTGTGATTTCTCCTGTCTAATCTTTTACCCTCCTTCCCCTCAGCCTCTCTGCTTTAGTGAACTTTCTCCTAGTTTTTTGAAGAAGTTCATCAATTCAAGCTTTTGTACATGGGATTTCCTAAACCTGAAATGTGCCTCCCGTTTTGTCCAAACAGACACGGGCTCCACTCTGCCCCCTGGCTCACACCTGCTTAACCTGTCAAGTCACATCTGAACCGTCACTATTAAGAGGGTCCTTCTCTGGCACCCTAATGTAATTGAGATCATCCTATTATTCTCTGTTCTAGAACTCCACACTTCCGACATTTCTCATTCCTGTCTAAGCTCTTGTGTGTTTGGTGTTGGGCCATCACTTTCACTGCTCTTTAAGCTCCCCCAGCGGAGTGGAGAGGTCTGTTTTCCCTCGTTTGGATTCCTAGAGGCAGCGCAGACCAGGCACAAGGTCAGCACTAAGGAAGGGTTCACAGGATGAACGCGGTGGGTGCTGTTTAAGGAACCGGTAAACGTGTGGGATGAGAGAAGGAGCAGAGTGTCTTTGGGGTGGAGGCTCCCAGGAGGAGGCGGCGCGGGCTGCGGTGCGGGGCGGATCCTCCTCCAGCTCCTGCTTGGAGGTCTCCAGAACAGGCTGGAGGTAGGGAGGGGGGTCCCAAAAGCCTGGGGATCAGACGTGGTTTTCCCGCCTGGTCCCCCAGGCCCCCTTTCGCCTCAGGAAGACAGAGGAGGAGCCCCTGGGCTGCAGGTGGTGGGCGTTGCGGCGGGGGCCGGTTAAGGTTCCCAGTGCCCGCACCCGGCCCAGGGAGCCCCGGATGGCGGCGTCACTGTCAGTGTCTTCTCAGGAGGCCGCCTGTGTGACTGGATCGTTCGTGTCCCCACAGCACGTTTCTTGGAGTACTCTACGTCTGAGTGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCTATAACCAAGAGGAGTACGTGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGATGAGGAGTACTGGAACAGCCAGAAGGACATCCTGGAAGACGAGCGGGCCGCGGTGGACACCTACTGCAGACACAACTACGGGGTTGTGGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCGCGGCGCGGGGCGGGGCCTGAGTCCCTGTGAGCTGGGAATCTGAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAGAGAGAGACAGAGAGACAGAGAGAGAGAGAGCGCCATCTGTGAGCATTTAGAATCCTCTCTATCCTGAGCAAGGAGTTCTGAGGGCACAGGTGTGTGTGTAGAGTGTGGATTTGTCTGTGTCTGTGAGGCTGTTGTGGGAGGGGAGGCAGGAGGGGGCTGCTTCTTATTCTTGGAGGACTCTGTGGGGAGGTGACAAGGGAGGTGGGTGCGGGCGGCTGGAGAGAGAGGTGACCTTGATTGTCTCGGGTCCTTAGAGATGCAGGGAAGGGAAATGTAAGGGGTGTGTGGTTGGGGTGAAGGTTTAGGGGAGGAGAGCTGAGGGGTAAGGAAGGTTTGGGATAATGTGAGGAGGCCAGTTCCAGACTGTCCCTGGCACACACCCTTCATGTAATCTCTGAAATAAAAGTGTGTGCTGTTTGTTTGTAAAAGCATTAGATTAATTTCTAGGGGAATTGAGGAGACCTCTGAGGCATCTCTGAAGCTTCTTTAGGTCTAAATTTCTTGCTAGTTTTTTGTTTTTTATTGTGTATATTTTTACATAGTAGAAATGACTGTGAAACTAACTTTTTGAATTAAAGTTTTAACACAGTTACTATTTTATTATAATGCTAATAGTTTTCTAGTAGTTACATATTATTCTTTTATATATAATAGTTGTGACACAACTTACCTCACTTTCCCCTTTGTTGACCTTTATTATGACATTCACCAAAATTTGAAAATGTATGTTTCTGGTTAATTTTTAATTTATATTTTTTTCATTTATAATTCTTTTGAATTATTTTGACCTATTTATTGGCCAGTTTTAATAACTGCTGTAAGAATTCCCTATTGTATTTGGTAGGGAATGGACAATGATCTACTGCCTAATATCTCGAGGGCTTAGTATTTTTCTCAGTGACTTTGTGGGTTCTTTGTACTGTGAGATTATTAACACTTTATTGATATTTGATTCAGCATTTGCTCCAGTTTGTGGTTTGTATGTTGATTTTGAAAATTCTTTTCCATGTTAAGAATTTGAACATTTTTATATAATAAAATATGTTGCAAAATTTTTATTAATGATTTACAATCCATCTTAAATCTGCCATTTTGTGGTATTGTTGTCTCCAGGTTTCTCCTTACTTCTAAAAAAAATTGCATTTATTGAGAGTCTGCTAGTGTTAGGGATTTTCCTGGGCATAAGCACCCCAAGTGACGAGTCCCAGACACTGCCTTAATCCAAATGTGATTCTGGAAAGAAAAATCATTTTACAATGATAGGCCTAATAATAATTAAGCTTGTGTTGCATGGGAGATGCATTGATCAGCTAAATGTAAATATAAGAACTTTCAAAACTAAAATGACGTTCCTTAATCCTTCTCTCTGCTTTATGACTCATGCTTTTCTGGGAAAGTAAAAATTTGGAGAATCATTTCTGTCTGTCCCACCTTCCCAGGGGCAGAACCATTTCTGTGGTGTTCTAAGGTGTGAGTGCATGGCGGTAGTATTCCTAAAAATTCATATTCGGTTTCGTCATGTACCCAACTCTGTCCCGTTATCTATCAACATTGTTTTAAATCATATATTTCTGTCAAGGTGTACAAGGATGATAAATAGGTGCCAAGTGGAGCACCCAAGTGTGATGAGCCCCCTCACAGTGGAATGGAGTGTGAAGCTTTATGACCTCATAAATTGAAGGTTATCTTCAGTCATTGTTTTATATATTTTACATGCATTAATCCTCATATAATCCCAAGAGGTAAATTAGTATAATTATCCTTCATTATAGGTGACAAAGTTGAGACACAGAAGAATCGAACTCTTAAGGCAGACCTTGGATTTGAACCAGGCAACCTGGCTCAGATATCAGTTTTAATTACTACACTCTGTACTTTCAAAGATTTGTAAACACTTTGACAATGCATGACAATTTCAAGCTATGAAGAAACAAACACAATTTTTCACAATATCTCTCAAATCTAATAGGTCCTCACTATCAAGATTAAGTTCCAGGCTGATGACACTGTAAGGCCACATGGCCAGCTGTGCTGGAGGCCTGGTCAAGGTCAGAGCCTGGGTTTGCAGAGAAGCAGACAAACAGCCAAACAAGGAGACTTACTCTGTCTTCATGACTCATTCCCTCTACCTTTTTTCTCCTAGTCCATCCTAAGGTGACTGTGTATCCTTCAAAGACCCAGCCCCTGCAGCACCACAACCTCCTGGTCTGTTCTGTGAGTGGTTTCTATCCAGGCAGCATTGAAGTCAGGTGGTTCCGGAATGGCCAGGAAGAGAAGACTGGGGTGGTGTCCACAGGCCTGATCCACAATGGAGACTGGACCTTCCAGACCCTGGTGATGCTGGAAACAGTTCCTCGGAGTGGAGAGGTTTACACCTGCCAAGTGGAGCACCCAAGCGTGACAAGCCCTCTCACAGTGGAATGGAGTGAGCAGCTTTCTGACTTCATAAATTTCTCACCCACCAAGAAGGGGACTGTGCTCATCCCTGAGTGTCAGGTTTCTCCTCTCCGACATCCTATTTTCATTTGCTCCATGTTCTCATCTCCATCAGCACAGGTCACTGGGGGTAGCCCTGTAGGTGTTTCTAGAAACACCTGTACCTCCTGGAGAAGCAGTCTCGCCTGCCAGGCAGGAGAGGCTGTCCCTCTTTTGAACCTCCCCATGATGTCACAGGTCAGGGTCACCCACCCTCCCCGGGCTCCAGGCACTGCCTCTGGGTCTGAGACTGAGTTTCTGGTGCTGTTGATCTGAGTTATTTGTTGTGATCTGGGAAGAGGAGAAGTGTAGGGGCCTTCCTGACATGAGGGGAGTCCAATCTCAGCTCTGCCTTTTATTAGCTCTGTCACTCTAGACAAACTACTTAGCCTCATTGAGTCTCAGGCTTTCTGTGGATCAGATGTTGAACTCTTGCCTTACATCAAGGCTGTAATATTTGAATGAGTTTGATGTCTGAACCTTGTAACTGTTCAGTGTGATTTGAAATCCTTTTTTTCTCCAGAAATGGCTAGTTATTTTAGTTCTTGTGGGGCAGACTTCTTCCCCATTTTCAAAGCTCTGAATCTTAGAGTCTCAATTAAAGAGGTTCAATTTGGAATAAACACTAAACCTGGCTTCCTCTCTCAGGAGCACGGTCTGAATCTGCACAGAGCAAGATGCTGAGTGGAGTCGGGGGCTTTGTGCTGGGCCTGCTCTTCCTTGGGGCCGGGCTGTTCATCTACTTCAGGAATCAGAAAGGTGAGGAGCCTTTGGTAGCTGGCTCTCTCCATAGGCTTTTCTGGAGGAGGAACTATGGCTTTGCTGAGGTTAGTTCTCAGTATATGAGTGGCCCTGAATAAAGCCTTTCTTTCCCCAAACGGCTCTAATGTCCTGCTAATCCAGAAATCATCAGTGCATGGTTACTATGTGAAAGCATAATAGCTTGTGGCCTGCAGAGACAAGAGGAAGGTTAACAAGTAGGGGTCCTTTGGTTTGAGATCTTGGAGCAGATTAAGGAAGAGCCACTAAGACTAATGGAATTACACTGGATCCTGTGACAGACACTTCACCCTTCATGGGTCACATGGTCTGTTTCTGCTCCTCTCTGCCCTGGCTGGTGTGGGTTGTAGTGACAGAGAACTCTCCGGTGGGAGATCTGGGGCTGGGACATTGTGTTGGAAGACAGATTTGCTTCCATAAATTTTAAGTGTATATATTTTCCTCTTTTTCCCAGGACACTCTGGACTTCAGCCAAGAGGTAATACCTTTTAATCCTCTTTTAGAAACAGATACGGTTTCCCTAGTGAGAGGTGAAGCCAGCTGGACTTCTGGGTCGGGTAGGGACTTGCAGAACTTTCCTGTCTTAGGAGAGGTTTCTAAATGCACCAATCAGTGCTCTGTAAAAACACACCAATTGGCACTCTGTGGCTAGATAGATGTTTGTAAAATGGACTAATCAGCACTCTGTAAAATGGAGCAATCCACACTCTGTAAAATGGACCAATCAATGCTCTTTAAAATGGACCAATCAGCAGGACATGGGCGGGGACAAATAAGGGAATACAAGCTGGCCACCCCAGCCAGCAGCAGCAACCCGCTCAGGTCGCCTTCCATGCTGTGGAAGCTTTGTTCTTTTGCTCTTCACAATAAATCTTGCTGTTGCTCACTCTTCGGGTCTGTGCCACCTTTAAGAGCTGTAACACTCACTGTGAAGATTCGCGGCTTCATTCTTGAAGTCAGCGAAACCACGAACCCACCGGAAGGAACAAACTCTGGACACACTAGAATTGATGGTAGAGGTGATAAGGCATGAGACAGAAATAATAGGAAAGACTTTGGATCCAAATTTCTGATCAGGCAATTTACACCAAAACTCCTCCTCTCCACTTAGAAAAGGCCTGTGCTCTGTGGGACTATTGGCTCTGGGAGACTCAGGAACTTGTTTTTCTTCTTCCTGCAGTGCTCTCATCTGAGTCCCTGAAAGAGAGGAAAAGAAACTGTTAGTAGAGTCAGGTTGAAAACAACACTCTCCTCTGTCTTTTGCAGGATTCCTGAGCTGAAGTGCAGATGACACATTCAAAGAAGAACTTTCTGCCCCAGCTTTGCAGGATGAAAAGCTTTCCCTCCTGGCTGT'
#     hlaLocus = submission1.submittedGene.geneLocus
#     print('I stored the gene name successfully:' + str(hlaLocus))
#     # roughFeatureSequence = 'aag\nCGTCGT\nccg\nGGCTGA\naat'
#     alleleCallWithGFE = fetchSequenceAlleleCallWithGFE(roughFeatureSequence, hlaLocus)
#     assert_true(len(alleleCallWithGFE) > 3)
#
#     # print 'ALLELE CALL RESULTS:'
#     # print alleleCallWithGFE
#
#     annotatedSequence = parseExons(roughFeatureSequence, alleleCallWithGFE)
#     #print ('ANNOTATED SEQUENCE:')
#     #print (annotatedSequence)
#     assert_true(len(annotatedSequence) > 3)
#
#     ipdGenerator = IpdSubGenerator()
#
#     ipdGenerator.sequenceAnnotation = identifyGenomicFeatures(annotatedSequence)
#
#
#     ipdSubmission = ipdGenerator.buildIPDSubmission()
#
#     print('IPD SUBMISSION:\n' + ipdSubmission)
#
#     assert_true(len(ipdSubmission) > 3)
#     assert_true(ipdSubmission is not None)
#
#
#     writeConfigurationFile()
#
#     # Print out the submission to a file, for easier validation.
#     outputFileObject = open('Test_IPD_Submission_HLA_C.txt', 'w')
#     submissionText = ipdSubmission
#     outputFileObject.write(submissionText)
#
#
#
# def testCreateIPDSubmissionFlatfileDRB1_wSnp():
#
#     print ('Test: Creating IPD Flatfile')
#     # assignConfigurationValue('nmdp_act_rest_address', 'http://act.b12x.org/type_align' )
#     assignConfigurationValue('nmdp_act_rest_address', 'http://localhost/type_align')
#     assert_true(True)
#     # >HLA-A*02:01:01:12 Full Length Allele. Will we allele call correctly?
#     # > HLA - DRB1 * 11:02:01 I will add a single SNP in exon 3
#     roughFeatureSequence = 'GAAAGACCTGAAAGATCACGGTGCCTTCATTTCAACTGTGAGACATGATGTAATTTTCACAAATCTACAACAGTAAGATATAGTGCAACAGGACCAGATTAAGGTCTCCTGGTTTGCAACCATGTCCCCTCCATCTCCTTTACTCCTGAACACACTCACTCCTGCAAACAGTTCTCTTGTCAAGTGGGAAATGAATGCTCTTACAAGGCTCAAACTTGTGAACACATCACTGACCAGCACAGAGCTAAAATAATTGGGGCTAAAAATACCGCCCCAATTAAAGTGTTTTACATGCAACTGGTTCAAACCTTTCAAGTACTAAAAACAATCCTGTAAAGAAGGAAATTCTGTTTCAGAAGAGGACCTTCATACAGCATCTCTGACCAGCAACTGATGATGCTATTGAACTCAGATGCTGATTCGTTCTCCAACACTAGATTACCCAATCCAGGAGCAAGGAAATCAGTAACTTACTCCCTATAACTTGGAATGTGGGTGGAGGGGTTCATAGTTCTCCCTGAGTGAGACTTGCCTGCTGCTCTGGCCCCTGGTCCTGTCCTGTTCTCCAGCATGGTGTGTCTGAGGCTCCCTGGAGGCTCCTGCATGGCAGTTCTGACAGTGACACTGATGGTGCTGAGCTCCCCACTGGCTTTGGCTGGGGACACCAGACGTAAGTGCACATTGTGGGTGCTGAGCTACTATGGGGTGGGGAAAATAGGGAGTTTTGTTAACATTGTGCCCAGGCCATGTACCTTAAGAAATTGTGACGTATTTTTCAGAGATTGCCCATCTTTATCATATGGATCCCAAATAATTTCCCCACCACAAAAGGAGCTTGGCTACTTGCCCACTCCATGAGACTTGTGTAAGGGGCCTCCATACAGGTCATTTCTACTCAAATCTCCACCAATAAAACCTTTGCATCACATGTCCTCAGGGTCTTTAGAGGATTTAGAAATAAGGATGCTAAAATAAATTCCCCATACAGCACTTGCCTTTATTATGTTGACTTATGTCAGACAAAGGAGGTTTTTTCTGAAAATTTTGTGGAAGTCAAGGGAATTTAAAGGGTCTCTCCTAAACGATCCTGGGTTATGTCACCCACAGGACCTTTGGTGTTGGCCCCTCTTCCTCATATGTGAGGATGGACCCAGTGGCCTCCCCATTATCTACTTTCTTTTCTTTCTGAACTCCAATGTTTATAAAGCCTGTACCCCTGTAGTGTATGTAGGTTGTCTGACAGAAGTTATACTTAGTGCTCTTTCTTTCTTGTGGGGAAAAATCCCTGGAACTGAAGCTGAGATCTTTAGTACTTGGAGTCACCTTACAGATACAGAGCATTTATGAGGTATTCTTTGGTGCCTAAAGAACTTAAGGCATCCTCTGAAAACCTGGCCCAGGTTAGTGTTTATTATGAATCTCTTTTAACCTTTCTATACTTGTTTCTCCTACATCTCCTAAGTGCTCCAACTAGACATGACAGAAGAGATTTAACTAACGTAGTATGAGTTATATAAAATTCTATTTTTGTAAGTCAAAAATAATCAAATATCAAAAATTTAATAATGTTCAAACTATATACTCTGTGTGGGGTTACCGAGACAATGTGGACATTGTTCACATCTCATAGGGCTGAAAGTCAATGGGCAAGTCCTGGGAACTCATTGTCTTACTGGGGTCTTGTCCTAAATTTCATAGGTTCACCCATCATGCCCTCAACTTTCCTTAATTAGCCATGTCTGCTTACCTCTTCCTCCAGTTTCCCTCTTATTTTTCCCCAGCTATGTTGTTATCATTTCCAGAAATCTCTAAAGCTTGCACAGATCCTTAGCACTATGAGATCCATTGAAAGAGATAGTTTTTTTCTTTTTGAGATAGGGCCTGGCTCTGTCACCCAGGCTGTAGCTCAGTGGTGCGATCGAGGCTCACTGCAACCTCTGCCTCCCACGCTCAAGTGATCCTCCCTCCTCAGGCTCCAGAGTAGCTGGGAATACAGGCAGGCAACCACGCCCAGCTAATTTTTGTAATTTTGGTAGAAATGAGATTTTGCCATGTTGCCCAGGCTGGTCTTAAACTGCTGGACTCAAGCAATCCTCCTGCCTTGGCCTCCCAACATGCTAGGATTATAGATGTGAGCCACTGTGCCCAGGCAAAAAGAGATGACTCTCAATAAAAAAAAGTCCTTTTTCTTAAATCACTGTTTCTTTATCTGTGAATTCTTCTTCCAACTAGAAGGAGGAGAAAGAAGTTTGCCTGTATTTCTCACCAGGAGGAGAAGGGGTCTAGTGTGACATCAGAATGAAAGAGTGCTGGAGTTTGAGCCCCTTCTTGCTTTCCAAGATCCCCACAGTTATCACTTCCCATACCCTGGTTTATTCATGTAAACCACACTTATTTTTCTTAGCAGCTACTGTGTACTCGGCTCCATTCTAGGTTCAGATCATTCTATTTGATTAAGACAGAGAGGGTCCCGACTCTCATGGAAGTTACACAACAATAGAGGAGACAGACACTAACCCAATAAGCATTTAACAAAGAATAAAATATTAGAGAGTCATAGTGCACTGAAGAAAAGACATCAGGTTTGTGAAGAAGAGAGACATGGATTCACCTACTTTAGTTCATATGTTTAGGGAGCTCTACCTGAGAAAGTGACATTCAGCTGAGACAACAAAATAAGTAGACAGTCATGAAGATCTAATGGACGAAAGCTCCAGAGAGACCAAATGGGGGGAAAGCCCTGGTGTGGGAAATTATGTGGAGAGAGAGAAAGACGGCTAGAGGGGCTGATGTATAGAAAGTAAGGAAATGGAGAGGCAGAAGATGAGGTAGGACACAGAGAGAAAGTCAGGAGCCTCATCATTATAGGCTCTGATGTCCACGGTAAAAAATTTGAATTTTATTTTATTTATTTATTTATTTATTTATTTATTTATTTTATTTATTTATTTATTTTTTTAAGATGGAGTCTCGCTTTGTTGCCCAGGCTGGAGTGCAGTGGCGTGATCTCAGCTCACTGCAAGCTCCACCTCCTGGGTTCATGCCATTCTCCTGCCTCAGCCTCCCTAGTAGCTGGGACTACAGGCACCTGCCACCACGCCTGGCTAATTTTTAGTATTTTTAGTAGTGAGGGCATTTTGTCATGTTAACCAGGGTGGTCTCCATCTCCTGACCTCGTGATCCACCTGCCTCAGCCTCCCAAAGTGCTGAGACTACAGGTGTGAGCCACCACGCCTGGCCTATTTTTTTTTTTTTGAGACGGAGTTTCGTTCTTGTTGCCCAGGCTGGAGTGCAATGGTGCAATCTCGGCTCACTGCAACCTTCGCCTCCCTGGTTCAAGTGATTCTCCTGCCTCAGCTTCCCAAGTATCTGAGATTACAGGCACCCACCACCATACCTGGCTAATTTTTATTTTTTTGTATTTTTAGTAGACATGGGGTATCACCATGTTGACCATGCTGGTCTCGAACTCCTGACCTCAGATAATCTGTCTGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCATGAGCCACCGCGGCCGGACTGAATTTTATTTAAATAGATACGAGAAGCTACTGTATGGTTACAAGGAGAGTCAATTTATATTCAATTTATTTTTATTTTATTTTATTTTTTTGGTGATGGGGTCTTGCTCTGTTGCCCAGGCTAGATTGTAGTGGCACAATCTCGGCTCACTGCAACCTCTGCCTTCTGGGTTCAGGCGATTCTCCTGACTCAGCCTCCAGAGTAGCTGGGACCACAGGTACATGCCACCACACCTGGCTAAGTTTTTGTATTTTTTAGTAGAGACAGGATTTCACCATGTTAGCCAGGATGGTCTCGATCTCCTGATCTCGTGATCCACCCACCTTGGCCTCCCAAAGTGCAGGGATTACAGGCCTGAGCCACTGCGCCCGGCCTATATTCAATTTTTAAAACTAATTCTAGCTACTCTGTGGGGATTGGAATGTTGGGGTTCACAAGTGGTCAGGAAGACTATTTAGGAGCACAGCAGGGAATTCTCCAGCGAAAACAGGCTTGTGGCTTCATGGAGTGCATTAGTGATAAAGACGGTGAAAAAGATAAAGTGGACAGACTTGGCATGTATTTTTCCTTAGCTTGTTAATGAATTACTGTAAAGGGGGTAGAACAATCAAGCTTATTCCTAAGGATTTTGTTTTGACAAATAAGTGGGTGGTAGTGTTGTTTATTGAGATAGGAAAAACTATGGGAGGAAATTATTTGAAGTGGGTGGTTGGAAATAAAAGTTTTGTTTAAATTTGAGATGATTTATTGACATTTATGTGGAGCAATCAGAAGGTCAATGGCATTTAAGAGACTCATGGTGAGGCTAGGGCTTCAAGTATTTATGTTGGCGGCATCAATACGTGTAGTGTGTTAAATTCCAGGGAGTGGAAGAGGATACATAGGGAGATGGATTGTGTGGAGAAAAAAGAACAGGGCACAGGCCAGCAAAGGGGGCTGAGAAAGAGCCCAGGGATGTTGGAGAAAAATCAAGAGAACATGATGCATGTAAGTCAAGGAAAATAGATTTTTTTCAAGGAGAAGGGAGAGGCCAATTGTGGTGAGTACCACTAAGCGGAGGGGGAAGTGAGAACGTGACAGAGAAGCAAGTGCTGGGTTTGGTGGAGTTGATATTTGCAGTCAGTGGAGTATCCAGGGAGGAAACTGGATTGGACAATTTGAAGAGCGAGTAGAAGTGAGGATGAGGTTAAGGTTGACTGTTTTGAGTAGAGAGGTTCAGGGAAGGACTGCACTCTGGGTTCAGGGAGCCAGCTGGATCAAAAGGAAAAGGCTAAAGAGGCTGAAGAGAAGCAGGAGGACCTGTGAACCAGAGATGCTCAGTCATTATTAGCGAGGAAATACTAGAAAGCCCCTGTGTGCAGTGATGACTACTCATGCAGAAGGTCACACAGCCAATATTTAACACAGCCAGTATTTCACACAGCTAATATTTATTAGTGACATAGAATATACCAGTTATTACTCTAGGTCATGAGAATGGAGTGATAAATAAAATGAATCCGGTCGCCATCAGTATATGCCATGTAACATTTTGCAGTGACTGTGTACCAGGCCTGTGAATTTCAGTATGCAATTTCAATAATGATCCTGCTGTATCTGTGGTGTTTAAAAACATATACATCTCTGGAATCTAAAATTGAGAGGTTATAAGTAAAACCCAGTATTACAAATTGAGTGCTGGAAATCAGATTGCAGTTTAAATCTGAGCATATAGAAAGTCCCTTTCTTCTATGTCAGCAGATGCCTTTTGTGTGAGGTTTAGCTGGACTGCATTATTAGACATAAACCAGTGTTTCTGCCCTATGTTTTCAGAATGACAATTCTTTATGAAACTCATAGAAGAACAGAAGACAACTGCAAAATCATGATGAAGATAGTAATTGCTTTAGAATTAAGGAATACAAAAAATAATGTGAGCTGTAGTTATAGGGATCATAAAAGTTTAAATGGGAATGTATTTGAGTATGTGATCAGTGCTAAGAAGAGTCATCATTTAATTTTACACTTAACAGTAATCTCGTGAGGATTACGCTATTATTAAATGCATTTGATAGATTACAAAAAGGCTTATGGTTGGTAAAAATTGACCCAAGTAGAAGAGATCATGTTTTTATTCAGGTTTTCTGATTCTAGAGTTTGAGAGTTTGTCCATCATTAGTGAGTAGTGACTATATTGTGTCTGAATTATTGACAGAATTTCTGATATTCATATGTACCAGGTTGTTTCTTAGAGTGGGGATAGAGATGCAAGGGCTGCTAGTTCCGATGTATTGGGGAAACTTTCATTCATTTTGCATTTATCATTTTAAAGTTCTGTATGTCTATAATGGTCATGTGTTGAAGAACACAAGGAAGTATTAAATCACTCCTTCTTCTAAGGTTTGACTAGCAAGTTGGGCTAGAGTTACCAAATAAAATACATGTTCCTAGTGAAATCTGAATTTCAGATACAAAACCATAATTTATTGAAAATCCAAATTTAACTGGGCATCCTCTGGTTTTATTTGCCACATCTGTCAACCCTAAGTGCGACACATGGACATGGATTACAGTGCTAACCATGCAAGCCACGGTGACAGCAACTTCACACATGTTTATTTTTAACTTTCTCTGTAAGAAAGTGCTTAGATAATTTAGGGATAAAAAGATAGACATTGCTTGATCCAGGGTGCACACCTCTCTGCCACCATTTCTAAAGGCAAAGGGAGATTTCTGCAGGTCTTGCTCACAGTCTGGGGAGCTGCTCATTTTTGTAAAGTGTCTGTATGAGAATGTCATTTTCTTGGTTTCCTCCTTTCCGAGGGGACTTGACTACAAAACCAAGAGTTCTGCCTCTGGCCAAGGCTGGTAATTTGATGCCTGCTAGTATTGTTGGGAGTGGGAGACTGAAAGAAATGAGTTAGTTGGGGCATTTAACGGGAATAAAATAGCTGTGGTTGTGACTCATTACTACAGATAATTAGTGGACCAGTGGCAGAGAAATTAAGAAAAAAGATGATGTGAATGATAAATGATATGATTAGTGACTGCTTGGTAAGGCAAGGAAATCATTAAATCTTGGTTCTCATCAAGTTCATTTTCTGGAAAGATAGCACTGTATTGGGAGCAGAATTCTACAAAACCTTCTTTTTACATAGGACCAAGATTTTCAACAAATATTTTTCAATGCAATTCTCAGCTGCTCCATAACTAATAGTAGCTTGTTCAACACAGATTTTTTCAGATGATTCACACCTGTGGTACTTACCCAGGGATGGTTCACCACCCCTCCCTTCCCTCTCATCATCGTTGGGGAACGGTGACAATGTTTGGAACAATTTTTGGTTGTCACAAACAGGGGTTTCTTCTGATATTCAATGAGTAGAAGCCAGGGACACTGCTAGAGAACCCACAATGTTCAGAACAGCCTCTGCCATCAACAAGGAATTATCTGGTCCAAAATGTCAATAGTGCTGAGGCTAAGAGCACTGGTTCACACTGTGCTCTTTCTGAAAATTCTAGACTCACATCTGTTATACACTCACCACACAGTTTAGTCTTTTATTTTTGCTTGTTTCATTATAAACAATTAGACAGTTGTATAAATTCAACCACTTTCTTGTTGAATCCATTTAGTCAATGCAAGCTCAATATTTTCATATTTATTTTTTGCCTTATGCAATATTTTTCAACATTTTCATGAGTTGTCGGTCATCACTATCTCTATTAACTTTCAACAACTTGCCCTTGTAAGTCACAAATAGTGCTGCTGCTGAAATTATTTCTCACTAACATGCCTCAGATTTCTGTAGTGATTCTACATTTGATATTATTCACAATGTAAAATGCTTCTATTTATTCATTTCGCTTTTACCCAAGGATTATTTTTAAGTTATTTTTGTCATTTTCACACTTCAAACATAAAGACAAAAACATCAAAAATATAGTGTTTTACATATGTGCATATTTTCACACATATATGTATGTATATTTATATGTATTGAAACTACAGAAGCACATGTCACCAATAAGAGCTCTGAGACACCTTCGACCACTTACCCTTATCAGATGAGTTGTGGAAACAAGTTTTTTTTAACTGAATTTCTGAGCTTTGTGGATTTAGAAATGCAAAGGAAGGTTTGTGGACATTCACAGGGATCATGATTTTATTCTCCTTAAAACTCTTCTGTACTTTCCAATTGTCCTTAGTATAAATCCAAAATCCTAACACCACCCAAGAGGCTTTTCAATACCTGGCTCCTGTGATTTCTCCTGTCTAATCTTTTACCCTCCTTCCCCTCAGCCTCTCTGCTTTAGTGAACTTTCTCCTAGTTTTTTGAAGAAGTTCATCAATTCAAGCTTTTGTACATGGGATTTCCTAAACCTGAAATGTGCCTCCCGTTTTGTCCAAACAGACACGGGCTCCACTCTGCCCCCTGGCTCACACCTGCTTAACCTGTCAAGTCACATCTGAACCGTCACTATTAAGAGGGTCCTTCTCTGGCACCCTAATGTAATTGAGATCATCCTATTATTCTCTGTTCTAGAACTCCACACTTCCGACATTTCTCATTCCTGTCTAAGCTCTTGTGTGTTTGGTGTTGGGCCATCACTTTCACTGCTCTTTAAGCTCCCCCAGCGGAGTGGAGAGGTCTGTTTTCCCTCGTTTGGATTCCTAGAGGCAGCGCAGACCAGGCACAAGGTCAGCACTAAGGAAGGGTTCACAGGATGAACGCGGTGGGTGCTGTTTAAGGAACCGGTAAACGTGTGGGATGAGAGAAGGAGCAGAGTGTCTTTGGGGTGGAGGCTCCCAGGAGGAGGCGGCGCGGGCTGCGGTGCGGGGCGGATCCTCCTCCAGCTCCTGCTTGGAGGTCTCCAGAACAGGCTGGAGGTAGGGAGGGGGGTCCCAAAAGCCTGGGGATCAGACGTGGTTTTCCCGCCTGGTCCCCCAGGCCCCCTTTCGCCTCAGGAAGACAGAGGAGGAGCCCCTGGGCTGCAGGTGGTGGGCGTTGCGGCGGGGGCCGGTTAAGGTTCCCAGTGCCCGCACCCGGCCCAGGGAGCCCCGGATGGCGGCGTCACTGTCAGTGTCTTCTCAGGAGGCCGCCTGTGTGACTGGATCGTTCGTGTCCCCACAGCACGTTTCTTGGAGTACTCTACGTCTGAGTGTCATTTCTTCAATGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCTATAACCAAGAGGAGTACGTGCGCTTCGACAGCGACGTGGGGGAGTTCCGGGCGGTGACGGAGCTGGGGCGGCCTGATGAGGAGTACTGGAACAGCCAGAAGGACATCCTGGAAGACGAGCGGGCCGCGGTGGACACCTACTGCAGACACAACTACGGGGTTGTGGAGAGCTTCACAGTGCAGCGGCGAGGTGAGCGCGGCGCGGGGCGGGGCCTGAGTCCCTGTGAGCTGGGAATCTGAGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGAGAGAGAGACAGAGAGACAGAGAGAGAGAGAGCGCCATCTGTGAGCATTTAGAATCCTCTCTATCCTGAGCAAGGAGTTCTGAGGGCACAGGTGTGTGTGTAGAGTGTGGATTTGTCTGTGTCTGTGAGGCTGTTGTGGGAGGGGAGGCAGGAGGGGGCTGCTTCTTATTCTTGGAGGACTCTGTGGGGAGGTGACAAGGGAGGTGGGTGCGGGCGGCTGGAGAGAGAGGTGACCTTGATTGTCTCGGGTCCTTAGAGATGCAGGGAAGGGAAATGTAAGGGGTGTGTGGTTGGGGTGAAGGTTTAGGGGAGGAGAGCTGAGGGGTAAGGAAGGTTTGGGATAATGTGAGGAGGCCAGTTCCAGACTGTCCCTGGCACACACCCTTCATGTAATCTCTGAAATAAAAGTGTGTGCTGTTTGTTTGTAAAAGCATTAGATTAATTTCTAGGGGAATTGAGGAGACCTCTGAGGCATCTCTGAAGCTTCTTTAGGTCTAAATTTCTTGCTAGTTTTTTGTTTTTTATTGTGTATATTTTTACATAGTAGAAATGACTGTGAAACTAACTTTTTGAATTAAAGTTTTAACACAGTTACTATTTTATTATAATGCTAATAGTTTTCTAGTAGTTACATATTATTCTTTTATATATAATAGTTGTGACACAACTTACCTCACTTTCCCCTTTGTTGACCTTTATTATGACATTCACCAAAATTTGAAAATGTATGTTTCTGGTTAATTTTTAATTTATATTTTTTTCATTTATAATTCTTTTGAATTATTTTGACCTATTTATTGGCCAGTTTTAATAACTGCTGTAAGAATTCCCTATTGTATTTGGTAGGGAATGGACAATGATCTACTGCCTAATATCTCGAGGGCTTAGTATTTTTCTCAGTGACTTTGTGGGTTCTTTGTACTGTGAGATTATTAACACTTTATTGATATTTGATTCAGCATTTGCTCCAGTTTGTGGTTTGTATGTTGATTTTGAAAATTCTTTTCCATGTTAAGAATTTGAACATTTTTATATAATAAAATATGTTGCAAAATTTTTATTAATGATTTACAATCCATCTTAAATCTGCCATTTTGTGGTATTGTTGTCTCCAGGTTTCTCCTTACTTCTAAAAAAAATTGCATTTATTGAGAGTCTGCTAGTGTTAGGGATTTTCCTGGGCATAAGCACCCCAAGTGACGAGTCCCAGACACTGCCTTAATCCAAATGTGATTCTGGAAAGAAAAATCATTTTACAATGATAGGCCTAATAATAATTAAGCTTGTGTTGCATGGGAGATGCATTGATCAGCTAAATGTAAATATAAGAACTTTCAAAACTAAAATGACGTTCCTTAATCCTTCTCTCTGCTTTATGACTCATGCTTTTCTGGGAAAGTAAAAATTTGGAGAATCATTTCTGTCTGTCCCACCTTCCCAGGGGCAGAACCATTTCTGTGGTGTTCTAAGGTGTGAGTGCATGGCGGTAGTATTCCTAAAAATTCATATTCGGTTTCGTCATGTACCCAACTCTGTCCCGTTATCTATCAACATTGTTTTAAATCATATATTTCTGTCAAGGTGTACAAGGATGATAAATAGGTGCCAAGTGGAGCACCCAAGTGTGATGAGCCCCCTCACAGTGGAATGGAGTGTGAAGCTTTATGACCTCATAAATTGAAGGTTATCTTCAGTCATTGTTTTATATATTTTACATGCATTAATCCTCATATAATCCCAAGAGGTAAATTAGTATAATTATCCTTCATTATAGGTGACAAAGTTGAGACACAGAAGAATCGAACTCTTAAGGCAGACCTTGGATTTGAACCAGGCAACCTGGCTCAGATATCAGTTTTAATTACTACACTCTGTACTTTCAAAGATTTGTAAACACTTTGACAATGCATGACAATTTCAAGCTATGAAGAAACAAACACAATTTTTCACAATATCTCTCAAATCTAATAGGTCCTCACTATCAAGATTAAGTTCCAGGCTGATGACACTGTAAGGCCACATGGCCAGCTGTGCTGGAGGCCTGGTCAAGGTCAGAGCCTGGGTTTGCAGAGAAGCAGACAAACAGCCAAACAAGGAGACTTACTCTGTCTTCATGACTCATTCCCTCTACCTTTTTTCTCCTAGTCCATCCTAAGGTGACTGTGTATCCTTCAAAGACCCAGCCCGTGCAGCACCACAACCTCCTGGTCTGTTCTGTGAGTGGTTTCTATCCAGGCAGCATTGAAGTCAGGTGGTTCCGGAATGGCCAGGAAGAGAAGACTGGGGTGGTGTCCACAGGCCTGATCCACAATGGAGACTGGACCTTCCAGACCCTGGTGATGCTGGAAACAGTTCCTCGGAGTGGAGAGGTTTACACCTGCCAAGTGGAGCACCCAAGCGTGACAAGCCCTCTCACAGTGGAATGGAGTGAGCAGCTTTCTGACTTCATAAATTTCTCACCCACCAAGAAGGGGACTGTGCTCATCCCTGAGTGTCAGGTTTCTCCTCTCCGACATCCTATTTTCATTTGCTCCATGTTCTCATCTCCATCAGCACAGGTCACTGGGGGTAGCCCTGTAGGTGTTTCTAGAAACACCTGTACCTCCTGGAGAAGCAGTCTCGCCTGCCAGGCAGGAGAGGCTGTCCCTCTTTTGAACCTCCCCATGATGTCACAGGTCAGGGTCACCCACCCTCCCCGGGCTCCAGGCACTGCCTCTGGGTCTGAGACTGAGTTTCTGGTGCTGTTGATCTGAGTTATTTGTTGTGATCTGGGAAGAGGAGAAGTGTAGGGGCCTTCCTGACATGAGGGGAGTCCAATCTCAGCTCTGCCTTTTATTAGCTCTGTCACTCTAGACAAACTACTTAGCCTCATTGAGTCTCAGGCTTTCTGTGGATCAGATGTTGAACTCTTGCCTTACATCAAGGCTGTAATATTTGAATGAGTTTGATGTCTGAACCTTGTAACTGTTCAGTGTGATTTGAAATCCTTTTTTTCTCCAGAAATGGCTAGTTATTTTAGTTCTTGTGGGGCAGACTTCTTCCCCATTTTCAAAGCTCTGAATCTTAGAGTCTCAATTAAAGAGGTTCAATTTGGAATAAACACTAAACCTGGCTTCCTCTCTCAGGAGCACGGTCTGAATCTGCACAGAGCAAGATGCTGAGTGGAGTCGGGGGCTTTGTGCTGGGCCTGCTCTTCCTTGGGGCCGGGCTGTTCATCTACTTCAGGAATCAGAAAGGTGAGGAGCCTTTGGTAGCTGGCTCTCTCCATAGGCTTTTCTGGAGGAGGAACTATGGCTTTGCTGAGGTTAGTTCTCAGTATATGAGTGGCCCTGAATAAAGCCTTTCTTTCCCCAAACGGCTCTAATGTCCTGCTAATCCAGAAATCATCAGTGCATGGTTACTATGTGAAAGCATAATAGCTTGTGGCCTGCAGAGACAAGAGGAAGGTTAACAAGTAGGGGTCCTTTGGTTTGAGATCTTGGAGCAGATTAAGGAAGAGCCACTAAGACTAATGGAATTACACTGGATCCTGTGACAGACACTTCACCCTTCATGGGTCACATGGTCTGTTTCTGCTCCTCTCTGCCCTGGCTGGTGTGGGTTGTAGTGACAGAGAACTCTCCGGTGGGAGATCTGGGGCTGGGACATTGTGTTGGAAGACAGATTTGCTTCCATAAATTTTAAGTGTATATATTTTCCTCTTTTTCCCAGGACACTCTGGACTTCAGCCAAGAGGTAATACCTTTTAATCCTCTTTTAGAAACAGATACGGTTTCCCTAGTGAGAGGTGAAGCCAGCTGGACTTCTGGGTCGGGTAGGGACTTGCAGAACTTTCCTGTCTTAGGAGAGGTTTCTAAATGCACCAATCAGTGCTCTGTAAAAACACACCAATTGGCACTCTGTGGCTAGATAGATGTTTGTAAAATGGACTAATCAGCACTCTGTAAAATGGAGCAATCCACACTCTGTAAAATGGACCAATCAATGCTCTTTAAAATGGACCAATCAGCAGGACATGGGCGGGGACAAATAAGGGAATACAAGCTGGCCACCCCAGCCAGCAGCAGCAACCCGCTCAGGTCGCCTTCCATGCTGTGGAAGCTTTGTTCTTTTGCTCTTCACAATAAATCTTGCTGTTGCTCACTCTTCGGGTCTGTGCCACCTTTAAGAGCTGTAACACTCACTGTGAAGATTCGCGGCTTCATTCTTGAAGTCAGCGAAACCACGAACCCACCGGAAGGAACAAACTCTGGACACACTAGAATTGATGGTAGAGGTGATAAGGCATGAGACAGAAATAATAGGAAAGACTTTGGATCCAAATTTCTGATCAGGCAATTTACACCAAAACTCCTCCTCTCCACTTAGAAAAGGCCTGTGCTCTGTGGGACTATTGGCTCTGGGAGACTCAGGAACTTGTTTTTCTTCTTCCTGCAGTGCTCTCATCTGAGTCCCTGAAAGAGAGGAAAAGAAACTGTTAGTAGAGTCAGGTTGAAAACAACACTCTCCTCTGTCTTTTGCAGGATTCCTGAGCTGAAGTGCAGATGACACATTCAAAGAAGAACTTTCTGCCCCAGCTTTGCAGGATGAAAAGCTTTCCCTCCTGGCTGT'
#     hlaLocus = 'HLA-DRB1'
#     # roughFeatureSequence = 'aag\nCGTCGT\nccg\nGGCTGA\naat'
#     alleleCallWithGFE = fetchSequenceAlleleCallWithGFE(roughFeatureSequence, hlaLocus)
#     assert_true(len(alleleCallWithGFE) > 3)
#
#     # print 'ALLELE CALL RESULTS:'
#     # print alleleCallWithGFE
#
#     annotatedSequence = parseExons(roughFeatureSequence, alleleCallWithGFE)
#     # print ('ANNOTATED SEQUENCE:')
#     # print (annotatedSequence)
#     assert_true(len(annotatedSequence) > 3)
#
#     ipdGenerator = IpdSubGenerator()
#
#     ipdGenerator.sequenceAnnotation = identifyGenomicFeatures(annotatedSequence)
#
#     assignConfigurationValue('ipd_submission_identifier', 'HWS100234567')
#     assignConfigurationValue('ipd_submission_version', '1')
#     # I think I dont need this ena release date.
#     # assignConfigurationValue('ena_release_date', '09/08/2017')
#     assignConfigurationValue('allele_name', 'HLA-DRB1_NEW_1')
#     # I need a method to do this...I am adding it now.
#     # assignConfigurationValue('closest_allele_written_description', 'TEMP_IPD_IDENTIFIER')
#     assignConfigurationValue('ena_sequence_accession', 'ENA_Seq_Acc')
#     assignConfigurationValue('ipd_submitter_id', 'IPD_Submitter_ID')
#     assignConfigurationValue('ipd_submitter_name', 'Ben Bioinformaticist')
#     assignConfigurationValue('ipd_alt_contact', 'IPD_Alt_Contact')
#     assignConfigurationValue('ipd_submitter_email', 'Ben@BioinformaticsResearch.com')
#     assignConfigurationValue('cell_id', 'Cell_Identifier123')
#     assignConfigurationValue('ethnic_origin', 'African American')
#     assignConfigurationValue('sex', 'Unknown')
#     assignConfigurationValue('consanguineous', 'Unknown')
#     assignConfigurationValue('homozygous', 'Unknown')
#     assignConfigurationValue('lab_of_origin', 'Maastricht University Medical Center')
#     assignConfigurationValue('lab_contact', 'Marspel JR Spilanus')
#     assignConfigurationValue('alleles', 'HLA-A*01:01:01:01;HLA-B*40:72:01')
#     assignConfigurationValue('material_availability', 'No Material Available')
#     assignConfigurationValue('cell_bank', 'Not Available')
#     assignConfigurationValue('primary_sequencing_methodology', 'Direct sequencing of PCR product from DNA (SBT)')
#     assignConfigurationValue('secondary_sequencing_methodology', 'Direct sequencing of PCR product from DNA (SBT)')
#     assignConfigurationValue('primer_type', 'Both allele and locus specific')
#     assignConfigurationValue('primers',
#                              '03PID03     CCCAAAGGGTTTCCCGGGAAATTT 3UT 3015-3042;04PID04     AAAGGGTTTCCCGGGAAATTTCCC 5UT 4015-4042')
#     assignConfigurationValue('sequenced_in_isolation', 'Yes')
#     assignConfigurationValue('no_of_reactions', '3')
#     assignConfigurationValue('method_comments', 'The method we used was MinION sequencing!')
#     # assignConfigurationValue('closest_known_allele', 'TEMP_IPD_IDENTIFIER')
#     # I'm specifically escapting semicolons here.
#     # Ben's unnoficial dirty trick: Convert any semi colon on the input to ;. Just store those 3 characters instead of a semicolon.
#     assignConfigurationValue('citations',
#                              '1-3503;PUBMED; 9349617.;Laforet M, Froelich N, Parissiadis A, Pfeiffer B, Schell A, Faller B, Woehl-Jaegle ML, Cazenave JP, Tongio MM;A nucleotide insertion in exon 4 is responsible for the absence of expression of an HLA-A*01 allele;Tissue Antigens 50:347 - 50(1997).;1-3503;PUBMED; 9349617.;Laforet M, Froelich N, Parissiadis A, Pfeiffer B, Schell A, Faller B, Woehl-Jaegle ML, Cazenave JP, Tongio MM;A nucleotide insertion in exon 4 is responsible for the absence of expression of an HLA-A*01 allele;Tissue Antigens 50:347 - 50(1997).')
#
#     ipdSubmission = ipdGenerator.buildIPDSubmission()
#
#     print('IPD SUBMISSION:\n' + ipdSubmission)
#
#     assert_true(len(ipdSubmission) > 3)
#     assert_true(ipdSubmission is not None)
#
#     # Maybe I don't want to safe my configuration with all these test values but maybe I don't care :D
#     writeConfigurationFile()
#
#     # Print out the submission to a file, for easier validation.
#     outputFileObject = open('Test_IPD_Submission_HLA_DRB1_SNP.txt', 'w')
#     submissionText = ipdSubmission
#     outputFileObject.write(submissionText)
#
# def testCreateENASubmissionFlatfile():
#     print ('Test: Create an ENA SubmissionFlatfile')
#     assert_true(True)
#
#     roughFeatureSequence = 'aag\nCGTCGT\nccg\nGGCTGA\naat'
#
#     assignConfigurationValue('sample_id', 'Donor_12345')
#     assignConfigurationValue('gene', 'HLA-C')
#     assignConfigurationValue('class', '1')
#     assignConfigurationValue('allele_name', 'Allele:01:02')
#
#     allGen = EnaSubGenerator()
#     # roughFeatureSequence = self.featureInputGuiObject.get('1.0', 'end')
#
#     allGen.sequenceAnnotation = identifyGenomicFeatures(roughFeatureSequence)
#
#     enaSubmission = allGen.buildENASubmission()
#
#     assert_true(len(enaSubmission) > 3)
#     assert_true(enaSubmission is not None)
#
#     #Print out the submission to a file, for easier validation.
#     outputFileObject = open('Test_ENA_Submission.txt', 'w')
#     submissionText = enaSubmission
#     outputFileObject.write(submissionText)
#
# def testAnnotateSimple():
#     # This test demonstrates a minimal snippet to annotate a raw sequence.
#     import gfe_client
#     #A*01:01:01:01
#     sequenceRaw = 'CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTGGCTCTCAGGGTCTCAGGCCCCGAAGGCGGTGTATGGATTGGGGAGTCCCAGCCTTGGGGATTCCCCAACTCCGCAGTTTCTTTTCTCCCTCTCCCAACCTACGTAGGGTCCTTCATCCTGGATACTCACGACGCGGACCCAGTTCTCACTCCCATTGGGTGTCGGGTTTCCAGAGAAGCCAATCAGTGTCGTCGCGGTCGCTGTTCTAAAGTCCGCACGCACCCACCGGGACTCAGATTCTCCCCAGACGCCGAGGATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGTGAGTGCGGGGTCGGGAGGGAAACCGCCTCTGCGGGGAGAAGCAAGGGGCCCTCCTGGCGGGGGCGCAGGACCGGGGGAGCCGCGCCGGGAGGAGGGTCGGGCAGGTCTCAGCCACTGCTCGCCCCCAGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTGAGTGACCCCGGCCCGGGGCGCAGGTCACGACCCCTCATCCCCCACGGACGGGCCAGGTCGCCCACAGTCTCCGGGTCCGAGATCCACCCCGAAGCCGCGGGACTCCGAGACCCTTGTCCCGGGAGAGGCCCAGGCGCCTTTACCCGGTTTCATTTTCAGTTTAGGCCAAAAATCCCCCCGGGTTGGTCGGGGCGGGGCGGGGCTCGGGGGACTGGGCTGACCGCGGGGTCGGGGCCAGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCCGCGGGTACCGGCAGGACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGATCACCAAGCGCAAGTGGGAGGCGGTCCATGCGGCGGAGCAGCGGAGAGTCTACCTGGAGGGCCGGTGCGTGGACGGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGGGTACCAGGGGCCACGGGGCGCCTCCCTGATCGCCTATAGATCTCCCGGGCTGGCCTCCCACAAGGAGGGGAGACAATTGGGACCAACACTAGAATATCACCCTCCCTCTGGTCCTGAGGGAGAGGAATCCTCCTGGGTTTCCAGATCCTGTACCAGAGAGTGACTCTGAGGTTCCGCCCTGCTCTCTGACACAATTAAGGGATAAAATCTCTGAAGGAGTGACGGGAAGACGATCCCTCGAATACTGATGAGTGGTTCCCTTTGACACCGGCAGCAGCCTTGGGCCCGTGACTTTTCCTCTCAGGCCTTGTTCTCTGCTTCACACTCAATGTGTGTGGGGGTCTGAGTCCAGCACTTCTGAGTCTCTCAGCCTCCACTCAGGTCAGGACCAGAAGTCGCTGTTCCCTTCTCAGGGAATAGAAGATTATCCCAGGTGCCTGTGTCCAGGCTGGTGTCTGGGTTCTGTGCTCTCTTCCCCATCCCGGGTGTCCTGTCCATTCTCAAGATGGCCACATGCGTGCTGGTGGAGTGTCCCATGACAGATGCAAAATGCCTGAATTTTCTGACTCTTCCCGTCAGACCCCCCCAAGACACATATGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTCTGCCCAAGCCCCTCACCCTGAGATGGGGTAAGGAGGGAGATGGGGGTGTCATGTCTCTTAGGGAAAGCAGGAGCCTCTCTGGAGACCTTTAGCAGGGTCAGGGCCCCTCACCTTCCCCTCTTTTCCCAGAGCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCCGTGATGTGGAGGAGGAAGAGCTCAGGTGGAGAAGGGGTGAAGGGTGGGGTCTGAGATTTCTTGTCTCACTGAGGGTTCCAAGCCCCAGCTAGAAATGTGCCCTGTCTCATTACTGGGAAGCACCTTCCACAATCATGGGCCGACCCAGCCTGGGCCCTGTGTGCCAGCACTTACTCTTTTGTAAAGCACCTGTTAAAATGAAGGACAGATTTATCACCTTGATTACGGCGGTGATGGGACCTGATCCCAGCAGTCACAAGTCACAGGGGAAGGTCCCTGAGGACAGACCTCAGGAGGGCTATTGGTCCAGGACCCACACCTGCTTTCTTCATGTTTCCTGATCCCGCCCTGGGTCTGCAGTCACACATTTCTGGAAACTTCTCTGGGGTCCAAGACTAGGAGGTTCCTCTAGGACCTTAAGGCCCTGGCTCCTTTCTGGTATCTCACAGGACATTTTCTTCCCACAGATAGAAAAGGAGGGAGTTACACTCAGGCTGCAAGTAAGTATGAAGGAGGCTGATGCCTGAGGTCCTTGGGATATTGTGTTTGGGAGCCCATGGGGGAGCTCACCCACCCCACAATTCCTCCTCTAGCCACATCTTCTGTGGGATCTGACCAGGTTCTGTTTTTGTTCTACCCCAGGCAGTGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGTAAAGGTGAGAGCTTGGAGGGCCTGATGTGTGTTGGGTGTTGGGTGGAACAGTGGACACAGCTGTGCTATGGGGTTTCTTTGCGTTGGATGTATTGAGCATGCGATGGGCTGTTTAAGGTGTGACCCCTCACTGTGATGGATATGAATTTGTTCATGAATATTTTTTTCTATAGTGTGAGACAGCTGCCTTGTGTGGGACTGAGAGGCAAGAGTTGTTCCTGCCCTTCCCTTTGTGACTTGAAGAACCCTGACTTTGTTTCTGCAAAGGCACCTGCATGTGTCTGTGTTCGTGTAGGCATAATGTGAGGAGGTGGGGAGAGCACCCCACCCCCATGTCCACCATGACCCTCTTCCCACGCTGACCTGTGCTCCCTCTCCAATCATCTTTCCTGTTCCAGAGAGGTGGGGCTGAGGTGTCTCCATCTCTGTCTCAACTTCATGGTGCACTGAGCTGTAACTTCTTCCTTCCCTATTAAAA'
#     config = gfe_client.Configuration()
#     config.host = 'http://act.b12x.org'
#
#
#     api = gfe_client.ApiClient(configuration=config)
#     ann_api = gfe_client.TypeSeqApi(api_client=api)
#     responseText = ann_api.typeseq_get(sequence=sequenceRaw, imgthla_version="3.31.0", locus='HLA-A')
#
#     # Convert the python objects to more proper json.
#     annotation = responseText.to_dict()
#     jsonResponse = dumps(annotation, indent=4)
#
#     #print('RESPONSE:\n' + str(jsonResponse))
#     assert(len(str(jsonResponse)) > 5)
#
#
#
# def testAnnotateRoughSequence():
#     print ('Test: Annotate a rough sequence using the ACT APIs')
#
#     #>HLA-A*02:01:01:12 Full Length Allele. Will we allele call correctly?
#     #roughFeatureSequence = 'CCAGTTCTCACTCCCATTGGGTGTCGGGTTTCCAGAGAAGCCAATCAGTGTCGTCGCGGTCGCGGTTCTAAAGTCCGCACGCACCCACCGGGACTCAGATTCTCCCCAGACGCCGAGGATGGCCGTCATGGCGCCCCGAACCCTCGTCCTGCTACTCTCGGGGGCTCTGGCCCTGACCCAGACCTGGGCGGGTGAGTGCGGGGTCGGGAGGGAAACGGCCTCTGTGGGGAGAAGCAACGGGCCCGCCTGGCGGGGGCGCAGGACCCGGGAAGCCGCGCCGGGAGGAGGGTCGGGCGGGTCTCAGCCACTCCTCGTCCCCAGGCTCTCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGTCCGGAGTATTGGGACGGGGAGACACGGAAAGTGAAGGCCCACTCACAGACTCACCGAGTGGACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCGGTGAGTGACCCCGGCCCGGGGCGCAGGTCACGACCTCTCATCCCCCACGGACGGGCCAGGTCGCCCACAGTCTCCGGGTCCGAGATCCGCCCCGAAGCCGCGGGACCCCGAGACCCTTGCCCCGGGAGAGGCCCAGGCGCCTTTACCCGGTTTCATTTTCAGTTTAGGCCAAAAATCCCCCCAGGTTGGTCGGGGCGGGGCGGGGCTCGGGGGACCGGGCTGACCGCGGGGTCCGGGCCAGGTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGGGTACCAGGGGCCACGGGGCGCCTCCCTGATCGCCTGTAGATCTCCCGGGCTGGCCTCCCACAAGGAGGGGAGACAATTGGGACCAACACTAGAATATCGCCCTCCCTCTGGTCCTGAGGGAGAGGAATCCTCCTGGGTTTCCAGATCCTGTACCAGAGAGTGACTCTGAGGTTCCGCCCTGCTCTCTGACACAATTAAGGGATAAAATCTCTGAAGGAATGACGGGAAGACGATCCCTCGAATACTGATGAGTGGTTCCCTTTGACACACACAGGCAGCAGCCTTGGGCCCGTGACTTTTCCTCTCAGGCCTTGTTCTCTGCTTCACACTCAATGTGTGTGGGGGTCTGAGTCCAGCACTTCTGAGTCCTTCAGCCTCCACTCAGGTCAGGACCAGAAGTCGCTGTTCCCTCTTCAGGGACTAGAATTTTCCACGGAATAGGAGATTATCCCAGGTGCCTGTGTCCAGGCTGGTGTCTGGGTTCTGTGCTCCCTTCCCCATCCCAGGTGTCCTGTCCATTCTCAAGATAGCCACATGTGTGCTGGAGGAGTGTCCCATGACAGATGCAAAATGCCTGAATGATCTGACTCTTCCTGACAGACGCCCCCAAAACGCATATGACTCACCACGCTGTCTCTGACCATGAAGCCACCCTGAGGTGCTGGGCCCTGAGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGACAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTTTGCCCAAGCCCCTCACCCTGAGATGGGGTAAGGAGGGAGACGGGGGTGTCATGTCTTTTAGGGAAAGCAGGAGCCTCTCTGACCTTTAGCAGGGTCAGGGCCCCTCACCTTCCCCTCTTTTCCCAGAGCCGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCTTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCTGTGATGTGGAGGAGGAAGAGCTCAGGTGGGGAAGGGGTGAAGGGTGGGTCTGAGATTTCTTGTCTCACTGAGGGTTCCAAGACCCAGGTAGAAGTGTGCCCTGCCTCGTTACTGGGAAGCACCACCCACAATTATGGGCCTACCCAGCCTGGGCCCTGTGTGCCAGCACTTACTCTTTTGTAAAGCACCTGTTAAAATGAAGGACAGATTTATCACCTTGATTACAGCGGTGATGGGACCTGATCCCAGCAGTCACAAGTCACAGGGGAAGGTCCCTGAGGACCTTCAGGAGGGCGGTTGGTCCAGGACCCACACCTGCTTTCTTCATGTTTCCTGATCCCGCCCTGGGTCTGCAGTCACACATTTCTGGAAACTTCTCTGAGGTCCAAGACTTGGAGGTTCCTCTAGGACCTTAAGGCCCTGACTCCTTTCTGGTATCTCACAGGACATTTTCTTCCCACAGATAGAAAAGGAGGGAGCTACTCTCAGGCTGCAAGTAAGTATGAAGGAGGCTGATGCCTGAGGTCCTTGGGATATTGTGTTTGGGAGCCCATGGGGGAGCTCACCCACCCCACAATTCCTCCTCTAGCCACATCTTCTGTGGGATCTGACCAGGTTCTGTTTTTGTTCTACCCCAGGCAGTGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGTAAAGGTGAGAGCCTGGAGGGCCTGATGTGTGTTGGGTGTTGGGCGGAACAGTGGACACAGCTGTGCTGTGGGGTTTCTTTCCATTGGATGTATTGAGCATGCGATGGGCTGTTTAAAGTGTGACCCCTCACTGTGACAGATACGAATTTGTTCATGAATATTTTTTTCTATAGTGTGAGACAGCTGCCTTGTGTGGGACTGAGAGGCAAGAGTTGTTCCTGCCCTTCCCTTTGTGACTTGAAGAACCCTGACTTTGTTTCTGCAAAGGCACCTGCATGTGTCTGTGTTCGTGTAGGCATAATGTGAGGAGGTGGGGAGACCACCCCACCCCCATGTCCACCATGACCCTCTTCCCACGCTGACCTGTGCTCCCTCCCCAATCATCTTTCCTGTTCCAGAGAGGTGGGGCTGAGGTGTCTCCATCTCTGTCTCAACTTCATGGTGCACTGAGCTGTAACTTCTTCCTTCCCTATTAAAA'
#
#     import gfe_client
#     # This one is HLA-A 01:01:01, should work fine.
#     roughFeatureSequence = 'CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTGGCTCTCAGGGTCTCAGGCCCCGAAGGCGGTGTATGGATTGGGGAGTCCCAGCCTTGGGGATTCCCCAACTCCGCAGTTTCTTTTCTCCCTCTCCCAACCTACGTAGGGTCCTTCATCCTGGATACTCACGACGCGGACCCAGTTCTCACTCCCATTGGGTGTCGGGTTTCCAGAGAAGCCAATCAGTGTCGTCGCGGTCGCTGTTCTAAAGTCCGCACGCACCCACCGGGACTCAGATTCTCCCCAGACGCCGAGGATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGTGAGTGCGGGGTCGGGAGGGAAACCGCCTCTGCGGGGAGAAGCAAGGGGCCCTCCTGGCGGGGGCGCAGGACCGGGGGAGCCGCGCCGGGAGGAGGGTCGGGCAGGTCTCAGCCACTGCTCGCCCCCAGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTGAGTGACCCCGGCCCGGGGCGCAGGTCACGACCCCTCATCCCCCACGGACGGGCCAGGTCGCCCACAGTCTCCGGGTCCGAGATCCACCCCGAAGCCGCGGGACTCCGAGACCCTTGTCCCGGGAGAGGCCCAGGCGCCTTTACCCGGTTTCATTTTCAGTTTAGGCCAAAAATCCCCCCGGGTTGGTCGGGGCGGGGCGGGGCTCGGGGGACTGGGCTGACCGCGGGGTCGGGGCCAGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCCGCGGGTACCGGCAGGACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGATCACCAAGCGCAAGTGGGAGGCGGTCCATGCGGCGGAGCAGCGGAGAGTCTACCTGGAGGGCCGGTGCGTGGACGGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGGGTACCAGGGGCCACGGGGCGCCTCCCTGATCGCCTATAGATCTCCCGGGCTGGCCTCCCACAAGGAGGGGAGACAATTGGGACCAACACTAGAATATCACCCTCCCTCTGGTCCTGAGGGAGAGGAATCCTCCTGGGTTTCCAGATCCTGTACCAGAGAGTGACTCTGAGGTTCCGCCCTGCTCTCTGACACAATTAAGGGATAAAATCTCTGAAGGAGTGACGGGAAGACGATCCCTCGAATACTGATGAGTGGTTCCCTTTGACACCGGCAGCAGCCTTGGGCCCGTGACTTTTCCTCTCAGGCCTTGTTCTCTGCTTCACACTCAATGTGTGTGGGGGTCTGAGTCCAGCACTTCTGAGTCTCTCAGCCTCCACTCAGGTCAGGACCAGAAGTCGCTGTTCCCTTCTCAGGGAATAGAAGATTATCCCAGGTGCCTGTGTCCAGGCTGGTGTCTGGGTTCTGTGCTCTCTTCCCCATCCCGGGTGTCCTGTCCATTCTCAAGATGGCCACATGCGTGCTGGTGGAGTGTCCCATGACAGATGCAAAATGCCTGAATTTTCTGACTCTTCCCGTCAGACCCCCCCAAGACACATATGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTCTGCCCAAGCCCCTCACCCTGAGATGGGGTAAGGAGGGAGATGGGGGTGTCATGTCTCTTAGGGAAAGCAGGAGCCTCTCTGGAGACCTTTAGCAGGGTCAGGGCCCCTCACCTTCCCCTCTTTTCCCAGAGCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCCGTGATGTGGAGGAGGAAGAGCTCAGGTGGAGAAGGGGTGAAGGGTGGGGTCTGAGATTTCTTGTCTCACTGAGGGTTCCAAGCCCCAGCTAGAAATGTGCCCTGTCTCATTACTGGGAAGCACCTTCCACAATCATGGGCCGACCCAGCCTGGGCCCTGTGTGCCAGCACTTACTCTTTTGTAAAGCACCTGTTAAAATGAAGGACAGATTTATCACCTTGATTACGGCGGTGATGGGACCTGATCCCAGCAGTCACAAGTCACAGGGGAAGGTCCCTGAGGACAGACCTCAGGAGGGCTATTGGTCCAGGACCCACACCTGCTTTCTTCATGTTTCCTGATCCCGCCCTGGGTCTGCAGTCACACATTTCTGGAAACTTCTCTGGGGTCCAAGACTAGGAGGTTCCTCTAGGACCTTAAGGCCCTGGCTCCTTTCTGGTATCTCACAGGACATTTTCTTCCCACAGATAGAAAAGGAGGGAGTTACACTCAGGCTGCAAGTAAGTATGAAGGAGGCTGATGCCTGAGGTCCTTGGGATATTGTGTTTGGGAGCCCATGGGGGAGCTCACCCACCCCACAATTCCTCCTCTAGCCACATCTTCTGTGGGATCTGACCAGGTTCTGTTTTTGTTCTACCCCAGGCAGTGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGTAAAGGTGAGAGCTTGGAGGGCCTGATGTGTGTTGGGTGTTGGGTGGAACAGTGGACACAGCTGTGCTATGGGGTTTCTTTGCGTTGGATGTATTGAGCATGCGATGGGCTGTTTAAGGTGTGACCCCTCACTGTGATGGATATGAATTTGTTCATGAATATTTTTTTCTATAGTGTGAGACAGCTGCCTTGTGTGGGACTGAGAGGCAAGAGTTGTTCCTGCCCTTCCCTTTGTGACTTGAAGAACCCTGACTTTGTTTCTGCAAAGGCACCTGCATGTGTCTGTGTTCGTGTAGGCATAATGTGAGGAGGTGGGGAGAGCACCCCACCCCCATGTCCACCATGACCCTCTTCCCACGCTGACCTGTGCTCCCTCTCCAATCATCTTTCCTGTTCCAGAGAGGTGGGGCTGAGGTGTCTCCATCTCTGTCTCAACTTCATGGTGCACTGAGCTGTAACTTCTTCCTTCCCTATTAAAA'
#
#     #Fairly sure I do not need the HLA locus, but there it is.
#     hlaLocus = 'HLA-A'
#     alleleCallWithGFE = fetchSequenceAlleleCallWithGFE(roughFeatureSequence, hlaLocus)
#
#
#     print('I got an allele call with gfe, of the type:' + str(type(alleleCallWithGFE)))
#     #print('the sequence from the Annotation object:' + str(alleleCallWithGFE))
#
#
#     #assert_true(len(alleleCallWithGFE) > 3)
#     #print 'ALLELE CALL RESULTS:'
#     #print alleleCallWithGFE
#
#     annotatedSequence = parseExons(roughFeatureSequence, alleleCallWithGFE)
#     #print ('ANNOTATED SEQUENCE:')
#     #print (annotatedSequence)
#
#     annotatedSingleLine = annotatedSequence.replace('\n','')
#     cleanedAnnotatedSequence = annotatedSingleLine.upper()
#
#     print ('RAW      :' + roughFeatureSequence + '\n')
#     print ('ANNOTATED:' + annotatedSingleLine + '\n')
#     print ('CLEANED  :' + cleanedAnnotatedSequence + '\n')
#
#     assert_equal(roughFeatureSequence,cleanedAnnotatedSequence)
#     assert_true(len(annotatedSequence) > 3)
#
#
#
#
# def testAnnotateNOVELSequence():
#     print ('Test: Annotate a NOVEL SEQUENCE USING GFE')
#     #assignConfigurationValue('nmdp_act_rest_address', 'http://act.b12x.org/type_align' )
#     #assignConfigurationValue('nmdp_act_rest_address', 'http://localhost/type_align' )
#     assert_true(True)
#     #> HLA-A*02:01:01:12 with single SNP.
#     roughFeatureSequence = 'CCAGTTCTCACTCCCATTGGGTGTCGGGTTTCCAGAGAAGCCAATCAGTGTCGTCGCGGTCGCGGTTCTAAAGTCCGCACGCACCCACCGGGACTCAGATTCTCCCCAGACGCCGAGGATGGCCGTCATGGCGCCCCGAACCCTCGTCCTGCTACTCTCGGGGGCTCTGGCCCTGACCCAGACCTGGGCGGGTGAGTGCGGGGTCGGGAGGGAAACGGCCTCTGTGGGGAGAAGCAACGGGCCCGCCTGGCGGGGGCGCAGGACCCGGGAAGCCGCGCCGGGAGGAGGGTCGGGCGGGTCTCAGCCACTCCTCGTCCCCAGGCTCTCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGTCCGGAGTATTGGGACGGGGAGACACGGAAAGTGAAGGCCCACTCACAGACTCACCGAGTGGACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGCCGGTGAGTGACCCCGGCCCGGGGCGCAGGTCACGACCTCTCATCCCCCACGGACGGGCCAGGTCGCCCACAGTCTCCGGGTCCGAGATCCGCCCCGAAGCCGCGGGACCCCGAGACCCTTGCCCCGGGAGAGGCCCAGGCGCCTTTACCCGGTTTCATTTTCAGTTTAGGCCAAAAATCCCCCCAGGTTGGTCGGGGCGGGGCGGGGCTCGGGGGACCGGGCTGACCGCGGGGTCCGGGCCAGGTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGGGTACCAGGGGCCACGGGGCGCCTCCCTGATCGCCTGTAGATCTCCCGGGCTGGCCTCCCACAAGGAGGGGAGACAATTGGGACCAACACTAGAATATCGCCCTCCCTCTGGTCCTGAGGGAGAGGAATCCTCCTGGGTTTCCAGATCCTGTACCAGAGAGTGACTCTGAGGTTCCGCCCTGCTCTCTGACACAATTAAGGGATAAAATCTCTGAAGGAATGACGGGAAGACGATCCCTCGAATACTGATGAGTGGTTCCCTTTGACACACACAGGCAGCAGCCTTGGGCCCGTGACTTTTCCTCTCAGGCCTTGTTCTCTGCTTCACACTCAATGTGTGTGGGGGTCTGAGTCCAGCACTTCTGAGTCCTTCAGCCTCCACTCAGGTCAGGACCAGAAGTCGCTGTTCCCTCTTCAGGGACTAGAATTTTCCACGGAATAGGAGATTATCCCAGGTGCCTGTGTCCAGGCTGGTGTCTGGGTTCTGTGCTCCCTTCCCCATCCCAGGTGTCCTGTCCATTCTCAAGATAGCCACATGTGTGCTGGAGGAGTGTCCCATGACAGACGCAAAATGCCTGAATGATCTGACTCTTCCTGACAGACGCCCCCAAAACGCATATGACTCACCACGCTGTCTCTGACCATGAAGCCACCCTGAGGTGCTGGGCCCTGAGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGACAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTTTGCCCAAGCCCCTCACCCTGAGATGGGGTAAGGAGGGAGACGGGGGTGTCATGTCTTTTAGGGAAAGCAGGAGCCTCTCTGACCTTTAGCAGGGTCAGGGCCCCTCACCTTCCCCTCTTTTCCCAGAGCCGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCTTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCTGTGATGTGGAGGAGGAAGAGCTCAGGTGGGGAAGGGGTGAAGGGTGGGTCTGAGATTTCTTGTCTCACTGAGGGTTCCAAGACCCAGGTAGAAGTGTGCCCTGCCTCGTTACTGGGAAGCACCACCCACAATTATGGGCCTACCCAGCCTGGGCCCTGTGTGCCAGCACTTACTCTTTTGTAAAGCACCTGTTAAAATGAAGGACAGATTTATCACCTTGATTACAGCGGTGATGGGACCTGATCCCAGCAGTCACAAGTCACAGGGGAAGGTCCCTGAGGACCTTCAGGAGGGCGGTTGGTCCAGGACCCACACCTGCTTTCTTCATGTTTCCTGATCCCGCCCTGGGTCTGCAGTCACACATTTCTGGAAACTTCTCTGAGGTCCAAGACTTGGAGGTTCCTCTAGGACCTTAAGGCCCTGACTCCTTTCTGGTATCTCACAGGACATTTTCTTCCCACAGATAGAAAAGGAGGGAGCTACTCTCAGGCTGCAAGTAAGTATGAAGGAGGCTGATGCCTGAGGTCCTTGGGATATTGTGTTTGGGAGCCCATGGGGGAGCTCACCCACCCCACAATTCCTCCTCTAGCCACATCTTCTGTGGGATCTGACCAGGTTCTGTTTTTGTTCTACCCCAGGCAGTGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGTAAAGGTGAGAGCCTGGAGGGCCTGATGTGTGTTGGGTGTTGGGCGGAACAGTGGACACAGCTGTGCTGTGGGGTTTCTTTCCATTGGATGTATTGAGCATGCGATGGGCTGTTTAAAGTGTGACCCCTCACTGTGACAGATACGAATTTGTTCATGAATATTTTTTTCTATAGTGTGAGACAGCTGCCTTGTGTGGGACTGAGAGGCAAGAGTTGTTCCTGCCCTTCCCTTTGTGACTTGAAGAACCCTGACTTTGTTTCTGCAAAGGCACCTGCATGTGTCTGTGTTCGTGTAGGCATAATGTGAGGAGGTGGGGAGACCACCCCACCCCCATGTCCACCATGACCCTCTTCCCACGCTGACCTGTGCTCCCTCCCCAATCATCTTTCCTGTTCCAGAGAGGTGGGGCTGAGGTGTCTCCATCTCTGTCTCAACTTCATGGTGCACTGAGCTGTAACTTCTTCCTTCCCTATTAAAA'
#     hlaLocus = 'HLA-A'
#     #roughFeatureSequence = 'aag\nCGTCGT\nccg\nGGCTGA\naat'
#     alleleCallWithGFE = fetchSequenceAlleleCallWithGFE(roughFeatureSequence, hlaLocus)
#     #assert_true(len(alleleCallWithGFE) > 3)
#
#     #print 'ALLELE CALL RESULTS:'
#     #print alleleCallWithGFE
#
#     annotatedSequence = parseExons(roughFeatureSequence, alleleCallWithGFE)
#     #print ('ANNOTATED SEQUENCE:')
#     #print (annotatedSequence)
#
#     annotatedSingleLine = annotatedSequence.replace('\n','')
#     cleanedAnnotatedSequence = annotatedSingleLine.upper()
#
#     print ('RAW      :' + roughFeatureSequence + '\n')
#     print ('ANNOTATED:' + annotatedSingleLine + '\n')
#     print ('CLEANED  :' + cleanedAnnotatedSequence + '\n')
#     assert_true(len(annotatedSequence) > 3)
    


    