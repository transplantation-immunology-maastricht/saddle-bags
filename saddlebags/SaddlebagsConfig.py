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
# You should have received a copy of the GNU Lesser General Publsic License
# along with saddle-bags. If not, see <http://www.gnu.org/licenses/>.

# TODO: tidy up these imports.
# import sys
from sys import exc_info
#
# from datetime import datetime
#
# import csv
#
# try:
#     from sys import _MEIPASS
# except Exception:
#     print('No MEIPASS Directory. This is not running from a compiled EXE file. No problem.')
#
# from os import makedirs, remove, rmdir, name
from os.path import join, isfile
    #, expanduser, abspath, isdir, split

# from tkinter import messagebox, simpledialog
#
from xml.etree import ElementTree as ET
from xml.dom import minidom as MD
#
# from pycurl import Curl
#
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
#
# from io import StringIO, BytesIO
# from urllib.parse import urlencode
#
# from json import loads, dumps # these method names are terrible IMO.
#
# from zipfile import ZipFile
from saddlebags.AlleleSubCommon import getSaddlebagsDirectory, createOutputFile, showInfoBox
from saddlebags.Logging import initializeLog
from saddlebags.AlleleSubmission import SubmissionBatch, AlleleSubmission

import logging

def clearGlobalVariables():
    # Clear my configuration, fearlessly and without hesitation.
    global globalVariables
    globalVariables = {}

def initializeGlobalVariables():
    # I'm storing global variables in a dictionary for now.
    global globalVariables

    if not ("globalVariables" in globals()):
        globalVariables={}

        # I'm removing this here, and putting it in the loadConfigurationFile method.
        # I think this is safe to delete.
        #submissionBatch = SubmissionBatch(True)
        #assignConfigurationValue('submission_batch', submissionBatch)

def assignConfigurationValue(configurationKey, configurationValue):
    # assignConfigurationValue will overwrite config value without question.
    initializeGlobalVariables()

    # Lists or Strings are handled well by the configuration serializer.
    if(type(configurationValue) is list or type(configurationValue) is str ):
        globalVariables[configurationKey] = serializeConfigValue(configurationValue)
    else:
        globalVariables[configurationKey] = configurationValue

    globalVariables[configurationKey] = configurationValue
    logging.debug ('Just stored configuration key ' + configurationKey + ' which is ' + str(configurationValue) + ' of type ' + str(type(configurationValue)))

def assignIfNotExists(configurationKey, configurationValue):
    # Use this assigner if we want to declare important, new configuration values.
    # Using this method, we will not overwrite custom values
    # But we will provide critical new config values.
    initializeGlobalVariables()
    if configurationKey not in globalVariables.keys():
        assignConfigurationValue(configurationKey, configurationValue)

def getConfigurationValue(configurationKey):
    if configurationKey in globalVariables.keys():

        configurationValue = globalVariables[configurationKey]

        logging.debug ('Retrieving configuration key ' + configurationKey + ' which is ' + str(configurationValue) + ' of type ' + str(type(configurationValue)))

        if (type(configurationValue) is str):
            return deserializeConfigValue(configurationValue)
        else:
            return configurationValue
    else:
        logging.warning ('Configuration Key Not Found:' + configurationKey)
        #raise KeyError('Key Not Found:' + configurationKey)
        return None

def loadFromCSV(csvFileName):
    # Read submission data from a .csv file.
    logging.debug ('loading data from this csv file:' + csvFileName)

    #TODO: If it's a zip file, I should be able to look for a .csv file in the root of the .zip.

    # Load our submission batch. Does it already exist? It should.
    submissionBatch = getConfigurationValue('submission_batch')
    if(submissionBatch == None):
        logging.warning('Loading from CSV file. There was no batch of submissions, so I had to create an empty batch')

        submissionBatch = SubmissionBatch()

        # Assign some default information about this batch of submissions
        submissionBatch.ipdSubmitterId = ''
        submissionBatch.ipdSubmitterName = ''
        submissionBatch.ipdAltContact = ''
        submissionBatch.ipdSubmitterEmail = ''
        submissionBatch.labOfOrigin = ''
        submissionBatch.labContact = ''

    # Open the CSV and read the header.
    csvFile = open(csvFileName, 'r')
    csvInputReader = csv.reader(csvFile)

    header = next(csvInputReader)
    #logging.debug("This is the header row:" + str(header))

    # Convert the header names to uppercase to allow the input of upper or lowercase header names.
    header = [item.upper() for item in header]

    # Check for the existence of each required field in the header.
    # Store the indices of the column headers in a dictionary. So I can lookup the fields in any order, input file is more flexible.
    requiredFields = ['CLASS','CELLBANK','CELLID','CITATIONS','CLOSESTALLELEWRITTENDESCRIPTION','CONSANGUINEOUS','ETHNICORIGIN','GENELOCUS','HOMOZYGOUS','IPDSUBMISSIONIDENTIFIER','IPDSUBMISSIONVERSION','ENASEQUENCEACCESSION','LOCALALLELENAME','MATERIALAVAILABILITY','METHODCOMMENTS','NUMOFREACTIONS','PRIMARYSEQUENCINGMETHODOLOGY','PRIMERS','PRIMERTYPE','SECONDARYSEQUENCINGMETHODOLOGY','SEQUENCEDINISOLATION','SEQUENCINGDIRECTION','SEX','TYPEDALLELES','SEQUENCE']
    requiredFieldIndices = {}
    try:
        for requiredField in requiredFields:
            #print('assigning index of ' + str(requiredField) + '=' + str(header.index(requiredField)))
            requiredFieldIndices[requiredField] = header.index(requiredField)

    except:
        logging.error('Error when reading from the CSV file (' + str(csvFileName) + '). Perhaps you are missing a required field:' + str(requiredField))

    for submissionCSVRow in csvInputReader:
        #logging.debug("loading this row from csv:" + str(submissionCSVRow))


        #print ('Im creating a submission, the required field indices look like this:' + str(requiredFieldIndices))

        # Create a submission Object
        # Get each column of data and store it in the submission.
        submission = AlleleSubmission()
        # TODO: What if the .csv file is not annotated? I believe identifyFeaturesFromFormattedSequence expects the annotated sequence.
        submission.submittedAllele.rawSequence = submissionCSVRow[requiredFieldIndices['SEQUENCE']]
        submission.submittedAllele.identifyFeaturesFromFormattedSequence()
        submission.submittedAllele.geneLocus = submissionCSVRow[requiredFieldIndices['GENELOCUS']]
        submission.submittedAllele.hlaClass = submissionCSVRow[requiredFieldIndices['CLASS']]
        submission.localAlleleName = submissionCSVRow[requiredFieldIndices['LOCALALLELENAME']]
        submission.closestAlleleWrittenDescription = submissionCSVRow[requiredFieldIndices['CLOSESTALLELEWRITTENDESCRIPTION']]
        submission.ipdSubmissionIdentifier = submissionCSVRow[requiredFieldIndices['IPDSUBMISSIONIDENTIFIER']]
        submission.ipdSubmissionVersion = submissionCSVRow[requiredFieldIndices['IPDSUBMISSIONVERSION']]
        submission.enaAccessionIdentifier = submissionCSVRow[requiredFieldIndices['ENASEQUENCEACCESSION']]
        submission.cellId = submissionCSVRow[requiredFieldIndices['CELLID']]
        submission.ethnicOrigin = submissionCSVRow[requiredFieldIndices['ETHNICORIGIN']]
        submission.sex = submissionCSVRow[requiredFieldIndices['SEX']]
        submission.consanguineous = submissionCSVRow[requiredFieldIndices['CONSANGUINEOUS']]
        submission.homozygous = submissionCSVRow[requiredFieldIndices['HOMOZYGOUS']]
        submission.typedAlleles = parseTypedAlleleInput(submissionCSVRow[requiredFieldIndices['TYPEDALLELES']])
        submission.materialAvailability = submissionCSVRow[requiredFieldIndices['MATERIALAVAILABILITY']]
        submission.cellBank = submissionCSVRow[requiredFieldIndices['CELLBANK']]
        submission.primarySequencingMethodology = submissionCSVRow[requiredFieldIndices['PRIMARYSEQUENCINGMETHODOLOGY']]
        submission.secondarySequencingMethodology = submissionCSVRow[requiredFieldIndices['SECONDARYSEQUENCINGMETHODOLOGY']]
        submission.primerType = submissionCSVRow[requiredFieldIndices['PRIMERTYPE']]
        submission.primers = submissionCSVRow[requiredFieldIndices['PRIMERS']]
        submission.sequencedInIsolation = submissionCSVRow[requiredFieldIndices['SEQUENCEDINISOLATION']]
        submission.sequencingDirection = submissionCSVRow[requiredFieldIndices['SEQUENCINGDIRECTION']]
        submission.numOfReactions = submissionCSVRow[requiredFieldIndices['NUMOFREACTIONS']]
        submission.methodComments = submissionCSVRow[requiredFieldIndices['METHODCOMMENTS']]
        submission.citations = submissionCSVRow[requiredFieldIndices['CITATIONS']]
        submissionBatch.submissionBatch.append(submission)

    csvFile.close()

def assignConfigName():
    # Join together the working directory, a subfolder called "saddlebags", and the config name.
    assignConfigurationValue('config_file_location', join(getSaddlebagsDirectory(), 'Saddlebags.Config.xml'))

def serializeConfigValue(configListObject):
    # Encode semicolons within a config value.
    # Store lists as a string separated by a semicolon.
    if(configListObject is None):
        logging.warning('Serializing a config list object that is None!')
    else:
        #logging.debug('serializing a config value:' + str(configListObject))
        pass

    if (type(configListObject) == str):
        return configListObject.replace(';', '@@@')
    elif (type(configListObject) == list):
        serializedString = ''
        # TODO: what if the list elements are not strings? Just don't do that. Lol.
        # TODO: shouldn't there be built in methods for serializing? probably.
        # Encode semicolons in the list elements.
        for configString in configListObject:
            serializedString = serializedString + str(configString).replace(';', '@@@') + ';'
        serializedString = serializedString[:-1]
        return
    elif (configListObject is None):
        return ''
    else:
        # TODO: Yeah I did have one random problem when the type was an int. I can just always use strings or maybe handle it smarter.

        raise (Exception('Unknown configuration type, can not serialize:' + str(type(configListObject))))

def deserializeConfigValue(serializedConfigString):
    # Split strings containing semicolon into a list.
    # Decode @@@ back into a semicolon.
    if(serializedConfigString is None):
        logging.warning('Deserializing a config string that is None!')
    else:
        #logging.debug('deserializing a config value: ' + str(serializedConfigString))
        pass

    if (';' in serializedConfigString):
        configList = serializedConfigString.split(';')
        # Maybe a for loop is not necessary but I want to make sure i can change the string in the list.
        for i in range(0, len(configList)):
            configList[i] = configList[i].replace('@@@', ';')
        return configList
    else:
        # TODO What if it's not actually a string? Could be an arbitrary object...
        return serializedConfigString.replace('@@@', ';')

def writeConfigurationFile():
    assignConfigName()
    logging.debug('Writing a config file to:\n' + getConfigurationValue('config_file_location'))

    # Root node stores "normal" configuration keys, stuff related to software
    # and not necessarily HLA or sequence submission.
    root = ET.Element("config")

    # Loop through configuration keys in the dictionary and write em out.
    for key in globalVariables.keys():
        # "normal" configuration keys, stuff related to software and not necessarily HLA or sequence submission.
        # Some config values I don't want to store. I can add more to this list if i want.
        # Don't store passwords. Passwords should be attached to the submission batch, no need to worry about it here.
        # I handle the submission batch manually.
        if (key not in [
            # 'ena_password'
            # , 'ipd_password'
            # ,
            'submission_batch'
        ]):
            # getConfigurationValue will handle serializing and encoding semicolons.
            ET.SubElement(root, key).text = getConfigurationValue(key)

    # Add keys for "each" batch of submissions.
    # TODO: May want to add functionality for multiple batches later. Put this in a loop.
    # TODO: Batches of Batches, I don't think this is necessary. 1 batch is fine for now.
    submissionBatch = getConfigurationValue('submission_batch')

    # If the config is not already initiated, this can be None. Make a new one.
    if (submissionBatch is None):
        submissionBatch = SubmissionBatch()

    # Create a node object, most of this stuff can be parameters on the node.
    # Dont write any passwords to the config.
    submissionBatchElement = ET.SubElement(root, 'submission_batch')
    submissionBatchElement.set('enausername', serializeConfigValue(submissionBatch.enaUserName))
    submissionBatchElement.set('ipdsubmitterid', serializeConfigValue(submissionBatch.ipdSubmitterId))
    submissionBatchElement.set('ipdsubmittername', serializeConfigValue(submissionBatch.ipdSubmitterName))
    submissionBatchElement.set('ipdaltcontact', serializeConfigValue(submissionBatch.ipdAltContact))
    submissionBatchElement.set('ipdsubmitteremail', serializeConfigValue(submissionBatch.ipdSubmitterEmail))
    submissionBatchElement.set('laboforigin', serializeConfigValue(submissionBatch.labOfOrigin))
    submissionBatchElement.set('labcontact', serializeConfigValue(submissionBatch.labContact))
    submissionBatchElement.set('choosestudy', serializeConfigValue(submissionBatch.chooseStudy))
    submissionBatchElement.set('studyaccession', serializeConfigValue(submissionBatch.studyAccession))
    submissionBatchElement.set('studyid', serializeConfigValue(submissionBatch.studyId))
    submissionBatchElement.set('studyshorttitle', serializeConfigValue(submissionBatch.studyShortTitle))
    submissionBatchElement.set('studyabstract', serializeConfigValue(submissionBatch.studyAbstract))

    # Keys for each submission.
    for hlaSubmission in submissionBatch.submissionBatch:
        submissionElement = ET.SubElement(submissionBatchElement, 'submission')

        # Most of this stuff is attributes. Store the Sequence as the text of this element.
        submissionElement.text = hlaSubmission.submittedAllele.getAnnotatedSequence(includeLineBreaks=False) # Problem: This is returning Nothing.
        submissionElement.set('genelocus', serializeConfigValue(hlaSubmission.submittedAllele.geneLocus))
        submissionElement.set('class', serializeConfigValue(hlaSubmission.submittedAllele.hlaClass))
        submissionElement.set('localallelename', serializeConfigValue(hlaSubmission.localAlleleName))
        submissionElement.set('closestallelewrittendescription', serializeConfigValue(hlaSubmission.closestAlleleWrittenDescription))
        submissionElement.set('ipdsubmissionidentifier', serializeConfigValue(hlaSubmission.ipdSubmissionIdentifier))
        submissionElement.set('ipdsubmissionversion', serializeConfigValue(hlaSubmission.ipdSubmissionVersion))
        submissionElement.set('enaaccessionidentifier', serializeConfigValue(hlaSubmission.enaAccessionIdentifier))
        submissionElement.set('cellid', serializeConfigValue(hlaSubmission.cellId))
        submissionElement.set('ethnicorigin', serializeConfigValue(hlaSubmission.ethnicOrigin))
        submissionElement.set('sex', serializeConfigValue(hlaSubmission.sex))
        submissionElement.set('consanguineous', serializeConfigValue(hlaSubmission.consanguineous))
        submissionElement.set('homozygous', serializeConfigValue(hlaSubmission.homozygous))

        # TypedAlleles is special. It's a dictionary, where the keys are an HLA locus.
        typedAlleleText = ""
        if hlaSubmission.typedAlleles is not None:
            for loci in sorted(hlaSubmission.typedAlleles.keys()):
                typedAlleleText += str(loci) + '*' + hlaSubmission.typedAlleles[loci] + ';'
        submissionElement.set('typedalleles', serializeConfigValue(typedAlleleText))
        submissionElement.set('materialavailability', serializeConfigValue(hlaSubmission.materialAvailability))
        submissionElement.set('cellbank', serializeConfigValue(hlaSubmission.cellBank))
        submissionElement.set('primarysequencingmethodology',
                              serializeConfigValue(hlaSubmission.primarySequencingMethodology))
        submissionElement.set('secondarysequencingmethodology',
                              serializeConfigValue(hlaSubmission.secondarySequencingMethodology))
        submissionElement.set('primertype', serializeConfigValue(hlaSubmission.primerType))
        submissionElement.set('primers', serializeConfigValue(hlaSubmission.primers))
        submissionElement.set('sequencedinisolation', serializeConfigValue(hlaSubmission.sequencedInIsolation))
        submissionElement.set('sequencingdirection', serializeConfigValue(hlaSubmission.sequencingDirection))
        submissionElement.set('numofreactions', serializeConfigValue(hlaSubmission.numOfReactions))
        submissionElement.set('methodcomments', serializeConfigValue(hlaSubmission.methodComments))
        submissionElement.set('citations', serializeConfigValue(hlaSubmission.citations))

    xmlText = ET.tostring(root, encoding='utf8', method='xml')
    prettyXmlText = MD.parseString(xmlText).toprettyxml()

    xmlOutput = createOutputFile(globalVariables['config_file_location'])
    xmlOutput.write(prettyXmlText)
    xmlOutput.close()

def loadConfigurationFile():
    # TODO: should I clear my configuration first? I have a method to purge my globals.
    # I don't know right now, but probably not.
    assignConfigName()

    try:

        if not isfile(globalVariables['config_file_location']):
            logging.info ('The config file does not exist yet. I will not load it:\n' + globalVariables['config_file_location'])
        else:
            logging.info ('The config file already exists, I will load it:\n' + globalVariables['config_file_location'])

            tree = ET.parse(globalVariables['config_file_location'])
            root = tree.getroot()

            for child in root:
                #logging.debug('The child tag is:' + child.tag)

                # If the child node is a submission batch
                if(child.tag == 'submission_batch'):

                    if(child):
                        # If the submission batch has children nodes, start with an empty batch.
                        submissionBatch = SubmissionBatch(False)
                    else:
                        # Otherwise we want to start with a single empty submission in the batch.
                        submissionBatch = SubmissionBatch(True)

                    # Assign some information about this batch of submissions.
                    submissionBatch.enaUserName = deserializeConfigValue(child.attrib['enausername'])
                    submissionBatch.studyAccession = deserializeConfigValue(child.attrib['studyaccession'])
                    submissionBatch.chooseStudy = deserializeConfigValue(child.attrib['choosestudy'])
                    submissionBatch.ipdSubmitterId = deserializeConfigValue(child.attrib['ipdsubmitterid'])
                    submissionBatch.ipdSubmitterName = deserializeConfigValue(child.attrib['ipdsubmittername'])
                    submissionBatch.ipdAltContact = deserializeConfigValue(child.attrib['ipdaltcontact'])
                    submissionBatch.ipdSubmitterEmail = deserializeConfigValue(child.attrib['ipdsubmitteremail'])
                    submissionBatch.labOfOrigin = deserializeConfigValue(child.attrib['laboforigin'])
                    submissionBatch.labContact = deserializeConfigValue(child.attrib['labcontact'])
                    submissionBatch.studyId = deserializeConfigValue(child.attrib['studyid'])
                    submissionBatch.studyShortTitle = deserializeConfigValue(child.attrib['studyshorttitle'])
                    submissionBatch.studyAbstract = deserializeConfigValue(child.attrib['studyabstract'])

                    # Loop the children, they are submission objects. Load up their information.
                    for submissionChild in child:
                        #logging.debug('The submission child tag is:' + submissionChild.tag)
                        #logging.debug('This submission has the text:' + submissionChild.text)
                        # Add a few submissions to this batch.
                        # Submission # 1
                        submission = AlleleSubmission()
                        submission.submittedAllele.rawSequence = submissionChild.text
                        submission.submittedAllele.identifyFeaturesFromFormattedSequence()
                        submission.submittedAllele.geneLocus = deserializeConfigValue(submissionChild.attrib['genelocus'])
                        submission.localAlleleName = deserializeConfigValue(submissionChild.attrib['localallelename'])
                        submission.submittedAllele.hlaClass = deserializeConfigValue(submissionChild.attrib['class'])
                        submission.closestAlleleWrittenDescription = deserializeConfigValue(submissionChild.attrib['closestallelewrittendescription'])
                        submission.ipdSubmissionIdentifier = deserializeConfigValue(submissionChild.attrib['ipdsubmissionidentifier'])
                        submission.ipdSubmissionVersion = deserializeConfigValue(submissionChild.attrib['ipdsubmissionversion'])
                        submission.enaAccessionIdentifier = deserializeConfigValue(submissionChild.attrib['enaaccessionidentifier'])
                        submission.cellId = deserializeConfigValue(submissionChild.attrib['cellid'])
                        submission.ethnicOrigin = deserializeConfigValue(submissionChild.attrib['ethnicorigin'])
                        submission.sex = deserializeConfigValue(submissionChild.attrib['sex'])
                        submission.consanguineous = deserializeConfigValue(submissionChild.attrib['consanguineous'])
                        submission.homozygous = deserializeConfigValue(submissionChild.attrib['homozygous'])
                        #print ('I am about to read and store my typed alleles.')
                        childElementText = submissionChild.attrib['typedalleles']
                        #print ('element text:' + childElementText)
                        deserializedText = deserializeConfigValue(childElementText)
                        #print ('deserialized text:' + deserializedText)
                        parsedObject = parseTypedAlleleInput(deserializedText)
                        #print('parsedObject:' + str(parsedObject))
                        submission.typedAlleles = parsedObject
                        #print ('Success.')
                        submission.materialAvailability = deserializeConfigValue(submissionChild.attrib['materialavailability'])
                        submission.cellBank = deserializeConfigValue(submissionChild.attrib['cellbank'])
                        submission.primarySequencingMethodology = deserializeConfigValue(submissionChild.attrib['primarysequencingmethodology'])
                        submission.secondarySequencingMethodology = deserializeConfigValue(submissionChild.attrib['secondarysequencingmethodology'])
                        submission.primerType = deserializeConfigValue(submissionChild.attrib['primertype'])
                        submission.primers = deserializeConfigValue(submissionChild.attrib['primers'])
                        submission.sequencedInIsolation = deserializeConfigValue(submissionChild.attrib['sequencedinisolation'])
                        submission.sequencingDirection = deserializeConfigValue(submissionChild.attrib['sequencingdirection'])
                        submission.numOfReactions = deserializeConfigValue(submissionChild.attrib['numofreactions'])
                        submission.methodComments = deserializeConfigValue(submissionChild.attrib['methodcomments'])
                        submission.citations = deserializeConfigValue(submissionChild.attrib['citations'])
                        submissionBatch.submissionBatch.append(submission)

                    # Store my submission batch in the global variables.
                    assignConfigurationValue('submission_batch', submissionBatch)

                    logging.debug('Just loaded the config file and stored the submission batch. Batch is length:' + str(len(submissionBatch.submissionBatch)))

                else:
                    # Any arbitrary configuration value, just store it.
                    assignConfigurationValue(child.tag, child.text)

        # Here is where I assign the common/critical configuration values
        # I do this if the config file already existed, or if it didnt.
        # test_submission indicates if we should use the "test" values.
        # I think I'll use this value for both ENA and IPD submissions, if it applies.
        assignIfNotExists('test_submission', '1')
        # Log levels are defined in the Saddlebags config, and passed into the python logging module.
        assignIfNotExists('logging', 'DEBUG')
        # Placeholder for proxy configuration. Not needed right now but maybe for the future.
        #assignIfNotExists('proxy', None)
        assignIfNotExists('ena_rest_address_test', 'https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/')
        assignIfNotExists('ena_rest_address_prod', 'https://www.ebi.ac.uk/ena/submit/drop-box/submit/')
        assignIfNotExists('nmdp_act_rest_address', 'http://act.b12x.org/annotate')
        assignIfNotExists('webin_jar_location','webin-cli.jar')
        assignIfNotExists('submission_batch', SubmissionBatch(True))

        writeConfigurationFile()

        # Last step is to initialize the log files. Why is this the last step? initializing log should be first
        # but I need some config values before starting the log.
        initializeLog()
    except:
        logging.error('Error when loading configuration file:' + str(globalVariables['config_file_location']) + '.\nTry deleting your configuration file and reload Saddlebags.\n' + str(exc_info()[1]))
        showInfoBox('Error Loading Configuration','Error when loading configuration file:' + str(globalVariables['config_file_location']) + '.\nTry deleting your configuration file and reload Saddlebags.\n' + str(exc_info()[1]))

        # TODO: Should I just delete the config file with any exception? Probably not.
        # TODO: I should ask permission, and clear the config and create a new one.
        # Messagebox:'New config?'
        # If yes?
            # New config
        # If no
            # Close program, they must fix config manually.

def parseTypedAlleleInput(alleleInputString):
    # Typed alleles are special. It's a dictionary, where the key is the Locus (HLA-A) and the value is a String, with a list of typings separated by a comma.
    # The input is string of standard nomenclature HLA alleles, separated by a semicolon. (HLA-A*02:01;HLA-A*03:02:14)
    # submission.typedAlleles = submissionCSVRow[requiredFieldIndices['TYPEDALLELES']]
    typedAlleleDictionary = {}

    try:
        #Semicolon separates the different typing
        typedAlleleTokens = alleleInputString.split(';')
        for typedAlleleToken in typedAlleleTokens:
            # * separates the locus and allele.
            if('*' in typedAlleleToken ):
                loci, allele = typedAlleleToken.split('*')
                if loci in typedAlleleDictionary.keys():
                    typedAlleleDictionary[loci] = str(typedAlleleDictionary[loci]) + ',' + allele
                else:
                    typedAlleleDictionary[loci] = allele

        return typedAlleleDictionary
    except:
        logging.error('Error when parsing HLA typing string: ' + str(alleleInputString))
        return {}