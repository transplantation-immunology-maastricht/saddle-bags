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

"""

from os import makedirs
from os.path import join, isdir, split, expanduser



from pycurl import Curl, FORM_FILE, HTTPHEADER, SSL_VERIFYHOST, SSL_VERIFYPEER

from io import BytesIO

from hashlib import md5

from xml.etree.ElementTree import fromstring
from ftplib import FTP

from tkinter import messagebox

from saddlebags.AlleleSubCommon import getConfigurationValue, createOutputFile, assignConfigurationValue, getSaddlebagsDirectory
from saddlebags.EnaSubXml import createProjectXML, createProjectSubmissionXML, createAnalysisSubmissionXML , createAnalysisXML
"""

import logging

from os.path import join, isdir
from os import makedirs
from subprocess import call

from sys import exc_info
from datetime import datetime
from shutil import copyfileobj
from gzip import open as gzipOpen

from saddlebags.AlleleSubCommon import getSaddlebagsDirectory, showYesNoBox, showInfoBox, getInfoBox, createOutputFile, getConfigurationValue, assignConfigurationValue
from saddlebags.EnaSubXml import createProjectXML, createProjectSubmissionXML
from saddlebags.EnaSubRest import performProjectSubmission, interpretAnalysisSubmissionResults
from saddlebags.EnaSubJar import findJarFile

# In this file we submit to EMBL/ENA using the webin .jar file.
# ENA Submission manual can be found here:
# https://ena-docs.readthedocs.io/en/latest/general-guide/webin-cli.html


def checkENAPrerequisites():
    # TODO: Check for prerequeisites for submission. The .jar file should exist and .java should be installed. Anything else?
    return True


def performBatchEnaSubmission(submissionBatch):
    logging.info('Submitting this batch of alleles:' + str(submissionBatch))

    for submission in submissionBatch.submissionBatch:
        performFullEnaSubmission(submission, submissionBatch)

def performFullEnaSubmission(submission, submissionBatch):
    logging.info('Performing an EMBL/ENA Submission.')

    submissionText = submission.enaSubmissionText

    # This includes a "seconds" measure, should be pretty unique.
    # TODO: I could add milliseconds to make this more unique.
    # TODO: I'm worried about batch submissions, is it a problem if the submission files have the same name? I think it's not a problem.
    dateTimeNow = '{:%Y_%m_%d_%H_%M_%S_%f}'.format(datetime.now())

    # TODO: If the submissionText is None, generate a submission now.
    if(submissionText is None or len(submissionText) < 5):
        logging.error('Oh no, there is no submission text. I probably forgot to generate a submission. Fix this bug please!')
        #generate a submission text now......


    if (checkENAPrerequisites):

        useTestServers = (int(getConfigurationValue('test_submission')) == 1)
        sequenceName = submission.localAlleleName
        # Are you sure? Test or Live?
        if useTestServers:
            logging.info('Using Test ENA Server.' + '\n')
            result = showYesNoBox("Submit to TEST / DEMO environment",
                "You are about to submit " + str(sequenceName) + " to the\n\nTEST / DEMO ENA environment.\n\nAre You Sure?")
        else:
            logging.info('Using Production ENA Server.' + '\n')
            result = showYesNoBox("Submit to LIVE / PROD environment",
                "You are about to submit " + str(sequenceName) + " to the\n\nLIVE / PROD ENA environment.\n\nAre You Sure?")
        if result:
            pass
        else:
            logging.error('Submission aborted by the user, because we do not want to submit to the Test/Live server')
            return

        # set some parameters real quick.
        jarFileLocation = findJarFile()
        logging.info('I will use this webin cli jar file:' + str(jarFileLocation))
        saddlebagsDirectory = getSaddlebagsDirectory()
        workingDirectory = join(saddlebagsDirectory, 'submission_temp')
        logging.info('I\'m working in this directory:' + str(workingDirectory))

        enaUserName = submissionBatch.enaUserName
        enaPassword = submissionBatch.enaPassword

        # Check the credentials, do they look okay?
        if (enaUserName is None or len(enaUserName) < 1):
            logging.warning('Missing ENA Username.' + '\n')
            enaUserName = getInfoBox("ENA Username Please", "You must provide ENA Username for submission.")
            submissionBatch.enaUserName = enaUserName
        else:
            logging.info('ENA Username ok.' + '\n')

        if (enaPassword is None or len(enaPassword) < 1):
            logging.warning('Missing ENA Password.' + '\n')
            enaPassword = getInfoBox("ENA Password Please", "You must provide ENA Password for user " + enaUserName + " for submission.\nSaddlebags will not store your password anywhere.")
            submissionBatch.enaPassword = enaPassword
        else:
            logging.info('ENA Password look ok.' + '\n')

        # Submission is divided into 3 stages (https://ena-docs.readthedocs.io/en/latest/general-guide/webin-cli.html)
        # Stage 1 - Register Study (Registering a Sample is not necessary for submitting a sequence)
        registerStudy(submissionBatch, workingDirectory, dateTimeNow)

        # Stage 2 - Prepare Files
        prepareSubmissionFiles(submission, submissionBatch, workingDirectory, dateTimeNow)

        # Stage 3 - Validate and Submit Files
        validateAndSubmit(submission, submissionBatch, workingDirectory, dateTimeNow)

        # Report Results
        # Send a summary of what was submitted.
        # Delete the files.
        # reportAndCleanup()
        # popup message with results.
        # Maybe: delete working directory, but the files first. Maybe delete output files?

    else:
        # TODO: Handle this better, tell the user somehow.
        logging.error('ENA Submission Requirements are not met. Details: (TODO)')

def registerStudy(submissionBatch, workingDirectory, dateTimeNow):
    logging.info('Registering a Project/Study')
    # effectively, study = project

    # The configuration value comes from a radio button in the configuration GUI. existing study = 1, new study = 2
    # Study = Project, because ENA is always sensible.
    newProject = (str(submissionBatch.chooseProject) == '2')


    if newProject:

        # Abstract = Description. ENA should have only one name for identical things, that's confusing.
        #popup confirmation for new project.
        if(showYesNoBox('Create New Project', 'Are you sure you want me to create a new Project to store your submission(s)?'
            + '\nID:\n' + str(submissionBatch.projectId)
            + '\nShort Title:\n' + str(submissionBatch.projectShortTitle)
            + '\nDescription:\n' + str(submissionBatch.projectAbstract))):


            # Generate Project and Project Submission XML Files
            try:
                projectFileName = join(workingDirectory, 'project.xml')
                #createProjectXML(fullXmlFilePath, projectID, projectShortTitle, projectAbstract):
                projectText = createProjectXML(projectFileName, submissionBatch.projectId, submissionBatch.projectShortTitle, submissionBatch.projectAbstract)

                logging.debug('Project Text:\n' + projectText)

                projectSubmissionFileName = join(workingDirectory, 'project_submission.xml')
                projectSubmissionText = createProjectSubmissionXML(projectSubmissionFileName
                    , 'proj_sub_' + dateTimeNow
                    , 'project.xml')

                logging.debug('Project Submission Text:\n' + projectSubmissionText)

            except Exception:
                logging.error('Cannot Create Project Submission XML')
                logging.error(exc_info())
                showInfoBox('Cannot Create Project Submission XML',
                    'Sorry, I failed to create a project XML file\nand I cannot continue.\n'
                    + str(exc_info()[1]))
                logging.error('Failure to create project submission file:' + str(exc_info()[1]) + '\n')
                return

            logging.info('Project/Submission XML files were created.\n')

            # Use REST to submit this project
            try:
                # Return value should be a tuple:
                # (Success, ProjectAccession, Messages[])
                (projectSubmissionSuccess, projectAccessionNumber, projectErrorMessages) = performProjectSubmission(projectSubmissionFileName, projectFileName, submissionBatch)

                if (projectSubmissionSuccess):
                    # Great. The project was created successfully.
                    # Lets use this new study accession moving forward. Any future submissions in this batch are submitted to the same project.
                    submissionBatch.studyAccession = projectAccessionNumber
                    submissionBatch.chooseProject = "1"

                    logging.info('New study has been uploaded, accession:' + str(submissionBatch.studyAccession) + '\n' + 'Subsequent submissions will use this project number.')
                    showInfoBox('Success!','New study has been uploaded, accession:' + str(
                        submissionBatch.studyAccession) + '\n' + 'Subsequent submissions will use this project number.')
                else:
                    messageText = ('There was a problem in the Project Submission.\n'
                        + 'I cannot continue.\n'
                        + 'These messages were reported by ENA:\n')
                    for errorMessage in projectErrorMessages:
                        messageText += ('\n' + errorMessage + '\n')
                    showInfoBox('Cannot Submit Project XML via REST', messageText)
                    logging.error('Failure to submit project submission file:' + str(exc_info()[1]) + '\n' + messageText + '\n')
                    return

            except Exception:
                logging.error('Cannot Submit Project XML')
                logging.error(exc_info()[1])
                showInfoBox('Cannot Submit Project XML',
                    'Sorry, I failed to submit the project XML file\nand I cannot continue.\n'
                    + str(exc_info()[1]))
                logging.error('Failure to upload project submission file:' + str(exc_info()[1]) + '\n')
                return


        else:
            logging.error('User asked me to not create a new project. I try to get a new project id from them.')
            submissionBatch.studyAccession = getInfoBox('A Project ID Please.','In that case, please provide a project ID I can use. You can find it on the ENA webin website after you login.\nIt should look like this:PRJEB12345')
            if(submissionBatch.studyAccession is not None and len(submissionBatch.studyAccession) > 3):
                submissionBatch.chooseProject = '1'
                logging.info('User gave me a new project number. I will use this one:' + str(submissionBatch.studyAccession))


    # existing project, we will use the supplied accession#. Easy.
    else:
        logging.info('Using existing study accession:' + str(submissionBatch.studyAccession) + '\n')
        pass

    assignConfigurationValue('submission_batch',submissionBatch)


def prepareSubmissionFiles(submission, submissionBatch, workingDirectory, dateTimeNow):
    logging.info('Preparing Submission Files')


    submissionShortFileName = 'HLA_Submission_' + dateTimeNow + '.txt'
    submissionFileName = join(workingDirectory, submissionShortFileName)

    zippedShortFileName = submissionShortFileName + '.gz'
    zippedFileName = join(workingDirectory, zippedShortFileName)

    manifestShortFileName = 'manifest_' + dateTimeNow + '.txt'
    manifestFileName = join(workingDirectory, manifestShortFileName)

    if not isdir(workingDirectory):
        makedirs(workingDirectory)

    # Create the submission file
    try:
        outputFileObject = open(submissionFileName, 'w')
        outputFileObject.write(submission.enaSubmissionText)
        logging.debug('Submission Text:\n' + submission.enaSubmissionText)
        outputFileObject.close()

    except Exception:
        logging.error('Cannot Write Submission Flatfile')
        logging.error(exc_info())
        showInfoBox('Cannot Write Submission Flatfile',
            'Sorry, I failed to create the submission file:\n'
            + str(submission.enaSubmissionText)
            + '\n and I cannot continue.\nMaybe this is a '
            + 'permissions issue, are these folders read only?\n'
            + str(exc_info()[1]))
        logging.error('Failure to create submission file:' + str(exc_info()[1]) + '\n')
        return

    logging.info('Submission file was created:\n' + str(submissionFileName) + '\n')

    # gzip the submission file.  Make a gz file.
    try:
        with open(submissionFileName, 'rb') as fileIn, gzipOpen(zippedFileName, 'wb') as fileOut:
            copyfileobj(fileIn, fileOut)

    except Exception:
        logging.error('Cannot Compress Submission File')
        logging.error(exc_info())
        showInfoBox('Cannot Compress Submission File',
            'Sorry, I failed to compress the submission file:\n'
            + str(zippedFileName)
            + '\n and I cannot continue.\n'
            + str(exc_info()[1]))
        logging.error('Failure to create zip file:' + str(exc_info()[1]) + '\n')
        return

    logging.info('Zip file was created:\n' + str(zippedFileName) + '\n')

    # Create the Manifest file, which looks like this:
    # NAME    Novel_HLA_Allele_A
    # STUDY   PRJEB22887
    # FLATFILE    ENA.HLA.Submission_A.txt.gz
    manifestFile = createOutputFile(manifestFileName)
    # I'm using the "allele_name" as the submission name in the manifest file. Is that okay?
    # Maybe it should be the analysis alias...which I am probably deleting.
    # TODO: For batch submissions, i shouldn't be able to get the configuration value.. Refactor for batch submission.
    manifestFile.write('NAME\t' + str(submission.localAlleleName) + '\n')
    manifestFile.write('STUDY\t' + str(submissionBatch.studyAccession) + '\n')
    manifestFile.write('FLATFILE\t' + str(zippedFileName) + '\n')
    manifestFile.close()

def validateAndSubmit(submission, submissionBatch, workingDirectory, dateTimeNow):
    logging.info('Validating and Submitting Files.')

    # Find the webin cli application
    webinJarLocation = findJarFile()
    logging.info('Webin Jar file found here:' + str(webinJarLocation))

    # Construct the cli Parameters (for both validation and submission)
    # <outputDir> can be specified using the -outputDir option, the
    # <context> is specified using the -context option, and the
    # <name> is a submitter provided unique name specified in the manifest file.
    # manifest file is specified using the -manifest <filename>
    # -userName=USER, -password=PASSWORD]
    outputDir = join(workingDirectory, 'SubmissionOutput')
    if not isdir(outputDir):
        makedirs(outputDir)
    manifestShortFileName = 'manifest_' + dateTimeNow + '.txt'
    manifestFileName = join(workingDirectory, manifestShortFileName)

    # TODO: they list an option to use a proxy. Maybe I need to use a proxy at some point, look at the ENA webin instructions. https://ena-docs.readthedocs.io/en/latest/general-guide/webin-cli.html
    # TODO: I'm not super familiar with subprocess.call. I create an array of parameters. Will it handle when my paths have a space? I hope it's smart enough to do that.
    validateCommand = [
        'java', '-jar', str(webinJarLocation)
        , '-validate'
        , '-outputDir', outputDir
        , '-context', 'sequence'
        , '-manifest', manifestFileName
        , '-userName', submissionBatch.enaUserName
        , '-password', submissionBatch.enaPassword

    ]
    submitCommand = [
        'java', '-jar', str(webinJarLocation)
        , '-submit'
        , '-outputDir', outputDir
        , '-context', 'sequence'
        , '-manifest', manifestFileName
        , '-userName', submissionBatch.enaUserName
        , '-password', submissionBatch.enaPassword
    ]

    if(int(getConfigurationValue('test_submission')) == 1):
        validateCommand.append('-test')
        submitCommand.append('-test')

    # TODO: This puts a password in the log file. Is that okay?
    logging.debug('validate webin command:' + str(validateCommand))
    logging.debug('submit webin command:' + str(submitCommand))

    # Validate Sequence
    #jarSubmissionResults = call(validateCommand)
    jarSubmissionResults = call(submitCommand)

    #if(str(jarSubmissionResults) != '1'):
    #    logging.error('Error executing the .jar file to submit sequence to ENA.')
    #    showInfoBox('Error executing jar file', 'Error executing the .jar file to submit sequence to ENA. Not sure what the problem is.')

    # Parse the Result Files. Might need an analysis accession but not sure.
    # Manifest Result file lists some errors, such as incorrect password. Tell the user to check this file if there are problems.
    # manifestResultFile = '/home/ben/saddlebags/submission_temp/SubmissionOutput/manifest_2019_06_27_17_08_50_626831.txt.report'
    # CLI report file has just a log of the commandline tool. It also lists the analysis accession #
    #cliReportFile = '/home/ben/saddlebags/submission_temp/SubmissionOutput/webin-cli.report'
    # the analysis reciept / result file has the analysis accession, submission accession, and messages.
    #analysisResultFile= '/home/ben/saddlebags/submission_temp/SubmissionOutput/sequence/HLA-DRA_MUMC_1/submit/receipt.xml'

    analysisResultFileLocation = join(join(join(join(outputDir, 'sequence'), submission.localAlleleName), 'submit'), 'receipt.xml')
    analysisResultFile =  open(analysisResultFileLocation, 'r')
    analysisResultText = analysisResultFile.read()

    (analysisSubmissionSuccess, analysisAccessionNumber, analysisErrorMessages) = interpretAnalysisSubmissionResults(analysisResultText)
    if (analysisSubmissionSuccess):
        # Great. The analysis was created successfully.
        showInfoBox('Successful Submission.','Successful submission. ' + str(submission.localAlleleName) + ' has analysis Accession number is:' + str(analysisAccessionNumber))
        pass
    else:
        messageText = ('There was a problem in the Analysis Submission.\n'
            + 'I cannot continue.\n'
            + 'These messages were reported by ENA:\n')
        for errorMessage in analysisErrorMessages:
            messageText += ('\n' + errorMessage + '\n')
        showInfoBox('Cannot Submit Analysis XML via REST', messageText)
        logging.error('Failure to submit analysis submission file:' + str(exc_info()[1]) + '\n')
        return







