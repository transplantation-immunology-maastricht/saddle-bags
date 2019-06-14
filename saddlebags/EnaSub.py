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

from sys import exc_info

from pycurl import Curl, FORM_FILE, HTTPHEADER, SSL_VERIFYHOST, SSL_VERIFYPEER

from io import BytesIO

from hashlib import md5
from datetime import datetime
from xml.etree.ElementTree import fromstring
from shutil import copyfileobj
from gzip import open as gzipOpen
from ftplib import FTP

from tkinter import messagebox

from saddlebags.AlleleSubCommon import getConfigurationValue, createOutputFile, assignConfigurationValue, getSaddlebagsDirectory
from saddlebags.EnaSubXml import createProjectXML, createProjectSubmissionXML, createAnalysisSubmissionXML , createAnalysisXML
"""

import logging

from os.path import join

from saddlebags.AlleleSubCommon import getSaddlebagsDirectory, showQuestionBox, showInfoBox, createOutputFile, getConfigurationValue

from saddlebags.EnaSubXml import createProjectXML, createProjectSubmissionXML
from saddlebags.EnaSubRest import performProjectSubmission

# In this file we submit to EMBL/ENA using the webin .jar file.
# ENA Submission manual can be found here:
# https://ena-docs.readthedocs.io/en/latest/general-guide/webin-cli.html


def checkENAPrerequisites():
    # TODO: Check for prerequeisites for submission. The .jar file should exist and .java should be installed. Anything else?
    return True


def performBatchEnaSubmission(submissionBatch):
    logging.info('Submitting this batch of alleles:' + str(submissionBatch))

    for submission in submissionBatch():
        performFullEnaSubmission(submission)

def performFullEnaSubmission(submission):
    logging.info('Performing an EMBL/ENA Submission.')

    submissionText = submission.enaSubmissionText
    # TODO: If the submissionText is None, generate a submission now.



    if (checkENAPrerequisites):

        useTestServers = (int(getConfigurationValue('test_submission')) == 1)
        # Are you sure? Test or Live?
        if useTestServers:
            logging.info('Using Test ENA Server.' + '\n')
            result = showQuestionBox("Submit to TEST / DEMO environment",
                "You are about to submit a sequence to the\n\nTEST / DEMO ENA environment.\n\nAre You Sure?")
        else:
            logging.info('Using Production ENA Server.' + '\n')
            result = showQuestionBox("Submit to LIVE / PROD environment",
                "You are about to submit a sequence to the\n\nLIVE / PROD ENA environment.\n\nAre You Sure?")
        if result == 'yes':
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
        enaUsername = getConfigurationValue('ena_username')
        enaPassword = getConfigurationValue('ena_password')

        # Check the credentials, do they look okay?
        if (enaUsername is None or len(enaUsername) < 1 or enaPassword is None or len(enaPassword) < 1):
            showInfoBox('Missing Login Credentials',
                'You must provide ENA username and password.\n'
                'Please use the "Submission Options" button.')
            logging.error('Missing ENA Username or Password.' + '\n')
            return
        else:
            logging.info('ENA Username and Password look ok.' + '\n')

        # Submission is divided into 3 stages (https://ena-docs.readthedocs.io/en/latest/general-guide/webin-cli.html)
        # Stage 1 - Register Study (Registering a Sample is not necessary for submitting a sequence)
        registerStudy(workingDirectory)

        # Stage 2 - Prepare Files
        prepareSubmissionFiles(submissionText, workingDirectory)

        # Stage 3 - Validate and Submit Files
        validateAndSubmit()

        # Report Results
        # Send a summary of what was submitted.
        # Delete the files.
        # reportAndCleanup()
        # popup message with results.
        # Maybe: delete working directory, but the files first. Maybe delete output files?

    else:
        # TODO: Handle this better, tell the user somehow.
        logging.error('ENA Submission Requirements are not met. Details: (TODO)')

def registerStudy(workingDirectory):
    logging.info('Registering a Project/Study')
    # effectively, study = project

    # The configuration value comes from a radio button in the configuration GUI. existing study = 1, new study = 2
    newProject = (getConfigurationValue('choose_project') == '2')
    if newProject:
        # Generate Project and Project Submission XML Files
        try:
            projectFileName = join(workingDirectory, 'project.xml')
            projectText = createProjectXML(projectFileName)

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
            (projectSubmissionSuccess, projectAccessionNumber, projectErrorMessages) = performProjectSubmission(projectSubmissionFileName, projectFileName)

            if (projectSubmissionSuccess):
                # Great. The project was created successfully.
                # Lets use this new study accession moving forward. Any future submissions in this batch are submitted to the same project.
                assignConfigurationValue('study_accession', projectAccessionNumber)
                assignConfigurationValue('choose_project', '1')

                logging.info('New study has been uploaded, accession:' + str(getConfigurationValue('study_accession')) + '\n')
                logging.info('Subsequent submissions will use this project number.')
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

    # existing project, we will use the supplied accession#. Easy.
    else:
        logging.info('Using existing study accession:' + str(getConfigurationValue('study_accession')) + '\n')
        # projectAccessionNumber = getConfigurationValue('study_accession')
        pass


def prepareSubmissionFiles(submissionText, workingDirectory):
    logging.info('Preparing Submission Files')

    # This includes a "seconds" measure, should be pretty unique.
    # TODO: I could add milliseconds to make this more unique.
    # TODO: I'm worried about batch submissions, is it a problem if the submission files have the same name? I think it's not a problem.
    dateTimeNow = '{:%Y_%m_%d_%H_%M_%S}'.format(datetime.now())

    submissionShortFileName = 'HLA_Submission_' + dateTimeNow + '.txt'
    submissionFileName = join(workingDirectory, submissionShortFileName)

    zippedShortFileName = submissionShortFileName + '.gz'
    zippedFileName = join(workingDirectory, zippedShortFileName)

    manifestShortFileName = 'manifest.txt'
    manifestFileName = join(workingDirectory, manifestShortFileName)

    # Create the submission file
    try:
        outputFileObject = open(submissionFileName, 'w')
        outputFileObject.write(submissionText)
        logging.debug('Submission Text:\n' + submissionText)
        outputFileObject.close()

    except Exception:
        logging.error('Cannot Write Submission Flatfile')
        logging.error(exc_info())
        showInfoBox('Cannot Write Submission Flatfile',
            'Sorry, I failed to create the submission file:\n'
            + str(submissionText)
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
    # TODO: For batch submissions, i shouldn't be able to get the configuration name. Refactor for batch submission.
    manifestFile.write('NAME\t' + str(getConfigurationValue('allele_name')) + '\n')
    manifestFile.write('STUDY\t' + str(getConfigurationValue('study_accession')) + '\n')
    manifestFile.write('FLATFILE\t' + str(zippedFileName) + '\n')
    manifestFile.close()


def validateAndSubmit(workingDirectory):
    logging.info('Validating and Submitting Files.')

    # Find the webin cli application

    # Construct the cli Parameters (for both validation and submission)

    # Validate Sequence

    # Submit Sequence

    # Parse the Result Files. Might need an analysis accession but not sure.






