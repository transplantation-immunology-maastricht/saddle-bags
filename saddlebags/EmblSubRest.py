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

import sys

from AlleleSubCommon import getConfigurationValue, createOutputFile
import pycurl
import StringIO

import datetime
import hashlib
import ftplib
import gzip
import shutil

from AlleleSubCommon import assignConfigurationValue
from EmblSubXml import createProjectXML, createProjectSubmissionXML, createAnalysisSubmissionXML , createAnalysisXML

from os import makedirs
from os.path import join, isdir, split, expanduser

import xml.etree.ElementTree as ET

import tkMessageBox
#import xml.dom.minidom

# Here we have methods to perform REST interactions necessary for EMBL submission.

def performProjectSubmission(submissionFileName, projectFileName, restLog):
    POST_DATA = [('SUBMISSION', (pycurl.FORM_FILE, submissionFileName)), 
        ('PROJECT', (pycurl.FORM_FILE, projectFileName))]
    
    responseText = performSubmission(submissionFileName, POST_DATA, restLog)    
    return interpretProjectSubmissionResults(responseText, restLog)

def performAnalysisSubmission(submissionFileName, analysisFileName, restLog):
    POST_DATA = [('SUBMISSION', (pycurl.FORM_FILE, submissionFileName)), 
        ('ANALYSIS', (pycurl.FORM_FILE, analysisFileName))]
    
    responseText = performSubmission(submissionFileName, POST_DATA, restLog)
    
    return interpretAnalysisSubmissionResults(responseText, restLog)
    
def performSubmission(submissionFileName, POST_DATA, restLog):
    
    restLog.write('Performing submission of ' + submissionFileName + '\n')
    restLog.write('POST Data:\n' + str(POST_DATA) + '\n')
    
    
    if (str(getConfigurationValue('test_submission')) == '0'):
        print 'THIS IS A LIVE SUBMISSION AT EMBL.'
        requestURL = str(getConfigurationValue('embl_rest_address_prod')) + '?auth=ENA%20' + str(getConfigurationValue('embl_username')) + '%20' + str(getConfigurationValue('embl_password'))
    else:
        print 'THIS IS A TEST SUBMISSION AT EMBL.'
        requestURL = str(getConfigurationValue('embl_rest_address_test')) + '?auth=ENA%20' + str(getConfigurationValue('embl_username')) + '%20' + str(getConfigurationValue('embl_password'))
    
    curlResponseBuffer = StringIO.StringIO()
    curlObject = pycurl.Curl()
    curlObject.setopt(curlObject.URL, requestURL)
    curlObject.setopt(curlObject.POST, 1)
    curlObject.setopt(curlObject.HTTPPOST, POST_DATA)
    curlObject.setopt(curlObject.USERAGENT, 'Curl')
    curlObject.setopt(curlObject.WRITEFUNCTION, curlResponseBuffer.write)
    curlObject.setopt(pycurl.HTTPHEADER, ['Accept:application/xml'])
    # Insecure.  Any security experts want to make this better?
    curlObject.setopt(pycurl.SSL_VERIFYHOST, 0)
    curlObject.setopt(pycurl.SSL_VERIFYPEER, 0)
    curlObject.perform()
    curlObject.close()

    responseText = curlResponseBuffer.getvalue()
    
    # write XML to file. 
    projectSubResultsFileName = submissionFileName.replace('.xml','_results.xml')
    resultsFile = createOutputFile(projectSubResultsFileName)
    resultsFile.write(responseText)
    resultsFile.close()
    
    return responseText
        
def interpretProjectSubmissionResults(responseText, restLog):
    
    restLog.write('Parsing Project Submission Results:\n' + str(responseText) + '\n')
     
       
    # Open XML to report results:
    root = ET.fromstring(responseText)  
    submissionSuccess = (root.attrib['success'] == 'true')
    
    projectAccession = None
    messages = []    

    for child in root:
        if(child.tag == 'PROJECT'):
            if ('accession' in child.attrib.keys()):
                projectAccession = child.attrib['accession']
            else:
                projectAccession = None
            #print('I found a project node.')
        elif(child.tag == 'MESSAGES'):
            print('I found some messages.')
            for messageNode in child:
                #print (messageNode.tag + ':' + messageNode.text)
                messages.append(messageNode.tag + ':' + messageNode.text)
        else:
            # Don't care about the other nodes
            pass
    
    # Return value should be a tuple:
    # (Success, ProjectAccession, Messages[])      
    return (submissionSuccess,projectAccession,messages)

def interpretAnalysisSubmissionResults(responseText, restLog):
    
    restLog.write('Parsing Analysis Submission Results:\n' + str(responseText) + '\n')
    
        
    root = ET.fromstring(responseText)  
    submissionSuccess = (root.attrib['success'] == 'true')
    
    analysisAccession = None
    messages = []    

    for child in root:
        if(child.tag == 'ANALYSIS'):
            if ('accession' in child.attrib.keys()):
                analysisAccession = child.attrib['accession']
            else:
                analysisAccession = None
            #print('I found a project node.')
        elif(child.tag == 'MESSAGES'):
            print('I found some messages.')
            for messageNode in child:
                #print (messageNode.tag + ':' + messageNode.text)
                messages.append(messageNode.tag + ':' + messageNode.text)
        else:
            # Don't care about the other nodes
            pass
    
    # Return value should be a tuple:
    # (Success, ProjectAccession, Messages[])      
    return (submissionSuccess,analysisAccession,messages)

        
def writeMd5(inputFileName, outputFileName):
    hash_md5 = hashlib.md5()
    with open(inputFileName, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    hashValue= hash_md5.hexdigest()
    
    outputFile = createOutputFile(outputFileName)
    # The Ubuntu md5sum program seems to write a single checksum and filename with 2 spaces between
    # I don't know why 2 spaces, but I'll roll with it.
    outputFile.write(str(hashValue) + '  ' + str(split(inputFileName)[1]))
    outputFile.close()
    
    return hashValue



# TODO: This method is long, can i clean it up at all?

def performFullSubmission(submissionText):
    print('Uploading Submission to EMBL')
    
    # Determine a working directory. Folder underneath executable called temp.
    try:
        workingDirectory = join(expanduser("~"), 'temp_upload_directory')
        print('I can work in this directory:' + workingDirectory)
        
        if not isdir(workingDirectory):
            print('Making Directory:' + workingDirectory)
            makedirs(workingDirectory)
    except Exception:
        print 'Cannot Initialize Working Directory'
        print sys.exc_info()[1]
        tkMessageBox.showinfo('Working Directory Error', 
            'Sorry, I failed to create this working directory:\n'
            + str(workingDirectory)
            + '\n and I cannot continue.\nMaybe this is a '
            + 'permissions issue, are these folders read only?\n' 
            +  str(sys.exc_info()[1]))
        return
    
    restLog = createOutputFile(join(workingDirectory, 'Submission_Log.txt'))
    
    
    
    # TODO: Make a REST log.
    # For each step report success or failure.  Same as popup messages.
    # Looks like I did thisalready, probably delete this TODO
    
    

    emblUsername = getConfigurationValue('embl_username')
    emblPassword = getConfigurationValue('embl_password')
    if(emblUsername is None 
        or len(emblUsername) < 1
        or emblPassword is None 
        or len(emblPassword) < 1):
        tkMessageBox.showinfo('Missing Login Credentials', 
            'You must provide EMBL username and password.\n'
            'Please use the "Submission Options" button.')
        restLog.write('Missing EMBL Username or Password.' + '\n')
        return
    else:
        restLog.write('EMBL Username and Password exist.' + '\n')
       

    useTestServers = (int(getConfigurationValue('test_submission')) == 1)
    # Are you sure?
    if useTestServers:
        restLog.write('Using Test EMBL Server.' + '\n')
        result = tkMessageBox.askquestion("Submit to TEST / DEMO environment", "You are about to submit a sequence to the\n\nTEST / DEMO EMBL environment.\n\nAre You Sure?", icon='warning')
    else:
        restLog.write('Using Production EMBL Server.' + '\n')
        result = tkMessageBox.askquestion("Submit to LIVE / PROD environment", "You are about to submit a sequence to the\n\nLIVE / PROD EMBL environment.\n\nAre You Sure?", icon='warning')

    if result == 'yes':
        pass
    else:
        return
    
    # TODO: Existing project? Maybe I should check if the study/project exists, before I get started
    
    
    # Give my submission a filename. SOmething with a datetime stamp
    try:
        # This includes a "seconds" measure, should be pretty unique.
        dateTimeNow = '{:%Y_%m_%d_%H_%M_%S}'.format(datetime.datetime.now())
        submissionShortFileName = 'HLA_Submission_' + dateTimeNow + '.txt'
        submissionFileName = join(workingDirectory, submissionShortFileName)
        zippedShortFileName = submissionShortFileName + '.gz'
        zippedFileName = join(workingDirectory, zippedShortFileName)
        md5FileName = zippedFileName + '.md5'
  
        #submissionText = self.submOutputGuiObject.get('1.0', 'end')         
        
        outputFileObject = open(submissionFileName, 'w') 
        outputFileObject.write(submissionText) 
        restLog.write('Submission Text:\n' + submissionText)
        outputFileObject.close()        
    
    except Exception:
        print 'Cannot Write Submission Flatfile'
        print sys.exc_info()[1]
        tkMessageBox.showinfo('Cannot Write Submission Flatfile', 
            'Sorry, I failed to create the submission file:\n'
            + str(submissionText)
            + '\n and I cannot continue.\nMaybe this is a '
            + 'permissions issue, are these folders read only?\n' 
            +  str(sys.exc_info()[1]))
        restLog.write('Failure to create submission file:' + str(sys.exc_info()[1]) + '\n')
        return
    
    restLog.write('Submission file was created:\n' + str(submissionFileName) + '\n')
    
    # gzip the submission file.  Make a gz file.
    try:
        #zippedFileName = submissionFileName + '.gz'
        
        with open(submissionFileName, 'rb') as fileIn, gzip.open(zippedFileName, 'wb') as fileOut:
            shutil.copyfileobj(fileIn, fileOut)
    
    except Exception:
        print 'Cannot Compress Submission File'
        print sys.exc_info()[1]
        tkMessageBox.showinfo('Cannot Compress Submission File', 
            'Sorry, I failed to compress the submission file:\n'
            + str(zippedFileName)
            + '\n and I cannot continue.\n' 
            +  str(sys.exc_info()[1]))
        restLog.write('Failure to create zip file:' + str(sys.exc_info()[1]) + '\n')
        return
    
    restLog.write('Zip file was created:\n' + str(zippedFileName) + '\n')
    
    # Calculate an MD5SUM
    try:
        #md5FileName = zippedFileName + '.md5'
        md5HashValue = writeMd5(zippedFileName,md5FileName)
        
    except Exception:
        print 'Cannot Calculate MD5'
        print sys.exc_info()[1]
        tkMessageBox.showinfo('Cannot Calculate an Md5 checksum', 
            'Sorry, I failed to calculate an md5 checksum\nand I cannot continue.\n' 
            +  str(sys.exc_info()[1]))
        restLog.write('Failure to create zip file:' + str(sys.exc_info()[1]) + '\n')
        return
    
    restLog.write('md5 file was created:\n' + str(md5FileName) + '\n')

    # Use FTP  to send the file to EMBL
    try:
        if useTestServers:
            ftpServerAddress = getConfigurationValue('embl_ftp_upload_site_test')      
        else:
            ftpServerAddress = getConfigurationValue('embl_ftp_upload_site_prod')   
        
        #print ('attempting to open ftp connection')
        ftp = ftplib.FTP(ftpServerAddress)
        ftp.login(getConfigurationValue('embl_username'), getConfigurationValue('embl_password'))
        ftp.storbinary('STOR ' + '/' + split(zippedFileName)[1], open(zippedFileName, 'rb'), 1024)
        ftp.storbinary('STOR ' + '/' + split(md5FileName)[1], open(md5FileName, 'rb'), 1024)
        ftp.close()
        # is that it?  Easy.

    except Exception:
        print 'Cannot Upload to FTP site'
        print sys.exc_info()[1]
        tkMessageBox.showinfo('Cannot Upload to FTP site', 
            'Sorry, I failed to upload your submission files to the EMBL FTP site\nand I cannot continue.\n' 
            +  str(sys.exc_info()[1]))
        restLog.write('Failure to upload to FTP site:' + str(sys.exc_info()[1]) + '\n')
        return
    
    restLog.write('Submission and MD5 successfully uploaded.\n')
    
    # Handle the new project
    # effectively, study = project 
    # existing study = 1, new study = 2
    newProject = (getConfigurationValue('choose_project') == '2')
    if newProject:
        
        # Generate Project and Project Submission XML Files
        try:
            projectFileName = join(workingDirectory, 'project.xml')
            projectText = createProjectXML(projectFileName)
            
            projectSubmissionFileName = join(workingDirectory, 'project_submission.xml')
            projectSubmissionText = createProjectSubmissionXML(projectSubmissionFileName
                ,'proj_sub_' + dateTimeNow
                ,'project.xml')
            
            #print('I made this project text:\n' + projectText)
            #print('I made this project submission text:\n' + projectSubmissionText)
            
        except Exception:
            print 'Cannot Create Project Submission XML'
            print sys.exc_info()[1]
            tkMessageBox.showinfo('Cannot Create Project Submission XML', 
                'Sorry, I failed to create a project XML file\nand I cannot continue.\n' 
                +  str(sys.exc_info()[1]))
            restLog.write('Failure to create project submission file:' + str(sys.exc_info()[1]) + '\n')
            return
        
        restLog.write('Project Submission XML files were created.\n')
        restLog.write('Project Text:\n' + projectText)
        restLog.write('Project Submission Text:\n' + projectText)
        
                    
        # Use REST to submit this project
        try:
            # Return value should be a tuple:
            # (Success, ProjectAccession, Messages[])   
            (projectSubmissionSuccess, projectAccessionNumber, projectErrorMessages) = performProjectSubmission(projectSubmissionFileName,projectFileName, restLog)
            
            if(projectSubmissionSuccess):
                # Great. The project was created successfully. 
                # Lets use this new study accession moving forward.
                assignConfigurationValue('study_accession', projectAccessionNumber)
                assignConfigurationValue('choose_project','1')
                pass
            else:
                messageText = ('There was a problem in the Project Submission.\n' 
                    + 'I cannot continue.\n'
                    + 'These messages were reported by EMBL:\n')
                for errorMessage in projectErrorMessages:
                    messageText += ('\n' + errorMessage + '\n')                    
                tkMessageBox.showinfo('Cannot Submit Project XML via REST', messageText)
                restLog.write('Failure to submit project submission file:' + str(sys.exc_info()[1]) + '\n' + messageText + '\n')
                return
            
        except Exception:
            print 'Cannot Submit Project XML'
            print sys.exc_info()[1]
            tkMessageBox.showinfo('Cannot Submit Project XML', 
                'Sorry, I failed to submit the project XML file\nand I cannot continue.\n' 
                +  str(sys.exc_info()[1]))
            restLog.write('Failure to upload project submission file:' + str(sys.exc_info()[1]) + '\n')
            return
        
        restLog.write('New study has been uploaded, accession:' + str(getConfigurationValue('study_accession')) + '\n')
           
    # existing project, we will use the supplied accession #    
    else: 
        restLog.write('Using existing study accession:' + str(getConfigurationValue('study_accession')) + '\n')
        # projectAccessionNumber = getConfigurationValue('study_accession')
        pass
    
    # Generate Analysis and Analysis Submission xmls
    try:
        analysisFileName = join(workingDirectory, 'analysis.xml')
        analysisText = createAnalysisXML(analysisFileName, md5HashValue, zippedShortFileName)
        
        analysisSubmissionFileName = join(workingDirectory, 'analysis_submission.xml')
        analysisSubmissionText = createAnalysisSubmissionXML(analysisSubmissionFileName
            ,'analysis_sub_' + dateTimeNow
            ,'analysis.xml')
        
    except Exception:
        print 'Cannot Create Analysis Submission XML'
        print sys.exc_info()[1]
        tkMessageBox.showinfo('Cannot Create Analysis Submission XML', 
            'Sorry, I failed to create a Analysis XML file\nand I cannot continue.\n' 
            +  str(sys.exc_info()[1]))
        restLog.write('Failure to create analysis submission file:' + str(sys.exc_info()[1]) + '\n')
        return
    
    restLog.write('Analysis Submission XML files were created.\n')
                
    # Use REST to submit this analysis
    try:
        # Return value should be a tuple:
        # (Success, analysisAccessionNumber, Messages[])   
        (analysisSubmissionSuccess, analysisAccessionNumber, analysisErrorMessages) = performAnalysisSubmission(analysisSubmissionFileName,analysisFileName, restLog)
        
        if(analysisSubmissionSuccess):
            # Great. The analysis was created successfully. 
            pass
        else:
            messageText = ('There was a problem in the Analysis Submission.\n' 
                + 'I cannot continue.\n'
                + 'These messages were reported by EMBL:\n')
            for errorMessage in analysisErrorMessages:
                messageText += ('\n' + errorMessage + '\n')                    
            tkMessageBox.showinfo('Cannot Submit Analysis XML via REST', messageText)
            restLog.write('Failure to submit analysis submission file:' + str(sys.exc_info()[1]) + '\n')
            return
        
    except Exception:
        print 'Cannot Submit Analysis XML'
        print sys.exc_info()[1]
        tkMessageBox.showinfo('Cannot Submit Analysis XML via REST', 
            'Sorry, I failed to submit the analysis XML file\nand I cannot continue.\n' 
            +  str(sys.exc_info()[1]))
        return

    restLog.write('New analysis has been Uploaded, accession:' + str(analysisAccessionNumber) + '\n')

    restLog.close()

    # Popup message with Results
    tkMessageBox.showinfo('Success uploading submission to EMBL.', 
        'The sequence and analysis was uploaded to EMBL ENA Successfully.\n\n' 
        + 'For your reference:\n\n'
        + 'You can use this Project/Study accession\nnumber on future submissions:\n'
        + 'Study Accession:' + str(getConfigurationValue('study_accession') + '\n\n')
        + 'Use the Analysis Accession number if you\ncontact EMBL regarding this\nsequence submission:\n'
        + 'Analysis Accession:' + str(analysisAccessionNumber) + '\n\n'
        + 'Find your submission files here:\n'
        + workingDirectory + '\n\n'
        + 'If EMBL successfully validates your sequence, you will\n'
        + 'receive an email with an EMBL Sequence accession number.\n'
        + 'This accession number is necessary for IMGT/HLA submission.\n'
        + 'Contact EMBL Support with your\nAnalysis Accession # if it has been\nmore than 48 hours since submission.\n'

        )