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

#import sys

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

from saddlebags.AlleleSubCommon import getConfigurationValue, createOutputFile, assignConfigurationValue, logEvent
from saddlebags.EmblSubXml import createProjectXML, createProjectSubmissionXML, createAnalysisSubmissionXML , createAnalysisXML



# Here we have methods to perform REST interactions necessary for EMBL submission.

def performProjectSubmission(submissionFileName, projectFileName):
    POST_DATA = [('SUBMISSION', (FORM_FILE, submissionFileName)), 
        ('PROJECT', (FORM_FILE, projectFileName))]
    
    responseText = performSubmission(submissionFileName, POST_DATA)    
    return interpretProjectSubmissionResults(responseText)

def performAnalysisSubmission(submissionFileName, analysisFileName):
    POST_DATA = [('SUBMISSION', (FORM_FILE, submissionFileName)), 
        ('ANALYSIS', (FORM_FILE, analysisFileName))]
    
    responseText = performSubmission(submissionFileName, POST_DATA)
    
    return interpretAnalysisSubmissionResults(responseText)
    
def performSubmission(submissionFileName, POST_DATA):
    
    logEvent('Performing submission of ' + submissionFileName + '\n')
    logEvent('POST Data:\n' + str(POST_DATA) + '\n')
    
    
    if (str(getConfigurationValue('test_submission')) == '0'):
        logEvent ('THIS IS A LIVE SUBMISSION AT EMBL.')
        requestURL = str(getConfigurationValue('embl_rest_address_prod')) + '?auth=ENA%20' + str(getConfigurationValue('embl_username')) + '%20' + str(getConfigurationValue('embl_password'))
    else:
        logEvent ('THIS IS A TEST SUBMISSION AT EMBL.')
        requestURL = str(getConfigurationValue('embl_rest_address_test')) + '?auth=ENA%20' + str(getConfigurationValue('embl_username')) + '%20' + str(getConfigurationValue('embl_password'))
    
    # Problem: StringIO Doesn't work with pycurl in python 3.6. Must replace this with a BytesIO.
    #curlResponseBuffer = StringIO()
    curlResponseBuffer = BytesIO()
    
    
    try:
        curlObject = Curl()
        curlObject.setopt(curlObject.URL, requestURL)
        curlObject.setopt(curlObject.POST, 1)
        curlObject.setopt(curlObject.HTTPPOST, POST_DATA)
        curlObject.setopt(curlObject.USERAGENT, 'Curl')

        curlObject.setopt(curlObject.WRITEFUNCTION, curlResponseBuffer.write)
        
        curlObject.setopt(HTTPHEADER, ['Accept:application/xml'])
        # Insecure.  Any security experts want to make this better?
        curlObject.setopt(SSL_VERIFYHOST, 0)
        curlObject.setopt(SSL_VERIFYPEER, 0)
        curlObject.perform()
        curlObject.close()
    except Exception:
        logEvent ('Exception when performing CURL:\n')
        #logEvent (str(exc_info()))
        logEvent('Exception when performing CURL.\n')
        logEvent('URL:' + str(requestURL))        
        
        raise
    
    responseText = curlResponseBuffer.getvalue()
    
    #logEvent ('the type of the responseText is:' + str(type(responseText)))
    #logEvent ('after it becomes a string:' + str(type(str(responseText))))
    
    # write XML to file. 
    projectSubResultsFileName = submissionFileName.replace('.xml','_results.xml')
    resultsFile = createOutputFile(projectSubResultsFileName)
    resultsFile.write(str(responseText))
    resultsFile.close()
    
    return responseText
        
def interpretProjectSubmissionResults(responseText):
    
    logEvent('Parsing Project Submission Results:\n' + str(responseText) + '\n')
     
       
    # Open XML to report results:
    root = fromstring(responseText)  
    submissionSuccess = (root.attrib['success'] == 'true')
    
    projectAccession = None
    messages = []    

    for child in root:
        if(child.tag == 'PROJECT'):
            if ('accession' in child.attrib.keys()):
                projectAccession = child.attrib['accession']
            else:
                projectAccession = None
            #logEvent('I found a project node.')
        elif(child.tag == 'MESSAGES'):
            logEvent('I found some messages.')
            for messageNode in child:
                #logEvent (messageNode.tag + ':' + messageNode.text)
                messages.append(messageNode.tag + ':' + messageNode.text)
        else:
            # Don't care about the other nodes
            pass
    
    # Return value should be a tuple:
    # (Success, ProjectAccession, Messages[])      
    return (submissionSuccess,projectAccession,messages)

def interpretAnalysisSubmissionResults(responseText):
    
    logEvent('Parsing Analysis Submission Results:\n' + str(responseText) + '\n')
    
        
    root = fromstring(responseText)  
    submissionSuccess = (root.attrib['success'] == 'true')
    
    analysisAccession = None
    messages = []    

    for child in root:
        if(child.tag == 'ANALYSIS'):
            if ('accession' in child.attrib.keys()):
                analysisAccession = child.attrib['accession']
            else:
                analysisAccession = None
            #logEvent('I found a project node.')
        elif(child.tag == 'MESSAGES'):
            logEvent('I found some messages.')
            for messageNode in child:
                #logEvent (messageNode.tag + ':' + messageNode.text)
                messages.append(messageNode.tag + ':' + messageNode.text)
        else:
            # Don't care about the other nodes
            pass
    
    # Return value should be a tuple:
    # (Success, ProjectAccession, Messages[])      
    return (submissionSuccess,analysisAccession,messages)

        
def writeMd5(inputFileName, outputFileName):
    hash_md5 = md5()
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
    logEvent('Uploading Submission to EMBL')
    
    # Determine a working directory. Folder underneath executable called temp.
    try:
        workingDirectory = join(expanduser("~"), 'temp_saddlebags')
        logEvent('I can work in this directory:' + workingDirectory)
        
        if not isdir(workingDirectory):
            logEvent('Making Directory:' + workingDirectory)
            makedirs(workingDirectory)
    except Exception:
        logEvent ('Cannot Initialize Working Directory')
        logEvent (exc_info())
        messagebox.showinfo('Working Directory Error', 
            'Sorry, I failed to create this working directory:\n'
            + str(workingDirectory)
            + '\n and I cannot continue.\nMaybe this is a '
            + 'permissions issue, are these folders read only?\n' 
            +  str(exc_info()[1]))
        return
    
    # I'm deprecating this log file. 
    # Use the logging i put into AlleleSubCommon instead.
    # logEvent.
    #restLog = createOutputFile(join(workingDirectory, 'Submission_Log.txt'))

    emblUsername = getConfigurationValue('embl_username')
    emblPassword = getConfigurationValue('embl_password')
    if(emblUsername is None 
        or len(emblUsername) < 1
        or emblPassword is None 
        or len(emblPassword) < 1):
        messagebox.showinfo('Missing Login Credentials', 
            'You must provide EMBL username and password.\n'
            'Please use the "Submission Options" button.')
        logEvent('Missing EMBL Username or Password.' + '\n')
        return
    else:
        logEvent('EMBL Username and Password exist.' + '\n')
       

    useTestServers = (int(getConfigurationValue('test_submission')) == 1)
    # Are you sure?
    if useTestServers:
        logEvent('Using Test EMBL Server.' + '\n')
        result = messagebox.askquestion("Submit to TEST / DEMO environment", "You are about to submit a sequence to the\n\nTEST / DEMO EMBL environment.\n\nAre You Sure?", icon='warning')
    else:
        logEvent('Using Production EMBL Server.' + '\n')
        result = messagebox.askquestion("Submit to LIVE / PROD environment", "You are about to submit a sequence to the\n\nLIVE / PROD EMBL environment.\n\nAre You Sure?", icon='warning')

    if result == 'yes':
        pass
    else:
        return
    
    # TODO: Existing project? Maybe I should check if the study/project exists, before I get started
    
    
    # Give my submission a filename. SOmething with a datetime stamp
    try:
        # This includes a "seconds" measure, should be pretty unique.
        dateTimeNow = '{:%Y_%m_%d_%H_%M_%S}'.format(datetime.now())
        submissionShortFileName = 'HLA_Submission_' + dateTimeNow + '.txt'
        submissionFileName = join(workingDirectory, submissionShortFileName)
        zippedShortFileName = submissionShortFileName + '.gz'
        zippedFileName = join(workingDirectory, zippedShortFileName)
        md5FileName = zippedFileName + '.md5'
  
        #submissionText = self.submOutputGuiObject.get('1.0', 'end')         
        
        outputFileObject = open(submissionFileName, 'w') 
        outputFileObject.write(submissionText) 
        logEvent('Submission Text:\n' + submissionText)
        outputFileObject.close()        
    
    except Exception:
        logEvent ('Cannot Write Submission Flatfile')
        logEvent (exc_info())
        messagebox.showinfo('Cannot Write Submission Flatfile', 
            'Sorry, I failed to create the submission file:\n'
            + str(submissionText)
            + '\n and I cannot continue.\nMaybe this is a '
            + 'permissions issue, are these folders read only?\n' 
            +  str(exc_info()[1]))
        logEvent('Failure to create submission file:' + str(exc_info()[1]) + '\n')
        return
    
    logEvent('Submission file was created:\n' + str(submissionFileName) + '\n')
    
    # gzip the submission file.  Make a gz file.
    try:
        #zippedFileName = submissionFileName + '.gz'
        
        with open(submissionFileName, 'rb') as fileIn, gzipOpen(zippedFileName, 'wb') as fileOut:
            copyfileobj(fileIn, fileOut)
    
    except Exception:
        logEvent ('Cannot Compress Submission File')
        logEvent (exc_info())
        messagebox.showinfo('Cannot Compress Submission File', 
            'Sorry, I failed to compress the submission file:\n'
            + str(zippedFileName)
            + '\n and I cannot continue.\n' 
            +  str(exc_info()[1]))
        logEvent('Failure to create zip file:' + str(exc_info()[1]) + '\n')
        return
    
    logEvent('Zip file was created:\n' + str(zippedFileName) + '\n')
    
    # Calculate an MD5SUM
    try:
        #md5FileName = zippedFileName + '.md5'
        md5HashValue = writeMd5(zippedFileName,md5FileName)
        
    except Exception:
        logEvent ('Cannot Calculate MD5')
        logEvent (exc_info())
        messagebox.showinfo('Cannot Calculate an Md5 checksum', 
            'Sorry, I failed to calculate an md5 checksum\nand I cannot continue.\n' 
            +  str(exc_info()[1]))
        logEvent('Failure to create zip file:' + str(exc_info()[1]) + '\n')
        return
    
    logEvent('md5 file was created:\n' + str(md5FileName) + '\n')

    # Use FTP  to send the file to EMBL
    try:
        if useTestServers:
            ftpServerAddress = getConfigurationValue('embl_ftp_upload_site_test')      
        else:
            ftpServerAddress = getConfigurationValue('embl_ftp_upload_site_prod')   
        
        #logEvent ('attempting to open ftp connection')
        ftp = FTP(ftpServerAddress)
        ftp.login(getConfigurationValue('embl_username'), getConfigurationValue('embl_password'))
        ftp.storbinary('STOR ' + '/' + split(zippedFileName)[1], open(zippedFileName, 'rb'), 1024)
        ftp.storbinary('STOR ' + '/' + split(md5FileName)[1], open(md5FileName, 'rb'), 1024)
        ftp.close()
        # is that it?  Easy.

    except Exception:
        logEvent ('Cannot Upload to FTP site')
        logEvent (exc_info())
        messagebox.showinfo('Cannot Upload to FTP site', 
            'Sorry, I failed to upload your submission files to the EMBL FTP site\nand I cannot continue.\n' 
            +  str(exc_info()[1]))
        logEvent('Failure to upload to FTP site:' + str(exc_info()[1]) + '\n')
        return
    
    logEvent('Submission and MD5 successfully uploaded.\n')
    
    # Handle the new project
    # effectively, study = project 
    # existing study = 1, new study = 2
    newProject = (getConfigurationValue('choose_project') == '2')
    if newProject:
        
        # Generate Project and Project Submission XML Files
        try:
            projectFileName = join(workingDirectory, 'project.xml')
            projectText = createProjectXML(projectFileName)
            
            #logEvent ('My project text looks like this:\n' + projectText)
            
            projectSubmissionFileName = join(workingDirectory, 'project_submission.xml')
            projectSubmissionText = createProjectSubmissionXML(projectSubmissionFileName
                ,'proj_sub_' + dateTimeNow
                ,'project.xml')
            
            #logEvent('I made this project text:\n' + projectText)
            #logEvent('I made this project submission text:\n' + projectSubmissionText)
            
        except Exception:
            logEvent ('Cannot Create Project Submission XML')
            logEvent (exc_info())
            messagebox.showinfo('Cannot Create Project Submission XML', 
                'Sorry, I failed to create a project XML file\nand I cannot continue.\n' 
                +  str(exc_info()[1]))
            logEvent('Failure to create project submission file:' + str(exc_info()[1]) + '\n')
            return
        
        logEvent('Project Submission XML files were created.\n')
        logEvent('Project Text:\n' + projectText)
        logEvent('Project Submission Text:\n' + projectText)
        
                    
        # Use REST to submit this project
        try:
            # Return value should be a tuple:
            # (Success, ProjectAccession, Messages[])   
            (projectSubmissionSuccess, projectAccessionNumber, projectErrorMessages) = performProjectSubmission(projectSubmissionFileName,projectFileName)
            
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
                messagebox.showinfo('Cannot Submit Project XML via REST', messageText)
                logEvent('Failure to submit project submission file:' + str(exc_info()[1]) + '\n' + messageText + '\n')
                return
            
        except Exception:
            logEvent ('Cannot Submit Project XML')
            logEvent (exc_info())
            messagebox.showinfo('Cannot Submit Project XML', 
                'Sorry, I failed to submit the project XML file\nand I cannot continue.\n' 
                +  str(exc_info()[1]))
            logEvent('Failure to upload project submission file:' + str(exc_info()[1]) + '\n')
            return
        
        logEvent('New study has been uploaded, accession:' + str(getConfigurationValue('study_accession')) + '\n')
           
    # existing project, we will use the supplied accession #    
    else: 
        logEvent('Using existing study accession:' + str(getConfigurationValue('study_accession')) + '\n')
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
        logEvent('Cannot Create Analysis Submission XML')
        logEvent (exc_info())
        messagebox.showinfo('Cannot Create Analysis Submission XML', 
            'Sorry, I failed to create a Analysis XML file\nand I cannot continue.\n' 
            +  str(exc_info()[1]))
        logEvent('Failure to create analysis submission file:' + str(exc_info()[1]) + '\n')
        return
    
    logEvent('Analysis Submission XML files were created.\n')
                
    # Use REST to submit this analysis
    try:
        # Return value should be a tuple:
        # (Success, analysisAccessionNumber, Messages[])   
        (analysisSubmissionSuccess, analysisAccessionNumber, analysisErrorMessages) = performAnalysisSubmission(analysisSubmissionFileName,analysisFileName)
        
        if(analysisSubmissionSuccess):
            # Great. The analysis was created successfully. 
            pass
        else:
            messageText = ('There was a problem in the Analysis Submission.\n' 
                + 'I cannot continue.\n'
                + 'These messages were reported by EMBL:\n')
            for errorMessage in analysisErrorMessages:
                messageText += ('\n' + errorMessage + '\n')                    
            messagebox.showinfo('Cannot Submit Analysis XML via REST', messageText)
            logEvent('Failure to submit analysis submission file:' + str(exc_info()[1]) + '\n')
            return
        
    except Exception:
        logEvent ('Cannot Submit Analysis XML')
        logEvent (exc_info())
        messagebox.showinfo('Cannot Submit Analysis XML via REST', 
            'Sorry, I failed to submit the analysis XML file\nand I cannot continue.\n' 
            +  str(exc_info()[1]))
        return

    logEvent('New analysis has been Uploaded, accession:' + str(analysisAccessionNumber) + '\n')


    # Popup message with Results
    messagebox.showinfo('Success uploading submission to EMBL.', 
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