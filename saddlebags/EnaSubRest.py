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

from pycurl import Curl, FORM_FILE, HTTPHEADER, SSL_VERIFYHOST, SSL_VERIFYPEER
from io import BytesIO
from xml.etree.ElementTree import fromstring
import logging
from saddlebags.AlleleSubCommon import createOutputFile
#getConfigurationValue,

# Here we have methods to perform REST interactions necessary for ENA submission.
# Some REST functionality was removed by ENA, but I must still perform a few things.
# I need to "Register Study",  Study = Project.
# Analysis is no longer needed, those are created by the Webin CLI tool automatically.


def performProjectSubmission(submissionFileName, projectFileName, submissionBatch):
    POST_DATA = [('SUBMISSION', (FORM_FILE, submissionFileName)), 
        ('PROJECT', (FORM_FILE, projectFileName))]
    
    responseText = performSubmission(submissionFileName, POST_DATA, submissionBatch.enaUserName, submissionBatch.enaPassword)
    return interpretProjectSubmissionResults(responseText)

def performSubmission(submissionFileName, POST_DATA, enaUserName, enaPassword):
    
    logging.info('Performing submission of ' + submissionFileName + '\n')
    logging.info('POST Data:\n' + str(POST_DATA) + '\n')
    
    
    if (str(getConfigurationValue('test_submission')) == '0'):
        logging.info ('THIS IS A LIVE SUBMISSION AT ENA.')
        requestURL = str(getConfigurationValue('ena_rest_address_prod')) + '?auth=ENA%20' + str(enaUserName) + '%20' + str(enaPassword)
    else:
        logging.info ('THIS IS A TEST SUBMISSION AT ENA.')
        requestURL = str(getConfigurationValue('ena_rest_address_test')) + '?auth=ENA%20' + str(enaUserName) + '%20' + str(enaPassword)
    
    # Problem: StringIO Doesn't work with pycurl in python 3.6. Must replace this with a BytesIO.
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
        logging.error ('Exception when performing CURL:\n')
        #logging.error (str(exc_info()))
        logging.error('Exception when performing CURL.\n')
        logging.error('URL:' + str(requestURL))
        
        raise
    
    responseText = curlResponseBuffer.getvalue()
    
    #logging.info ('the type of the responseText is:' + str(type(responseText)))
    #logging.info ('after it becomes a string:' + str(type(str(responseText))))
    
    # write XML to file. 
    projectSubResultsFileName = submissionFileName.replace('.xml','_results.xml')
    resultsFile = createOutputFile(projectSubResultsFileName)
    resultsFile.write(str(responseText))
    resultsFile.close()
    
    return responseText
        
def interpretProjectSubmissionResults(responseText):
    
    logging.info('Parsing Project Submission Results:\n' + str(responseText) + '\n')
       
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
            #logging.info('I found a project node.')
        elif(child.tag == 'MESSAGES'):
            logging.debug('I found some messages.')
            for messageNode in child:
                #logging.debug (messageNode.tag + ':' + messageNode.text)
                messages.append(messageNode.tag + ':' + messageNode.text)
        else:
            # Don't care about the other nodes
            pass
    
    # Return value should be a tuple:
    # (Success, ProjectAccession, Messages[])      
    return (submissionSuccess,projectAccession,messages)


def interpretAnalysisSubmissionResults(responseText):
    logging.info('Parsing Analysis Submission Results:\n' + str(responseText) + '\n')

    root = fromstring(responseText)
    submissionSuccess = (root.attrib['success'] == 'true')

    analysisAccession = None
    messages = []

    for child in root:
        if (child.tag == 'ANALYSIS'):
            if ('accession' in child.attrib.keys()):
                analysisAccession = child.attrib['accession']
            else:
                analysisAccession = None
            # logging.info('I found a project node.')
        elif (child.tag == 'MESSAGES'):
            logging.info('I found some messages.')
            for messageNode in child:
                # logging.info (messageNode.tag + ':' + messageNode.text)
                messages.append(messageNode.tag + ':' + messageNode.text)
        else:
            # Don't care about the other nodes
            pass

    # Return value should be a tuple:
    # (Success, ProjectAccession, Messages[])
    return (submissionSuccess, analysisAccession, messages)


