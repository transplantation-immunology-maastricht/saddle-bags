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




from AlleleSubCommon import *
import pycurl
import StringIO

import xml.etree.ElementTree as ET
import xml.dom.minidom

# Here we have methods to perform REST interactions necessary for EMBL submission.

def performProjectSubmission(submissionFileName, projectFileName):
    POST_DATA = [('SUBMISSION', (pycurl.FORM_FILE, submissionFileName)), 
        ('PROJECT', (pycurl.FORM_FILE, projectFileName))]
    
    responseText = performSubmission(submissionFileName, POST_DATA)    
    return interpretProjectSubmissionResults(responseText)

def performAnalysisSubmission(submissionFileName, analysisFileName):
    POST_DATA = [('SUBMISSION', (pycurl.FORM_FILE, submissionFileName)), 
        ('ANALYSIS', (pycurl.FORM_FILE, analysisFileName))]
    
    responseText = performSubmission(submissionFileName, POST_DATA)
    
    return interpretAnalysisSubmissionResults(responseText)
    
def performSubmission(submissionFileName, POST_DATA):
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
        
def interpretProjectSubmissionResults(responseText):    
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

def interpretAnalysisSubmissionResults(responseText):    
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

