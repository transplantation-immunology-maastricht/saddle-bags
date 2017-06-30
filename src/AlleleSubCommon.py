# This file is part of EMBL-HLA-Submission.
#
# EMBL-HLA-Submission is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EMBL-HLA-Submission is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with EMBL-HLA-Submission. If not, see <http://www.gnu.org/licenses/>.


#from numpy.compat.setup import configuration

#SoftwareVersion = "Bhast Version 1.0"

import xml.etree.ElementTree as ET
import xml.dom.minidom

from os.path import isdir, split
from os import makedirs


import sys
from os.path import dirname, join, abspath, isfile

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput





# I'm storing global variables in a dictionary for now. 
def initializeGlobalVariables():    
    global globalVariables 
    
    if not ("globalVariables" in globals()):
        globalVariables={}
        
def assignConfigurationValue(configurationKey, configurationValue):
    initializeGlobalVariables()
    globalVariables[configurationKey] = configurationValue
    
def getConfigurationValue(configurationKey):
    if configurationKey in globalVariables.keys():
        return globalVariables[configurationKey]
    else:
        print ('Configuration Key Not Found:' + configurationKey)
        #raise KeyError('Key Not Found:' + configurationKey)
        return None
    

def assignConfigName():
    initializeGlobalVariables()

    # Find the directory the program is running from.
    # It is not straight-forward, because sometimes we run this program inside an .exe
    # pyinstaller puts the exe path in sys._MEIPASS
    # This is useful because we want a config file in the same directory.
    if getattr(sys, 'frozen', False):
        globalVariables['saddlebags_application_path'] = sys._MEIPASS
    else:
        globalVariables['saddlebags_application_path'] = dirname(abspath(__file__))
        
    # TODO: Store the directory someone saves in.  
    # I should assign the directory to a default value.
        
    print 'This application is running from the following directory:\n' + globalVariables['saddlebags_application_path']
    globalVariables['config_file_location'] = join(globalVariables['saddlebags_application_path'], 'Saddlebags.Config.xml')
    
    
def writeConfigurationFile():
    assignConfigName()
    print ('Writing a config file to:\n' + globalVariables['config_file_location'])
    
    root = ET.Element("config")
    
    for key in globalVariables.keys():
        # Some config values I don't want to store.
        # Add to this: Sequence Text, EMBL Submission Text, IMGT Submission Text
        if(key not in [
            'embl_password'
            ,'imgt_password'
            , 'saddlebags_application_path'
            , 'config_file_location'
            , 'sequence'
            ]):
            ET.SubElement(root, key).text = globalVariables[key]

    xmlText = ET.tostring(root, encoding='utf8', method='xml')
    prettyXmlText = xml.dom.minidom.parseString(xmlText).toprettyxml()
    
    xmlOutput = createOutputFile(globalVariables['config_file_location'])
    xmlOutput.write(prettyXmlText)
    xmlOutput.close()


def loadConfigurationFile():
    assignConfigName()
    
    if not isfile(globalVariables['config_file_location']):
        print ('The config file does not exist yet. I will not load it:\n' + globalVariables['config_file_location'])
        
        # Here is where I assign the common configuration values
        # test_submission indicates if we should use the "test" values.
        # I think I'll use this value for both EMBL and IMGT submissions, if it applies.
        assignConfigurationValue('test_submission', '1')
        
        # I'm storing FTP without the ftp:// identifier, because it is not necessary.
        # The test and prod ftp sites have the same address. This is intentional, embl doesn't have a test ftp
        assignConfigurationValue('embl_ftp_upload_site_test', 'webin.ebi.ac.uk')
        assignConfigurationValue('embl_ftp_upload_site_prod', 'webin.ebi.ac.uk')
        assignConfigurationValue('embl_rest_address_test', 'https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/')
        assignConfigurationValue('embl_rest_address_prod', 'https://www.ebi.ac.uk/ena/submit/drop-box/submit/')
      
        
        
    else:
        print ('The config file already exists, I will load it:\n' + globalVariables['config_file_location'])
        
        tree = ET.parse(globalVariables['config_file_location'])
        root = tree.getroot()
        
        for child in root:
            assignConfigurationValue(child.tag, child.text)
            
       