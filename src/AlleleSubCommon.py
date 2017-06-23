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

import xml.etree.cElementTree as ET


import sys
from os.path import dirname, join, abspath, isfile

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
    
    # TODO: Some variables should be hard-coded here, I think.
    # If they are undefined:
    # Define the TEST website.
    # Define the PRODUCTION website
    # Make an option pointing to the production website.
    # Change the MAIN INTERFACE so that we know if we're pointing at the TEST server.

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
    
    tree = ET.ElementTree(root)

    tree.write(globalVariables['config_file_location'])
    
    

def loadConfigurationFile():
    assignConfigName()
    
    if not isfile(globalVariables['config_file_location']):
        print ('The config file does not exist yet. I will not load it:\n' + globalVariables['config_file_location'])
        
        # Here is where I assign the common configuration values
        # TODO: Assign the commonly used config values
        
    else:
        print ('The config file already exists, I will load it:\n' + globalVariables['config_file_location'])
        
        tree = ET.parse(globalVariables['config_file_location'])
        root = tree.getroot()
        
        for child in root:
            #print 'childtag:' + str(child.tag)
            #print 'childattrib:' + str(child.attrib)
            #print 'childtext:' + str(child.text)
            assignConfigurationValue(child.tag, child.text)
            
        
    
    
    
    