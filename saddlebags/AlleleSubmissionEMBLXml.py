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

import xml.etree.ElementTree as ET
import xml.dom.minidom

# Here we have methods to create XML files necessary for EMBL submission.
# Schemas are defined on github.
# https://github.com/enasequence/schema

def writeToXml(fullXmlFilePath, xmlElementTree):
    xmlText = ET.tostring(xmlElementTree, encoding='utf8', method='xml')
    prettyXmlText = xml.dom.minidom.parseString(xmlText).toprettyxml()
    
    xmlOutput = createOutputFile(fullXmlFilePath)
    xmlOutput.write(prettyXmlText)
    xmlOutput.close
    
    return prettyXmlText

def getCenterName():
    # TODO: Should I use REST here? 
    # Probably not, center_name is not required in the xmls.
    return 'Center_Name'

def createProjectXML(fullXmlFilePath):
    # They are called "Project" in xml, but "Study" on the website.
    # Project = Study
    root = ET.Element('PROJECT_SET')

    projectID = getConfigurationValue('study_identifier')
    projectShortTitle = getConfigurationValue('study_short_title')
    projectAbstract = getConfigurationValue('study_abstract')
    
    projectElement = ET.SubElement(root, 'PROJECT')
    projectElement.set('alias', projectID)
    # Center Name is optional according to schemas.  Forget it. EMBL Knows our login info.
    #projectElement.set('center_name', getCenterName() )
    titleElement = ET.SubElement(projectElement, 'TITLE')
    titleElement.text = projectShortTitle
    descriptionElement = ET.SubElement(projectElement, 'DESCRIPTION')
    descriptionElement.text = projectAbstract
    submissionProjectElement = ET.SubElement(projectElement, 'SUBMISSION_PROJECT')
    sequencingProjectElement = ET.SubElement(submissionProjectElement, 'SEQUENCING_PROJECT')

    return writeToXml(fullXmlFilePath, root)
    
def createProjectSubmissionXML(fullXmlFilePath, submissionAlias, shortProjectFileName):    
    root = ET.Element('SUBMISSION')
    root.set('alias', submissionAlias)
    # Center Name is optional according to schemas.  Forget it.
    #root.set('center_name', getCenterName() )
    actionsElement = ET.SubElement(root, 'ACTIONS')    
    actionElement = ET.SubElement(actionsElement, 'ACTION')
    addElement = ET.SubElement(actionElement, 'ADD')  
    addElement.set('source',shortProjectFileName)
    addElement.set('schema','project')

    return writeToXml(fullXmlFilePath, root)

def createAnalysisXML(fullXmlFilePath, checksumValue, flatfileZipFileName):
    # An analysis xml is just a wrapper for a sequence submission. 
    root = ET.Element('ANALYSIS_SET')

    # TODO: I haven't created these three analysis configuration values yet.
    # Probably need to add this to the GUI, or somehow generate them automagically.   
    analysisElement = ET.SubElement(root, 'ANALYSIS')
    analysisElement.set('alias', getConfigurationValue('analysis_alias'))
    
    titleElement = ET.SubElement(analysisElement, 'TITLE')
    titleElement.text = (getConfigurationValue('analysis_title'))
    
    descriptionElement = ET.SubElement(analysisElement, 'DESCRIPTION')
    descriptionElement.text = (getConfigurationValue('analysis_description'))
    
    studyRefElement = ET.SubElement(analysisElement, 'STUDY_REF')
    studyRefElement.set('accession', getConfigurationValue('study_accession'))
    
    analysisTypeElement = ET.SubElement(analysisElement, 'ANALYSIS_TYPE')
    sequenceFlatfileElement = ET.SubElement(analysisTypeElement, 'SEQUENCE_FLATFILE')
    
    filesElement = ET.SubElement(analysisElement, 'FILES')
    
    fileElement = ET.SubElement(filesElement, 'FILE')
    fileElement.set('checksum', checksumValue)
    fileElement.set('checksum_method', 'MD5')
    fileElement.set('filename', flatfileZipFileName)
    fileElement.set('filetype', 'flatfile')

    return writeToXml(fullXmlFilePath, root)
    
def createAnalysisSubmissionXML(fullXmlFilePath, submissionAlias, shortAnalysisFileName):    
    root = ET.Element('SUBMISSION')
    
    root.set('alias', submissionAlias)
    actionsElement = ET.SubElement(root, 'ACTIONS')    
    actionElement = ET.SubElement(actionsElement, 'ACTION')
    addElement = ET.SubElement(actionElement, 'ADD')  
    addElement.set('source',shortAnalysisFileName)
    addElement.set('schema','analysis')

    return writeToXml(fullXmlFilePath, root)
    

