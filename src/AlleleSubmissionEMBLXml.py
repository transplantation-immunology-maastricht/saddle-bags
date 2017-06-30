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


#import os

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

def createProjectXML(fullXmlFilePath, projectName, projectDescription, projectAbstract):
    root = ET.Element('PROJECT_SET')
    # TODO: How do I get the center name?  
    # Maybe I should get this via REST?
    
    # Ok this is really confusing.
    # According to the docs, http://ena-docs.readthedocs.io/en/latest/prog_01.html
    # "alias" attribute on the project node contains projectName
    # "title" node contains project description
    # "description" node is the project abstract 
    # EMBL should be more consistent in their terminology.
    projectElement = ET.SubElement(root, 'PROJECT')
    projectElement.set('alias', projectName)
    projectElement.set('center_name', 'Maastricht University Medical Center' )
    titleElement = ET.SubElement(projectElement, 'TITLE')
    titleElement.text = projectDescription
    descriptionElement = ET.SubElement(projectElement, 'DESCRIPTION')
    descriptionElement.text = projectAbstract
    submissionProjectElement = ET.SubElement(projectElement, 'SUBMISSION_PROJECT')
    sequencingProjectElement = ET.SubElement(submissionProjectElement, 'SEQUENCING_PROJECT')

    return writeToXml(fullXmlFilePath, root)
    
def createProjectSubmissionXML(submissionAlias, fullXmlFilePath):    
    root = ET.Element('SUBMISSION')
    # TODO: How do I get the center name?  
    # Maybe I should get this via REST?
    root.set('alias', submissionAlias)
    root.set('center_name', 'Maastricht University Medical Center' )
    actionsElement = ET.SubElement(root, 'ACTIONS')    
    actionElement = ET.SubElement(actionsElement, 'ACTION')
    addElement = ET.SubElement(actionElement, 'ADD')  
    addElement.set('source','project.xml')
    addElement.set('schema','project')

    return writeToXml(fullXmlFilePath, root)
    

