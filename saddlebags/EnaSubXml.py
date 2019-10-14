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



from saddlebags.AlleleSubCommon import createOutputFile, getConfigurationValue

# XML packages use very generic method names (like parse and tostring) so I will alias this reference.
from xml.etree import ElementTree as ET
from xml.dom import minidom as MD

# Here we have methods to create XML files necessary for EMBL submission.
# Schemas are defined on github.
# https://github.com/enasequence/schema
# I must still register a Study/Project using XML and REST.

def writeToXml(fullXmlFilePath, xmlElementTree):
    xmlText = ET.tostring(xmlElementTree, encoding='utf8', method='xml')
    prettyXmlText = MD.parseString(xmlText).toprettyxml()
    
    xmlOutput = createOutputFile(fullXmlFilePath)
    xmlOutput.write(prettyXmlText)
    xmlOutput.close
    
    return prettyXmlText

def createProjectXML(fullXmlFilePath, projectID, projectShortTitle, projectAbstract):
    # They are called "Project" in xml, but "Study" on the website.
    # Project = Study
    root = ET.Element('PROJECT_SET')
    projectElement = ET.SubElement(root, 'PROJECT')
    projectElement.set('alias', projectID)
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
    actionsElement = ET.SubElement(root, 'ACTIONS')
    actionElement = ET.SubElement(actionsElement, 'ACTION')
    addElement = ET.SubElement(actionElement, 'ADD')  
    addElement.set('source',shortProjectFileName)
    addElement.set('schema','project')

    return writeToXml(fullXmlFilePath, root)

