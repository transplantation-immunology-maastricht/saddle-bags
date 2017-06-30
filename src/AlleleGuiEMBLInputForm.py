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

import os

import Tkinter, Tkconstants, tkFileDialog, tkMessageBox
from Tkinter import *

from AlleleSubCommon import *

class AlleleGuiEMBLInputForm(Tkinter.Frame):
        
    # Initialize the GUI
    def __init__(self, root):
        Tkinter.Frame.__init__(self, root)
        root.title("Choose EMBL Submission Options")
        self.parent = root

        button_opt = {'fill': Tkconstants.BOTH, 'padx': 35, 'pady': 5}
        
        # To define the exit behavior.  Save and exit.
        self.parent.protocol('WM_DELETE_WINDOW', self.saveOptions)
  
        # This window should not be resizeable. I guess.
        self.parent.resizable(width=False, height=False)
        
        self.instructionsFrame = Tkinter.Frame(self)  
        self.instructionText = Tkinter.StringVar()       
        self.instructionText.set('\nThese options are required for an EMBL allele submission.\n'
            + 'Login Credentials will not be stored, but they will be sent to EMBL via\n'
            + 'secure https connection.\n')        
        Tkinter.Label(self.instructionsFrame, width=85, height=6, textvariable=self.instructionText).pack()
        self.instructionsFrame.pack()
        
        #Standard Inputs widths for the form elements
        formInputWidth = 30
        labelInputWidth = 30
                
                
        # Make a frame to contain the Test/Production radio buttons.
        self.testProductionFrame = Tkinter.Frame(self)
        
        self.testProductionInstrText = Tkinter.StringVar()
        self.testProductionInstrText.set('\nBy default, you submit to the EMBL test servers,\n'
            + 'where submissions are regularly deleted.\n'
            + 'change this option if you want to submit to the live EMBL environment.\n'
            )
        self.alleleInstrLabel = Tkinter.Label(self.testProductionFrame, width=70, height=5, textvariable=self.testProductionInstrText).pack()#.grid(row=2, column=0)
 
        # 1 = Test.  0 = Production/live server
        self.chooseTestServersIntVar = IntVar()
        self.chooseTestServersIntVar.set(int(getConfigurationValue('test_submission')))
 
        Radiobutton(self.testProductionFrame, text="Submit to EMBL TEST / DEMO environment.", variable=self.chooseTestServersIntVar, value=1).pack()
        Radiobutton(self.testProductionFrame, text="Submit to EMBL LIVE / PROD environment.", variable=self.chooseTestServersIntVar, value=0).pack()
        
        self.testProductionFrame.pack()
     
        # Make a frame to contain the input variables
        self.submissionDetailsInputFrame = Tkinter.Frame(self)

        self.usernameInstrText = Tkinter.StringVar()
        self.usernameInstrText.set('EMBL Username:')
        self.usernameInstrLabel = Tkinter.Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.usernameInstrText).grid(row=0, column=0)
        self.inputUsername = Tkinter.StringVar()
        self.inputUsernameEntry = Tkinter.Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputUsername).grid(row=0, column=1)

        self.passwordInstrText = Tkinter.StringVar()
        self.passwordInstrText.set('EMBL Password:')
        self.passwordInstrLabel = Tkinter.Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.passwordInstrText).grid(row=1, column=0)
        self.inputPassword = Tkinter.StringVar()
        self.inputPasswordEntry = Tkinter.Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputPassword, show="*").grid(row=1, column=1)
  
        self.sampleIDInstrText = Tkinter.StringVar()
        self.sampleIDInstrText.set('Sample ID:')
        self.sampleIDinstrLabel = Tkinter.Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.sampleIDInstrText).grid(row=2, column=0)
        self.inputSampleID = Tkinter.StringVar()
        self.inputSampleIDEntry = Tkinter.Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputSampleID).grid(row=2, column=1)

        self.geneInstrStringVar = Tkinter.StringVar()
        self.geneInstrStringVar.set('Gene:')
        self.geneInstrLabel = Tkinter.Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.geneInstrStringVar).grid(row=3, column=0)
        self.inputGene = Tkinter.StringVar()        
        self.inputGeneEntry = Tkinter.Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputGene).grid(row=3, column=1)

        self.chooseClassIntVar = IntVar()
        self.chooseClassIntVar.set(1)
        Radiobutton(self.submissionDetailsInputFrame, text="HLA Class I ", variable=self.chooseClassIntVar, value=1).grid(row=4, column=0)
        Radiobutton(self.submissionDetailsInputFrame, text="HLA Class II", variable=self.chooseClassIntVar, value=2).grid(row=4, column=1)

        self.alleleInstrText = Tkinter.StringVar()
        self.alleleInstrText.set('Allele Local Name:')
        self.alleleInstrLabel = Tkinter.Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.alleleInstrText).grid(row=5, column=0)
        self.inputAllele = Tkinter.StringVar() 
        self.inputAlleleEntry = Tkinter.Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputAllele).grid(row=5, column=1)

        self.submissionDetailsInputFrame.pack()
        
        # A Frame for specifing the details of the Study / Project
        self.projectDetailsFrame = Tkinter.Frame(self)
        
        self.alleleInstrText = Tkinter.StringVar()
        self.alleleInstrText.set('\nEMBL requires that submissions are assigned to a Study/Project.\n'
            + 'Will you provide an existing EMBL study accession #?\n'
            + '(ex. \'PRJEB01234\')\n'
            + 'Or will you specify a new study?\n')
        self.alleleInstrLabel = Tkinter.Label(self.projectDetailsFrame, width=70, height=6, textvariable=self.alleleInstrText).pack()#.grid(row=2, column=0)
 
        self.chooseProjectIntVar = IntVar()
        self.chooseProjectIntVar.set(2)
 
        # A frame for the "new study" radio button
        self.existingProjectFrame = Tkinter.Frame(self.projectDetailsFrame)
        Radiobutton(self.existingProjectFrame, text="Use this study accession:", variable=self.chooseProjectIntVar, value=1).grid(row=0,column=0)
        self.inputStudyAccession = Tkinter.StringVar()
        self.inputStudyNameEntry = Tkinter.Entry(self.existingProjectFrame, width=formInputWidth, textvariable=self.inputStudyAccession).grid(row=0, column=1)
        self.existingProjectFrame.pack()
        
        
        # Filler Label
        Tkinter.Label(self.projectDetailsFrame, width=labelInputWidth, height=1, text=' ').pack()
        
        # This radio button is on the project details frame, but not 
        # on one of it's sub-frames (existingProjectFrame or newProjectFrame)  
        # That's so i can pack it, and not use a grid
        Radiobutton(self.projectDetailsFrame, text="Create a new study with this information:", variable=self.chooseProjectIntVar, value=2).pack()

        self.newProjectFrame = Tkinter.Frame(self.projectDetailsFrame)            

        self.studyNameInstrText = Tkinter.StringVar()
        self.studyNameInstrText.set('Study Name:')
        self.studyNameInstrLabel = Tkinter.Label(self.newProjectFrame, width=labelInputWidth, height=1, textvariable=self.studyNameInstrText).grid(row=0, column=0)
        self.inputStudyName = Tkinter.StringVar()
        self.inputStudyNameEntry = Tkinter.Entry(self.newProjectFrame, width=formInputWidth, textvariable=self.inputStudyName).grid(row=0, column=1)

        self.studyShortDescriptionInstrText = Tkinter.StringVar()
        self.studyShortDescriptionInstrText.set('Short Description:')
        self.studyShortDescriptionInstrLabel = Tkinter.Label(self.newProjectFrame, width=labelInputWidth, height=1, textvariable=self.studyShortDescriptionInstrText).grid(row=1, column=0)
        self.inputStudyShortDescription = Tkinter.StringVar()
        self.inputStudyShortDescriptionEntry = Tkinter.Entry(self.newProjectFrame, width=formInputWidth, textvariable=self.inputStudyShortDescription).grid(row=1, column=1)

        self.studyAbstractInstrText = Tkinter.StringVar()
        self.studyAbstractInstrText.set('Study Abstract:')
        self.studyAbstractInstrLabel = Tkinter.Label(self.newProjectFrame, width=labelInputWidth, height=1, textvariable=self.studyAbstractInstrText).grid(row=2, column=0)
        self.inputStudyAbstract = Tkinter.StringVar()
        self.inputStudyAbstractEntry = Tkinter.Entry(self.newProjectFrame, width=formInputWidth, textvariable=self.inputStudyAbstract).grid(row=2, column=1)
        
        self.newProjectFrame.pack()
        
        self.projectDetailsFrame.pack()

        # Make a frame for the save options button.
        self.saveOptionsFrame = Tkinter.Frame(self)
        Tkinter.Button(self.saveOptionsFrame, text='Save Options', command=self.saveOptions).grid(row=0, column=0)
        self.saveOptionsFrame.pack()
        
        self.loadOptions()

    # submissionOptions is a dictionary, passed by the parent.
    def loadOptions(self):
        if getConfigurationValue('embl_username') is not None:
            self.inputUsername.set(getConfigurationValue('embl_username'))
            
        if getConfigurationValue('embl_password') is not None:
            self.inputPassword.set(getConfigurationValue('embl_password'))
            
        if getConfigurationValue('sample_id') is not None:
            self.inputSampleID.set(getConfigurationValue('sample_id'))
            
        if getConfigurationValue('gene') is not None:
            self.inputGene.set(getConfigurationValue('gene'))
            
        if getConfigurationValue('class') is not None:
            if (str(getConfigurationValue('class')) == '1'):
                self.chooseClassIntVar.set(1)
            elif (str(getConfigurationValue('class')) == '2'):
                self.chooseClassIntVar.set(2)
            else:
                raise Exception('Error loading EMBL submission options. Invalid class:' + str(getConfigurationValue('class')))
    
        if getConfigurationValue('allele_name') is not None:
            self.inputAllele.set(getConfigurationValue('allele_name'))
            
        if getConfigurationValue('choose_project') is not None:
            if (str(getConfigurationValue('choose_project')) == '1'):
                self.chooseProjectIntVar.set(1)
            elif (str(getConfigurationValue('choose_project')) == '2'):
                self.chooseProjectIntVar.set(2)
            else:
                raise Exception('Error loading EMBL submission options. Invalid Project choice:' + str(getConfigurationValue('choose_project')))
            
        if getConfigurationValue('study_accession') is not None:
            self.inputStudyAccession.set(getConfigurationValue('study_accession'))
            
        if getConfigurationValue('study_name') is not None:
            self.inputStudyName.set(getConfigurationValue('study_name'))
            
        if getConfigurationValue('study_description') is not None:
            self.inputStudyShortDescription.set(getConfigurationValue('study_description'))
            
        if getConfigurationValue('study_abstract') is not None:
            self.inputStudyAbstract.set(getConfigurationValue('study_abstract'))
            
        if getConfigurationValue('test_submission') is not None:
            # 1 = Test.  0 = Production/live server
            self.chooseTestServersIntVar.set(int(getConfigurationValue('test_submission')))
 
        
            

    def saveOptions(self):
        # TODO: Save the options to our configuration dictionary
        # Close the window
        if (self.checkOptions()):
            print ('Saving Options....')
            
            assignConfigurationValue('embl_username', self.inputUsername.get())
            # I store this password so I can use it in the submission
            # I don't ever want to save the password. Make sure it isn't being saved in the config, in AlleleSubCommon.py
            assignConfigurationValue('embl_password', self.inputPassword.get())
            assignConfigurationValue('sample_id', self.inputSampleID.get())
            assignConfigurationValue('gene', self.inputGene.get())
            assignConfigurationValue('class', str(self.chooseClassIntVar.get()))             
            assignConfigurationValue('allele_name', self.inputAllele.get())
            assignConfigurationValue('choose_project', str(self.chooseProjectIntVar.get()))    
            assignConfigurationValue('study_accession', self.inputStudyAccession.get())                    
            assignConfigurationValue('study_name', self.inputStudyName.get())            
            assignConfigurationValue('study_description', self.inputStudyShortDescription.get())
            assignConfigurationValue('study_abstract', self.inputStudyAbstract.get())
            assignConfigurationValue('test_submission', str(self.chooseTestServersIntVar.get()))
            
            self.parent.destroy() 
            
        else:
            #print('Not ready to save, you are missing options.')
            pass
        
    def checkOptions(self):
        # TODO this method
        print ('Checking options.')

        #chooseProjectIntVar
        
        # Don't check the EMBL Username
        # Don't check the EMBL Password
        
        if (not self.inputSampleID.get()):
            tkMessageBox.showwarning('Missing Form Value',
                'You are missing a Sample ID. Please try again.')
            return False
        
        if (not self.inputGene.get()):
            tkMessageBox.showwarning('Missing Form Value',
                'You are missing a Gene. Please try again.')
            return False
        
        # Don't check the class boolean
        
        if (not self.inputAllele.get()):
            tkMessageBox.showwarning('Missing Form Value',
                'You are missing an Allele Name. Please try again.')
            return False
        
        if (str(self.chooseProjectIntVar.get()) == '1'):
            # Use Existing Project
            if (not self.inputStudyAccession.get()):
                tkMessageBox.showwarning('Missing Form Value',
                    'You are missing a Study Accession number. Please try again.')
                return False
            
        elif(str(self.chooseProjectIntVar.get()) == '2'):
            # Use New Project
            if (not self.inputStudyName.get()):
                tkMessageBox.showwarning('Missing Form Value',
                    'You are missing a Study Name. Please try again.')
                return False
            
            if (not self.inputStudyShortDescription.get()):
                tkMessageBox.showwarning('Missing Form Value',
                    'You are missing a Study Description. Please try again.')
                return False
            
            
            if (not self.inputStudyAbstract.get()):
                tkMessageBox.showwarning('Missing Form Value',
                    'You are missing a Study Accession number. Please try again.')
                return False
            
        else:
            raise Exception ('Unknown value of self.chooseProjectIntVar. I expect 1 or 2. Observed:' + str(self.chooseProjectIntVar)) 

        # All options look good, right?
        return True
    
    
    def closeWindow(self):
        #writeConfigurationFile()

        self.parent.destroy()        
    