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

import Tkinter, Tkconstants, tkMessageBox
from Tkinter import Radiobutton, IntVar

from AlleleSubCommon import getConfigurationValue, assignConfigurationValue

class EmblSubOptionsForm(Tkinter.Frame):
        
    # Initialize the GUI
    def __init__(self, root):
        Tkinter.Frame.__init__(self, root)
        root.title("Choose EMBL Submission Options")
        self.parent = root

        #button_opt = {'fill': Tkconstants.BOTH, 'padx': 35, 'pady': 5}
        
        # To define the exit behavior.  Save and exit.
        self.parent.protocol('WM_DELETE_WINDOW', self.saveOptions)
        
        # Define the return behavior.  Same as "close window" etc
        root.bind('<Return>', self.returnFunction)
  
        # This window should not be resizeable. I guess.
        self.parent.resizable(width=False, height=False)
        
                #Standard Inputs widths for the form elements
        formInputWidth = 30
        labelInputWidth = 30
        
        self.instructionsFrame = Tkinter.Frame(self)  
        self.instructionText = Tkinter.StringVar()       
        self.instructionText.set('\nThese options are required for an EMBL allele submission.\n')        
        Tkinter.Label(self.instructionsFrame, width=85, height=3, textvariable=self.instructionText).pack()
        self.instructionsFrame.pack()
        
        self.submissionDetailsInputFrame2 = Tkinter.Frame(self)
  
        self.sampleIDInstrText = Tkinter.StringVar()
        self.sampleIDInstrText.set('Sample ID:')
        self.sampleIDinstrLabel = Tkinter.Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.sampleIDInstrText).grid(row=0, column=0)
        self.inputSampleID = Tkinter.StringVar()
        self.inputSampleIDEntry = Tkinter.Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputSampleID).grid(row=0, column=1)

        self.geneInstrStringVar = Tkinter.StringVar()
        self.geneInstrStringVar.set('Gene:')
        self.geneInstrLabel = Tkinter.Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.geneInstrStringVar).grid(row=1, column=0)
        self.inputGene = Tkinter.StringVar()        
        self.inputGeneEntry = Tkinter.Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputGene).grid(row=1, column=1)

        self.chooseClassIntVar = IntVar()
        self.chooseClassIntVar.set(1)
        Radiobutton(self.submissionDetailsInputFrame2, text="HLA Class I ", variable=self.chooseClassIntVar, value=1).grid(row=2, column=0)
        Radiobutton(self.submissionDetailsInputFrame2, text="HLA Class II", variable=self.chooseClassIntVar, value=2).grid(row=2, column=1)

        self.alleleInstrText = Tkinter.StringVar()
        self.alleleInstrText.set('Allele Local Name:')
        self.alleleInstrLabel = Tkinter.Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.alleleInstrText).grid(row=3, column=0)
        self.inputAllele = Tkinter.StringVar() 
        self.inputAlleleEntry = Tkinter.Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputAllele).grid(row=3, column=1)

        self.submissionDetailsInputFrame2.pack()
                
                
        # Make a frame to contain the Test/Production radio buttons.
        self.testProductionFrame = Tkinter.Frame(self)
        
        self.testProductionInstrText = Tkinter.StringVar()
        self.testProductionInstrText.set('\nBy default, you submit to the EMBL test servers,\n'
            + 'where submissions are regularly deleted.\n'
            + 'change this option if you want to submit to the live EMBL environment.\n'
            + 'Login Credentials will not be stored, but they will be sent\n'
            + 'to EMBL via secure https connection.\n'
            )
        self.alleleInstrLabel = Tkinter.Label(self.testProductionFrame, width=70, height=7, textvariable=self.testProductionInstrText).pack()#.grid(row=2, column=0)
 
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
  
        self.submissionDetailsInputFrame.pack()
        
        
        # Frame to specify Analysis Information
        self.newAnalysisFrame = Tkinter.Frame(self)            

        self.analysisAliasInstrText = Tkinter.StringVar()
        self.analysisAliasInstrText.set('Analysis Alias:')
        self.analysisAliasInstrLabel = Tkinter.Label(self.newAnalysisFrame, width=labelInputWidth, height=1, textvariable=self.analysisAliasInstrText).grid(row=0, column=0)
        self.inputAnalysisAlias = Tkinter.StringVar()
        self.inputStudyIdEntry = Tkinter.Entry(self.newAnalysisFrame, width=formInputWidth, textvariable=self.inputAnalysisAlias).grid(row=0, column=1)

        self.analysisTitleInstrText = Tkinter.StringVar()
        self.analysisTitleInstrText.set('Analysis Title:')
        self.analysisTitleInstrLabel = Tkinter.Label(self.newAnalysisFrame, width=labelInputWidth, height=1, textvariable=self.analysisTitleInstrText).grid(row=1, column=0)
        self.inputAnalysisTitle = Tkinter.StringVar()
        self.inputAnalysisTitleEntry = Tkinter.Entry(self.newAnalysisFrame, width=formInputWidth, textvariable=self.inputAnalysisTitle).grid(row=1, column=1)

        self.analysisDescriptionInstrText = Tkinter.StringVar()
        self.analysisDescriptionInstrText.set('Analysis Description:')
        self.analysisDescriptionInstrLabel = Tkinter.Label(self.newAnalysisFrame, width=labelInputWidth, height=1, textvariable=self.analysisDescriptionInstrText).grid(row=2, column=0)
        self.inputAnalysisDescription = Tkinter.StringVar()
        self.inputAnalysisDescriptionEntry = Tkinter.Entry(self.newAnalysisFrame, width=formInputWidth, textvariable=self.inputAnalysisDescription).grid(row=2, column=1)
        
        self.newAnalysisFrame.pack()
        
        
        
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
        self.inputStudyIdEntry = Tkinter.Entry(self.existingProjectFrame, width=formInputWidth, textvariable=self.inputStudyAccession).grid(row=0, column=1)
        self.existingProjectFrame.pack()
        
        
        # Filler Label
        Tkinter.Label(self.projectDetailsFrame, width=labelInputWidth, height=1, text=' ').pack()
        
        # This radio button is on the project details frame, but not 
        # on one of it's sub-frames (existingProjectFrame or newProjectFrame)  
        # That's so i can pack it, and not use a grid
        Radiobutton(self.projectDetailsFrame, text="Create a new study with this information:", variable=self.chooseProjectIntVar, value=2).pack()

        self.newProjectFrame = Tkinter.Frame(self.projectDetailsFrame)            

        self.studyIdInstrText = Tkinter.StringVar()
        self.studyIdInstrText.set('Short Study Identifier:')
        self.studyIdInstrLabel = Tkinter.Label(self.newProjectFrame, width=labelInputWidth, height=1, textvariable=self.studyIdInstrText).grid(row=0, column=0)
        self.inputStudyId = Tkinter.StringVar()
        self.inputStudyIdEntry = Tkinter.Entry(self.newProjectFrame, width=formInputWidth, textvariable=self.inputStudyId).grid(row=0, column=1)

        self.studyShortTitleInstrText = Tkinter.StringVar()
        self.studyShortTitleInstrText.set('Descriptive Study Title:')
        self.studyShortTitleInstrLabel = Tkinter.Label(self.newProjectFrame, width=labelInputWidth, height=1, textvariable=self.studyShortTitleInstrText).grid(row=1, column=0)
        self.inputStudyShortTitle = Tkinter.StringVar()
        self.inputStudyShortTitleEntry = Tkinter.Entry(self.newProjectFrame, width=formInputWidth, textvariable=self.inputStudyShortTitle).grid(row=1, column=1)

        self.studyAbstractInstrText = Tkinter.StringVar()
        self.studyAbstractInstrText.set('Study Description / Abstract:')
        self.studyAbstractInstrLabel = Tkinter.Label(self.newProjectFrame, width=labelInputWidth, height=1, textvariable=self.studyAbstractInstrText).grid(row=2, column=0)
        self.inputStudyAbstract = Tkinter.StringVar()
        self.inputStudyAbstractEntry = Tkinter.Entry(self.newProjectFrame, width=formInputWidth, textvariable=self.inputStudyAbstract).grid(row=2, column=1)
        
        self.newProjectFrame.pack()
        
        self.projectDetailsFrame.pack()

        # Make a frame for the save options button.
        self.saveOptionsFrame = Tkinter.Frame(self)
        Tkinter.Button(self.saveOptionsFrame, text='Save Options', command=self.saveOptions).grid(row=0, column=0)
        self.saveOptionsFrame.pack()
        
        # TODO: Should there be a cancel button, to close this window without saving?
        
        self.loadOptions()

    # I needed a function for the return keypress to latch onto.
    # It is just a wrapper for the saveOptions method.
    def returnFunction(self, event):
        self.saveOptions()
       

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
            
        if getConfigurationValue('study_identifier') is not None:
            self.inputStudyId.set(getConfigurationValue('study_identifier'))
            
        if getConfigurationValue('study_short_title') is not None:
            self.inputStudyShortTitle.set(getConfigurationValue('study_short_title'))
            
        if getConfigurationValue('study_abstract') is not None:
            self.inputStudyAbstract.set(getConfigurationValue('study_abstract'))
            
        if getConfigurationValue('test_submission') is not None:
            # 1 = Test.  0 = Production/live server
            self.chooseTestServersIntVar.set(int(getConfigurationValue('test_submission')))
            
        if getConfigurationValue('analysis_alias') is not None:
            self.inputAnalysisAlias.set(getConfigurationValue('analysis_alias'))
        if getConfigurationValue('analysis_title') is not None:
            self.inputAnalysisTitle.set(getConfigurationValue('analysis_title'))
        if getConfigurationValue('analysis_description') is not None:
            self.inputAnalysisDescription.set(getConfigurationValue('analysis_description'))


    def saveOptions(self):
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
            assignConfigurationValue('study_identifier', self.inputStudyId.get())            
            assignConfigurationValue('study_short_title', self.inputStudyShortTitle.get())
            assignConfigurationValue('study_abstract', self.inputStudyAbstract.get())
            assignConfigurationValue('test_submission', str(self.chooseTestServersIntVar.get()))
            assignConfigurationValue('analysis_alias', str(self.inputAnalysisAlias.get()))
            assignConfigurationValue('analysis_title', str(self.inputAnalysisTitle.get()))
            assignConfigurationValue('analysis_description', str(self.inputAnalysisDescription.get()))
            
            self.parent.destroy() 
            
        else:
            #print('Not ready to save, you are missing options.')
            pass
        
    def checkOptions(self):
        
        # TODO: THis is really annoying. I should allow them to close this screen without providing all the information
        # Im sure a user will complain about this.
        
        #print ('Checking options.')
       
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
            if (not self.inputStudyId.get()):
                tkMessageBox.showwarning('Missing Form Value',
                    'You are missing a Study Name. Please try again.')
                return False
            
            if (not self.inputStudyShortTitle.get()):
                tkMessageBox.showwarning('Missing Form Value',
                    'You are missing a Study Description. Please try again.')
                return False
            
            
            if (not self.inputStudyAbstract.get()):
                tkMessageBox.showwarning('Missing Form Value',
                    'You are missing a Study Accession number. Please try again.')
                return False
            
        else:
            raise Exception ('Unknown value of self.chooseProjectIntVar. I expect 1 or 2. Observed:' + str(self.chooseProjectIntVar))
        
        
        if (not self.inputAnalysisAlias.get()):
            tkMessageBox.showwarning('Missing Form Value',
                'You are missing an Analysis Alias. Please try again.')
            return False        
        
        if (not self.inputAnalysisTitle.get()):
            tkMessageBox.showwarning('Missing Form Value',
                'You are missing an Analysis Title. Please try again.')
            return False        
        
        if (not self.inputAnalysisDescription.get()):
            tkMessageBox.showwarning('Missing Form Value',
                'You are missing an Analysis Description. Please try again.')
            return False
        
    

        # All options look good, right?
        
        
        return True
    
    
    def closeWindow(self):
        self.parent.destroy()        
    