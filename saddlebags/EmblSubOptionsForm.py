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

#import Tkinter, Tkconstants, tkMessageBox
#from Tkinter import Radiobutton, IntVar

from tkinter import Frame, StringVar, IntVar, Label, Entry, Radiobutton, Button, messagebox

from saddlebags.AlleleSubCommon import getConfigurationValue, assignConfigurationValue, assignIcon, logEvent
from saddlebags.ScrolledWindow import VerticalScrolledFrame

class EmblSubOptionsForm(VerticalScrolledFrame):
        
    # Initialize the GUI
    def __init__(self, root):
        
        #Frame.__init__(self, root)
        VerticalScrolledFrame.__init__(self, root)
        
        root.title("Choose EMBL-ENA Submission Options")
        self.parent = root
        
        # Assign the icon of this sub-window.
        assignIcon(self.parent)

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
        
        self.instructionsFrame = Frame(self.interior)  
        self.instructionText = StringVar()       
        self.instructionText.set('\nProvide this required sequence metadata, please.\n' 
                + 'The Allele Local Name is provided by the submitting laboratory.\n'
                + 'You may name it similar to the closest known allele,\n'
                + 'if you wish.')        
        Label(self.instructionsFrame, width=85, height=6, textvariable=self.instructionText).pack()
        self.instructionsFrame.pack()
        
        self.submissionDetailsInputFrame2 = Frame(self.interior)
  
        self.sampleIDInstrText = StringVar()
        self.sampleIDInstrText.set('Sample ID:')
        self.sampleIDinstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.sampleIDInstrText).grid(row=0, column=0)
        self.inputSampleID = StringVar()
        self.inputSampleIDEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputSampleID).grid(row=0, column=1)

        self.geneInstrStringVar = StringVar()
        self.geneInstrStringVar.set('Gene:')
        self.geneInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.geneInstrStringVar).grid(row=1, column=0)
        self.inputGene = StringVar()        
        self.inputGeneEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputGene).grid(row=1, column=1)

        self.chooseClassIntVar = IntVar()
        self.chooseClassIntVar.set(1)
        Radiobutton(self.submissionDetailsInputFrame2, text="HLA Class I ", variable=self.chooseClassIntVar, value=1).grid(row=2, column=0)
        Radiobutton(self.submissionDetailsInputFrame2, text="HLA Class II", variable=self.chooseClassIntVar, value=2).grid(row=2, column=1)

        self.alleleInstrText = StringVar()
        self.alleleInstrText.set('Allele Local Name:')
        self.alleleInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.alleleInstrText).grid(row=3, column=0)
        self.inputAllele = StringVar() 
        self.inputAlleleEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputAllele).grid(row=3, column=1)

        self.submissionDetailsInputFrame2.pack()
                
                
        # Make a frame to contain the Test/Production radio buttons.
        self.testProductionFrame = Frame(self.interior)
        
        self.testProductionInstrText = StringVar()
        self.testProductionInstrText.set('\nBy default, you submit to the EMBL-ENA test servers,\n'
            + 'where submissions are regularly deleted.\n'
            + 'change this option if you want to submit to the live environment.\n'
            + 'Login Credentials will not be stored, but they will be sent\n'
            + 'to EMBL-ENA via REST during allele submission.\n'
            )
        self.alleleInstrLabel = Label(self.testProductionFrame, width=70, height=7, textvariable=self.testProductionInstrText).pack()#.grid(row=2, column=0)
 
        # 1 = Test.  0 = Production/live server
        self.chooseTestServersIntVar = IntVar()
        self.chooseTestServersIntVar.set(int(getConfigurationValue('test_submission')))
 
        Radiobutton(self.testProductionFrame, text="Submit to EMBL-ENA TEST / DEMO environment.", variable=self.chooseTestServersIntVar, value=1).pack()
        Radiobutton(self.testProductionFrame, text="Submit to EMBL-ENA LIVE / PROD environment.", variable=self.chooseTestServersIntVar, value=0).pack()
        
        self.testProductionFrame.pack()
     
        # Make a frame to contain the input variables
        self.submissionDetailsInputFrame = Frame(self.interior)

        self.usernameInstrText = StringVar()
        self.usernameInstrText.set('EMBL Username:')
        self.usernameInstrLabel = Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.usernameInstrText).grid(row=0, column=0)
        self.inputUsername = StringVar()
        self.inputUsernameEntry = Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputUsername).grid(row=0, column=1)

        self.passwordInstrText = StringVar()
        self.passwordInstrText.set('EMBL Password:')
        self.passwordInstrLabel = Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.passwordInstrText).grid(row=1, column=0)
        self.inputPassword = StringVar()
        self.inputPasswordEntry = Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputPassword, show="*").grid(row=1, column=1)
  
        self.submissionDetailsInputFrame.pack()
  
        # This frame just has instructions for how to use the Analysis options
        self.analysisOptionsInputFrame = Frame(self.interior)  
  
        # TODO: Better instructions for how a user will use this?
        self.analysisMetadataInstrText = StringVar()
        self.analysisMetadataInstrText.set('\nEMBL-ENA will store the sequence in an "Analysis" object.\n'
                + 'You must specify an analysis Alias, Title, and\n'
                + 'Description for every submitted sequence.\n'
                + 'The Alias must be unique.\n'
            )
        self.analysisMetadataLabel = Label(self.analysisOptionsInputFrame, width=70, height=5, textvariable=self.analysisMetadataInstrText).pack()#.grid(row=2, column=0)
        self.analysisOptionsInputFrame.pack()
       
        
        
        # Frame to specify Analysis Information
        self.newAnalysisFrame = Frame(self.interior)    
        
        self.analysisAliasInstrText = StringVar()
        self.analysisAliasInstrText.set('Analysis Alias:')
        self.analysisAliasInstrLabel = Label(self.newAnalysisFrame, width=labelInputWidth, height=1, textvariable=self.analysisAliasInstrText).grid(row=0, column=0)
        self.inputAnalysisAlias = StringVar()
        self.inputStudyIdEntry = Entry(self.newAnalysisFrame, width=formInputWidth, textvariable=self.inputAnalysisAlias).grid(row=0, column=1)

        self.analysisTitleInstrText = StringVar()
        self.analysisTitleInstrText.set('Analysis Title:')
        self.analysisTitleInstrLabel = Label(self.newAnalysisFrame, width=labelInputWidth, height=1, textvariable=self.analysisTitleInstrText).grid(row=1, column=0)
        self.inputAnalysisTitle = StringVar()
        self.inputAnalysisTitleEntry = Entry(self.newAnalysisFrame, width=formInputWidth, textvariable=self.inputAnalysisTitle).grid(row=1, column=1)

        self.analysisDescriptionInstrText = StringVar()
        self.analysisDescriptionInstrText.set('Analysis Description:')
        self.analysisDescriptionInstrLabel = Label(self.newAnalysisFrame, width=labelInputWidth, height=1, textvariable=self.analysisDescriptionInstrText).grid(row=2, column=0)
        self.inputAnalysisDescription = StringVar()
        self.inputAnalysisDescriptionEntry = Entry(self.newAnalysisFrame, width=formInputWidth, textvariable=self.inputAnalysisDescription).grid(row=2, column=1)
        
        self.newAnalysisFrame.pack()
       
        # A Frame for specifing the details of the Study / Project
        self.projectDetailsFrame = Frame(self.interior)
        
        self.alleleInstrText = StringVar()
        self.alleleInstrText.set('\nEMBL-ENA requires that submissions are assigned to a Study/Project.\n'
            + 'You must either:\n\n'
            + '1.) Provide an existing EMBL-ENA study/project accession #.\n'
            + 'This was provided in a previous submission, or you can see a list of\n'
            + 'your projects by logging into EMBL Webin. (ex. \'PRJEB01234\')\n'
            + '\n'
            + '2.) Specify metadata for a new study at EMBL-ENA.\n'
            + 'Provide an Identifier, Title, and short Abstract for the new Study,\n'
            + 'and I will create it automatically.\n'
            + 'Study Identifier must be Unique.'
            )
        self.alleleInstrLabel = Label(self.projectDetailsFrame, width=70, height=12, textvariable=self.alleleInstrText).pack()#.grid(row=2, column=0)
 
        self.chooseProjectIntVar = IntVar()
        self.chooseProjectIntVar.set(2)
 
        # A frame for the "new study" radio button
        self.existingProjectFrame = Frame(self.projectDetailsFrame)
        Radiobutton(self.existingProjectFrame, text="Use this study accession:", variable=self.chooseProjectIntVar, value=1).grid(row=0,column=0)
        self.inputStudyAccession = StringVar()
        self.inputStudyIdEntry = Entry(self.existingProjectFrame, width=formInputWidth, textvariable=self.inputStudyAccession).grid(row=0, column=1)
        self.existingProjectFrame.pack()
        
        
        # Filler Label
        Label(self.projectDetailsFrame, width=labelInputWidth, height=1, text=' ').pack()
        
        # This radio button is on the project details frame, but not 
        # on one of it's sub-frames (existingProjectFrame or newProjectFrame)  
        # That's so i can pack it, and not use a grid
        Radiobutton(self.projectDetailsFrame, text="Create a new study with this information:", variable=self.chooseProjectIntVar, value=2).pack()

        self.newProjectFrame = Frame(self.projectDetailsFrame)            

        self.studyIdInstrText = StringVar()
        self.studyIdInstrText.set('Short Study Identifier:')
        self.studyIdInstrLabel = Label(self.newProjectFrame, width=labelInputWidth, height=1, textvariable=self.studyIdInstrText).grid(row=0, column=0)
        self.inputStudyId = StringVar()
        self.inputStudyIdEntry = Entry(self.newProjectFrame, width=formInputWidth, textvariable=self.inputStudyId).grid(row=0, column=1)

        self.studyShortTitleInstrText = StringVar()
        self.studyShortTitleInstrText.set('Descriptive Study Title:')
        self.studyShortTitleInstrLabel = Label(self.newProjectFrame, width=labelInputWidth, height=1, textvariable=self.studyShortTitleInstrText).grid(row=1, column=0)
        self.inputStudyShortTitle = StringVar()
        self.inputStudyShortTitleEntry = Entry(self.newProjectFrame, width=formInputWidth, textvariable=self.inputStudyShortTitle).grid(row=1, column=1)

        self.studyAbstractInstrText = StringVar()
        self.studyAbstractInstrText.set('Study Description / Abstract:')
        self.studyAbstractInstrLabel = Label(self.newProjectFrame, width=labelInputWidth, height=1, textvariable=self.studyAbstractInstrText).grid(row=2, column=0)
        self.inputStudyAbstract = StringVar()
        self.inputStudyAbstractEntry = Entry(self.newProjectFrame, width=formInputWidth, textvariable=self.inputStudyAbstract).grid(row=2, column=1)
        
        self.newProjectFrame.pack()
        
        self.projectDetailsFrame.pack()

        # Make a frame for the save options button.
        self.saveOptionsFrame = Frame(self.interior)
        Button(self.saveOptionsFrame, text='Save Options', command=self.saveOptions).grid(row=0, column=0)
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
                raise Exception('Error loading EMBL-ENA submission options. Invalid Project choice:' + str(getConfigurationValue('choose_project')))
            
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
            logEvent ('Saving Options....')
            
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
            #logEvent('Not ready to save, you are missing options.')
            pass
        
    def checkOptions(self):
        
        # TODO: THis is really annoying. I should allow them to close this screen without providing all the information
        # Im sure a user will complain about this.
        
        #logEvent ('Checking options.')
       
        # Don't check the EMBL Username
        # Don't check the EMBL Password
        
        # TODO: I should be more strict with this enforcement.  I got rid of the False returns because it was too darn annoying.
        
        if (not self.inputSampleID.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing a Sample ID. Please try again.')
            #return False
            return True
        
        if (not self.inputGene.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing a Gene. Please try again.')
            #return False
            return True
        
        # Don't check the class boolean
        
        if (not self.inputAllele.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing an Allele Name. Please try again.')
            #return False
            return True
        
        if (str(self.chooseProjectIntVar.get()) == '1'):
            # Use Existing Project
            if (not self.inputStudyAccession.get()):
                messagebox.showwarning('Missing Form Value',
                    'You are missing a Study Accession number. Please try again.')
                #return False
                return True
            
        elif(str(self.chooseProjectIntVar.get()) == '2'):
            # Use New Project
            if (not self.inputStudyId.get()):
                messagebox.showwarning('Missing Form Value',
                    'You are missing a Study Name. Please try again.')
                #return False
                return True
            
            if (not self.inputStudyShortTitle.get()):
                messagebox.showwarning('Missing Form Value',
                    'You are missing a Study Description. Please try again.')
                #return False
                return True
            
            
            if (not self.inputStudyAbstract.get()):
                messagebox.showwarning('Missing Form Value',
                    'You are missing a Study Accession number. Please try again.')
                #return False
                return True
            
        else:
            raise Exception ('Unknown value of self.chooseProjectIntVar. I expect 1 or 2. Observed:' + str(self.chooseProjectIntVar))
        
        
        if (not self.inputAnalysisAlias.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing an Analysis Alias. Please try again.')
            #return False
            return True        
        
        if (not self.inputAnalysisTitle.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing an Analysis Title. Please try again.')
            #return False
            return True        
        
        if (not self.inputAnalysisDescription.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing an Analysis Description. Please try again.')
            #return False
            return True
        
    

        # All options look good, right?
        
        
        return True
    
    
    def closeWindow(self):
        self.parent.destroy()        
    