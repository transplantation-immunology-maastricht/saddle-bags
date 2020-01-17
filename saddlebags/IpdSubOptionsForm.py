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

from tkinter import Frame, Label, StringVar, IntVar, Entry, Radiobutton, messagebox, Button

from saddlebags.SaddlebagsConfig import getConfigurationValue, assignConfigurationValue
from saddlebags.ScrolledWindow import VerticalScrolledFrame

import logging

class IpdSubOptionsForm(VerticalScrolledFrame):
        
    # Initialize the GUI
    def __init__(self, root):
        
        
        VerticalScrolledFrame.__init__(self, root)
        #Frame.__init__(self, root)
        #super(500, 500)
        root.title("Choose IPD-IMGT/HLA Submission Options")
        self.parent = root

        #button_opt = {'fill': Tkconstants.BOTH, 'padx': 35, 'pady': 5}
                        
        # This window should not be resizeable. I guess.
        self.parent.resizable(width=False, height=False)
        #self.parent.resizable(width=True, height=True)
        
        # To define the exit behavior.  Save and exit.
        self.parent.protocol('WM_DELETE_WINDOW', self.saveOptions)
        
        # Define the return behavior.  Same as "close window" etc
        root.bind('<Return>', self.returnFunction)
        
        self.instructionsFrame = Frame(self.interior)  
        self.instructionText = StringVar()       
        self.instructionText.set('\nThese options are required for an IPD allele submission.\n'
            + 'Login Credentials will not be stored, but they will be sent to IPD via\n'
            + 'secure https connection.\n')        
        Label(self.instructionsFrame, width=85, height=6, textvariable=self.instructionText).pack()
        self.instructionsFrame.pack()
        
        #Standard Inputs widths for the form elements
        formInputWidth = 35
        labelInputWidth = 35
                
        # Make a frame to contain the input variables
        # self.interior is defined in the ScrolledWindow class
        self.submissionDetailsInputFrame = Frame(self.interior)

        self.usernameInstrText = StringVar()
        self.usernameInstrText.set('IPD Username:')
        self.usernameInstrLabel = Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.usernameInstrText).grid(row=0, column=0)
        self.inputUsername = StringVar()
        self.inputUsernameEntry = Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputUsername).grid(row=0, column=1)

        self.passwordInstrText = StringVar()
        self.passwordInstrText.set('IPD Password:')
        self.passwordInstrLabel = Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.passwordInstrText).grid(row=1, column=0)
        self.inputPassword = StringVar()
        self.inputPasswordEntry = Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputPassword, show="*").grid(row=1, column=1)
        
        
        # TODO: Submitter / Laboratory ID.  
        # This is on the IPD form.
        #Do I know this infromation? Do I need to tell user how to get it?
  
        self.sampleIDInstrText = StringVar()
        self.sampleIDInstrText.set('Cell/Sample ID:')
        self.sampleIDinstrLabel = Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.sampleIDInstrText).grid(row=2, column=0)
        self.inputSampleID = StringVar()
        self.inputSampleIDEntry = Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputSampleID).grid(row=2, column=1)

        self.geneInstrStringVar = StringVar()
        self.geneInstrStringVar.set('Gene:')
        self.geneInstrLabel = Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.geneInstrStringVar).grid(row=3, column=0)
        self.inputGene = StringVar()        
        self.inputGeneEntry = Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputGene).grid(row=3, column=1)

        self.chooseClassIntVar = IntVar()
        self.chooseClassIntVar.set(1)
        Radiobutton(self.submissionDetailsInputFrame, text="HLA Class I ", variable=self.chooseClassIntVar, value=1).grid(row=4, column=0)
        Radiobutton(self.submissionDetailsInputFrame, text="HLA Class II", variable=self.chooseClassIntVar, value=2).grid(row=4, column=1)

        self.alleleInstrText = StringVar()
        self.alleleInstrText.set('Allele Local Name:')
        self.alleleInstrLabel = Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.alleleInstrText).grid(row=5, column=0)
        self.inputAllele = StringVar() 
        self.inputAlleleEntry = Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputAllele).grid(row=5, column=1)
        
        
        
        
        # New form stuff
        # Gotta add this to the load/save config nonsense below.
        
        
        # TODO: Can I just load an ENA accession? I think that is possible.  Easier than filling it in here
        
        
        # TODO: When ENA Sequence Accession # Is provided, I can probably lookup an annotated sequence.
        # Should I put a button next to this field 
        # Button: "Lookup This ENA Sequence Accession #"
        # If it is found, then i already know the sequence with exon boundaries.
        
        
        # TODO: Do I need to specify if it is ENA / Genbank / The other one?  Probably not.
        # I can require an ENA code and disregard Genbank.
        # Radio Buttons?  
        # ENA / Genbank Accession #
        # No, this tool is for ENA submission. But this is a question for James Robinson.
        # Should i choose between which intermediate databse they use?
        self.enaAccInstrText = StringVar()
        self.enaAccInstrText.set('ENA Sequence Accession #:')
        self.enaAccInstrLabel = Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.enaAccInstrText).grid(row=6, column=0)
        self.inputEnaAcc = StringVar()
        self.inputEnaAccEntry = Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputEnaAcc).grid(row=6, column=1)

        
        # Release Date
        self.releaseDateInstrText = StringVar()
        self.releaseDateInstrText.set('IPD Release Date:')
        self.releaseDateInstrLabel = Label(self.submissionDetailsInputFrame, width=labelInputWidth, height=1, textvariable=self.releaseDateInstrText).grid(row=7, column=0)
        self.inputReleaseDate = StringVar() 
        self.inputReleaseDateEntry = Entry(self.submissionDetailsInputFrame, width=formInputWidth, textvariable=self.inputReleaseDate).grid(row=7, column=1)
        
        # Reference Details
        # Is this allele in a published paper or not?
        # 0=unpublished, 1=published
        self.publishedReferenceIntVar = IntVar()
        self.publishedReferenceIntVar.set(0)
        
        self.submissionDetailsInputFrame.pack()
 
        
        self.unpublishedReferenceFrame = Frame(self.interior)   
        
        self.referenceInstrText = StringVar()
        self.referenceInstrText.set('\nPlease provide some information about a\npublished paper relevant to this sequence.\n')
        self.referenceInstrLabel = Label(self.unpublishedReferenceFrame, width=70, height=4, textvariable=self.referenceInstrText).pack()#.grid(row=2, column=0)
             
        Radiobutton(self.unpublishedReferenceFrame, text="No Published Reference.", variable=self.publishedReferenceIntVar, value=0).pack()
        self.unpublishedReferenceFrame.pack()

        self.publishedReferenceFrame = Frame(self.interior)
        
        # Radio Button: Published
        Radiobutton(self.unpublishedReferenceFrame, text="Use This Reference:", variable=self.publishedReferenceIntVar, value=1).pack()
        
        # Reference Title
        self.referenceTitleInstrText = StringVar()
        self.referenceTitleInstrText.set('Reference Title:')
        self.referenceTitleInstrLabel = Label(self.publishedReferenceFrame, width=labelInputWidth, height=1, textvariable=self.referenceTitleInstrText).grid(row=1, column=0)
        self.inputReferenceTitle = StringVar() 
        self.inputReferenceTitleEntry = Entry(self.publishedReferenceFrame, width=formInputWidth, textvariable=self.inputReferenceTitle).grid(row=1, column=1)
        
        # Authors
        self.referenceAuthorsInstrText = StringVar()
        self.referenceAuthorsInstrText.set('Reference Authors:')
        self.referenceAuthorsInstrLabel = Label(self.publishedReferenceFrame, width=labelInputWidth, height=1, textvariable=self.referenceAuthorsInstrText).grid(row=2, column=0)
        self.inputReferenceAuthors = StringVar() 
        self.inputReferenceAuthorsEntry = Entry(self.publishedReferenceFrame, width=formInputWidth, textvariable=self.inputReferenceAuthors).grid(row=2, column=1)
        
        # Journal
        self.referenceJournalInstrText = StringVar()
        self.referenceJournalInstrText.set('Reference Journal:')
        self.referenceJournalInstrLabel = Label(self.publishedReferenceFrame, width=labelInputWidth, height=1, textvariable=self.referenceJournalInstrText).grid(row=3, column=0)
        self.inputReferenceJournal = StringVar() 
        self.inputReferenceJournalEntry = Entry(self.publishedReferenceFrame, width=formInputWidth, textvariable=self.inputReferenceJournal).grid(row=3, column=1)

        self.publishedReferenceFrame.pack()
               
        # Make a frame to contain the input variables.
        # I had to make 2 of them to organize my gui, maybe I can name this better.
        self.submissionDetailsInputFrame2 = Frame(self.interior)
            
        # /alignment -> defined by IPD sequence alignment service
        # In this case, it is the closest known allele.
        self.closestAlleleInstrText = StringVar()
        self.closestAlleleInstrText.set('Closest Known HLA Allele:')
        self.closestAlleleInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.closestAlleleInstrText).grid(row=1, column=0)
        self.inputClosestAllele = StringVar() 
        self.inputClosestAlleleEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputClosestAllele).grid(row=1, column=1)

        # Written Description
        # Looks like this is a description of how the sequence differes from closest knnown allele
        self.closestAlleleWrittenDescriptionInstrText = StringVar()
        self.closestAlleleWrittenDescriptionInstrText.set('Differences from Closest Allele:')
        self.closestAlleleWrittenDescriptionInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.closestAlleleWrittenDescriptionInstrText).grid(row=2, column=0)
        self.inputClosestAlleleWrittenDescription = StringVar() 
        self.inputClosestAlleleWrittenDescriptionEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputClosestAlleleWrittenDescription).grid(row=2, column=1)

        
        # DONOR INFORMATION
        
        # Cell ID (cellnum)
        # Wait, is this the same as the sample ID? Should I move the sample ID field down here?
        # No. I am disregarding this sample ID.
        
        # Ethnic Origin - Text
        self.ethnicOriginInstrText = StringVar()
        self.ethnicOriginInstrText.set('Ethnic Origin:')
        self.ethnicOriginInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.ethnicOriginInstrText).grid(row=3, column=0)
        self.inputEthnicOrigin = StringVar() 
        self.inputEthnicOriginEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputEthnicOrigin).grid(row=3, column=1)

        # Sex - Text
        self.sexInstrText = StringVar()
        self.sexInstrText.set('Sex:')
        self.sexInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.sexInstrText).grid(row=4, column=0)
        self.inputSex = StringVar() 
        self.inputSexEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputSex).grid(row=4, column=1)

        # TODO Make a boolean
        # Consanguineous (T/F)
        self.consanguineousInstrText = StringVar()
        self.consanguineousInstrText.set('Sample is Consanguineous:')
        self.consanguineousInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.consanguineousInstrText).grid(row=5, column=0)
        self.inputConsanguineous = StringVar() 
        self.inputConsanguineousEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputConsanguineous).grid(row=5, column=1)
        
        # TODO Make a boolean
        # Homozygous (T/F)
        # TODO: Accepted values are 'Yes', 'No', 'Unknown'
        # Make dropdown for this, or radio buttons.
        self.homozygousInstrText = StringVar()
        self.homozygousInstrText.set('Sample is Homozygous:')
        self.homozygousInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.homozygousInstrText).grid(row=6, column=0)
        self.inputHomozygous = StringVar() 
        self.inputHomozygousEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputHomozygous).grid(row=6, column=1)

        # Lab of Origin (text)
        self.labOriginInstrText = StringVar()
        self.labOriginInstrText.set('Lab of Origin:')
        self.labOriginInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.labOriginInstrText).grid(row=7, column=0)
        self.inputLabOrigin = StringVar() 
        self.inputLabOriginEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputLabOrigin).grid(row=7, column=1)
        
        # Lab Contact
        self.labContactInstrText = StringVar()
        self.labContactInstrText.set('Lab Contact:')
        self.labContactInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.labContactInstrText).grid(row=8, column=0)
        self.inputLabContact = StringVar() 
        self.inputLabContactEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputLabContact).grid(row=8, column=1)
              
        # Cell Availability
        # Material Available (T/F)
        self.materialAvailableInstrText = StringVar()
        self.materialAvailableInstrText.set('Material Availability:')
        self.materialAvailableInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.materialAvailableInstrText).grid(row=9, column=0)
        self.inputMaterialAvailable = StringVar() 
        self.inputMaterialAvailableEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputMaterialAvailable).grid(row=9, column=1)
        
        # Cell Bank (Text)
        self.cellBankInstrText = StringVar()
        self.cellBankInstrText.set('Cell Bank:')
        self.cellBankInstrLabel = Label(self.submissionDetailsInputFrame2, width=labelInputWidth, height=1, textvariable=self.cellBankInstrText).grid(row=10, column=0)
        self.inputCellBank = StringVar() 
        self.inputCellBankEntry = Entry(self.submissionDetailsInputFrame2, width=formInputWidth, textvariable=self.inputCellBank).grid(row=10, column=1)
        
        # Cell Workshop Details
        # I think Cell Workshop Details is just a header. there isn't new information here, just a header. 
        # TODO: Compare this with the IPD Submission website, im not missing something?
        
        self.submissionDetailsInputFrame2.pack()
        
        # numbering these input frames are not informative. Oh well. This frame has HLA Allele calls on it.
        self.submissionDetailsInputFrame3 = Frame(self.interior)

         
        # Alternative HLA DNA Typing
        # Dropdown Box with another Entry Field?
        # I need:
        # Label
        self.sourceHLAInstrText = StringVar()
        self.sourceHLAInstrText.set('Source HLA Types (Sequenced):')
        self.sourceHLAInstrLabel = Label(self.submissionDetailsInputFrame3, width=labelInputWidth, height=1, textvariable=self.sourceHLAInstrText).grid(row=1, column=0)
        # Combo Box, with source_hla dictionary keys. Sorted.
        # Text input, with the gene specified.
        # "Clear" button. clear out all allele calls.
        # Configuration should be assigned whenever text changes.
        # I Think i need a new panel for this. Yeah.
        
        
        # Source Serology Typing
        # Maybe the same as DNA typing?
        # Ignore for now.
        
        # Sequencing Methods
            
        # Primers
        # This is probably a Dropdown with Entry field also.
        
        # TODO: Comments.  Where does this stuff go?  This is details about the lab of origin. I haven't tried specifying this one yet, ask James how to do it.
        # Comments
        
        
        
        
        self.submissionDetailsInputFrame3.pack()


        

        # Make a frame for the save options button.
        self.saveOptionsFrame = Frame(self.interior)
        Button(self.saveOptionsFrame, text='Save Options', command=self.saveOptions).grid(row=0, column=0)
        self.saveOptionsFrame.pack()
        
        self.loadOptions()
        
    # I needed a function for the return keypress to latch onto.
    # It is just a wrapper for the saveOptions method.
    def returnFunction(self, event):
        self.saveOptions()

    # submissionOptions is a dictionary, passed by the parent.
    def loadOptions(self):
        if getConfigurationValue('ipd_username') is not None:
            self.inputUsername.set(getConfigurationValue('ipd_username'))
            
        if getConfigurationValue('ipd_password') is not None:
            self.inputPassword.set(getConfigurationValue('ipd_password'))
            
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
                raise Exception('Error loading IPD submission options. Invalid class:' + str(getConfigurationValue('class')))
    
        if getConfigurationValue('allele_name') is not None:
            self.inputAllele.set(getConfigurationValue('allele_name'))
            
        if getConfigurationValue('ena_sequence_accession') is not None:
            self.inputEnaAcc.set(getConfigurationValue('ena_sequence_accession'))
        
        if getConfigurationValue('ena_release_date') is not None:
            self.inputReleaseDate.set(getConfigurationValue('ena_release_date'))
   
        # 0=unpublished, 1=published
        #print('1Setting is_published value to:' + getConfigurationValue('is_published'))
        if (getConfigurationValue('is_published') is None or getConfigurationValue('is_published') == 'None'):
            self.publishedReferenceIntVar.set(0)
            #print('2Setting is_published value to:' + getConfigurationValue('is_published'))

        else:
            self.publishedReferenceIntVar.set(getConfigurationValue('is_published'))


        if getConfigurationValue('reference_title') is not None:
            self.inputReferenceTitle.set(getConfigurationValue('reference_title'))            
        if getConfigurationValue('reference_authors') is not None:
            self.inputReferenceAuthors.set(getConfigurationValue('reference_authors'))            
        if getConfigurationValue('reference_journal') is not None:
            self.inputReferenceJournal.set(getConfigurationValue('reference_journal'))
     
        if getConfigurationValue('reference_journal') is not None:
            self.inputReferenceJournal.set(getConfigurationValue('reference_journal'))

        if getConfigurationValue('closest_known_allele') is not None:
            self.inputClosestAllele.set(getConfigurationValue('closest_known_allele'))
        if getConfigurationValue('closest_allele_written_description') is not None:
            self.inputClosestAlleleWrittenDescription.set(getConfigurationValue('closest_allele_written_description'))
          
        if getConfigurationValue('ethnic_origin') is not None:
            self.inputEthnicOrigin.set(getConfigurationValue('ethnic_origin'))          
        if getConfigurationValue('sex') is not None:
            self.inputSex.set(getConfigurationValue('sex'))              
        if getConfigurationValue('consanguineous') is not None:
            self.inputConsanguineous.set(getConfigurationValue('consanguineous'))          
        if getConfigurationValue('homozygous') is not None:
            self.inputHomozygous.set(getConfigurationValue('homozygous'))  

        if getConfigurationValue('lab_of_origin') is not None:
            self.inputLabOrigin.set(getConfigurationValue('lab_of_origin'))          
        if getConfigurationValue('lab_contact') is not None:
            self.inputLabContact.set(getConfigurationValue('lab_contact')) 
            
        if getConfigurationValue('material_availability') is not None:
            self.inputMaterialAvailable.set(getConfigurationValue('material_availability'))          
        if getConfigurationValue('cell_bank') is not None:
            self.inputCellBank.set(getConfigurationValue('cell_bank')) 
            
        # TODO: 
        # Load options for HLA allele calls.
        # Clear the combo-box, and fill it with the keys of my hla allele call dictionary.  
            
    def saveOptions(self):
        # Close the window
        # TODO: If i force the user to fill in all the options, this form is really obnoxious.
        # Instead they should be allowed to close it and I will still warn them.
        # I can re-think this plan if people are trying to submit bad data.
        
        #if (self.checkOptions()):
        
        #Don't force user to fill in all the options:
        self.checkOptions()
        if(True):
            
            
            logging.info ('Saving Options....')
            
            assignConfigurationValue('ipd_username', self.inputUsername.get())
            # I store this password so I can use it in the submission
            # I don't ever want to save the password. Make sure it isn't being saved in the config, in AlleleSubCommon.py
            assignConfigurationValue('ipd_password', self.inputPassword.get())
            assignConfigurationValue('sample_id', self.inputSampleID.get())
            assignConfigurationValue('gene', self.inputGene.get())
            assignConfigurationValue('class', str(self.chooseClassIntVar.get()))             
            assignConfigurationValue('allele_name', self.inputAllele.get())

            assignConfigurationValue('ena_sequence_accession', self.inputEnaAcc.get())
            assignConfigurationValue('ena_release_date', self.inputReleaseDate.get())
            
            assignConfigurationValue('is_published', str(self.publishedReferenceIntVar.get()))
            #print('Saving is_published configuration as :' + str(self.publishedReferenceIntVar.get()))
            
            assignConfigurationValue('reference_title',self.inputReferenceTitle.get())    
            assignConfigurationValue('reference_authors',self.inputReferenceAuthors.get())    
            assignConfigurationValue('reference_journal',self.inputReferenceJournal.get())            
           
            assignConfigurationValue('closest_known_allele', self.inputClosestAllele.get())
            assignConfigurationValue('closest_allele_written_description', self.inputClosestAlleleWrittenDescription.get())
      
            assignConfigurationValue('ethnic_origin', self.inputEthnicOrigin.get())
            assignConfigurationValue('sex', self.inputSex.get())
            
            # TODO: Accepted values are 'Yes', 'No', 'Unknown'
            # Make dropdown for these
            assignConfigurationValue('consanguineous', self.inputConsanguineous.get())
            assignConfigurationValue('homozygous', self.inputHomozygous.get())
            
            assignConfigurationValue('lab_of_origin', self.inputLabOrigin.get())
            assignConfigurationValue('lab_contact', self.inputLabContact.get())
            
            assignConfigurationValue('material_availability', self.inputMaterialAvailable.get())
            assignConfigurationValue('cell_bank', self.inputCellBank.get())
   
            # I have saved hla calls in a dictionary. They should have been saved individually.
   
            self.parent.destroy() 
            
            
        else:
            #logging.info('Not ready to save, you are missing options.')
            pass
        
    def checkOptions(self):
        # TODO this method
        logging.info ('Checking options.')

        # Don't check the IPD Username
        # Don't check the IPD Password
        
        if (not self.inputSampleID.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing a Sample ID. Please try again.')
            return False        
        if (not self.inputGene.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing a Gene. Please try again.')
            return False
        if (not self.inputAllele.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing an Allele Name. Please try again.')
            return False
        
        if (not self.inputEnaAcc.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing an ENA Accession Number. Please try again.')
            return False
        if (not self.inputReleaseDate.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing an IPD Submission Release Date. Please try again.')
            return False

        # This is NOne if nothing is selected.
        if (self.publishedReferenceIntVar.get() == 0):
        # unpublished, nothing to check
            pass
        else:
            if ((not self.inputReferenceTitle.get())
                or (not self.inputReferenceAuthors.get())
                or (not self.inputReferenceJournal.get())
                ):
                messagebox.showwarning('Missing Form Value',
                    'You are must supply information about the published Reference. Please try again.')
                return False

        if (not self.inputClosestAllele.get()):
            messagebox.showwarning('Missing Form Value',
                'You are missing the closest known reference allele to this sequence. Please provide this information.')
            return False
        if (not self.inputEthnicOrigin.get()):
            messagebox.showwarning('Missing Form Value',
                'Please provide a description of an ethnic origin for this sample.')
            return False
        if (not self.inputSex.get()):
            messagebox.showwarning('Missing Form Value',
                'Please identify the sex for this sample.')
            return False
        
        # TODO: Accepted values are 'Yes', 'No', 'Unknown' I think
        if (not self.inputConsanguineous.get()):
            messagebox.showwarning('Missing Form Value',
                'Please indicate if the sample is consanguineous or not.')
            return False
        if (not self.inputHomozygous.get()):
            messagebox.showwarning('Missing Form Value',
                'Please indicate if the sample is homozygous or not.')
            return False
        
        
        if (not self.inputLabOrigin.get()):
            messagebox.showwarning('Missing Form Value',
                'Please provide the name of the submitting laboratory.')
            return False
        if (not self.inputLabContact.get()):
            messagebox.showwarning('Missing Form Value',
                'Please provide the name of the laboratory contact.')
            return False

        if (not self.inputMaterialAvailable.get()):
            messagebox.showwarning('Missing Form Value',
                'Please indicate if the cell material is available.')
            return False
        if (not self.inputCellBank.get()):
            messagebox.showwarning('Missing Form Value',
                'Please provide the name of the cell bank where the sample can be found.')
            return False

        # TODO Validate the HLA ALlele calls. I won't do IMGT/HLA validation, I will leave that validation up to IMGT/HLA
        # Validate A, B, DRB1. The rest, I don't care.

        # All options look good, right?
        return True
    
    
    def closeWindow(self):
        #writeConfigurationFile()

        self.parent.destroy()        
    