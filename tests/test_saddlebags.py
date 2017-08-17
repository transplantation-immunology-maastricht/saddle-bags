#from hf_cypher.cypherQuery import CypherQuery
from nose.tools import assert_equal, assert_true

#import inspect
import os
#import json

from saddlebags.AlleleSubCommon import assignConfigurationValue, annotateRoughInputSequence
from saddlebags.SubmissionGeneratorEMBL import SubmissionGeneratorEMBL


#json_file = os.path.join(os.path.dirname(__file__), '')
#file_data = open(json_file).read()
#expected = json.loads(file_data)


def testCreateEMBLSubmissionFlatfile():
    print ('Test: Create an EMBL SubmissionFlatfile')
    assert_true(True)
        
    roughFeatureSequence = 'aag\nCGTCGT\nccg\nGGCTGA\naat'
    
    assignConfigurationValue('sample_id', 'Donor_12345')
    assignConfigurationValue('gene','HLA-C')
    assignConfigurationValue('class','1')
    assignConfigurationValue("allele_name",'Allele:01:02')

    allGen = SubmissionGeneratorEMBL()
    #roughFeatureSequence = self.featureInputGuiObject.get('1.0', 'end')

    allGen.sequenceAnnotation = annotateRoughInputSequence(roughFeatureSequence)

    enaSubmission = allGen.buildENASubmission()
    
    assert_true(len(enaSubmission) > 3)
    assert_true(enaSubmission is not None)
                
                
def testAnnotateRoughSequence():
    print ('Test: Annotate a rough sequence using GFE')
    assert_true(True)
    