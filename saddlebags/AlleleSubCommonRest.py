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

from saddlebags.AlleleSubCommon import getConfigurationValue, cleanSequence

#from io import StringIO
from io import BytesIO
from pycurl import Curl, HTTPHEADER, SSL_VERIFYHOST, SSL_VERIFYPEER


def fetchSequenceAlleleCallWithGFE(rawSequence, locus):
    #roughFeatureSequence
    
    #print ('annotating this sequence: ' + str(rawSequence) + '\nat locus: ' + str(locus))
    #curl -X GET --header 'Accept: application/json' 'http://act.b12x.org/act?locus=HLA-A&sequence=GAT'

    cleanedSequence = cleanSequence(rawSequence.upper())
    
    #print('cleanedSequence:' + cleanedSequence)
    
    # TODO: Wasn't I worried about passing large strings in the url adresses pycurl? 
    # I thought that caused problems...maybe just in windows?
    # If the request is too long, the problem is that HLA is too big.
    requestURL = (str(getConfigurationValue('nmdp_act_rest_address')) 
        + '?locus=' + str(locus)
        + '&sequence=' + str(cleanedSequence))
    
    #print ('request URL:' + requestURL)
    
    
    #curlResponseBuffer = StringIO()
    curlResponseBuffer = BytesIO()    
    
    curlObject = Curl()
    curlObject.setopt(curlObject.URL, requestURL)
    #curlObject.setopt(curlObject.GET, 1)
    #curlObject.setopt(curlObject.HTTPPOST, POST_DATA)
    curlObject.setopt(curlObject.USERAGENT, 'Curl')
    curlObject.setopt(curlObject.WRITEFUNCTION, curlResponseBuffer.write)
    curlObject.setopt(HTTPHEADER, ['Accept:application/json'])
    # TODO:
    # Insecure.  Any security experts want to make this better?
    curlObject.setopt(SSL_VERIFYHOST, 0)
    curlObject.setopt(SSL_VERIFYPEER, 0)
    curlObject.perform()
    curlObject.close()

    responseText = curlResponseBuffer.getvalue().decode('utf-8')
        
    #print ('Received Response Text:\n'
    #    + str(responseText))
    
    return responseText
    

