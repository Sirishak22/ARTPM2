###################################################
# Authentication class for UMLS. Modeled after UMLS API Documentation.
# Lowell Milliken
###################################################
import requests
import lxml.html as lh
from lxml.html import fromstring

uri = 'https://utslogin.nlm.nih.gov/cas/v1/api-key'
head = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain', 'User-Agent': 'python'}


class Authentication:

    def __init__(self, apikeyfile):
        self.apikeyfile = apikeyfile

    def get_ticket_granting_ticket(self):
        with open(self.apikeyfile, 'r') as keyfile:
            apikey = keyfile.readline().strip()
            print(apikey)
        params = {'apikey': apikey}
        #print(apikey)
        #print("No error yet")
        request = requests.post(uri, data=params, headers=head)
        print(request)
        response = fromstring(request.text)
        #print("error here")
        return response.xpath('//form/@action')[0]

    def get_service_ticket(self, tgt):
        params = {'service': 'http://umlsks.nlm.nih.gov'}
        head = {'Content-type': 'application/x-www-form-urlencoded', 'Accept': 'text/plain', 'User-Agent': 'python'}
        request = requests.post(tgt, data=params, headers=head)
        return request.text
