#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 13:20:35 2019

@author: lowellmilliken
"""

"""
totalwordmatches=0
base="Liposarcoma CDK4 Amplification 38 year old male GERD"
abstract="We reported a 36-year-old woman with metastatic liposarcoma originating in the retroperitoneum, which responded well to adjuvant chemotherapy. The primary tumor was removed by surgery. Two months later, the patient developed metastasis to the brain, and to the lung four months later. Metastatic liposarcomas to the brain generally are extremely rare. The patient was treated with combination chemotherapy using cyclophosphamide, vincristine, adriamycin, and dacarbazine (CYVADIC). After she was examined, the former two drugs were alternated with vindesine and ifosfamide, and another regimen with cisplatin and etoposide was given after a three-week interval. As a result, both of the metastases totally disappeared. No recurrent lesion has been noted for two years. Although the role of chemotherapy for liposarcoma has not been well defined and little data support its use in an adjuvant setting, this combination chemotherapy seemed to be effective for advanced liposarcoma."
for querywords in base.split():
            if querywords.lower() in abstract.lower().split():
                totalwordmatches = totalwordmatches + abstract.lower().split().count(querywords.lower())
print(totalwordmatches)
print(abstract.split())
publications=['Adaptive Clinical Trial', 'Address', 'Autobiography', 'Bibliography', 'Biography', 'Case Reports', 'Classical Article', 'Clinical Conference', 'Clinical Study', 'Clinical Trial', 'Clinical Trial, Phase I', 'Clinical Trial, Phase II', 'Clinical Trial, Phase III', 'Clinical Trial, Phase IV', 'Clinical Trial Protocol', 'Clinical Trial, Veterinary', 'Collected Works', 'Comparative Study', 'Congress', 'Consensus Development Conference', 'Consensus Development Conference, NIH', 'Controlled Clinical Trial', 'Dataset', 'Dictionary', 'Directory', 'Duplicate Publication', 'Editorial', 'English Abstract', 'Equivalence Trial', 'Evaluation Studies', 'Expression of Concern', 'Festschrift', 'Government Document', 'Guideline', 'Historical Article', 'Interactive Tutorial', 'Interview', 'Introductory Journal Article', 'Journal Article', 'Lecture', 'Legal Case', 'Legislation', 'Letter', 'Meta-Analysis', 'Multicenter Study', 'News', 'Newspaper Article', 'Observational Study', 'Observational Study, Veterinary', 'Overall', 'Patient Education Handout', 'Periodical Index', 'Personal Narrative', 'Portrait', 'Practice Guideline', 'Pragmatic Clinical Trial', 'Publication Components', 'Publication Formats', 'Publication Type Category', 'Randomized Controlled Trial', 'Research Support, American Recovery and Reinvestment Act', 'Research Support, N.I.H., Extramural', 'Research Support, N.I.H., Intramural', "Research Support, Non-U.S. Gov't Research Support, U.S. Gov't, Non-P.H.S.", "Research Support, U.S. Gov't, P.H.S.", 'Review', 'Scientific Integrity Review', 'Study Characteristics', 'Support of ResearchSystematic Review', 'Technical Report', 'Twin Study', 'Validation Studies', 'Video-Audio Media', 'Webcasts']
print(len(publications))
"""
"""
import scholarly
search_query = scholarly.search_pubs_query("A case of metastatic liposarcoma originating in the retroperitoneum successfully treated with combination chemotherapy")
#print(search_query)
        # feature 29
val=next(search_query).citedby
print(val)
"""
pubtype="Journal Article,Research Support, Non-U.S. Gov't,Research Support, U.S. Gov't, P.H.S."
i=0
pubtypes=[]
start=i
for char in pubtype:
   if char==',' and pubtype[i+1]!=" ":
       pubtypes.append(pubtype[start:i])
       start=i+1
       i=i+1
   else:
       i=i+1
print(pubtypes)