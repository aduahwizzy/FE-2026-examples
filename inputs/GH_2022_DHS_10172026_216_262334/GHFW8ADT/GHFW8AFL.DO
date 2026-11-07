infile using "C:\USERS\21429\ICF\DATA ARCHIVE - DOCUMENTS\EXPORT MENU\GH8A\GHFW8AFL.dct"

label variable fw101    "Fieldworker code"
label variable fw000    "Country code and phase"
label variable fw102    "Fieldworker region of residence"
label variable fw103    "Fieldworker type of place of residence"
label variable fw104    "Fieldworker age"
label variable fw105    "Fieldworker sex"
label variable fw106    "Fieldworker marital status"
label variable fw107    "Fieldworker number of living children"
label variable fw108    "Fieldworker ever had a child who died"
label variable fw109    "Highest level of school attended by fieldworker"
label variable fw110    "Highest grade/form/year at that level completed by fieldworker"
label variable fw110a   "Received clinical, medical, or laboratory training or worked in healthcare"
label variable fw110b   "Current occupational category or qualification"
label variable fw111    "Fieldworker religion"
label variable fw112    "Fieldworker ethnicity"
label variable fw113a   "AKAN"
label variable fw113b   "GA/DANGME"
label variable fw113c   "EWE"
label variable fw113d   "DAGBANI"
label variable fw113e   "NZEMA"
label variable fw113f   "BRONG"
label variable fw113g   "Language G"
label variable fw113h   "Language H"
label variable fw113i   "Language I"
label variable fw113j   "Language J"
label variable fw113k   "Language K"
label variable fw113l   "Language L"
label variable fw113m   "Language M"
label variable fw113x   "Other language"
label variable fw113y   "No other language"
label variable fw114    "Fieldworker's mother tongue/native language"
label variable fw115a   "Has fieldworker ever worked on a DHS prior to this survey"
label variable fw115b   "Has fieldworker ever worked on an MIS prior to this survey"
label variable fw115c   "Has fieldworker ever worked on any other survey prior to this survey"
label variable fw116    "Was fieldworker working for implementing agency at the time employed for DHS"
label variable fw117    "Is fieldworker a permanent or temporary employee of agency"
label variable fw118    "Does fieldworker has any comments"

#delimit ;
label define FW102   
     1 "Western"
     2 "Central"
     3 "Greater accra"
     4 "Volta"
     5 "Eastern"
     6 "Ashanti"
     7 "Western north"
     8 "Ahafo"
     9 "Bono"
    10 "Bono east"
    11 "Oti"
    12 "Northern"
    13 "Savannah"
    14 "North east"
    15 "Upper east"
    16 "Upper west"
;
label define FW103   
     1 "City"
     2 "Town"
     3 "Rural"
;
label define FW105   
     1 "Male"
     2 "Female"
;
label define FW106   
     1 "Currently married"
     2 "Living with a man/woman"
     3 "Widowed"
     4 "Divorced"
     5 "Separated"
     6 "Never married or lived with a man/woman"
;
label define FW108   
     1 "Yes"
     2 "No"
;
label define FW109   
     1 "Primary"
     2 "Middle"
     3 "JSS/JHS"
     4 "Secondary"
     5 "SSS/SHS"
     6 "Higher"
;
label define FW110   
     0 "Less than one year at level"
;
label define FW110A  
     1 "Yes"
     2 "No"
;
label define FW110B  
     1 "MEDICAL DOCTOR"
     2 "ASSISTANT MEDICAL OFFICER"
     3 "CLINICAL OFFICER"
     4 "ASSISTANT CLINICAL OFFICER"
     5 "REGISTERED NURSE/MIDWIFE"
     6 "ENROLLED NURSE/MIDWIFE"
     7 "NURSE ASSISTANT/ATTENDANT"
     8 "LABORATORY SCIENTIST"
     9 "LABORATORY TECHNOLOGIST"
    10 "LABORATORY TECHNICIAN"
    11 "LABORATORY ASSISTANT"
    95 "NO TECHNICAL QUALIFICATION"
    96 "OTHER"
;
label define FW111   
     1 "Catholic"
     2 "Anglican"
     3 "Methodist"
     4 "Presbyterian"
     5 "Pentecostal/charismatic"
     6 "Other Christian"
     7 "Islam"
     8 "Traditional/spiritualist"
    95 "No religion"
    96 "Other"
;
label define FW112   
     1 "AKAN"
     2 "GA/DANGME"
     3 "EWE"
     4 "GUAN"
     5 "MOLE-DAGBANI"
     6 "GRUSI"
     7 "GURMA"
     8 "MANDE"
    96 "Other"
;
label define FW114   
     1 "AKAN"
     2 "GA/DANGME"
     3 "EWE"
     4 "DAGBANI"
     5 "NZEMA"
     6 "BRONG"
    96 "Other"
;
label define FW115A  
     1 "Yes"
     2 "No"
;
label define FW115B  
     1 "Yes"
     2 "No"
;
label define FW115C  
     1 "Yes"
     2 "No"
;
label define FW116   
     1 "Yes, implementing agency 1"
     2 "Yes, implementing agency 2"
     3 "No"
;
label define FW117   
     1 "Permanent"
     2 "Temporary"
;
label define FW118   
     1 "Yes"
     2 "No"
;

#delimit cr
label values fw102    FW102   
label values fw103    FW103   
label values fw105    FW105   
label values fw106    FW106   
label values fw108    FW108   
label values fw109    FW109   
label values fw110    FW110   
label values fw110a   FW110A  
label values fw110b   FW110B  
label values fw111    FW111   
label values fw112    FW112   
label values fw114    FW114   
label values fw115a   FW115A  
label values fw115b   FW115B  
label values fw115c   FW115C  
label values fw116    FW116   
label values fw117    FW117   
label values fw118    FW118   
