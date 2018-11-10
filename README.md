# NCBI oneoff-tools
## minimally tested uses at your own risk

### pmid2csv
This simple do-hicky does the following:
   * Read a simple text file containing pubmed ID's
   * Retrieve the xml data affiliated with the PMID's from NCBI's pubmed via the api
   * Parse the xml, writing a csv formated as  
        | PMID | TITLE | AuthorName | AffiliationInfo1 | AffiliationInfo2 | ...  
        one row per PMID per AUTHOR
        
     **The windows executable will run on an internet connected Window's 10.**  
     For most issues it will *not* exit gracefully  

#### To use the executable
open up a console, the windows command prompt will work
use the command  
`> pmid2csv.exe -i <pmids.txt> -o <results.csv> `


or for help
`pmid2csv.exe -h`

#### Notes:
not yet tested with api-keys
input should handle all whitespace seperated PMIDS, however so far it has only been tested on one per line
12345678 PMID12345678 PMID:12345678

they python script was writen for python 3 (3.5) and depends on Beautifulsoup4 (bs4) and lxml.
it will probably work with 3.4+ 
