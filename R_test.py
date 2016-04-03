### R_test.py
import sys
from pyper import R
import unique

if 'darwin' in sys.platform:
    try:
        path = unique.filepath("AltDatabase/tools/R/Mac/R")
        r = R(RCMD=path,use_numpy=True)
        #print '1'
    except Exception:
        #print '2'
        r = R(use_numpy=True)
else:
    #print '3'
    r = R(use_numpy=True)

print_out = r('library("Biobase")')
if "Error" in print_out:
    print 'Installing the R package "Biobase" in Config/R'
    print_out = r('source("http://bioconductor.org/biocLite.R"); biocLite("Biobase")')
    if "Error" in print_out: print 'unable to download the package "Biobase"'; forceError 
    print_out = r('library("Biobase")')
            
print_out = r('library("hopach")')
if "Error" in print_out:
    print 'Installing the R package "hopach" in Config/R'
    print_out = r('source("http://bioconductor.org/biocLite.R"); biocLite("hopach")')
    if "Error" in print_out: print 'unable to download the package "hopach"'; forceError