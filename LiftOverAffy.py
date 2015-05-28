###LiftOverAffy

import sys, os, string
import export
import unique

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    #data = string.replace(data,'"','')
    return data

def covertAffyFormatToBED(filename, ConversionDB=None):
    print 'processing:',filename
    parent = export.findParentDir(filename)
    if ConversionDB==None:
        output_file = 'simple_chr.bed'
    else:
        output_file = export.findFilename(filename)
        output_file = string.replace(output_file,'mm9','mm10')
    export_obj = export.ExportFile(parent+'/'+output_file)
    fn=filepath(filename); entry_count=0; readfiles = False
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if data[0]=='#': readfiles = False
        elif readfiles==False:
            readfiles = True
            if ConversionDB!=None:
               export_obj.write(line) ### Write header 
        else:
            try:
                t = string.split(data[1:-1],'","')
                probeset_id,chr,strand,start,stop = t[:5]
                int(start)
                if ConversionDB==None:
                    if 'chr' in chr:
                        export_obj.write(chr+'\t'+start+'\t'+stop+'\t'+probeset_id+'\n')
                else:
                    chr,start,stop = ConversionDB[probeset_id]
                    t = [probeset_id,chr,strand,start,stop] + t[5:]
                    values = '"'+string.join(t,'","')+'"\n'
                    export_obj.write(values)
                entry_count+=1
            except Exception:
                pass
    export_obj.close()
    print entry_count, 'entries saved to:',parent+'/'+output_file
    
def importConvertedBED(filename):
    print 'processing:',filename
    parent = export.findParentDir(filename)
    fn=filepath(filename); entry_count=0; newCoordinates={}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if data[0]!='#':
            try:
                t = string.split(data,'\t')
                chr,start,stop,probeset_id = t
                int(start)
                if 'chr' in chr:
                    entry_count+=1
                newCoordinates[probeset_id] = chr,start,stop
            except ZeroDivisionError:
                pass
    print entry_count, 'imported and saved.'
    return newCoordinates
     
filename = '/Users/saljh8/Downloads/MoGene-1_0-st-v1-1/MoGene-1_0-st-v1.na33.2.mm9.probeset.csv'
bed_output = '/Users/saljh8/Downloads/MoGene-1_0-st-v1-1/input.bed'
#covertAffyFormatToBED(filename);sys.exit()

newCoordinates = importConvertedBED(bed_output)
covertAffyFormatToBED(filename,ConversionDB=newCoordinates)