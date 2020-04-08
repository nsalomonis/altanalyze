import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import export
import unique
import traceback

""" Intersecting Coordinate Files """

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def importLookupTable(fn):
    """ Import a gRNA to valid tag lookup table """
    lookup_table = []
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        gRNA,tag = t
        lookup_table.append((gRNA,tag))
    return lookup_table
  
def importCountMatrix(fn,mask=False): 
    """ Import a count matrix """
    classification = {}
    firstRow = True
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if firstRow:
            headers = t[1:]
            firstRow = False
        else:
            barcode = t[0]
            values = map(int,t[1:])
            if mask:
                sum_counts = sum(values[2:])
            else:
                sum_counts = sum(values)
            def threshold(val):
                if val>0.3: return 1
                else: return 0

            if sum_counts>0:
                ratios = map(lambda x: (1.000*x)/sum_counts, values)
                if mask:
                    original_ratios = ratios
                    ratios = ratios[2:] ### mask the first two controls which are saturating
                else:
                    original_ratios = ratios
                hits = map(lambda x: threshold(x), ratios)
                hits = sum(hits)
                if sum_counts>20 and hits == 1:
                    index=0
                    for ratio in ratios:
                        if ratio>0.3:
                            header = headers[index]
                        index+=1
                    classification[barcode] = header
    print len(classification),fn
    return classification

def exportGuideToTags(lookup_table,gRNA_barcode,tag_barcode,output):
    export_object = open(output,'w')
    for barcode in gRNA_barcode:
        gRNA = gRNA_barcode[barcode]
        if barcode in tag_barcode:
            tag = tag_barcode[barcode]
            if (gRNA,tag) in lookup_table:
                uid = tag+'__'+gRNA
                export_object.write(barcode+'\t'+uid+'\t'+uid+'\n')
    export_object.close()
                        
if __name__ == '__main__':
    ################  Comand-line arguments ################
    import getopt

    if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
        print 'WARNING!!!! Too commands supplied.'
        
    else:
        options, remainder = getopt.getopt(sys.argv[1:],'', ['species=','gRNA=', 'tag=', 'lookup=','output='])
        #print sys.argv[1:]
        for opt, arg in options:
            if opt == '--gRNA':
                gRNA = arg
            elif opt == '--tag':
                tag = arg
            elif opt == '--lookup':
                lookup = arg
            elif opt == '--output':
                output = arg
                
    lookup_table = importLookupTable(lookup)
    gRNA_barcode = importCountMatrix(gRNA)
    tag_barcode = importCountMatrix(tag)
    exportGuideToTags(lookup_table,gRNA_barcode,tag_barcode,output)
