### cellHarmony_validate.py
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import traceback

import numpy
import scipy
import time
import unique
from stats_scripts import statistics
import sys, os, string
import export
import UI

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def rankExpressionRescueFromCellHarmony(organized_diff_ref, repair1_folds, repair2_folds, reference_fold_dir, repair_dir1, repair_dir2, mode='validate'):

    def importCellHarmonyDEGs(folder, repair=False):
        print folder
        files = os.listdir(folder)
        DEG_db = {}
        for file in files:
            filename = folder + '/' + file
            if '.txt' in file and 'GE.' in file and '_vs_cellHarmony-Reference' not in file:
                header = True
                count = 0
                header = True
                file = file[:-4]
                file = string.split(file[3:], '_')[0]
                if '_vs_cellHarmony-Reference' in file:
                    file = 'global'
                    continue
                for line in open(filename, 'rU').xreadlines():
                    if header:
                        header = False
                    else:
                        data = cleanUpLine(line)
                        t = string.split(data, '\t')
                        GeneID, SystemCode, LogFold, rawp, adjp, Symbol, avg_g2, avg_g1 = t
                        rawp = float(rawp)
                        adjp = float(adjp)
                        if float(LogFold) > 0:
                            direction = 'positive'
                        else:
                            direction = 'negative'
                        if repair:
                            if float(LogFold) > 0:
                                fold = math.pow(2, float(LogFold))
                            else:
                                fold = -1 / math.pow(2, float(LogFold))
                            if abs(fold) > 0 and rawp < 0.05:
                                try:
                                    DEG_db[Symbol].append([file, direction])
                                except:
                                    DEG_db[Symbol] = [[file, direction]]

                        else:
                            try:
                                DEG_db[Symbol].append([file, direction])
                            except:
                                DEG_db[Symbol] = [[file, direction]]
        return DEG_db

    ref_DEGs = importCellHarmonyDEGs(reference_fold_dir)
    repaired_DEGs = importCellHarmonyDEGs(repair_dir1)
    repaired2_DEGs = importCellHarmonyDEGs(repair_dir2)
    
    ### Check to see if all three datasets agree
    total_repaired_genes ={}
    total_repaired_by_comp={}
    if mode == 'validate':
        for gene in ref_DEGs:
            if gene in repaired_DEGs and gene in repaired2_DEGs:
                for (file1,direction1) in ref_DEGs[gene]:
                    for (file2,direction2) in repaired_DEGs[gene]:
                        for (file3,direction3) in repaired2_DEGs[gene]: 
                            if direction1 == direction2 and file1 == file2 and direction3 == direction1 and file3 == file1:
                                try: total_repaired_genes[gene].append([direction1,file1,file2])
                                except: total_repaired_genes[gene] = [[direction1,file1,file2]]
                                try: total_repaired_by_comp[file1].append(gene)
                                except: total_repaired_by_comp[file1] = [gene]
                                
        files = os.listdir(reference_fold_dir)
        for file in files:
            filename = reference_fold_dir + '/' + file
            oef = export.ExportFile(reference_fold_dir + '-Filtered/'+file)
            if '.txt' in file and 'GE.' in file and '_vs_cellHarmony-Reference' not in file:
                header = True
                file = file[:-4]
                file = string.split(file[3:], '_')[0]
                for line in open(filename, 'rU').xreadlines():
                    if header:
                        header = False
                        oef.write(line)
                    else:
                        data = cleanUpLine(line)
                        t = string.split(data, '\t')
                        GeneID, SystemCode, LogFold, rawp, adjp, Symbol, avg_g2, avg_g1 = t
                        if file in total_repaired_by_comp:
                            if GeneID in total_repaired_by_comp[file] or Symbol in total_repaired_by_comp[file]:
                                oef.write(line)
                oef.close()
                          
    print len(total_repaired_genes)
    
    """
    for gene in total_repaired_genes:
        print gene

    def importCellHarmonyPseudoBulkFolds(filename):
        fold_db = {}
        header = True
        for line in open(filename, 'rU').xreadlines():
            data = cleanUpLine(line)
            t = string.split(data, '\t')
            if header:
                fold_db['header'] = t[1:]
                header = False
            else:
                uid = t[0]
                folds = t[1:]
                fold_db[uid] = folds

        return fold_db

    repaired_fold_db = importCellHarmonyPseudoBulkFolds(repair1_folds)
    repaired2_fold_db = importCellHarmonyPseudoBulkFolds(repair2_folds)
    import collections
    ordered_ref_degs = collections.OrderedDict()
    ordered_cluster_genes = collections.OrderedDict()
    repair_verified = collections.OrderedDict()
    repair2_verified = collections.OrderedDict()
    cluster_ordered_ref_db = collections.OrderedDict()
    header = True
    """
    header = True
    print organized_diff_ref
    eo1 = export.ExportFile(organized_diff_ref[:-4] + '-Filtered.txt')
    for line in open(organized_diff_ref, 'rU').xreadlines():
        data = cleanUpLine(line)
        t = string.split(data, '\t')
        if header:
            eo1.write(line)
            ref_header = t
            header = False
        else:
            cluster, geneID = string.split(t[0], ':')
            if geneID in total_repaired_genes:
                eo1.write(line)
    eo1.close()
    
    """
    repaired_verified = {}
    verified = {}
    for geneID, ref_cluster in ordered_ref_degs:
        for cluster, ref_direction in ref_DEGs[geneID]:
            if geneID in repaired_DEGs:
                for repair_cluster, repair_direction in repaired_DEGs[geneID]:
                    if repair_cluster == cluster and ref_direction != repair_direction and ('Neu' in repair_cluster or 'global' in repair_cluster):
                        try:
                            repair_verified[repair_cluster].append(geneID)
                        except:
                            repair_verified[repair_cluster] = [geneID]
                        else:
                            print geneID + '\t' + repair_direction + '\t' + repair_cluster + '\tR412X-HMZ'
                            try:
                                verified[geneID].append('R412X-HMZ')
                            except:
                                verified[geneID] = ['R412X-HMZ']

            if geneID in repaired2_DEGs:
                for repair_cluster, repair_direction in repaired2_DEGs[geneID]:
                    if repair_cluster == cluster and ref_direction != repair_direction and ('Neu' in cluster or 'global' in cluster):
                        try:
                            repair2_verified[repair_cluster].append(geneID)
                        except:
                            repair2_verified[repair_cluster] = [geneID]
                        else:
                            print geneID + '\t' + repair_direction + '\t' + repair_cluster + '\t' + 'R412X-Irf8'
                            try:
                                verified[geneID].append('R412X-Irf8')
                            except:
                                verified[geneID] = ['R412X-Irf8']

    for gene in verified:
        verified[gene] = unique.unique(verified[gene])

    eo1 = export.ExportFile(organized_diff_ref[:-4] + '-Repair-Sorted.txt')
    eo2 = export.ExportFile(organized_diff_ref[:-4] + '-Repaired-Only.txt')
    header = ref_header + repaired_fold_db['header'] + repaired2_fold_db['header']
    eo1.write(string.join(header, '\t') + '\n')
    eo2.write(string.join(header, '\t') + '\n')
    print len(ordered_ref_degs)
    print len(repaired_fold_db)
    print len(repaired2_fold_db)
    print len(repair_verified)
    print len(repair2_verified)
    print len(verified)
    print len(ordered_ref_degs)
    prior_cluster = None
    added_genes = []
    for geneID, cluster in ordered_ref_degs:
        try:
            folds = ordered_ref_degs[(geneID, cluster)] + repaired_fold_db[geneID] + repaired2_fold_db[geneID]
        except:
            print '...Error in identifying match UID for:', geneID
            added_genes.append(geneID)
            continue
        else:
            if geneID not in verified:
                eo1.write(string.join(folds, '\t') + '\n')
            elif len(verified[geneID]) > 1:
                added_genes.append(geneID)
            elif 'R412X-HMZ' in verified[geneID]:
                added_genes.append(geneID)
            else:
                eo2.write(string.join(folds, '\t') + '\n')
                added_genes.append(geneID)

    eo1.close()
    eo2.close()
    """
if __name__ == '__main__':
    repair1_folds = None
    repair2_folds = None
    organized_diff_ref = '/Users/saljh8/Downloads/Transfer/cellHarmony_10-control-clusters/cellHarmony-4-samples/OrganizedDifferentials.txt'
    repair_dir1 = '/Users/saljh8/Downloads/Transfer/cellHarmony-A3/DifferentialExpression_Fold_1.1_rawp_0.05'
    repair_dir2 = '/Users/saljh8/Downloads/Transfer/cellHarmony-A3/DifferentialExpression_Fold_1.1_rawp_0.05'
    reference_fold_dir = '/Users/saljh8/Downloads/Transfer/cellHarmony_10-control-clusters/cellHarmony-4-samples/DifferentialExpression_Fold_1.2_adjp_0.05'

    print 'comparing cellHarmony outputs'
    rankExpressionRescueFromCellHarmony(organized_diff_ref, repair1_folds, repair2_folds, reference_fold_dir, repair_dir1, repair_dir2);sys.exit()
  