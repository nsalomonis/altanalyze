import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import scipy, numpy
import statistics
from visualization_scripts import clustering

def evaluateMultiLinRegulatoryStructure(all_genes_TPM,MarkerFinder,SignatureGenes,state=None,query=None):
    all_indexes, group_index, expressionData = loopThroughEachState(all_genes_TPM)

    if state!=None:
        states = [state] ### For example, we only want to look in annotated Multi-Lin's
    else:
        states = group_index
    state_scores=[]
    
    for state in states:
        print '\n',state, 'running now.'
        score = evaluateStateRegulatoryStructure(expressionData,all_indexes,group_index,MarkerFinder,SignatureGenes,state,query=query)
        state_scores.append([score,state])
        print state, score
    state_scores.sort()
    state_scores.reverse()
    print state_scores
    
def loopThroughEachState(all_genes_TPM):
    ### Import all genes with TPM values for all cells
    matrix, column_header, row_header, dataset_name, group_db = clustering.importData(all_genes_TPM)
    group_index={}
    all_indexes=[]
    for sampleName in group_db:
        ICGS_state = group_db[sampleName][0]
        try: group_index[ICGS_state].append(column_header.index(sampleName))
        except Exception: group_index[ICGS_state] = [column_header.index(sampleName)]
        all_indexes.append(column_header.index(sampleName))
    for ICGS_state in group_index:
        group_index[ICGS_state].sort()
    all_indexes.sort()
    expressionData = matrix, column_header, row_header, dataset_name, group_db
    return all_indexes, group_index, expressionData
    
def evaluateStateRegulatoryStructure(expressionData, all_indexes,group_index,MarkerFinder,SignatureGenes,state,query=None):
    """Predict multi-lineage cells and their associated coincident lineage-defining TFs"""

    useProbablityOfExpression=False
    ICGS_State_as_Row = False

    matrix, column_header, row_header, dataset_name, group_db = expressionData
    
    def importGeneLists(fn):
        genes={}
        for line in open(fn,'rU').xreadlines():
            data = clustering.cleanUpLine(line)
            gene,cluster = string.split(data,'\t')[0:2]
            genes[gene]=cluster
        return genes
    
    def importMarkerFinderHits(fn):
        genes={}
        genes_to_symbol={}
        ICGS_State_ranked={}
        skip=True
        for line in open(fn,'rU').xreadlines():
            data = clustering.cleanUpLine(line)
            if skip: skip=False
            else:
                try:
                    gene,symbol,rho,ICGS_State = string.split(data,'\t')
                except Exception:
                    gene,symbol,rho,rho_p,ICGS_State = string.split(data,'\t')
                genes_to_symbol[gene]=symbol
                #if ICGS_State!=state and float(rho)>0.0:
                if float(rho)>0.3:
                    try: ICGS_State_ranked[ICGS_State].append([float(rho),gene,symbol])
                    except Exception: ICGS_State_ranked[ICGS_State] = [[float(rho),gene,symbol]]

        for ICGS_State in ICGS_State_ranked:
            ICGS_State_ranked[ICGS_State].sort()
            ICGS_State_ranked[ICGS_State].reverse()
            #print ICGS_State, ICGS_State_ranked[ICGS_State][:50]
            for (rho,gene,symbol) in ICGS_State_ranked[ICGS_State][:50]:
                genes[gene]=rho,ICGS_State ### Retain all population specific genes (lax)
                genes[symbol]=rho,ICGS_State
                
        return genes, genes_to_symbol
    
    def importQueryDataset(fn):
        matrix, column_header, row_header, dataset_name, group_db = clustering.importData(fn)
        return matrix, column_header, row_header, dataset_name, group_db
    
    signatureGenes = importGeneLists(SignatureGenes)
    markerFinderGenes, genes_to_symbol = importMarkerFinderHits(MarkerFinder)
    #print len(signatureGenes),len(markerFinderGenes)

    ### Determine for each gene, its population frequency per cell state
    index=0
    expressedGenesPerState={}
    stateAssociatedMarkers={}
    
    def freqCutoff(x,cutoff):
        if x>cutoff: return 1 ### minimum expression cutoff
        else: return 0

    for row in matrix:
        ICGS_state_gene_frq={}
        gene = row_header[index]
        for ICGS_state in group_index:
            state_values = map(lambda i: row[i],group_index[ICGS_state])    
            def freqCheck(x):
                if x>1: return 1 ### minimum expression cutoff
                else: return 0
                
            expStateCells = sum(map(lambda x: freqCheck(x),state_values))
            statePercentage = (float(expStateCells)/len(group_index[ICGS_state]))
            ICGS_state_gene_frq[ICGS_state] = statePercentage
        
 
        datasets_values = map(lambda i: row[i],all_indexes)    
        all_cells_frq = sum(map(lambda x: freqCheck(x),datasets_values))/(len(datasets_values)*1.0)
        all_states_frq = map(lambda x: ICGS_state_gene_frq[x],ICGS_state_gene_frq)
        all_states_frq.sort() ### frequencies of all non-multilin states
        states_expressed = sum(map(lambda x: freqCutoff(x,0.5),all_states_frq))/(len(all_states_frq)*1.0)
        
        for State in ICGS_state_gene_frq:
            state_frq = ICGS_state_gene_frq[State]
            rank = all_states_frq.index(state_frq)
            if state_frq > 0.25 and rank>0: #and states_expressed<0.75 #and all_cells_frq>0.75
                if 'Rik' not in gene and 'Gm' not in gene and '-' not in gene:
                    if gene in markerFinderGenes:# and gene in markerFinderGenes:
                        if ICGS_State_as_Row:
                            ICGS_State = signatureGenes[gene]
                        if gene in markerFinderGenes:
                            if ICGS_State_as_Row == False:
                                rho, ICGS_State = markerFinderGenes[gene]
                            else:
                                rho, ICGS_Cell_State = markerFinderGenes[gene] #ICGS_Cell_State
                            #try: gene = genes_to_symbol[gene]
                            #except: gene = gene
                            score = int(rho*100*state_frq)*(float(rank)/len(all_states_frq))
                            try: expressedGenesPerState[ICGS_State].append((score,gene))
                            except Exception: expressedGenesPerState[ICGS_State]=[(score,gene)] #(rank*multilin_frq)
                            try: stateAssociatedMarkers[gene,ICGS_State].append(State)
                            except Exception: stateAssociatedMarkers[gene,ICGS_State] = [State]
        index+=1

    if query!=None:
        matrix, column_header, row_header, dataset_name, group_db = importQueryDataset(query)
    
    markers_to_exclude=[]
    expressedGenesPerState2={}
    for (gene,ICGS_State) in stateAssociatedMarkers:
        if len(stateAssociatedMarkers[(gene,ICGS_State)])<2: # or len(stateAssociatedMarkers[(gene,ICGS_State)])>len(ICGS_state_gene_frq)/2.0:
            markers_to_exclude.append(gene)
        else:
            print ICGS_State, gene, stateAssociatedMarkers[(gene,ICGS_State)]
    for ICGS_State in expressedGenesPerState:
        for (score,gene) in expressedGenesPerState[ICGS_State]:
            if gene not in markers_to_exclude:
                try: expressedGenesPerState2[ICGS_State].append((score,gene))
                except Exception: expressedGenesPerState2[ICGS_State] = [(score,gene)]
    expressedGenesPerState = expressedGenesPerState2
        
    createPseudoCell=True
    ### The expressedGenesPerState defines genes and modules co-expressed in the multi-Lin
    ### Next, find the cells that are most frequent in mulitple states
    representativeMarkers={}
    for ICGS_State in expressedGenesPerState:
        expressedGenesPerState[ICGS_State].sort()
        expressedGenesPerState[ICGS_State].reverse()
        if '1Multi' not in ICGS_State:
            markers = expressedGenesPerState[ICGS_State]#[:5]
            markers_unique = list(set(map(lambda x: x[1],list(markers))))
            print ICGS_State,":",string.join(markers_unique,', ')
            if createPseudoCell:
                for gene in markers:
                    def getBinary(x):
                        if x>1: return 1
                        else: return 0
                    if gene[1] in row_header: ### Only for query datasets
                        row_index = row_header.index(gene[1])
                        if useProbablityOfExpression:
                            pvalues = calculateGeneExpressProbilities(matrix[row_index]) ### probability of expression
                            values = pvalues
                        else:
                            binaryValues = map(lambda x: getBinary(x), matrix[row_index])
                            values = binaryValues
                            #values = matrix[row_index]
                        #if gene[1]=='S100a8': print binaryValues;sys.exit()
                        try: representativeMarkers[ICGS_State].append(values)
                        except Exception: representativeMarkers[ICGS_State] = [values]    
            else:
                representativeMarkers[ICGS_State]=markers[0][-1]
        #int(len(markers)*.25)>5:
        #print ICGS_State, markers
    #sys.exit()

    for ICGS_State in representativeMarkers:
        if createPseudoCell:
            signature_values = representativeMarkers[ICGS_State]
            if useProbablityOfExpression:
                signature_values = [numpy.sum(value) for value in zip(*signature_values)]
            else:
                signature_values = [float(numpy.mean(value)) for value in zip(*signature_values)]
            representativeMarkers[ICGS_State] = signature_values
        else:
            gene = representativeMarkers[ICGS_State]
            row_index = row_header.index(gene)
            gene_values = matrix[row_index]
            representativeMarkers[ICGS_State] = gene_values

    ### Determine for each gene, its population frequency per cell state
    expressedStatesPerCell={}
    multilin_probability={}
    import export
    print 'Writing results matrix to:',MarkerFinder[:-4]+'-cellStateScores.txt'
    eo = export.ExportFile(MarkerFinder[:-4]+'-cellStateScores.txt')
    eo.write(string.join(['UID']+column_header,'\t')+'\n')
    print 'a'
    print len(representativeMarkers)
    for ICGS_State in representativeMarkers:
        gene_values = representativeMarkers[ICGS_State]
        index=0
        scoreMatrix=[]
        HitsCount=0
        for cell in column_header:
            value = gene_values[index]
            """
            expressedLiklihood = '0'
            if (value<0.05 and useProbablityOfExpression==True) or (value==1 and useProbablityOfExpression==False):
                try: expressedStatesPerCell[cell].append(ICGS_State)
                except Exception: expressedStatesPerCell[cell] = [ICGS_State]
                expressedLiklihood = '1'
                HitsCount+=1
            if useProbablityOfExpression:
                try: multilin_probability[cell].append(value)
                except Exception: multilin_probability[cell] = [value]
            """
            index+=1
            
            HitsCount+=1
            scoreMatrix.append(str(value))
        if HitsCount>1:
            #print ICGS_State,HitsCount
            eo.write(string.join([ICGS_State]+scoreMatrix,'\t')+'\n')
    eo.close()
    sys.exit()
    def multiply(values):
        p = 1
        for i in values:
            if i>0:
                p = p*i
            else:
                p = p*1.e-16
        return p
    cell_mutlilin_ranking=[]   
    for cell in expressedStatesPerCell:
        #if 'Multi-Lin:Gmp.R3.10' in cell: sys.exit()
        if useProbablityOfExpression:
            p = numpy.mean(multilin_probability[cell]) ### mean state probability
        lineageCount = expressedStatesPerCell[cell]
        if useProbablityOfExpression:
            cell_mutlilin_ranking.append((p,len(lineageCount),cell))
        else:
            cell_mutlilin_ranking.append((len(lineageCount),cell))
    cell_mutlilin_ranking.sort()
    if useProbablityOfExpression == False:
        cell_mutlilin_ranking.reverse()
        
    scores = []
    state_scores={}
    cellsPerState={} ### Denominator for z-score analysis
    for cell in cell_mutlilin_ranking:
        score = cell[0]
        scores.append(score)
        cell_state = string.split(cell[-1],':')[0]
        try: cellsPerState[cell_state]+=1
        except Exception: cellsPerState[cell_state]=1
        try: state_scores[cell_state].append(float(score))
        except Exception: state_scores[cell_state] = [float(score)]

    scoreMean = numpy.mean(scores)
    scoreSD = numpy.std(scores)
    oneSD = scoreMean+scoreSD
    twoSD = scoreMean+scoreSD+scoreSD
    oneStandDeviationAway={}
    twoStandDeviationsAway={}
    oneStandDeviationAwayTotal=0
    twoStandDeviationsAwayTotal=0
    
    print 'Mean:',scoreMean
    print 'STDev:',scoreSD
    state_scores2=[]
    for cell_state in state_scores:
        state_scores2.append((numpy.mean(state_scores[cell_state]),cell_state))
    i=0
    for cell in cell_mutlilin_ranking:
        score,cellName = cell
        CellState,CellName = string.split(cellName,':')
        if score>=oneSD:
            try: oneStandDeviationAway[CellState]+=1
            except Exception: oneStandDeviationAway[CellState]=1
            oneStandDeviationAwayTotal+=1
            if score>=twoSD:
                try: twoStandDeviationsAway[CellState]+=1
                except Exception: twoStandDeviationsAway[CellState]=1
                twoStandDeviationsAwayTotal+=1
            print cell, string.join(expressedStatesPerCell[cell[-1]],'|')
        i+=1
    state_scores2
    state_scores2.sort()
    state_scores2.reverse()
    
    twoStandDeviationsAway = oneStandDeviationAway
    twoStandDeviationsAwayTotal = oneStandDeviationAwayTotal
    
    print '\n\n'
    import statistics
    zscores = []
    for CellState in twoStandDeviationsAway:
        #print CellState
        highMetaScoreCells = twoStandDeviationsAway[CellState]
        totalCellsPerState = cellsPerState[CellState]
        r = highMetaScoreCells
        n = twoStandDeviationsAwayTotal
        R = totalCellsPerState
        N = len(column_header)
        z = statistics.zscore(r,n,N,R)
        scores = [z, CellState,statistics.p_value(z)]
        zscores.append(scores)
        
    zscores.sort()
    zscores.reverse()
    for scores in zscores:
        scores = string.join(map(str,scores),'\t')
        print scores
        
    """
    for i in state_scores2:
        print str(i[0])+'\t'+str(i[1])"""
    sys.exit()
    return numpy.mean(state_scores)

def calculateGeneExpressProbilities(values, useZ=False):
    ### First calculate z-scores - scipy.stats.mstats.zscore for the entire matrix
    avg = numpy.mean(values)
    std = numpy.std(values)
    if std ==0:
        std = 0.1
    if useZ:
        values = map(lambda x: (x-avg)/std,values)
    else:
        values = map(lambda x: x*2,values)
    p_values = 1 - scipy.special.ndtr(values)
    return p_values
    
if __name__ == '__main__':
    #query_dataset = '/Users/saljh8/Desktop/Old Mac/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionInput/exp.GSE81682_HTSeq-cellHarmony-filtered.txt'
    all_tpm = '/Users/saljh8/Desktop/Old Mac/Desktop/demo/BoneMarrow/ExpressionInput/exp.BoneMarrow-scRNASeq.txt'
    markerfinder = '/Users/saljh8/Desktop/Old Mac/Desktop/demo/BoneMarrow/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt'
    signature_genes = '/Users/saljh8/Desktop/Old Mac/Desktop/Grimes/KashishNormalization/test/Panorama.txt'
    state = 'DC'

    all_tpm = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/ExpressionInput/exp.Guide3-cellHarmony-revised.txt'
    #markerfinder = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/updated.8.29.17/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt'
    signature_genes = '/Users/saljh8/Desktop/Old Mac/Desktop/Grimes/KashishNormalization/test/Panorama.txt'
    
    query_dataset = None
    query_dataset = '/Users/saljh8/Desktop/dataAnalysis/Collaborative/Grimes/All-Fluidigm/exp.NaturePan-PreGM-CD150-.txt'
    query_dataset = None
    """
    #all_tpm = '/Users/saljh8/Desktop/Old Mac/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionInput/MultiLin/exp.Gottgens_HarmonizeReference.txt'
    all_tpm = '/Users/saljh8/Desktop/Old Mac/Desktop/demo/Mm_Gottgens_3k-scRNASeq/MultiLin/ExpressionInput/exp.Gottgens_HarmonizeReference.txt'
    #signature_genes = '/Users/saljh8/Desktop/Old Mac/Desktop/demo/Mm_Gottgens_3k-scRNASeq/MultiLin/ExpressionInput/Gottgens_HarmonizeReference.txt'
    signature_genes = '/Users/saljh8/Desktop/Old Mac/Desktop/demo/Mm_Gottgens_3k-scRNASeq/MultiLin/Gottgens_HarmonizeReference.txt'
    #markerfinder = '/Users/saljh8/Desktop/Old Mac/Desktop/demo/Mm_Gottgens_3k-scRNASeq/ExpressionOutput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt'
    markerfinder = '/Users/saljh8/Desktop/Old Mac/Desktop/demo/Mm_Gottgens_3k-scRNASeq/MultiLin/ExpressionInput/MarkerFinder/AllGenes_correlations-ReplicateBased.txt'
    state = 'Eryth_Multi-Lin'
    """
    state = None
    import getopt
    options, remainder = getopt.getopt(sys.argv[1:],'', ['q=','expdir=','m=','ICGS=','state='])
    #print sys.argv[1:]
    for opt, arg in options:
        if opt == '--q': query_dataset=arg
        elif opt == '--expdir': all_tpm=arg
        elif opt == '--m': markerfinder=arg
        elif opt == '--ICGS': signature_genes=arg
        elif opt == '--state': state=arg
            
    #state = None
    #evaluateMultiLinRegulatoryStructure(all_tpm,markerfinder,signature_genes,state);sys.exit()
    

    evaluateMultiLinRegulatoryStructure(all_tpm,markerfinder,signature_genes,state,query = query_dataset);sys.exit()