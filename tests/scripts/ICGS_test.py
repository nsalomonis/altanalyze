"""
    Usuage: runICGStest(testType = string,inputData = string)
        
    Determines whether ICGS produces near identical log files to referenced AltAnalyze versions
    across operating systems, software versions and machines. Can be run using different options
    to optimize speed and comprehension for different input types. Compared to pre-computed logs.
    
    Parameters:
        testType = "fast", "complete"
        inputData = "text", "FASTQ", "BAM", "chromium", "excel"
    Returns:
        All ICGS results (ICGS directory), t-SNE plots (DataPlots directory), MarkerFinder
        (DataPlots/MarkerFinder directory) and log file (output root directory).
    
    Commandline example
    # ICGS on a tab-delimited text file
    python AltAnalyze.py --platform RNASeq --species Mm --column_method hopach --column_metric cosine --rho 0.4
    --ExpressionCutoff 1 --FoldDiff 4 --SamplesDiffering 4 --restrictBy protein_coding
    --excludeCellCycle conservative --removeOutliers no --row_method ward
    --expdir tests/demo_data/Fluidigim_TPM/input/BoneMarrow-scRNASeq.txt
    --output tests/demo_data/Fluidigim_TPM/output/
    --runICGS yes --genes "Gfi1 Irf8 Vwf Mmp9" --expname BoneMarrow-scRNASeq
    
    # ICGS on genome aligned BAM files
    python AltAnalyze.py --platform RNASeq --species Hs --column_method hopach --column_metric cosine --rho 0.4 
    --SamplesDiffering 3 --restrictBy protein_coding --excludeCellCycle no --removeOutliers no --row_method hopach --FoldDiff 4
    --bedDir tests/demo_data/BAM/input/ --expname test --ExpressionCutoff 1 
    --output /tests/demo_data/BAM/input/ --runICGS yes
    
    # ICGS on the 10X Chromium Platform
    python AltAnalyze.py --platform RNASeq --species Mm --column_method hopach --column_metric cosine --rho 0.3 --ExpressionCutoff 1
    --FoldDiff 4 --SamplesDiffering 4 --restrictBy None --excludeCellCycle no --removeOutliers yes --row_method ward
    --ChromiumSparseMatrix tests/demo_data/10X/input/mm10/matrix.mtx
    --output tests/demo_data/10X/output/ --runICGS yes --expname test
"""
import os,sys,string,inspect
### Import python modules from an upstream directory
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
parentdir = os.path.dirname(parentdir)
sys.path.insert(0,parentdir)

import subprocess
import export
import unique

def runICGStest(testType = "complete",inputData = "BAM"):
    
    import UI
    species_names = UI.getSpeciesInfo()
    
    additional=None
    genes = ''
    genes_option = "--genes"
    cluster_method = "ward"
    rho = "0.4"
    removeOutliers = "no"
    
    if inputData == "text":
        if 'Mm' not in species_names:
            return 'WARNING!!! Species Mm database not installed.'
        inputDataType = "Fluidigim_TPM"
        species = "Mm"
        SamplesDiffering = "4"
        excludeCellCycle = "conservative"
        restrictBy = "None"
        input_path = unique.filepath("AltDatabase/demo_data/"+inputDataType+"/input/BoneMarrow-scRNASeq.txt")
        expname = "BoneMarrow-scRNASeq"
        expdir = "--expdir"
        removeOutliers = "no"
    elif inputData == "BAM":
        if 'Mm' not in species_names:
            return 'WARNING!!! Species Mm database not installed.'
        species = "Mm"
        inputDataType = inputData
        SamplesDiffering = "3"
        excludeCellCycle = "no"
        restrictBy = "None"
        expdir = "--bedDir"
        expname = "test"
        rho = "0.6"
        input_path = unique.filepath("AltDatabase/demo_data/"+inputDataType+"/input")
    elif inputData == 'FASTQ':
        if 'Mm' not in species_names:
            return 'WARNING!!! Species Mm database not installed.'
        species = "Mm"
        inputDataType = inputData
        SamplesDiffering = "3"
        excludeCellCycle = "no"
        restrictBy = "None"
        rho = "0.6"
        expdir = "--fastq_dir"
        expname = "test"
        input_path = unique.filepath("AltDatabase/demo_data/"+inputDataType+"/input")
        custom_FASTA_path = unique.filepath("AltDatabase/demo_data/FASTA/Homo_sapiens.GRCh37.72.cdna.all.filtered.fa")
        additional = ["--runKallisto","True","--customFASTA",custom_FASTA_path]
        additional = []
    elif inputData == '10X':
        if 'Mm' not in species_names:
            return 'WARNING!!! Species Mm database not installed.'
        species = "Mm"
        inputDataType = inputData
        SamplesDiffering = "4"
        rho = "0.3"
        excludeCellCycle = "conservative"
        removeOutliers = "yes"
        restrictBy = "protein_coding"
        expdir = "--ChromiumSparseMatrix"
        expname = "test"
        input_path = unique.filepath("AltDatabase/demo_data/10X/input/mm10/matrix.mtx")
        
    if testType=="complete":
        cluster_method = "hopach"
    elif inputData == "text":
        ### Optionally, peed up the analysis by restricting to genes correlated to these guides
        genes = "Gfi1 Irf8 Vwf Mmp9"
        
    outputPath = unique.filepath("AltDatabase/demo_data/"+inputDataType+"/output/")
    input_root_dir = unique.filepath("AltDatabase/demo_data/"+inputDataType+"/input/")

    ### Create the output directory path if not created
    try: os.mkdir(outputPath)
    except Exception:
        ### Remove this directory
        print 'Removing previous output test directory...'
        export.deleteFolder(outputPath)
        try: os.mkdir(outputPath)
        except Exception: pass
    
    ### Commands replicate Olsson et al. 2016 Nature
    commands = ["python", "AltAnalyze.py", "--platform","RNASeq", "--species", species,
                "--column_method", cluster_method, "--column_metric", "cosine", "--rho", rho,
                "--ExpressionCutoff", "1", "--FoldDiff", "4", "--SamplesDiffering",SamplesDiffering,
                "--restrictBy", restrictBy, "--excludeCellCycle", excludeCellCycle, "--removeOutliers",
                removeOutliers, "--row_method", cluster_method, expdir, input_path, "--output", outputPath,
                "--runICGS", "yes", genes_option, genes, "--expname",expname]

    if additional != None: ### additional optional parameters
        commands+=additional
    
    print string.join(commands,' ')
    retcode = subprocess.call(commands)
    #sys.exit()
    reference_log, ref_exceptions = readLogFile(input_root_dir)
    test_log, test_exceptions = readLogFile(outputPath)
    overlap = list(set(reference_log).intersection(test_log))

    matching_content = float(len(overlap))/len(reference_log)
    percent_matching = int(matching_content*100)

    if percent_matching<90: status = "\nWARNING!!! "
    else: status = "\nSUCCESSFUL OUTPUT!!! "
    status+= str(percent_matching)+'% of reported ICGS outputs match to reference outputs. '
    status+= str(len(test_exceptions))+' errors encountered.'
    print status
    return status

def readLogFile(root_dir):
    log_output=[]
    exceptions=[]
    for file in os.listdir(root_dir):
        if 'AltAnalyze_report' in file and '.log' in file:
            log_file = root_dir+'/'+file
            log_contents = open(log_file, "rU")
            for line in log_contents:
                line = line.rstrip()
                if 'Pruning' not in line and 'seconds' not in line and 'AltAnalyze' not in line:
                    log_output.append(line)
                if 'Error' in line or 'xception' in line:
                    exceptions.append(line)
    return list(set(log_output)), exceptions

if __name__ == '__main__':
    runICGStest()