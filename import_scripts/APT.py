###APT
#Copyright 2012 J. David Gladstone Institutes, San Francisco California
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

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies
import os.path
import unique
import UI
import export; reload(export)
import time
import shutil
import traceback

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def osfilepath(filename):
    fn = filepath(filename)
    fn = string.replace(fn,'\\','/')
    return fn

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def annotateMetaProbesetGenes(summary_exp_file, expression_file, metaprobeset_file, species):
    metaprobeset_cv_file = string.replace(metaprobeset_file,species+'_',species+'_Conversion_')
    metaprobeset_cv_file = string.replace(metaprobeset_cv_file,'.mps','.txt')

    fn=filepath(metaprobeset_cv_file); uid_db={}
    for line in open(fn,'rU').xreadlines():
        data = UI.cleanUpLine(line)
        uid,ens_gene = string.split(data,'\t')
        uid_db[uid] = ens_gene

    export_data = export.ExportFile(expression_file)        
    fn=filepath(summary_exp_file); x=0
    for line in open(fn,'rU').xreadlines():
        if line[0] == '#': null=[]
        elif x == 0: export_data.write(line); x+=1
        else:
            data = cleanUpLine(line)
            t = string.split(data,'\t')
            uid = t[0]; ens_gene = uid_db[uid]
            export_data.write(string.join([ens_gene]+t[1:],'\t')+'\n')
    export_data.close()

def reformatResidualFile(residual_exp_file,residual_destination_file):
    ### Re-write the residuals file so it has a single combined unique ID (arbitrary gene ID + probe ID)
    print 'Re-formatting and moving the calculated residuals file...'
    export_data = export.ExportFile(residual_destination_file)   
    fn=filepath(residual_exp_file); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0 and data[0]=='#': null=[]
        elif x == 0:
            x+=1; t = string.split(data,'\t')
            export_data.write(string.join(['UID']+t[5:],'\t')+'\n')
        else:
            t = string.split(data,'\t')
            uid = t[0]+'-'+t[2] ### arbitrary numeric gene ID + probes ID
            export_data.write(string.join([uid]+t[5:],'\t')+'\n')
    export_data.close()
    os.remove(residual_exp_file)

def verifyFileLength(filename):
    count = 0
    try:
        fn=filepath(filename)
        for line in open(fn,'rU').xreadlines():
            count+=1
            if count>9: break
    except Exception: null=[]
    return count

def APTDebugger(output_dir):
    fatal_error = ''
    fn = filepath(output_dir+'/apt-probeset-summarize.log')
    for line in open(fn,'rU').xreadlines():
        if 'FATAL ERROR:' in line:
            fatal_error = line
    return fatal_error
        
def probesetSummarize(exp_file_location_db,analyze_metaprobesets,probeset_type,species,root):
    for dataset in exp_file_location_db: ### Instance of the Class ExpressionFileLocationData
        fl = exp_file_location_db[dataset]
        apt_dir =fl.APTLocation()
        array_type=fl.ArrayType()  
        pgf_file=fl.InputCDFFile()
        clf_file=fl.CLFFile()
        bgp_file=fl.BGPFile()
        xhyb_remove = fl.XHybRemoval()
        cel_dir=fl.CELFileDir() + '/cel_files.txt'
        expression_file = fl.ExpFile()
        stats_file = fl.StatsFile()
        output_dir = fl.OutputDir() + '/APT-output'
        cache_dir = output_dir + '/apt-probeset-summarize-cache'
        architecture = fl.Architecture() ### May over-ride the real architecture if a failure occurs
        get_probe_level_results = 'yes'
        
        if get_probe_level_results == 'yes': export_features = 'yes'
        if xhyb_remove == 'yes' and (array_type == 'gene' or array_type == 'junction'): xhyb_remove = 'no' ### This is set when the user mistakenly selects exon array, initially
        if analyze_metaprobesets == 'yes':
            export_features = 'true'
            metaprobeset_file = filepath('AltDatabase/'+species+'/'+array_type+'/'+species+'_'+array_type+'_'+probeset_type+'.mps')
            count = verifyFileLength(metaprobeset_file)
            if count<2:
                from build_scripts import ExonArray
                ExonArray.exportMetaProbesets(array_type,species) ### Export metaprobesets for this build
        import subprocess; import platform
        print 'Processor architecture set =',architecture,platform.machine()
        if '/bin' in apt_dir: apt_file = apt_dir +'/apt-probeset-summarize' ### if the user selects an APT directory
        elif os.name == 'nt':
            if '32bit' in architecture: apt_file = apt_dir + '/PC/32bit/apt-probeset-summarize'; plat = 'Windows'
            elif '64bit' in architecture: apt_file = apt_dir + '/PC/64bit/apt-probeset-summarize'; plat = 'Windows'
        elif 'darwin' in sys.platform: apt_file = apt_dir + '/Mac/apt-probeset-summarize'; plat = 'MacOSX'
        elif 'linux' in sys.platform:
            if '32bit' in platform.architecture(): apt_file = apt_dir + '/Linux/32bit/apt-probeset-summarize'; plat = 'linux32bit'
            elif '64bit' in platform.architecture(): apt_file = apt_dir + '/Linux/64bit/apt-probeset-summarize'; plat = 'linux64bit'
        apt_file = filepath(apt_file)
        apt_extract_file = string.replace(apt_file,'probeset-summarize','cel-extract')
        
        #print 'AltAnalyze has choosen APT for',plat
        print "Beginning probeset summarization of input CEL files with Affymetrix Power Tools (APT)..."
        if 'cdf' in pgf_file or 'CDF' in pgf_file:
            if xhyb_remove == 'yes' and array_type == 'AltMouse':
                kill_list_dir = osfilepath('AltDatabase/'+species+'/AltMouse/'+species+'_probes_to_remove.txt')
            else: kill_list_dir = osfilepath('AltDatabase/affymetrix/APT/probes_to_remove.txt')
            
            try:
                ### Below code attempts to calculate probe-level summarys and absent/present p-values
                ### for 3'arrays (may fail for arrays with missing missmatch probes - AltMouse)
                cdf_file = pgf_file; algorithm = 'rma'
                retcode = subprocess.call([
                    apt_file, "-d", cdf_file, "--kill-list", kill_list_dir, "-a", algorithm, "-o", output_dir,
                    "--cel-files", cel_dir, "-a", "pm-mm,mas5-detect.calls=1.pairs=1"])
                try:
                    extract_retcode = subprocess.call([
                        apt_extract_file, "-d", cdf_file, "--pm-with-mm-only", "-o", output_dir+'/probe.summary.txt',
                        "--cel-files", cel_dir, "-a"]) ### "quant-norm,pm-gcbg", "--report-background" -requires a BGP file
                except Exception,e:
                    #print traceback.format_exc()
                    retcode = False ### On some system there is a no file found error, even when the analysis completes correctly
                if retcode: status = 'failed'
                else:
                    status = 'run'
                    summary_exp_file = output_dir+'/'+algorithm+'.summary.txt'
                    export.customFileCopy(summary_exp_file, expression_file) ### Removes the # containing lines
                    #shutil.copyfile(summary_exp_file, expression_file)
                    os.remove(summary_exp_file)
                    
                    summary_stats_file = output_dir+'/pm-mm.mas5-detect.summary.txt'
                    try: shutil.copyfile(summary_stats_file, stats_file)
                    except Exception: None ### Occurs if dabg export failed
                    os.remove(summary_stats_file)
            except Exception:
                #print traceback.format_exc()
                try:
                    cdf_file = pgf_file; algorithm = 'rma'; pval = 'dabg'
                    retcode = subprocess.call([
                        apt_file, "-d", cdf_file, "--kill-list", kill_list_dir, "-a", algorithm, "-o", output_dir, "--cel-files", cel_dir]) # "-a", pval,
                    if retcode: status = 'failed'
                    else:
                        status = 'run'
                        summary_exp_file = output_dir+'/'+algorithm+'.summary.txt'
                        export.customFileCopy(summary_exp_file, expression_file) ### Removes the # containing lines
                        #shutil.copyfile(summary_exp_file, expression_file)
                        os.remove(summary_exp_file)
                except NameError:
                    status = 'failed'
                    #print traceback.format_exc()
        else:
            if xhyb_remove == 'yes':
                kill_list_dir = osfilepath('AltDatabase/'+species+'/exon/'+species+'_probes_to_remove.txt')
            else: kill_list_dir = osfilepath('AltDatabase/affymetrix/APT/probes_to_remove.txt')
            if 'Glue' in pgf_file:
                kill_list_dir = string.replace(pgf_file,'pgf','kil') ### Needed to run DABG without crashing
                #psr_dir = string.replace(pgf_file,'pgf','PSR.ps') ### used with -s
            try:
                algorithm = 'rma-sketch'; pval = 'dabg'
                if analyze_metaprobesets != 'yes':
                    retcode = subprocess.call([
                    apt_file, "-p", pgf_file, "-c", clf_file, "-b", bgp_file, "--kill-list", kill_list_dir,
                    "-a", algorithm, "-a", pval, "-o", output_dir, "--cel-files", cel_dir]) # "--chip-type", "hjay", "--chip-type", "HJAY" http://www.openbioinformatics.org/penncnv/penncnv_tutorial_affy_gw6.html
                    if retcode:
                        summary_exp_file = output_dir+'/'+pval+'.summary.txt'
                        try: os.remove(summary_exp_file)
                        except Exception: null=[] ### Occurs if dabg export failed
                        fatal_error = APTDebugger(output_dir)
                        if len(fatal_error)>0:
                            print fatal_error
                            print 'Skipping DABG p-value calculation to resolve (Bad library files -> contact Affymetrix support)'
                            retcode = subprocess.call([
                            apt_file, "-p", pgf_file, "-c", clf_file, "-b", bgp_file, "--kill-list", kill_list_dir,
                            "-a", algorithm, "-o", output_dir, "--cel-files", cel_dir]) ### Exclude DABG p-value - known issue for Glue junction array
                        else: bad_exit
                else:
                    retcode = subprocess.call([
                    apt_file, "-p", pgf_file, "-c", clf_file, "-b", bgp_file, "--kill-list", kill_list_dir, "-m", metaprobeset_file,
                    "-a", algorithm, "-a", pval, "-o", output_dir, "--cel-files", cel_dir, "--feat-details", export_features])
                    if retcode:
                        summary_exp_file = output_dir+'/'+pval+'.summary.txt'
                        try: os.remove(summary_exp_file)
                        except Exception: null=[] ### Occurs if dabg export failed
                        fatal_error = APTDebugger(output_dir)
                        if len(fatal_error)>0:
                            print fatal_error
                            print 'Skipping DABG p-value calculation to resolve (Bad library files -> contact Affymetrix support)'
                            retcode = subprocess.call([
                            apt_file, "-p", pgf_file, "-c", clf_file, "-b", bgp_file, "--kill-list", kill_list_dir, "-m", metaprobeset_file,
                            "-a", algorithm, "-o", output_dir, "--cel-files", cel_dir, "--feat-details", export_features]) ### Exclude DABG p-value - known issue for Glue junction array
                        else: bad_exit
                if retcode: status = 'failed'
                else:
                    status = 'run'
                    summary_exp_file = output_dir+'/'+algorithm+'.summary.txt'
                    #if analyze_metaprobesets == 'yes': annotateMetaProbesetGenes(summary_exp_file, expression_file, metaprobeset_file, species)
                    export.customFileCopy(summary_exp_file, expression_file) ### Removes the # containing lines
                    #shutil.copyfile(summary_exp_file, expression_file)
                    os.remove(summary_exp_file)

                    summary_exp_file = output_dir+'/'+pval+'.summary.txt'
                    #if analyze_metaprobesets == 'yes': annotateMetaProbesetGenes(summary_exp_file, stats_file, metaprobeset_file, species)
                    try:
                        shutil.copyfile(summary_exp_file, stats_file)
                        os.remove(summary_exp_file)
                    except Exception:
                        print traceback.format_exc()
                        null=[] ### Occurs if dabg export failed
                    
                    if analyze_metaprobesets == 'yes':
                        residual_destination_file = string.replace(expression_file,'exp.','residuals.')
                        residual_exp_file = output_dir+'/'+algorithm+'.residuals.txt' 
                        #shutil.copyfile(residual_exp_file, residual_destination_file);os.remove(residual_exp_file)
                        reformatResidualFile(residual_exp_file,residual_destination_file)
                        residual_dabg_file = output_dir+'/dabg.residuals.txt'; os.remove(residual_dabg_file)
            except NameError:
                status = 'failed'
                #print traceback.format_exc()
            
        cache_delete_status = export.deleteFolder(cache_dir)
        if status == 'failed':
            if architecture == '64bit' and platform.architecture()[0] == '64bit' and (os.name == 'nt' or 'linux' in sys.platform):
                print 'Warning! 64bit version of APT encountered an error, trying 32bit.'
                ### If the above doesn't work, try 32bit architecture instead of 64bit (assuming the problem is related to known transient 64bit build issues)
                for dataset in exp_file_location_db: ### Instance of the Class ExpressionFileLocationData
                    fl = exp_file_location_db[dataset]; fl.setArchitecture('32bit')
                probesetSummarize(exp_file_location_db,analyze_metaprobesets,probeset_type,species,root)            
            else:
                print_out = 'apt-probeset-summarize failed. See log and report file in the output folder under "ExpressionInput/APT-output" for more details.'
                try:
                    WarningWindow(print_out,'Exit')
                    root.destroy()
                except Exception:
                    print print_out; force_exit
        else:
            print 'CEL files successfully processed. See log and report file in the output folder under "ExpressionInput/APT-output" for more details.' 
            
