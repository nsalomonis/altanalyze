#!/usr/local/bin/python2.5 

###update
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
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

"""This module contains generic methods for downloading and decompressing files designated
online files and coordinating specific database build operations for all AltAnalyze supported
gene ID systems with other program modules. """

import os
import sys
import unique
import string
import export

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def zipDirectory(dir):
    #http://www.testingreflections.com/node/view/8173
    import zipfile
    dir = filepath(dir); zip_file = dir+'.zip'
    p = string.split(dir,'/'); top=p[-1]
    zip = zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            fullpath = os.path.join(root, f)
            archive_name = os.path.join(top+archive_root, f)
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file

def unzipFiles(filename,dir):
    import zipfile
    output_filepath = filepath(dir+filename)
    try:
        zfile = zipfile.ZipFile(output_filepath)
        for name in zfile.namelist():
            if name.endswith('/'):null=[] ### Don't need to export
            else: 
                try: outfile = export.ExportFile(dir+name)
                except Exception: outfile = export.ExportFile(dir+name[1:])
                outfile.write(zfile.read(name)); outfile.close()
        #print 'Zip extracted to:',output_filepath
        status = 'completed'
    except Exception, e:
        print e
        print 'WARNING!!!! The zip file',output_filepath,'does not appear to be a valid zip archive file or is currupt.'
        status = 'failed'
    return status

def download(url,dir,file_type):
    try: dp = download_protocol(url,dir,file_type); gz_filepath, status  = dp.getStatus()
    except Exception:
        gz_filepath='failed'; status = "Internet connection was not established. Re-establsih and try again."

    if status == 'remove':
        #print "\nRemoving zip file:",gz_filepath
        try: os.remove(gz_filepath); status = 'removed'
        except Exception: null=[] ### Not sure why this error occurs since the file is not open
    return gz_filepath, status

class download_protocol:
    def __init__(self,url,dir,file_type):
        """Copy the contents of a file from a given URL to a local file."""
        filename = url.split('/')[-1]
        if len(file_type) == 2: filename, file_type = file_type ### Added this feature for when a file has an invalid filename
        output_filepath_object = export.createExportFile(dir+filename,dir[:-1])
        output_filepath = filepath(dir+filename)

        print "Downloading the following file:",filename,' ',
        self.original_increment = 10
        self.increment = 0
        import urllib
        from urllib import urlretrieve
        try:
            try: webfile, msg = urlretrieve(url,output_filepath,reporthook=self.reporthookFunction)
            except IOError:
                if 'Binary' in traceback.format_exc(): #IOError: [Errno ftp error] 200 Switching to Binary mode.
                    ### https://bugs.python.org/issue1067702 - some machines the socket doesn't close and causes an error - reload to close the socket
                    reload(urllib)
                    webfile, msg = urlretrieve(url,output_filepath,reporthook=self.reporthookFunction)
                    reload(urllib)
        except:
            print 'Unknown URL error encountered...'; forceURLError
        print ''
        print "\nFile downloaded to:",output_filepath
        if '.zip' in filename:
            try: decompressZipStackOverflow(filename,dir); status = 'completed'
            except Exception:
                status = unzipFiles(filename,dir)
                if status == 'failed': print 'Zip Extraction failed'
            self.gz_filepath = filepath(output_filepath); self.status = 'remove'
            print "zip file extracted..."
        elif '.gz' in filename:
            self.gz_filepath = output_filepath
            if len(file_type)==0: extension = '.gz'
            else: extension = 'gz'
            decompressed_filepath = string.replace(self.gz_filepath,extension,file_type)
            ### Below code can be too memory intensive
            #file_size = os.path.getsize(output_filepath)
            #megabtyes = file_size/1000000.00
            #if megabtyes>5000: force_error ### force_error is an undefined variable which causes an exception
            import gzip; content = gzip.GzipFile(self.gz_filepath, 'rb')
            data = open(decompressed_filepath,'wb')
            #print "\nExtracting downloaded file:",self.gz_filepath
            import shutil; shutil.copyfileobj(content,data)
            self.status = 'remove'
        else: self.gz_filepath = ''; self.status = 'NA'
    def getStatus(self): return self.gz_filepath, self.status
    def reporthookFunction(self, blocks_read, block_size, total_size):
        if not blocks_read:
            print 'Connection opened. Downloading (be patient)'
        if total_size < 0:
            # Unknown size
            print 'Read %d blocks' % blocks_read
        else:
            amount_read = blocks_read * block_size
            percent_read = ((amount_read)*1.00/total_size)*100
            if percent_read>self.increment:
                #print '%d%% downloaded' % self.increment
                print '*',
                self.increment += self.original_increment
            #print 'Read %d blocks, or %d/%d' % (blocks_read, amount_read, (amount_read/total_size)*100.000)

def createExportFile(new_file,dir):
    try:
        fn=filepath(new_file); file_var = open(fn,'w')
    except IOError:
        #print "IOError", fn
        fn = filepath(dir)
        try:
            os.mkdir(fn) ###Re-Create directory if deleted
            #print fn, 'written'
        except OSError: createExportDir(new_file,dir) ###Occurs if the parent directory is also missing
        fn=filepath(new_file); file_var = open(fn,'w')
    return file_var

def downloadCurrentVersionUI(filename,secondary_dir,file_type,root):
    continue_analysis = downloadCurrentVersion(filename,secondary_dir,file_type)
    if continue_analysis == 'no':
        import UI
        root.destroy(); UI.getUserParameters('no'); sys.exit()
    root.destroy()
        
def downloadCurrentVersion(filename,secondary_dir,file_type):
    import UI
    file_location_defaults = UI.importDefaultFileLocations()

    uds = file_location_defaults['url'] ### Get the location of the download site from Config/default-files.csv
    for ud in uds: url_dir = ud.Location() ### Only one entry
    
    dir = export.findParentDir(filename)  
    filename = export.findFilename(filename)
    url = url_dir+secondary_dir+'/'+filename
    
    file,status = download(url,dir,file_type); continue_analysis = 'yes'
    if 'Internet' in status:
        print_out = "File:\n"+url+"\ncould not be found on server or internet connection is unavailable."
        try:
            UI.WarningWindow(print_out,'WARNING!!!')
            continue_analysis = 'no'
        except Exception:
            print url
            print 'cannot be downloaded';die
    elif status == 'remove':
        try: os.remove(file) ### Not sure why this works now and not before
        except Exception: status = status
    return continue_analysis

def decompressZipStackOverflow(zip_file,dir):
    zip_file = filepath(dir+zip_file)
    ###http://stackoverflow.com/questions/339053/how-do-you-unzip-very-large-files-in-python
    import zipfile
    import zlib
    src = open(zip_file,"rb")
    zf = zipfile.ZipFile(src)
    for m in zf.infolist():
        # Examine the header
        #print m.filename, m.header_offset, m.compress_size, repr(m.extra), repr(m.comment)
        src.seek(m.header_offset)
        src.read(30) # Good to use struct to unpack this.
        nm= src.read(len(m.filename))
        if len(m.extra) > 0: ex= src.read(len(m.extra))
        if len(m.comment) > 0: cm=src.read(len(m.comment)) 

        # Build a decompression object
        decomp= zlib.decompressobj(-15)

        # This can be done with a loop reading blocks
        out=open(filepath(dir+m.filename), "wb")
        result=decomp.decompress(src.read(m.compress_size))
        out.write(result); result=decomp.flush()
        out.write(result); out.close()
    zf.close()
    src.close()

if __name__ == '__main__':
    dp = download_protocol('http://may2009.archive.ensembl.org/biomart/martresults/136?file=martquery_1117221814_599.txt.gz','downloaded','')
    