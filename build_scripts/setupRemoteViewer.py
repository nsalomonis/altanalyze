#!/usr/local/bin/python2.6

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

_script = 'AltAnalyzeViewer.py'
_appName = "AltAnalyzeViewer"
_appVersion = '1.0.1'
_appDescription = "AltAnalyze is a freely available, open-source and cross-platform program that allows you to take RNASeq or "
_appDescription +="relatively raw microarray data (CEL files or normalized), identify predicted alternative splicing or alternative "
_appDescription +="promoter changes and view how these changes may affect protein sequence, domain composition, and microRNA targeting."
_authorName = 'Nathan Salomonis'
_authorEmail = 'nsalomonis@gmail.com'
_authorURL = 'http://www.altanalyze.org'
_appIcon = "Viewer.ico"

excludes = ["igraph","patsy","pandas","suds","lxml","cairo","cairo2","ImageTk","PIL","Pillow","mpmath",'pysam','Bio']
excludes = ["igraph","patsy","pandas","suds","lxml","cairo","cairo2","mpmath",'virtualenv','Tkinter','matplotlib.tests',"Pillow"]
#excludes = []
includes = ["wx"] #["suds", "mpmath", "numpy"]
includes = ["mpmath", "numpy",'pysam.TabProxies','pysam.ctabixproxies','dbhash','anydbm']
""" By default, suds will be installed in site-packages as a .egg file (zip compressed). Make a duplicate, change to .zip and extract
here to allow it to be recognized by py2exe (must be a directory) """


matplot_exclude = [] #['MSVCP90.dll']
scipy_exclude = [] #['libiomp5md.dll','libifcoremd.dll','libmmd.dll']

""" xml.sax.drivers2.drv_pyexpat is an XML parser needed by suds that py2app fails to include. Identified by looking at the line: parser_list+self.parsers in
/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/PyXML-0.8.4-py2.7-macosx-10.6-intel.egg/_xmlplus/sax/saxexts.py
check the py2app print out to see where this file is in the future

(reported issue - may or may not apply) For mac and igraph, core.so must be copied to a new location for py2app:
sudo mkdir /System/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/lib-dynload/igraph/
cp /Library/Python/2.6/site-packages/igraph/core.so /System/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/lib-dynload/igraph/

"""

if sys.platform.startswith("darwin"):
        ### Local version: /usr/local/bin/python2.6
        ### example command: python setup.py py2app
        from distutils.core import setup
        import py2app
        #import lxml
        includes+= ["pkg_resources","distutils","lxml.etree","lxml._elementpath",'pysam.TabProxies','pysam.ctabixproxies'] #"xml.sax.drivers2.drv_pyexpat"
        """
        resources = ['/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6/pyconfig.h']
        frameworks = ['/System/Library/Frameworks/Python.framework/Versions/2.6/include/python2.6/pyconfig.h']
        frameworks += ['/System/Library/Frameworks/Python.framework/Versions/2.6/Extras/lib/python/pkg_resources.py']
        frameworks += ['/System/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/distutils/util.py']
        frameworks += ['/System/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/distutils/sysconfig.py']
        import pkg_resources
        import distutils
        import distutils.sysconfig
        import distutils.util
        """

        options = {"py2app":
                    {"excludes": excludes,
                     "includes": includes,
                     #"frameworks": frameworks,
                     #"resources": resources,
                     "argv_emulation": True,
                     'arch': 'i386', ## for wx 
                     "iconfile": "Viewer.icns"}
        }
        setup(name=_appName,
                        app=[_script],
                        version=_appVersion,
                        description=_appDescription,
                        author=_authorName,
                        author_email=_authorEmail,
                        url=_authorURL,
                        options=options,
                        #data_files=data_files,
                        setup_requires=["py2app"]
        )

	import unique, shutil
	root_path = unique.filepath('')
	software_path = root_path+'/dist/AltAnalyzeViewer.app/Contents/Frameworks/Tcl.framework'
	shutil.rmtree(software_path)
	software_path = root_path+'/dist/AltAnalyzeViewer.app/Contents/Frameworks/Tk.framework'
	shutil.rmtree(software_path)
	software_path = root_path+'/dist/AltAnalyzeViewer.app/Contents/Resources/mpl-data/sample_data'
	shutil.rmtree(software_path)
	software_path = root_path+'/dist/AltAnalyzeViewer.app/Contents/Resources/lib/python2.7/matplotlib/tests'
	shutil.rmtree(software_path)

if sys.platform.startswith("win"):
        ### example command: python setup.py py2exe
        from distutils.core import setup
        import py2exe
        import suds
        import numpy
        import matplotlib
        import unique
        import lxml
        import sys
        import six ### relates to a date-time dependency in matplotlib
        import pysam
        import TabProxies
        import ctabix
        import csamtools
        import cvcf
        import dbhash
        import anydbm
        #sys.path.append(unique.filepath("Config\DLLs")) ### This is added, but DLLs still require addition to DLL python dir
        from distutils.filelist import findall
        import os
        
        data_files=matplotlib.get_py2exe_datafiles()
        
        matplotlibdatadir = matplotlib.get_data_path()
        matplotlibdata = findall(matplotlibdatadir)
        matplotlibdata_files = []
        for f in matplotlibdata:
            dirname = os.path.join('matplotlibdata', f[len(matplotlibdatadir)+1:])
            matplotlibdata_files.append((os.path.split(dirname)[0], [f]))

        windows=[{"script":_script,"icon_resources":[(1,_appIcon)]}]
        options={'py2exe':
                        {
                        "includes": 'suds',
			"includes": 'lxml',
			'includes': 'lxml.etree',
			'includes': 'lxml._elementpath',
                        "includes": 'matplotlib',
                        "includes": 'mpl_toolkits',
                        "includes": 'matplotlib.backends.backend_tkagg',
                        "dll_excludes": matplot_exclude+scipy_exclude,
                        }}
        setup(
                        windows = windows,
                        options = options,
                        version=_appVersion,
                        description=_appDescription,
                        author=_authorName,
                        author_email=_authorEmail,
                        url=_authorURL,
                        data_files=matplotlibdata_files+data_files,
        )

if sys.platform.startswith("2linux"):
	# bb_setup.py
	from bbfreeze import Freezer
	 
	f = Freezer(distdir="bb-binary")
	f.addScript("RemoteViewer.py")
	f()

if sys.platform.startswith("linux"):
        ### example command: python setup.py build
        includes = ['matplotlib','mpl_toolkits','matplotlib.backends.backend_tkagg']
        includefiles = []

        from cx_Freeze import setup, Executable
        ### use to get rid of library.zip and move into the executable, along with appendScriptToLibrary and appendScriptToExe
        #buildOptions = dict(create_shared_zip = False) 

	setup(
		name = _appName,
		version=_appVersion,
		description=_appDescription,
		author=_authorName,
		author_email=_authorEmail,
		url=_authorURL,
		#options = dict(build_exe = buildOptions),
		options = {"build_exe": {"includes":includes, "include_files": includefiles}},
		executables = [Executable(_script,
				#appendScriptToExe=True,
				#appendScriptToLibrary=False,
				#icon='goelite.ico',
				compress=True)],
	)

