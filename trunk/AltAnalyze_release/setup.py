from distutils.core import setup
import py2exe

setup(windows=[{"script":'AltAnalyze.py',"icon_resources":[(1,"AltAnalyze_W7.ico")]}])

#setup(console=[{"script":'AltAnalyze.py',"icon_resources":[(0,"altanalyze.ico")]}]) 