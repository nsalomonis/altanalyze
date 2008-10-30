from distutils.core import setup
import py2exe

setup(windows=[{"script":'AltAnalyze.py',"icon_resources":[(1,"altanalyze.ico")]}])