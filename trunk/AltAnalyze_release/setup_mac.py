import sys

_script = 'AltAnalyze.py'
_appVersion = '2.0.7'
_appDescription = "AltAnalyze is a freely available, open-source and cross-platform program that allows you to take RNASeq"
_appDescription +="or relatively raw microarray data (CEL files or normalized), identify predicted alternative splicing or"
_appDescription +="alternative promoter changes and view how these changes may affect protein sequence, domain composition,"
_appDescription +="and microRNA targeting."
_authorName = 'Nathan Salomonis'
_authorEmail = 'nsalomonis@gmail.com'
_authorURL = 'http://altanalyze.org'

excludes = ["matplotlib", "numpy","scipy","wxPython"]

if sys.platform.startswith("darwin"):
        from distutils.core import setup
        import py2app
        options = {"py2app":
            {"includes": ["mpmath"], "excludes": excludes,
            "iconfile": "altanalyze.icns"}}
        setup(name="AltAnalyze",
            app=[_script],
            version=_appVersion,
            description=_appDescription,
            author=_authorName,
            author_email=_authorEmail,
            url=_authorURL,
            options=options,
            #data_files=data_files,
            setup_requires=["py2app"])
"""
from distutils.core import setup
import py2app
setup(app=['AltAnalyze.py'])"""