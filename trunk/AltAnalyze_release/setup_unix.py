from cx_Freeze import setup, Executable

setup(
        name = "AltAnalyze",
        version = "2.0",
        description = "Conventional and Alternative Exon Microarray Analysis Suite",
        executables = [Executable("AltAnalyze.py")])
