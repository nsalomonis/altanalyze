from cx_Freeze import setup, Executable

setup(
        name = "AltAnalyze",
        version = "1.15",
        description = "Conventional and Alternative Exon Microarray Analysis Suite",
        executables = [Executable("AltAnalyze.py")])
