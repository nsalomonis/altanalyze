
import subprocess
retcode = subprocess.call(['python', 'AltAnalyze.py', "--","-i", indexFile, "-o", output_path,"--single","-l","200"]+p)