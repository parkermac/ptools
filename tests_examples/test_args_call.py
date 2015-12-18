"""
Test of calling a function that expects args.

RESULT: both work.
"""

# Use os.system:
import os
os.system('python ./test_args.py -t "howdy from os"')

# Use subprocess
import subprocess
run_cmd = ['python', './test_args.py',  '-t', 'hello from subprocess']

# As soon as you instantiate Popen the program runs.
proc = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# Then read the output and error messages.
out, err = proc.communicate()   
# And print them.    
print(out)

# you could also call the function this way
#subprocess.call(run_cmd)