"""
Code to test running multiple jobs using subprocess, and waiting for all to finish.

Following the excellent tutorial in:

https://medium.com/python-features/achieving-parallelism-in-python-20b973862310

"""

import subprocess
from time import time

# function defining the subprocess
def run_sleep(period):
    proc = subprocess.Popen(['sleep', str(period)])
    return proc

# run a number of jobs
start = time()
procs = []
for ii in range(10):
    proc = run_sleep(3)
    procs.append(proc)

# the proc.communicate() method will only return after a job is done
# so this loop effectively checks on all jobs sequentially
for proc in procs:
    proc.communicate()
    
# you can tell that this experiment worked because to total time
# was just over 3 sec, even though we had 10 jobs.
end = time()
print ('Finished in %.3F' % (end - start))
