"""
test parallel subprocess
https://stackoverflow.com/questions/26774781/
python-multiple-subprocess-with-a-pool-queue-recover-output-as-soon-as-one-finis
"""
import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

from multiprocessing.pool import ThreadPool
import subprocess

import multiprocessing
ncpu = multiprocessing.cpu_count()
print('Number of cpu available = ' + str(ncpu))

def work(sample):
    #my_tool_subprocess = subprocess.Popen('mytool {}'.format(sample), shell=True, stdout=subprocess.PIPE)
    
    if False:
        func = 'test_parallel_worker'
        cmd = Ldir['which_matlab']
        run_cmd = [cmd, "-nodisplay", "-r", func, "&"]
    else:
        func = 'test_parallel_worker.py'
        cmd = 'python'
        run_cmd = [cmd, func, "&"]
        
    my_tool_subprocess = subprocess.run(run_cmd, stdout=subprocess.PIPE, cwd='/Users/pm7/Documents/ptools/tests_examples')
    

num = None  # set to the number of workers you want (it defaults to the cpu count of your machine)
tp = ThreadPool(num)

from datetime import datetime
for sample in range(10):
    tp.apply_async(work, (sample,))

tp.close()
tp.join()