"""
Test code to try to figure out why subprocess hangs.

"""

import subprocess

if True:
    # This hung when we call from the command line, or from
    # a python IDE, but after I added the "&" as the last command item
    # it ran correctly from the command line and from the python IDE.
    cmd = "/Applications/MATLAB_R2015b.app/bin/matlab"
    func = "pr()"
    run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func,"> pr00.txt &"]
else:
    # this runs fine anywhere
    run_cmd = ['ls', '-lt']

proc = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

out, err = proc.communicate()

print(out.decode())