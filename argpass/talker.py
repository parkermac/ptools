"""
Test of using command line arguments, and calling a MATLAB function
with command line arguments
"""
import argparse
parser = argparse.ArgumentParser()
# positional arguments
parser.add_argument("gtag", type=str, help="e.g. cascadia1_base")
parser.add_argument("date_string", type=str, help="e.g. 2014.02.14")
parser.add_argument("backfill", type=str, help="0 for a forecast, 1 for a backfill")
# an optional flag
parser.add_argument("-o", "--opt_str", type=str, help="an optional string")
args = parser.parse_args()

# pass arguments to a matlab program
import subprocess
func = ("talk2me(\'" +
    args.gtag + "\',\'" +
    args.date_string + "\',\'" +
    args.backfill + "\',\'" +
    args.opt_str + "\')")
cmd = "/Applications/MATLAB_R2014b.app/bin/matlab"
run_cmd = [cmd, "-nojvm", "-nodisplay", "-r", func]
proc = subprocess.Popen(run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate() # "out" is the screen output of the matlab code

# check on the text file if desired
if False:
    a = open('blah.txt','r')
    for line in a:
        print line
    a.close()
