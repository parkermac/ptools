#!/bin/bash

# driver to run python with command line arguments

# define inputs
gtag='cascadia1_base'
date_string='2010.01.03'
backfill='1'
opt_string='hello'
echo Input = $gtag, $date_string, $backfill, $opt_string

# run program
source $HOME/.bashrc
python ./talker.py $gtag $date_string $backfill -o $opt_string

# print results
checkfile='blah.txt'
echo ''
echo 'Output'
while read line
do
  echo $line
done <$checkfile
echo ''