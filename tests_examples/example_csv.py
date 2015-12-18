"""
This is code to experiment with reading and writing information
such as csv files.

RESULTS: it seems easiest to write csv files using csv.writer, and
easiest to read them using pd.csv_read (especially if you have a header
row and an index column)
"""
# WRITING

# define the data
header = ('month','temp')
month_data = ['jan', 'feb', 'mar']
temp_data = [10, 20, 30]
# create a list of tuples, both strings and ints
data = zip(month_data,temp_data)

dir0 = '../../ptools_output/'

# (1) create a csv file using basic commands
fn1 = dir0 + 'test1.csv'
# open and write to the file
a1 = open(fn1,'wb')
a1.write(header[0] + ',' + header[1] + '\n')
for val in data:
    a_line = str(val[0]) + ',' + str(val[1]) + '\n'
    a1.write(a_line)
a1.close()
# RESULT: this works but is a lot of typing

# (2) try using the csv module
import csv
fn2 = dir0 + 'test2.csv'
with open(fn2,'wb') as a2:
    ww = csv.writer(a2)
    ww.writerow(header)
    for val in data:
        ww.writerow(val)
# RESULT: much cleaner


# (3) try using the csv module, but with no header
fn3 = dir0 + 'test3.csv'
with open(fn3,'wb') as a3:
    ww = csv.writer(a3)
    for val in data:
        ww.writerow(val)
# RESULT: also works fine
# but need to use header=None if reading with pd.read_csv
#
# and here is how you would read this into a matlab structure "info":
"""
fid = fopen('test3.csv','r');
C = textscan(fid,'%s%s','Delimiter',',');
fclose(fid);
items = C{1};
values = C{2};
info = struct();
for ii = 1:length(items)
    info.(items{ii}) = values{ii};
end
"""

# READING

# read the result into a pandas DataFrame
import pandas as pd
df1 = pd.read_csv(fn1, index_col='month')
df2 = pd.read_csv(fn2, index_col='month')
df3 = pd.read_csv(fn3, header=None)
# RESULT: this a super easy, but you have to be willing to 
# work with a DataFrame
# (a benefit of DataFrame is that it intuits type)
# EXAMPLES
x_point = df2.ix['jan', 'temp'] # returns 10 as an int [row, column]
x_data = df2.values # returns a numpy array of the data
x_index = df2.index # returns a pandas Index
x_column_series = df2['temp'] # returns that column as a pandas Series
x_row_series = df2.ix['jan'] # returns that row as a pandas Series
# and note that items in a Series can be accessed simply by number
xc_point = x_column_series[1]
xr_point = x_row_series[0]
# both return ints in this example

# or read results using csv
data_dict = dict()
data_list = []
with open(fn3) as a3:
    for val in csv.reader(a3):
        data_dict[val[0]] = int(val[1])
        data_list.append(val) # list items all strings
# this returns a dict and a list of lists
# and automatically closes the input file
# and the only disadvantage of a dict is that it is unordered.
# Another disadvantage of this method over a DataFrame is that it assumes
# there is just one values column.
