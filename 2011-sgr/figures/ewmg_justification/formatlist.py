from numpy import loadtxt

import os

olddata = open('linelist.original', 'r')


os.system('rm -f linelist.formatted')
newdata = open('linelist.formatted', 'w')

data = olddata.readlines()

for line in data:
    if not line.startswith('-'):

        line = line.split()
        newdata.write("  %4.3f    %s     %1.3f    %1.3f\n" % (float(line[0])*10., line[2][0:5]+'0', float(line[5]), float(line[1]),))

newdata.close()
olddata.close()

