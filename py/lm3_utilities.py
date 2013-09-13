# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 11:41:42 2013

@author: adamwolf
"""

import glob

mainfiledir = 'work/WOLF/'
newfiledir = 'work/WOLF/tmp/'

def append_files(keyword):
    
    mainfiles = glob.glob(mainfiledir+keyword+'*')
    newfiles = glob.glob(newfiledir+keyword+'*')
    for i in range(mainfiles):
        mainfile = open(mainfiles[i],"a")
        lines = open(newfiles[i]).readlines()
        mainfile.writelines(lines[3:-1])

if __name__ == "__main__":
    import sys
    append_files(int(sys.argv[1]))