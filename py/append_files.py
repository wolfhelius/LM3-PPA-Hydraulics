# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:56:57 2013

@author: adamwolf
"""

import glob

def append_files(rootdir, keyword):
    
    mainfiles = glob.glob(rootdir+'/'+keyword+'*')
    newfiles = glob.glob(rootdir+'/tmp/'+keyword+'*')
    for i in range(len(mainfiles)):
        mainfile = open(mainfiles[i],"a")
        lines = open(newfiles[i]).readlines()
        if (len(lines)==3): 
            mainfile.writelines(lines[-1])
        else:        
            mainfile.writelines(lines[2:-1])
        
if __name__ == "__main__":
    import sys
    append_files(sys.argv[1], sys.argv[2])