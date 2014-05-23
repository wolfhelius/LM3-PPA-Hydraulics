# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 12:56:18 2013

Helper files called by Cohortify

@author: adamwolf
"""

def READ_VEGN_INIT(filename, flag):
    
    f = open(filename)
    line = f.readline()
    line = f.readline()
    line = f.readline()
    line = line.rstrip()
    line = line.split(',')
    date = line[0].split(' ')
    data = line[1].split()    
    if flag:
        data = [int(d) for d in data]
    else:
        data = [float(d) for d in data]
        
    return (f, date, data)
    
def READ_VEGN(f, flag):

    line = f.readline()
    line = line.rstrip()
    line = line.split(',')
    date = line[0].split(' ')
    data = line[1].split()
    if flag:
        data = [int(d) for d in data]
    else:
        data = [float(d) for d in data]
    
    return (f, date, data)    
    
class fdata:
    def __init__(self, _s, _flag):    
        self.s = _s
        self.flag=_flag
        self.f=0
        self.data=[]
