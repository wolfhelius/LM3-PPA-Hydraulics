# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 11:57:00 2013

@author: adamwolf
"""

import json
from Cohort import Cohort
from cohort_utilities import *

prefix = "/Users/adamwolf/Research/LM3-PPA-r44/work/WOLF/vegn_cc_"
vegn_cc = ["id", "species", "age","bl","bliving","blv","br","bseed","bsw","bwood","crownarea","dbh","height","layer","nindivs","nsc"]

vegnfile = prefix+"id"

filelist = []
datalist = []
fileJSON = {}    

# s is suffix

for s in vegn_cc:
    vegnfile = prefix+s
    if ('id' == s) | ('species' == s):
        flag = True
    else:
        flag = False
    f, date, data = READ_VEGN_INIT(vegnfile, flag)
    filelist.append(f)
    datalist.append(data)
    fileJSON[s] = fdata(s, flag)
    fileJSON[s].f = f
    fileJSON[s].data = data    
        
ymd = date[0].split('-')
year = int(ymd[0])

# initialize objects
time = []

tile = {}

# get ccid row for current year
ccid = [int(i) for i in datalist[0]]

# initialize cohorts in tile
for i in range(len(ccid)):
    tile[ccid[i]] = Cohort(ccid[i])
    tile[ccid[i]].startyear = year

for i in range(len(ccid)):
    for s in vegn_cc:
        d = fileJSON[s].data
        if ('id' == s) | ('species' == s):
            tile[ccid[i]].__dict__[s] = d[i]
        else:
            d = fileJSON[s].data
            tile[ccid[i]].__dict__[s].append(d[i])

# continue for each available year
while True:
    try:
        # load a new set of data
        for s in vegn_cc:
            fileJSON[s].f, date, fileJSON[s].data = READ_VEGN(fileJSON[s].f, fileJSON[s].flag)
            ymd = date[0].split('-')
            year = int(ymd[0])
            ccid = fileJSON['id'].data
        # for each cohort append the proper data
        for i in range(len(ccid)):
            for s in vegn_cc:
                d = fileJSON[s].data
                
                # initialize new cohorts
                if not tile.has_key(ccid[i]):
                    tile[ccid[i]] = Cohort(ccid[i])
                    tile[ccid[i]].startyear = year
                
                # add id and species for new cohorts
                if ('id' == s) | ('species' == s):
                    if tile[ccid[i]].startyear == year:
                        tile[ccid[i]].__dict__[s] = d[i]
                    else:
                        continue
                else:
                    d = fileJSON[s].data
                    tile[ccid[i]].__dict__[s].append(d[i])
                        
    except:
        print "end year: ", year
        break


# final step to make into JSON
tileJSON = {}
for key in tile:
    tileJSON[key] = tile[key].__dict__
    
    

    
# dump object as json:    
# json.dumps(tile[2].__dict__)
  
# see if cohort already in tile:  
# tile.has_key(1)

# [json.dumps(tile[cohort].__dict__) for cohort in tile]  
    
