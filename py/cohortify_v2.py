# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 11:57:00 2013

@author: adamwolf
"""

import json
from Cohort import Cohort
from cohort_utilities import *

prefix = "/Users/adamwolf/Research/LM3-PPA-r44/work/WOLF/vegn_cc_"
vegn_cc = ["age","bl","bliving","blv","br","bseed","bsw","bwood","crownarea","dbh","height","layer","nindivs","nsc","species"]

vegnfile = prefix+"id"

f = open(vegnfile)
line = f.readline()
line = f.readline()
line = f.readline()
line = line.rstrip()
line = line.split(',')
date = line[0].split(' ')
data = line[1].split()
ymd = date[0].split('-')
year = int(ymd[0])

# initialize objects
time = []

tile = {}

# get ccid row for current year
ccid = [int(i) for i in data]

# initialize cohorts in tile
for i in range(len(ccid)):
    tile[ccid[i]] = Cohort(ccid[i])

tileJSON = {}
for key in tile:
    tileJSON[key] = tile[key].__dict__
    
# open each file in turn, taking each line    
#for suffix in vegn_cc:
suffix = "age"
vegnfile = prefix+suffix

g, data = READ_VEGN_INIT(vegnfile)
    
    
# add value to array    
# tile[2].age.append(12)
    
# dump object as json:    
# json.dumps(tile[2].__dict__)
  
# see if cohort already in tile:  
# tile.has_key(1)

# [json.dumps(tile[cohort].__dict__) for cohort in tile]  
    
