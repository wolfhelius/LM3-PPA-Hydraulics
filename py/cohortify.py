# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 11:57:00 2013

@author: adamwolf
"""

"""

"cohorts": [
    {
    "id" = 1
    "species" = 3
    "startyear" = 1956 # make this into a datetime
    "age" : []
    "bl" : []
    "bw" : []    
    }
    ]
    
"""

import json
import cohort

prefix = "/Users/adamwolf/Research/LM3-PPA-r44/work/WOLF/vegn_cc_"
vegn_cc = ["age","bl","bliving","blv","br","bseed","bsw","bwood","crownarea","dbh","height","id","layer","nindivs","nsc","species"]

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

tile = {"cohorts" : [],
        "ccidMax" : 0
        }

# get ccid row for current year
ccid = [int(i) for i in data]

# initialize cohorts in tile
for i in range(len(ccid)):
    tile['cohorts'].append({"ccid":ccid[i], 
        "age":[], "bl":[], "bliving":[], "blv":[], "br":[], 
        "bseed":[],"bsw":[],"bwood":[],"crownarea":[],"dbh":[],
        "height":[],"layer":[],"nindivs","nsc"})
    
# maintain record of member ccids
ccids = [cohort["ccid"] for cohort in tile["cohorts"]]

def add_keyval_to_cc(id, key, value):
    for cohort in tile["cohorts"]:
        if cohort["ccid"] == id:
            cohort.update({key: value})

#tile['cohorts'].append({"ccid":1, "species":2, "startyear":1950})
#tile['cohorts'].append({"ccid":2, "species":1, "startyear":1950})

# access as
# tile['cohorts'][0]['ccid']
# json.dumps(tile, sort_keys=True, indent=4, separators=(',', ': '))

# two equivalent ways to get data out of this structure:
"""
ccids = [cohort["ccid"] for cohort in tile["cohorts"]]
ccids2 = []
for cohort in tile["cohorts"]:
    ccids2.append(cohort["ccid"])
"""

"""
import json

data = json.loads(json_string)
counties = [item for item in data["features"] 
            if item["classifiers"][0]["subcategory"] == "County"]
            
# or

counties = [item["name"] for item in data["features"]
            if item["classifiers"][0]["subcategory"] ==  "County"]            

"""

"""
# check if item is in list:
some_list = ['abc-123', 'def-456', 'ghi-789', 'abc-456']
if any("abc" in s for s in some_list):
    # whatever

"""
