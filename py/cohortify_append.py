# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 11:57:00 2013

@author: adamwolf
"""

import json
from Cohort import Cohort
from cohort_utilities import *
from pymongo import MongoClient    
from bson.objectid import ObjectId

class JSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, ObjectId):
            return str(o)
        return json.JSONEncoder.default(self, o)

def cohortify_initialize(run_id, rootdir):

	prefix = rootdir+"/vegn_cc_"
	vegn_cc = ["id", "species", "age","bliving","br","bsw","bwood","crownarea","dbh","height","layer","nindivs"]

	vegnfile = prefix+"id"

	fileJSON = {}    

	# s is suffix
	for s in vegn_cc:
		vegnfile = prefix+s
		if ('id' == s) | ('species' == s):
			flag = True
		else:
			flag = False
		f, date, data = READ_VEGN_INIT(vegnfile, flag)
		fileJSON[s] = fdata(s, flag)
		fileJSON[s].f = f
		fileJSON[s].data = data    
	
	ymd = date[0].split('-')
	year = int(ymd[0])

	# retrieve JSON into mongoDB    
	from pymongo import MongoClient    
	# client = MongoClient('mongodb://localhost:27017/')    
	client = MongoClient('mongodb://lm3user:ezpass@emma.mongohq.com:10060/lm3')    
	tile = client.lm3.lm3_runs.find_one({'_id': ObjectId(run_id)})

	# get ccid row for current year
	ccid = [str(i) for i in fileJSON['id'].data]

	# initialize cohorts in tile
	# for i in range(len(ccid)):
	#     tile[ccid[i]] = Cohort(ccid[i])
	#     tile[ccid[i]].startyear = year

	# work on the first row of data
	for i in range(len(ccid)):
		for s in vegn_cc:
			d = fileJSON[s].data
			# initialize new cohorts
			if not tile.has_key(ccid[i]):
				tile[ccid[i]] = Cohort(ccid[i])
				tile[ccid[i]] = tile[ccid[i]].__dict__
				tile[ccid[i]]["startyear"] = year
	
			# add id and species for new cohorts
			if ('id' == s) | ('species' == s):
				if tile[ccid[i]]["startyear"] == year:
					tile[ccid[i]][s] = d[i]
				else:
					continue
			else:
				d = fileJSON[s].data
				tile[ccid[i]][s].append(d[i])

	# continue for each available year
	while True:
		try:
			# load a new set of data
			for s in vegn_cc:
				fileJSON[s].f, date, fileJSON[s].data = READ_VEGN(fileJSON[s].f, fileJSON[s].flag)
				ymd = date[0].split('-')
				year = int(ymd[0])
				#ccid = fileJSON['id'].data
				ccid = [str(i) for i in fileJSON['id'].data]
			# for each cohort append the proper data
			for i in range(len(ccid)):
				for s in vegn_cc:
					d = fileJSON[s].data
			
					# initialize new cohorts
					if not tile.has_key(ccid[i]):
						tile[ccid[i]] = Cohort(ccid[i])
						tile[ccid[i]] = tile[ccid[i]].__dict__
						tile[ccid[i]]["startyear"] = year
			
					# add id and species for new cohorts
					if ('id' == s) | ('species' == s):
						if tile[ccid[i]]["startyear"] == year:
							tile[ccid[i]][s] = d[i]
						else:
							continue
					else:
						d = fileJSON[s].data
						tile[ccid[i]][s].append(d[i])
					
		except:
			break


	# insert JSON into mongoDB    
	from pymongo import MongoClient    
	# client = MongoClient('mongodb://localhost:27017/')    
	# client.lm3.lm3_runs.save(tile)
	client = MongoClient('mongodb://lm3user:ezpass@emma.mongohq.com:10060/lm3')    
	client.lm3.lm3_runs.save(tile)


if __name__ == "__main__":
    import sys
    cohortify_initialize(sys.argv[1], sys.argv[2])