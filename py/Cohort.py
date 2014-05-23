# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 22:41:22 2013

@author: adamwolf
"""

"""
"cohorts": [
    {
    "id" = 1
    "species" = 3
    "startyear" = 1950 # make this into a datetime
    "age" : []
    "bl" : []
    "bw" : []    
    }
    ]

"""
    
class Cohort:
    def __init__(self, id):
        self.id = id
        self.startyear = 0
        self.species=0
        self.age=[]
        self.bliving=[]
        self.br=[]
        self.bseed=[]
        self.bsw=[]
        self.bwood=[]
        self.crownarea=[]
        self.dbh=[]
        self.height=[]
        self.layer=[]
        self.nindivs=[]
        self.nsc=[]
        self.bl=[]
        self.bl=[]
        