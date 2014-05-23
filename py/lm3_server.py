# Server for HFM data

import bottle
from pymongo import MongoClient
import json
from bson.objectid import ObjectId

lm3vis = bottle.Bottle()

class JSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, ObjectId):
            return str(o)
        return json.JSONEncoder.default(self, o)

@lm3vis.route('/')
def query():
	#client = MongoClient('mongodb://localhost:27017/')
	client = MongoClient('mongodb://lm3user:ezpass@emma.mongohq.com:10060/lm3')    
	tiles = {}
	for tile in client.lm3.lm3_runs.find():
		tile_id = str(tile['_id'])
		timestamp = tile['_id'].generation_time.isoformat()
		tiles.update({str(tile['_id']): tile['_id'].generation_time.isoformat()})
	return bottle.template('tile_list',tiles=tiles)

@lm3vis.route('/lm3_dashboard/<tile_id>')
def dashboard(tile_id):
	#client = MongoClient('mongodb://localhost:27017/')    
	client = MongoClient('mongodb://lm3user:ezpass@emma.mongohq.com:10060/lm3')    
	tile = client.lm3.lm3_runs.find_one({"_id": ObjectId(tile_id)})
	tile = JSONEncoder().encode(tile)

	@lm3vis.route('/tile.json')
	def tileJSON():
		return tile

	return open('lm3_dashboard')

@lm3vis.route('/tile/<tile_id>')
def tile(tile_id):
	#client = MongoClient('mongodb://localhost:27017/')   
	client = MongoClient('mongodb://lm3user:ezpass@emma.mongohq.com:10060/lm3')     
	tile = client.lm3.lm3_runs.find_one({"_id": ObjectId(tile_id)})
	tile = JSONEncoder().encode(tile)
	return tile

@lm3vis.route('/delete/<tile_id>')
def tile(tile_id):
	#client = MongoClient('mongodb://localhost:27017/')    
	client = MongoClient('mongodb://lm3user:ezpass@emma.mongohq.com:10060/lm3')    
	client.lm3.lm3_runs.remove({"_id": ObjectId(tile_id)})
	return '''<a href="http://localhost:8080/">Return</a>'''

@lm3vis.route('/lm3_dashboard.html')
def lm3db():
	#client = MongoClient('mongodb://localhost:27017/')  
	client = MongoClient('mongodb://lm3user:ezpass@emma.mongohq.com:10060/lm3')      
	tile = client.lm3.lm3_runs.find_one()
	tile = JSONEncoder().encode(tile)

	@lm3vis.route('/tile.json')
	def tileJSON():
		return tile
		
	return open('lm3_dashboard.html')

@lm3vis.route('/g/climate.html')
def climate():
	return open('g/climate.html')

@lm3vis.route('/g/climate4.csv')
def climate_data():
	return open('g/climate4.csv')

if __name__ == "__main__":
        bottle.run(lm3vis, host='localhost', port=8080, debug=True)
