from pymongo import MongoClient
import json
from bson.objectid import ObjectId
from datetime import datetime

client = MongoClient('mongodb://localhost:27017/')
tiles = {}
for tile in client.lm3.lm3_runs.find():
    tile_id = str(tile['_id'])
    timestamp = tile['_id'].generation_time.isoformat()
    tiles.update({str(tile['_id']): tile['_id'].generation_time.isoformat()})
