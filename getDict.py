#! /usr/bin/env python
import sys

def stripLine(item) :
	return(item.strip())

lines = map(stripLine,[line for line in open(sys.argv[1])])
words = []
for line in lines :
	words.extend(line.split())

for line in open(sys.argv[2]) : #reading the dictionary and 
	dict_items = line.split()
	if dict_items[0] in words :
		print line.strip()


#print items
	


