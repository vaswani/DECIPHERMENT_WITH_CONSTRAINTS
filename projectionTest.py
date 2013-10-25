#! /usr/bin/env python

import sys
from numpy import *
from math import *

def main() :
	temp = array([13.0,2.0,3.0,14.0,13.2,14.5,7.0])
	print ' temp is'
	print temp
	projection = projectOntoSimplex(temp) 
	print 'projection is '
	print projection

def projectOntoSimplex(v) :
	array_size = v.size
	mu_vector = array(v)
	mu_vector.sort()
	mu_vector = mu_vector[::-1]
	#print mu_vector
	#raw_input()
	sum_vector =  zeros(array_size)
	sum=0
	max = -999999999999999
	max_index = -1
	max_index_sum = 0.0
	for i,value in enumerate(mu_vector) :
		sum += value
		#print 'sum', sum
		#print 'value',value
		#print 'i',i
		#raw_input()
		#print 'blah1 ',sum-1
		#print 'blah',1/float(i+1) * (sum - 1)
		temp_rho = value - 1/float(i+1) * (sum - 1)	
		#print 'temp rho',temp_rho
		#raw_input()
		if temp_rho > 0 :
			max_index = i
			max_index_sum = sum
		#print 'current max index', max_index
		#print 'current max index sum', max_index_sum
		#print 'current max',max
		#raw_input()
		
	theta = 1/float(max_index+1) * (max_index_sum -1)

	#print 'theta is ',theta			
	final_vector = array([getMax(value-theta,0) for value in v])
	return(final_vector)

def getMax(value,theta) :
	return(max(value-theta,0))
main()
