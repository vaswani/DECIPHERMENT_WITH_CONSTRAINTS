#! /usr/bin/env python

import sys
import math
import random
from numpy import *
from scipy import *
import time

#implementing the linear time projection algorithm from John Dhuchi's paper. Efficient projections onto the L1-ball for learning in High Dimensions

global parameter_counter #keeps the size of the vector 

def linearProjectionOntoSimplex(v) :
	global parameter_counter
	u = array(range(0,parameter_counter))
	print u
	s = 0.0
	rho = 0.0
	current_u_size  = u.size
	while current_u_size != 0 :
		k = random.randint(0, current_u_size)		
		g = array([])
		l = array([])
		delta_s = 0.0
		delta_rho = 0.0
		index_for_k = 0
		num_elements_in_g = 0
		print 's in beginning of loop is ',s
		print 'rho in beginning of loop is ',rho

		for j in u :
			if (v[j] >= v[k]) :
				g = append(g,j)
				if j == k :
					index_for_k = num_elements_in_g
				delta_s += v[j]
				#delta_rho += 1
				num_elements_in_g += 1
			else :
				l = append(l,j)
		print 'delta s is ',delta_s
		delta_rho = num_elements_in_g
		condition1 = s + delta_s
		condition2 = (rho + delta_rho)*v[k]
		if condition1 - condition2 < 1 :
			print 'yeah'
			s += delta_s
			rho += delta_rho
			u = l
		else :
			u  = delete(g,index_for_k)
			#u.pop(index_for_k)
		
		current_u_size= u.size
		print 's in loop is ',s
		print 'rho in loop is ',rho

	#	print 'iteration'
	#	print u	
	print s 
	print rho
	theta = (s - 1)/rho
	final_vector = array([getMax(value-theta,0) for value in v])
	return(final_vector)

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
		temp_rho = value - 1/float(i+1) * (sum - 1)	
		if temp_rho > 0 :
			max_index = i
			max_index_sum = sum
	
	theta = 1/float(max_index+1) * (max_index_sum -1)

	#print 'theta is ',theta			
	final_vector = array([getMax(value-theta,0) for value in v])
	return(final_vector)


def getMax(value,theta) :
	return(max(value-theta,0))


	
def main() :
	global parameter_counter 
	v = array([0.1,0.2,0.3,0.4,0.5,0.1,0.6,0.7,0.8,1.1])#,2.3])#,5.6,1.3,5.0])
	parameter_counter = v.size
	start = time.clock()
	for i in range (0,100) :
	 	linear_projected_vector = linearProjectionOntoSimplex(v)
	elapsed_time = time.clock() - start
	print 'elapsed time is ', elapsed_time
	print 'linear projected vector is '
	print linear_projected_vector

	start = time.clock()

	for i in range (0,100) :
		log_n_projected_vector = projectOntoSimplex(v)
	elapsed_time = time.clock() - start
	print 'elapsed time is ', elapsed_time

	print 'log projected vector is '
	print linear_projected_vector

main()
