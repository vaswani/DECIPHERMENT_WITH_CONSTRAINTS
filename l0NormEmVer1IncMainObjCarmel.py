#!/usr/bin/env python
#   =================================================================
#   File: toyprob.py
#   =================================================================
#
#   =================================================================
#   Module: Subroutines that define the problem
#   =================================================================
#
#   Last update of any of the components of this module:
#
#   March 25, 2008.
#
#   Users are encouraged to download periodically updated versions of
#   this code at the TANGO home page:
#
#   www.ime.usp.br/~egbirgin/tango/
#
#   ******************************************************************
#   ******************************************************************

#the change in this one is going to be that we are going to include both tag bigrams and tags

from numpy import *
import algencan
import commands
import re
from math import e
import sys
import copy


probability_re = re.compile(r'Corpus probability=(.*) per-symbol')

#implementation of the em algorithm
#import numpy 
#import array
import latticeNode
#import emMethodsMDLVer2 as emMethods
import emMethodsGeneralNoPrint as emMethods
import copy
import sys


#global variables
current_optimization_params = '' #this will hold the parameters that are currently being optmized
fractional_counts_language = {}
fractional_counts_channel = {}
probabilities_channel = {}
probabilities_language = {}
current_probabilities = {}
current_fractional_counts = {}
#constraint_parameters = []
initial_parameter_args = {}
parameter_counter = 0
num_constraints = 0
constraint_tags_dict = {} #this will hold the tag for the constraint. the key is the constraint id and the value is the tag corresponding to this constraint
sigma = float(0)
alpha = float(0)
global_num_parameters = 0
init_option = ''
def inip():
	global parameter_counter
	global num_constraints
	global probabilities_language, probabilities_channel,initial_parameter_args,init_option
	"""This subroutine must set some problem data.a

	For achieving this objective YOU MUST MODIFY it according to your
	problem. See below where your modifications must be inserted.

	Parameters of the subroutine:

	On Entry:

	This subroutine has no input parameters.

	On Return

	n		integer,
			 number of variables,

	x		double precision x(n),
			 initial point,

	l		double precision l(n),
			 lower bounds on x,

	u		double precision u(n),
			 upper bounds on x,

	m		integer,
			 number of constraints (excluding the bounds),

	lambda   double precision lambda(m),
			 initial estimation of the Lagrange multipliers,

	equatn   logical equatn(m)
			 for each constraint j, set equatn(j) = .true. if it is an
			 equality constraint of the form c_j(x) = 0, and set
			 equatn(j) = .false. if it is an inequality constraint of
			 the form c_j(x) <= 0,

	linear   logical linear(m)
			 for each constraint j, set linear(j) = .true. if it is a
			 linear constraint, and set linear(j) = .false. if it is a
			 nonlinear constraint.
	"""

#   Number of variables

	#n = 2
	n = parameter_counter

#   Number of constraints (equalities plus inequalities)

	#m = 2
	m = num_constraints

#   Initial point

	#x = zeros(n)

	x = array([0.00000000001]*n)
	if init_option == 'current_prob' : #then we have to initize it 
		print 'initializing x to the current prob'
		for tag in initial_parameter_args.keys() :
			for next_tag in initial_parameter_args[tag].keys() :
				if initial_parameter_args[tag][next_tag] ==  'ONE' :
					continue
				else :
					parameter_number = initial_parameter_args[tag][next_tag]
					x[parameter_number] = probabilities_channel[tag][next_tag]	
		print x
	else :
		print 'initializing x to zero'
		print x	
#   Lower and upper bounds

	#l = array([-10.0,-1.0e20])a
	l = array([0.000000000001]*n)
	u = array([1.0]*n)
	#u = array([ 10.0, 1.0e20])

#   Lagrange multipliers approximation. Most users prefer to use the
#   null initial Lagrange multipliers estimates. However, if the
#   problem that you are solving is "slightly different" from a
#   previously solved problem of which you know the correct Lagrange
#   multipliers, we encourage you to set these multipliers as initial
#   estimates. Of course, in this case you are also encouraged to use
#   the solution of the previous problem as initial estimate of the
#   solution. Similarly, most users prefer to use rho = 10 as initial
#   penalty parameters. But in the case mentioned above (good
#   estimates of solution and Lagrange multipliers) larger values of
#   the penalty parameters (say, rho = 1000) may be more useful. More
#   warm-start procedures are being elaborated.

	lambda_ = zeros(m)

#   For each constraint i, set equatn[i] = 1. if it is an equality
#   constraint of the form c_i(x) = 0, and set equatn[i] = 0 if
#   it is an inequality constraint of the form c_i(x) <= 0.

	#equatn = [False, False]
	equatn = [1]*m

#   For each constraint i, set linear[i] = 1 if it is a linear
#   constraint, otherwise set linear[i] = 0.

	#linear = [False, True]
	linear = [1]*m
#   In this Python interface evalf, evalg, evalh, evalc, evaljac and
#   evalhc are present. evalfc, evalgjac, evalhl and evalhlp are not.

	coded = [True,  # evalf
			 True,  # evalg
			 True,  # evalh
			 True,  # evalc
			 True,  # evaljac
			 True,  # evalhc
			 False, # evalfc
			 False, # evalgjac
			 False, # evalhl
			 False] # evalhlp

#   Set checkder = 1 if you code some derivatives and you would
#   like them to be tested by finite differences. It is highly
#   recommended.

	#checkder = True
	checkder = False

	return n,x,l,u,m,lambda_,equatn,linear,coded,checkder

#   ******************************************************************
#   ******************************************************************

def evalf(x):
	little_infinity = 90000000
	global initial_parameter_args,current_fractional_counts,sigma,alpha
	#print 'the fractional counts is '
	#print current_fractional_counts
	"""This subroutine must compute the objective function.
	
	For achieving this objective YOU MUST MODIFY it according to your
	problem. See below where your modifications must be inserted.

	Parameters of the subroutine:

	On Entry:

	x		double precision x(n),
			 current point,

	On Return

	f		double precision,
			 objective function value at x,

	flag	 integer,
			 You must set it to any number different of 0 (zero) if
			 some error ocurred during the evaluation of the objective
			 function. (For example, trying to compute the square root
			 of a negative number, dividing by zero or a very small
			 number, etc.) If everything was o.k. you must set it
			 equal to zero.
	"""
	func_value = float(0.0)
	flag = 0
	#we will go through the initial parameter args 
	for tag in initial_parameter_args.keys() :
		for next_tag in initial_parameter_args[tag].keys() :
			if initial_parameter_args[tag][next_tag] == 'ONE' :
				break
			parameter_number = initial_parameter_args[tag][next_tag]
			if x[parameter_number] == float(0) :
				#flag = -1
				#print 'the flag became negavie'
				func_value += current_fractional_counts[tag][next_tag]  * -little_infinity + alpha * math.exp(-x[parameter_number]/sigma)
			else :
				#print 'sigma is %f'%sigma
				#print 'alpha is %f'%alpha
				func_value += current_fractional_counts[tag][next_tag]  * math.log(x[parameter_number]) + alpha * math.exp(-x[parameter_number]/sigma)

	#n = len(x)

	#f = x[n-1]

	#flag = 0
	#print 'the func value is %f' %-func_value
	#print 'sigma was %f' %sigma
	#print 'alpha was %f' %alpha
	#print 'the flag is %d'%flag
	#raw_input()
	return -func_value,flag

#   ******************************************************************
#   ******************************************************************

def evalg(x):
	little_infinity = 2000
	global initial_parameter_args,current_fractional_counts,alpha,sigma,parameter_counter,num_constraints
	g = array([0.0]*parameter_counter)
	"""
	This subroutine must compute the gradient vector of the objective
	function.

	For achieving these objective YOU MUST MODIFY it in the way specified
	below. However, if you decide to use numerical derivatives (we dont
	encourage this option at all!) you dont need to modify evalg.

	Parameters of the subroutine:

	On Entry:

	x		double precision x(n),
			 current point,

	On Return
https://www.google.com/accounts/ServiceLogin?service=mail&passive=true&rm=false&continue=http%3A%2F%2Fmail.google.com%2Fmail%2F%3Fui%3Dhtml%26zy%3Dl&bsv=zpwhtygjntrz&scc=1&ltmpl=default&ltmplcache=2&hl=en
	g		double precision g(n),
			 gradient vector of the objective function evaluated at x,

	flag	 integer,
			 You must set it to any number different of 0 (zero) if
			 some error ocurred during the evaluation of any component
			 of the gradient vector. (For example, trying to compute
			 the square root of a negative number, dividing by zero or
			 a very small number, etc.) If everything was o.k. you
			 must set it equal to zero.
	"""
	#flag = 0
	for tag in initial_parameter_args.keys() :
		for next_tag in initial_parameter_args[tag].keys() :
			if initial_parameter_args[tag][next_tag] == 'ONE' :
				break
			parameter_number = initial_parameter_args[tag][next_tag]
			if x[parameter_number] == float(0) :
				#print' the flag became negative'
				#flag = -1
				#break
				g[parameter_number]  = -1*(current_fractional_counts[tag][next_tag] * little_infinity - alpha*1/sigma)
			else :
				g[parameter_number]  = -1*(current_fractional_counts[tag][next_tag] /x[parameter_number] - alpha*math.exp(-x[parameter_number]/sigma)/sigma)
		#if flag== -1 :
		#	break
	#n = len(x)

	#g = zeros(n)

	#g[n-1] = 1.0

	flag = 0

	return g,flag

#   ******************************************************************
#   ******************************************************************

def evalh(x):
	global current_fractional_counts,initial_parameter_args,alpha,sigma
	"""This subroutine might compute the Hessian matrix of the objective
	function.

	For achieving this objective YOU MAY MODIFY it according to your
	problem. To modify this subroutine IS NOT MANDATORY. See below
	where your modifications must be inserted.

	Parameters of the subroutine:

	On Entry:

	x		double precision x(n),
			 current point,

	On Return

	hnnz	 integer,
			 number of perhaps-non-null elements of the computed
			 Hessian,

	hlin	 integer hlin(nnzh),
			 see below,

	hcol	 integer hcol(nnzh),
			 see below,

	hval	 double precision hval(nnzh),
			 the non-null value of the (hlin(k),hcol(k)) position
			 of the Hessian matrix of the objective function must
			 be saved at hval(k). Just the lower triangular part of
			 Hessian matrix must be computed,

	flag	 integer,
			 You must set it to any number different of 0 (zero) if
			 some error ocurred during the evaluation of the Hessian
			 matrix of the objective funtion. (For example, trying
			 to compute the square root of a negative number,
			 dividing by zero or a very small number, etc.) If
			 everything was o.k. you must set it equal to zero.
	"""
	global parameter_counter
	num_parameters = parameter_counter
	
	#hnnz = 0
	hnnz = num_parameters #* (num_parameters  +1) / 2
	hlin = zeros(hnnz, int)
	hcol = zeros(hnnz, int)
	hval = zeros(hnnz, float)
	
	counter = 0
	for tag in initial_parameter_args.keys() :
		for next_tag in initial_parameter_args[tag].keys() :
			if initial_parameter_args[tag][next_tag] == 'ONE' :
				continue
			else :
				parameter_number = initial_parameter_args[tag][next_tag]
				hessian_term = current_fractional_counts[tag][next_tag] / (x[parameter_number]*x[parameter_number]) - alpha * math.exp(-x[parameter_number]/sigma)/(sigma*sigma)
				hlin[parameter_number] = parameter_number
				hcol[parameter_number] = parameter_number
				hval[parameter_number] = hessian_term
	
	'''
	print 'the value of the counter is %d'%counter
	for i in range (0,hnnz) :
		if hlin[i] == hcol[i] :
			hval[i] = 1 
		else :
			hval[i] = 0
	print 'printing the hessian of the objective function'
	'''

	'''	
	for i in range (0,hnnz) :
		print 'printing counter'
		print i
		print hlin[i]
		print hcol[i]
		print hval[i]
	'''
	'''
	hlin = zeros(hnnz, int)
	hcol = zeros(hnnz, int)
	hval = zeros(hnnz, float)
	'''
	flag = 0

	#return hlin,hcol,hval,hnnz,flag
	return hlin,hcol,hval,hnnz,flag
	#return ()

#   ******************************************************************
#   ******************************************************************

def evalc(x,ind):
	global constraint_tags_dict, current_fractional_counts,initial_parameter_args,num_constraints
	"""This subroutine must compute the ind-th constraint.

	For achieving this objective YOU MUST MOFIFY it according to your
	problem. See below the places where your modifications must be
	inserted.

	Parameters of the subroutine:

	On Entry:

	x		double precision x(n),
			 current point,

	ind	  integer,
			 index of the constraint to be computed,

	On Return

	c		double precision,
			 i-th constraint evaluated at x,

	flag	 integer
			 You must set it to any number different of 0 (zero) if
			 some error ocurred during the evaluation of the
			 constraint. (For example, trying to compute the square
			 root of a negative number, dividing by zero or a very
			 small number, etc.) If everything was o.k. you must set
			 it equal to zero.
	"""
	
	'''
	n = len(x)

	if ind == 1:
		c = x[0] * x[0] + 1.0 - x[n-1]

		flag = 0

	elif ind == 2:
		c = 2.0 - x[0] - x[n-1]

		flag = 0

	else:
		c = 0.0

		flag = -1
	'''
	#first we have to get the tag for the constraint
	tag_corresponding_to_constraint = constraint_tags_dict[ind]
	#print 'the x vector is'
	#print x
	#print 'the constraint number is %d'%ind
	#print 'the tag corresponding to the constraint is %s'%tag_corresponding_to_constraint

	c = float(-1)
	for next_tag in initial_parameter_args[tag_corresponding_to_constraint].keys() :
		#print 'the next tag is %s'%next_tag
		parameter_number = initial_parameter_args[tag_corresponding_to_constraint][next_tag] 
		c += x[parameter_number]
	#raw_input()
	if ind > num_constraints :
		flag = -1
	else :
		flag = 0
	#print 'the value of the constraint was'
	#print c
	#raw_input()
	return c,flag

#   ******************************************************************
#   ******************************************************************

def evaljac(x,ind):
	global initial_parameter_args,constraint_tags_dict
	"""This subroutine must compute the gradient of the ind-th constraint.

	For achieving these objective YOU MUST MODIFY it in the way specified
	below.

	Parameters of the subroutine:

	On Entry:

	x		double precision x(n),
			 current point,

	ind	  integer,
			 index of the constraint whose gradient will be computed,

	On Return

	jcnnz	integer,
			 number of perhaps-non-null elements of the computed
			 gradient,

	jcvar	integer jcvar(jcnnz),
			 see below,

	jcval	double precision jcval(jcnnz),
			 the non-null value of the partial derivative of the i-th
			 constraint with respect to the jcvar(k)-th variable must
			 be saved at jcval(k).

	flag	 integer
			 You must set it to any number different of 0 (zero) if
			 some error ocurred during the evaluation of the
			 constraint. (For example, trying to compute the square
			 root of a negative number, dividing by zero or a very
			 small number, etc.) If everything was o.k. you must set
			 it equal to zero.
	"""

	n = len(x)
	
	
	#jcnnz = 2
	jcnnz = n

	jcvar = zeros(jcnnz, int)
	jcval = zeros(jcnnz, float)
	
	'''
	if ind == 1:
		jcvar[0] = 0
		jcval[0] = 2.0 * x[0]

		jcvar[1] = n - 1
		jcval[1] = - 1.0

		flag = 0

	elif ind == 2:
		jcvar[0] =  0
		jcval[0] = -1.0

		jcvar[1] = n - 1
		jcval[1] = - 1.0

		flag = 0

	else:
		flag = -1
	'''
	flag  = 0
	constraint_tag = constraint_tags_dict[ind] #getting the tag corresponding to the constraint id
	#now, we will go over the tag bigrams in initial parameter args which correspond to this constraint and then put a one in those place and leave the rest as zeros 
	for next_tag in initial_parameter_args[constraint_tag].keys() :
		parameter_number = initial_parameter_args[constraint_tag][next_tag]
		#jcvar[parameter_number] = parameter_number
		jcval[parameter_number] =  1.0

	#the ith position corresponds to the ith variable
	for i in range (0,len(jcvar)) :
		jcvar[i] = i
	#print 'teh constraint number is %d'%ind
	#print 'the number of variables is %d'%len(jcvar)
	#print 'the number of values is %d'%len(jcval)	
	return jcvar,jcval,jcnnz,flag

#   ******************************************************************
#   ******************************************************************

def evalhc(x,ind):
	global num_constraints
	#print 'the constraint number is %d'%ind
	#print 'the number of constraints is %d'%num_constraints
	"""This subroutine might compute the Hessian matrix of the ind-th
	constraint.

	For achieving this objective YOU MAY MODIFY it according to your
	problem. To modify this subroutine IS NOT MANDATORY. See below
	where your modifications must be inserted.

	Parameters of the subroutine:

	On Entry:

	x		double precision x(n),
			 current point,

	ind	  integer,
			 index of the constraint whose Hessian will be computed,

	On Return

	hcnnz	integer,
			 number of perhaps-non-null elements of the computed
			 Hessian,

	hclin	integer hclin(hcnnz),
			 see below,

	hccol	integer hccol(hcnnz),
			 see below,

	hcval	double precision hcval(hcnnz),
			 the non-null value of the (hclin(k),hccol(k)) position
			 of the Hessian matrix of the ind-th constraint must
			 be saved at hcval(k). Just the lower triangular part of
			 Hessian matrix must be computed,

	flag	 integer,
			 You must set it to any number different of 0 (zero) if
			 some error ocurred during the evaluation of the Hessian
			 matrix of the ind-th constraint. (For example, trying
			 to compute the square root of a negative number,
			 dividing by zero or a very small number, etc.) If
			 everything was o.k. you must set it equal to zero.
	"""

	#n = len(x)

	#hcnnz = 1
	#because the number of constraints will be zero
	hcnnz = 0	
	'''
	hcnnz = 0

	hclin = zeros(hcnnz, int)
	hccol = zeros(hcnnz, int)
	hcval = zeros(hcnnz, float)
	'''
	global parameter_counter
	num_parameters = parameter_counter
	#print 'the length of x is %d'%parameter_counter
	#print 'x array is'
	#for i in x :
	#	print i
	#hnnz = 0
	hclin = zeros(hcnnz, int)
	hccol = zeros(hcnnz, int)
	hcval = zeros(hcnnz, float)

	'''	
	counter = 0
	for i in range (0,num_parameters) :
		for j in range (0,i+1) :
			hclin[counter]  = i
			hccol[counter] = j
			counter += 1
	print 'the value of the counter in the hessian of the constraints is %d'%counter
	print 'the value of hcnnz was %d'%hcnnz

	'''

	'''
	if ind == 1:
		hclin[0] = 0
		hccol[0] = 0
		hcval[0] = 2.0

		flag = 0

	elif ind == 2:
		hcnnz = 0

		flag = 0

	else:
		flag = -1
	'''
	flag = 0
	return hclin,hccol,hcval,hcnnz,flag

#   ******************************************************************
#   ******************************************************************

def evalfc(x,m):

	f = 0.0

	c = zeros(m)

	flag = -1

	return f,c,flag

#   ******************************************************************
#   ******************************************************************

def evalgjac(x,m):

	n = len(x)

	g = zeros(n)

	jcnnz = 0

	jcfun = zeros(jcnnz, int)
	jcvar = zeros(jcnnz, int)
	jcval = zeros(jcnnz, float)

	flag = -1

	return g,jcfun,jcvar,jcval,jcnnz,flag

#   ******************************************************************
#   ******************************************************************

def evalhl(x,m,lambda_,scalef,scalec):

	hlnnz = 0

	hllin = zeros(hlnnz, int)
	hlcol = zeros(hlnnz, int)
	hlval = zeros(hlnnz, float)

	flag = -1

	return hllin,hlcol,hlval,hlnnz,flag

#   ******************************************************************
#   ******************************************************************

def evalhlp(x,m,lambda_,scalef,scalec,p,goth):

	n = len(x)

	hp = zeros(n)

	flag = -1

	return hp,goth,flag

#   ******************************************************************
#   ******************************************************************

def endp(x,l,u,m,lambda_,equatn,linear):
	global initial_parameter_args, current_optimization_params
	"""
	This subroutine can be used to do some extra job.

	This subroutine can be used to do some extra job after the solver
	has found the solution, like some extra statistics, or to save the
	solution in some special format or to draw some graphical
	representation of the solution. If the information given by the
	solver is enough for you then leave the body of this subroutine
	empty.

	Parameters of the subroutine:

	The parameters of this subroutine are the same parameters of
	subroutine inip. But in this subroutine there are not output
	parameter. All the parameters are input parameters.
	"""
	print 'the value of x was '
	print x
	print 'the value of lambda was '
	print lambda_
	#also going to check if they sum to 1 or not
	'''
	for tag in initial_parameter_args.keys() :
		prob_mass_sum = 0.0
		if len(initial_parameter_args[tag].keys()) == 1 :
			continue
		for next_tag in initial_parameter_args[tag].keys() :
			parameter_number = 	initial_parameter_args[tag][next_tag]
			prob_mass_sum += x[parameter_number] 
			#print x[parameter_number]
		#if not prob_mass_sum == float(1) :
		#	print 'the sum of the probabilities corresponding to the tag %s was not one'%tag
		#	print 'it was %f'%prob_mass_suma
		print 'probability sum is %f'%prob_mass_sum
			#exit(0)
	'''
	global probabilities_language,probabilities_channel
	if current_optimization_params == 'tag_bigrams' :
		print 'we are replacing bigram probs'
		assignProbs(x,probabilities_language)
	#	checkZeros(probabilities_language)
	else :
		print 'we are replacing channel probs'
		assignProbs(x,probabilities_channel)
		#checkZeros(probabilities_channel)
	pass

#initializing the parameter vectors
def createParameters(probabilities,fractional_counts,free_parameters,alpha,sigma) :
	#first make a list of all the parametersa
	internal_parameter_counter = int(0)
	internal_num_constraints = int(0)
	global initial_parameter_args, constraint_tags_dict
	initial_parameter_args.clear()
	constraint_tags_dict.clear()
	#constraint_parameters.clear()
	#so, I have to take the current probabilities and then create the constraint parameters and the initial parameters
	#you should also calculate the current number of zero parameters and pass it to the argmax function
	current_num_zero_parameters = int(0)
	for tag in probabilities.keys() :
		temp_dict_initial_args = {}
		#temp_dict_constraint_params = {}
		for next_tag in probabilities[tag] :	
			#print probabilities[tag][next_tag]
			if not probabilities[tag][next_tag] == 0 :
			#	print 'this was true'
				temp_dict_initial_args[next_tag] = 0
				#temp_dict_constraint_params[next_tag] = 0
			else :
				current_num_zero_parameters += 1
		initial_parameter_args[tag] = temp_dict_initial_args
		#now, i can go over the initial parameter args and number them now itself
		if len(initial_parameter_args[tag]) == 1 :
			for next_tag in initial_parameter_args[tag] :
				initial_parameter_args[tag][next_tag] = 'ONE'
		else :
			internal_num_constraints += 1
			constraint_tags_dict[internal_num_constraints] = tag 
			#print 'the constraint id is %d'%internal_num_constraints	
			#print 'the tag is %s'%tag
			for next_tag in initial_parameter_args[tag] :
				initial_parameter_args[tag][next_tag] = internal_parameter_counter
				internal_parameter_counter += 1
	
	#print 'the current number of zero parameters is %d'%current_num_zero_parameters	
	

	#for tag in initial_parameter_args.keys() :
	#	if len(initial_parameter_args[tag].keys()) == 1 :
	#		for next_tag in initial_parameter_args[tag].keys() :
	#			del(initial_parameter_args[tag])#[next_tag])
	#		continue
	#	for next_tag in initial_parameter_args[tag].keys() :
	#		if constraint_parameters[tag][next_tag] == 1 :
	#			del(initial_parameter_args[tag][next_tag])
	#			continue
	#		initial_parameter_args[tag][next_tag] = parameter_counter
	#		parameter_counter += 1

	global parameter_counter, num_constraints
	parameter_counter = internal_parameter_counter
	num_constraints = internal_num_constraints

			

def main() :
	global fractional_counts_language ,fractional_counts_channel,probabilities_channel,sigma,alpha,current_fractional_counts,current_optimization_params,init_option

	num_iterations = int(sys.argv[1])
	alpha = float(sys.argv[2])
	sigma = float(sys.argv[3])

#	dictionary = emMethods.createDictionary('DICT2')#.small')
#	word_list = emMethods.readWordList('TEXT.linear')
	#word_list_five = emMethods.readWordList('TEXT.5.linear')
#	gold_tag_sequence = emMethods.readWordList('GOLD.linear')
	gold_cipher = emMethods.readCipherFile('cipher.gold.noq')
	print gold_cipher
	#dictionary = emMethods.createDictionary('complete.dict.new-formatted')#.small')
	#word_lines= emMethods.readWordLines('test.words.new-formatted')
	cipher_letter_dict = emMethods.getUniqCipherLetters('cipher.data.noq')
	#word_list_five = emMethods.readWordList('TEXT.3.linear')
	#plaintext = map(chr, range(97, 123))
	plaintext = []
	for k in range(65, 91):
		plaintext.append(chr(k))
	print plaintext
	print 'the number of unique cipher letter is %d'%len(cipher_letter_dict.keys())
	print cipher_letter_dict
	#gold_tag_sequence = emMethods.readWordList('test.tags.new-formatted.linear')
 		
	free_parameters_channel = {}
	free_parameters_language = {}
	print 'starting to create parameters'
	total_language_parameters = 0
	total_channel_parameters = 0
	#for line in cipher_lines :
		#print 'created parameters for a line'
		#(language_parameters,channel_parameters) = emMethods.getFreeParametersBigram(line,dictionary,free_parameters_language,free_parameters_channel)
	emMethods.getFreeCipherParametersChannel(cipher_letter_dict,plaintext,free_parameters_channel)
	temp = {'_':0.0}
	free_parameters_channel['_'] = temp
	#print free_parameters_channel
	#sys.exit()
		#print language_parameters
		#print channel_parameters

		#total_language_parameters += language_parameters
		#total_channel_parameters += channel_parameters
	#print 'total language parameters is %d' %(total_language_parameters)
	#print 'total channel parameters is %d' %(total_channel_parameters)
	#now, we will build all the lattices, and create a special start node and end node for every sentence
	start_node_end_node_list = []

#	print len(word_list)
#	num_taggings = emMethods.getNumTaggings(word_list,dictionary)
#	print 'num_taggings '
#	print type(num_taggings)
	#print num_taggings
	fractional_counts_channel = copy.deepcopy(free_parameters_channel)
	probabilities_channel = copy.deepcopy(free_parameters_channel)
	emMethods.initUniformProbs(probabilities_channel)
#	emMethods.initUniformProbs(probabilities_language,probabilities_channel)
	emMethods.writeFst('cipher.fst',probabilities_channel)


	run_training = r'./carmel.static --train-cascade -M 0 -m -HJ cipher.data cipher.wfsa cipher.fst'
#		skel_size += len(col)
	#running the EM iterations
	for i in range (0,num_iterations) :
		'''
		print 'checking the initial zeros inlanguage'
		checkZeros(probabilities_language)
	
		print 'checking the initial zeros in channel'
		checkZeros(probabilities_channel)
		'''
		#best_tag_sequence = emMethods.viterbiSearch(start,end,probabilities_channel,probabilities_language,lattice_skeleton)
		#emMethods.calcAccuracy(gold_tag_sequence,best_tag_sequence)
		#raw_input()
		#this will create the parameters
		total_corpus_probability = 0.0
		(status,output) = commands.getstatusoutput(run_training)	
		print 'we just ran the training'
		prob_match = probability_re.search(output)
		if prob_match == None :
			print'we should have found a probability'
		else :
			print 'the probability is %s'%prob_match.group(1)
		total_corpus_probability = float(prob_match.group(1)[2:len(prob_match.group(1))])
		print 'reading channel fractional counts'
		emMethods.readCarmelFractionalCounts('cipher.fst.trained',fractional_counts_channel,'channel')
		print 'read the fst'
			
		print' the probability of the corpus was %f' %total_corpus_probability
		print 'we are now checking the accuracies'
		noe_command = 'cat tagging.fsa | sed \'s/*e*//g\' > tagging.fsa.noe'
		(status,output) = commands.getstatusoutput(noe_command)
		print 'we wrote the noe fsa'
		viterbi_command = r'cat cipher.data | ./carmel.static -srbk -QEWI 1 cipher.wfsa.noe cipher.fst > decipherment_output'
		(status,output) = commands.getstatusoutput(viterbi_command)
		#tagged_sequence = emMethods.readTaggingOutput('tagging_output')	
		deciphered_sequence = emMethods.readCipherFile('decipherment_output')
		accuracy = emMethods.calcAccuracy(gold_cipher,deciphered_sequence)

		print 'The accuracy was %s and the objective function value was %s'%(str(accuracy),str(evaluateObjectiveFuncValue(total_corpus_probability,probabilities_channel,alpha,sigma)))


		#first optimizing the channel
		current_fractional_counts = fractional_counts_channel
		createParameters(probabilities_channel,current_fractional_counts,free_parameters_channel,alpha,sigma)
		temp_channel_probs = dict(probabilities_channel)	
		init_option = 'current_prob'
		current_optimization_params = 'channel'
		algencan.solvers(evalf,evalg,evalh,evalc,evaljac,evalhc,evalfc,evalgjac,evalhl,evalhlp,inip,endp)
	
		
		
		#I should check if the objective function is increasing 
		channel_probs_after_init_current_prob = copy.deepcopy(probabilities_channel)
		#obj_val1 = evaluateObjectiveFuncValue(total_corpus_probability,probabilities_language,probabilities_channel,alpha,sigma)
		total_obj_val1 = evaluateObjectiveFuncValue(total_corpus_probability,probabilities_channel,alpha,sigma)

		obj_val1 = evaluateOptimizationFunction(initial_parameter_args,probabilities_channel,fractional_counts_channel,alpha,sigma) 
		print 'the function value was obj 1 %f'%obj_val1
   		#emMethods.clearAlphaBeta(lattice_skeleton)


		init_option = 'zeros'
		current_optimization_params = 'channel'
		algencan.solvers(evalf,evalg,evalh,evalc,evaljac,evalhc,evalfc,evalgjac,evalhl,evalhlp,inip,endp)
	
		channel_probs_after_init_zeros = copy.deepcopy(probabilities_channel)
		#obj_val2 = evaluateObjectiveFuncValue(total_corpus_probability,probabilities_language,probabilities_channel,alpha,sigma)a
		
		total_obj_val2 = evaluateObjectiveFuncValue(total_corpus_probability,probabilities_channel,alpha,sigma)
		obj_val2 = evaluateOptimizationFunction(initial_parameter_args,probabilities_channel,fractional_counts_channel,alpha,sigma) 
		print 'the function value was obj 2 %f'%obj_val2
   		#emMethods.clearAlphaBeta(lattice_skeleton)

		if (total_obj_val1 >= total_obj_val2) :
			#init_option = 'current_prob'
			#current_optimization_params = 'tag_bigrams'
			if (obj_val1 < obj_val2) :
				print 'the final objective function value was opposite'

			probabilities_channel = copy.deepcopy(channel_probs_after_init_current_prob)
			#algencan.solvers(evalf,evalg,evalh,evalc,evaljac,evalhc,evalfc,evalgjac,evalhl,evalhlp,inip,endp)
	

			print 'the final objective function value was obj 1 %f'%total_obj_val1
		
		else :
			#init_option = 'zeros'
			#current_optimization_params = 'tag_bigrams'
			probabilities_channel = copy.deepcopy(channel_probs_after_init_zeros)
	
	
			if (obj_val2 < obj_val1) :
				print 'the final objective function value was opposite'

			print 'the final objective function value was obj 2 %f'%total_obj_val2
			#raw_input()

#		emMethods.reEstimateProbabilities(probabilities_channel,probabilities_language,fractional_counts_channel,fractional_counts_language)
		#emMethods.reEstimateProbabilities(probabilities_channel,probabilities_language,fractional_counts_channel,fractional_counts_language)
		#print 'writing the fsa'
		#now writing the fsa back again
		#emMethods.writeFsa('tagging.fsa',probabilities_language)
		emMethods.writeFst('cipher.fst',probabilities_channel)
		#print 'checking the zeros in tag bigram model'
		#checkZeros(probabilities_language)
	
		print 'checking the zeros in channel model'
		checkZeros(probabilities_channel)

		#fractional_counts_language.clear()
		#fractional_counts_channel.clear()
		#fractional_counts_language = copy.deepcopy(free_parameters_language)
		fractional_counts_channel = copy.deepcopy(free_parameters_channel)
		#probabilities_language = copy.deepcopy(free_parameters_language)
		#probabilities_channel = copy.deepcopy(free_parameters_channel)
		
		#sys.exit()
	#raw_input()

def emIteration(start_node,end_node,probabilities_channel,probabilities_language,fractional_counts_channel,fractional_counts_language) :

	emMethods.forwardPass(start_node,end_node,probabilities_channel,probabilities_language)
	emMethods.backwardPass(start_node,end_node,probabilities_channel,probabilities_language)
	emMethods.getFractionalCounts(start_node,end_node,probabilities_channel,probabilities_language,fractional_counts_channel,fractional_counts_language)

'''
def emIteration(start,end,probabilities_channel,probabilities_language,fractional_counts_channel,fractional_counts_language,lattice_skeleton) :

	emMethods.forwardPass(start,probabilities_channel,probabilities_language,lattice_skeleton)
#	raw_input()
	emMethods.backwardPass(end,probabilities_channel,probabilities_language,lattice_skeleton)
	emMethods.getFracCounts(start,end,probabilities_channel,probabilities_language,fractional_counts_channel,fractional_counts_language,lattice_skeleton) 
'''

def assignProbs(x,probabilities) :
	global initial_parameter_args
	for tag in initial_parameter_args.keys() :
		for next_tag in initial_parameter_args[tag].keys() :
			if initial_parameter_args[tag][next_tag] == 'ONE' :
				break
			else :
				parameter_number = 	initial_parameter_args[tag][next_tag]
				probability = x[parameter_number]
				probabilities[tag][next_tag] = probability

def checkZeros(probabilities) :
	global initial_parameter_args
	num_zeros = 0
	non_zeros = 0
	for tag in probabilities.keys() :
		for next_tag in probabilities[tag].keys() :
			if probabilities[tag][next_tag] <= 10E-4 :
				num_zeros += 1	
			else :
				non_zeros += 1
	
	print 'num zeros is %d' %num_zeros
	print 'num non zeros is %d'%non_zeros	



def	evaluateObjectiveFuncValue(corpus_prob,probabilities_channel,alpha,sigma) :
	#corpus_prob = 
	l0_norm_terms  =float (0)
	'''
	for tag in probabilities_language.keys() :
		for next_tag in probabilities_language[tag].keys() :
			l0_norm_terms +=  alpha * math.exp(-probabilities_language[tag][next_tag]/sigma)
	'''
	for tag in probabilities_channel.keys() :
		for next_tag in probabilities_channel[tag].keys() :
			l0_norm_terms +=  alpha * math.exp(-probabilities_channel[tag][next_tag]/sigma)
	
	objective_function_val = l0_norm_terms + corpus_prob
	print 'the objective func value %f'%objective_function_val
	return objective_function_val
	#raw_input()

def	evaluateOptimizationFunction(initial_parameter_args,probabilities_language,fractional_counts_language,alpha,sigma) :
	#corpus_prob = 
	func_val = float(0)
	for tag in initial_parameter_args.keys() :
		for next_tag in initial_parameter_args[tag].keys() :
			if initial_parameter_args[tag][next_tag] == 'ONE' :
				continue
			else :
				func_val += fractional_counts_language[tag][next_tag]  * math.log(probabilities_language[tag][next_tag])+ alpha * math.exp(-probabilities_language[tag][next_tag]/sigma)

	print 'the optimization func value %f'%func_val
	return func_val


#checking if the sum to one constraints have been satisfied or not
def checkConstraints (probabilities_language,probabilities_channel) :
	print 'we are checking the constraints for the tag bigram parameters'
	for tag in probabilities_language.keys() :
		sum = float(0)
		for next_tag in probabilities_language[tag].keys() :
			sum += probabilities_language[tag][next_tag]
		print 'the sum for the constraints corresponding to tag %s was %f'%(tag,sum)
	print 'we are checking the constraints for the word given tag parameters'
	for tag in probabilities_channel.keys() :
		sum = float(0)
		for word in probabilities_channel[tag].keys() :
			sum += probabilities_channel[tag][word]
		print 'the sum for the constraints corresponding to tag %s was %f'%(tag,sum)

def createRandomPoint(probabilities_channel) :
	for tag in probabilities_channel.keys() :
                if len(probabilities_channel[tag].keys()) == 1 :
                        for word in probabilities_channel[tag].keys() :
                                probabilities_channel[tag][word] = 1
                else :
                        random_number_dict= {}
                        random_number_sum = 0
                        for word in probabilities_channel[tag].keys() :
                                random_number = random.randint(1,100000)
                                random_number_sum += random_number
                                random_number_dict[word] = random_number

                        for word in probabilities_channel[tag].keys() :
                                probabilities_channel[tag][word] = float(random_number_dict[word])/random_number_sum


main()
