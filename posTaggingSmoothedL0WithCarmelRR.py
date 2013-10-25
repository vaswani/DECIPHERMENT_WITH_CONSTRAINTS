#!/usr/bin/env python

#the change in this one is going to be that we are going to include both tag bigrams and tags
#I have to change the optimization here. I will only optimze every set of conditional probabilities on its own. So, there will be 44 optimizations in each iteration
#i dump theta

from numpy import *
#from scipy import *
#import algencan
import commands
import re
from math import e
import sys
import copy
import time
from optparse import OptionParser


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
current_initial_parameter_args = {}
parameter_counter = 0
num_constraints = 1
constraint_tags_dict = {} #this will hold the tag for the constraint. the key is the constraint id and the value is the tag corresponding to this constraint
beta = float(0)
alpha = float(0)
global_num_parameters = 0
init_option = ''
current_optimization_tag = '' #this will hold the of the constraint for which we are doing the optimization



#initializing the parameter vectors. Here, each parameters is given a number which is the index in the parameter vector that algencan optimizes. The constraints are also given a number . We need these constraint numbers for algencan
#the new technique here does the optimization for each constraint separately. Therefore,it makes sense to prepare these in advance and then send them off to algencan when they are needed
#here, we will create the indexing of the parameters needed by algencan so that we can scale it. Instead of throwing the entire problem to algencan, we will #also, I will call this just once because the only code that is affected by being called in different places is the one that checks if the probability is zero. Since we will not make any probabilities zero in algencan, we just need to call this once
def createParametersForScaling(probabilities) :
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
		for next_tag,prob in sorted(probabilities[tag].iteritems()) :	
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
			internal_parameter_counter = 0
			internal_num_constraints += 1
			#constraint_tags_dict[internal_num_constraints] = tag 
			#print 'the constraint id is %d'%internal_num_constraints	
			#print 'the tag is %s'%tag
			for next_tag in initial_parameter_args[tag] :
				initial_parameter_args[tag][next_tag] = internal_parameter_counter
				internal_parameter_counter += 1


def main() :
	global fractional_counts_language ,fractional_counts_channel,probabilities_channel,probabilities_language,beta,alpha,current_fractional_counts,current_optimization_params,init_option,initial_parameter_args,constraint_tags_dict,parameter_counter,current_optimization_tag,slack_option
	
	#adding parser options


	parser = OptionParser()
	parser.add_option("--num_iter", action="store", type="int", dest="num_iterations",default=100,help="number of iterations you would like to run em+smoothed l0. Default is 50")
	parser.add_option("--alpha", action="store", type="float", dest="alpha",default = 0.0,help="alpha,the weight of the smoothed l0 penalty. Defaut is 0 ")
	parser.add_option("--beta", action="store", type="float", dest="beta",default = 0.5,help="beta, the smoothness of the l0 prior, smaller the sigma, closer the approximation to true L0. Default is 0.5")
	parser.add_option("--slack", action="store_true", dest="slack_option",default = False,help="if you want to project on the simplex with slack.")
	parser.add_option("--noslack", action="store_false", dest="slack_option",default = False, help="if you want to project on the simplex with no slack. This is the regular projection approach")
	parser.add_option("--num_pgd_iterations", action="store", type="int", dest="num_pgd_iterations",default = 50,help="Number of Projected Gradient Descent to run. Default is 100")
	parser.add_option("--eta", action="store", type="float", dest="eta",default = 0.5,help="Eta, the constant step size for PGD using armijo line search. Default is 0.1")
	parser.add_option("--armijo_beta", action="store", type="float", dest="armijo_beta",default = 0.5,help="Set value for Armijo beta, the beta used in in armijo line search. Default value is 0.2")
	parser.add_option("--armijo_sigma", action="store", type="float", dest="armijo_sigma",default = 0.5,help="Set value for Armijo sigma, the sigma used in in armijo line search. Lower bound is 0.0001")
	parser.add_option("--lower_bound", action="store", type="float", dest="lower_bound",default = 0.000001,help="Set value for the lower bound on the probability.Default is 10E-6")
	(options, args) = parser.parse_args()

	#print options
	#print args
	
	#getting the values from optparse
	num_iterations = options.num_iterations
	alpha = options.alpha
	beta = options.beta
	slack_option = options.slack_option
	num_pgd_iterations = options.num_pgd_iterations
	eta = options.eta
	armijo_beta = options.armijo_beta
	armijo_sigma = options.armijo_sigma
	lower_bound = options.lower_bound

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

	createRandomPoint(probabilities_channel)
	#emMethods.initUniformProbs(probabilities_channel)
#	emMethods.initUniformProbs(probabilities_language,probabilities_channel)
	emMethods.writeFst('cipher.fst',probabilities_channel)

	final_probabilities_channel = copy.deepcopy(free_parameters_channel)

	run_training = r'./carmel --train-cascade -M 0 -m -HJ cipher.data cipher.wfsa cipher.fst'

	#running the EM iterations
	#we are creating the indexes for algencan . Notice that here, the probabilities language is already uniform and therefore none of them will be zero
	createParametersForScaling(probabilities_channel)
	start_time = time.clock()
	print 'start time was ',start_time
	#fractional_counts_dump_file = open('fractional.params','w')
	#probabilities_dump_file = open('probs.params','w')
	#optimized_probabilities_dump_file = open('probs.optimized.params','w')
	for i in range (0,num_iterations) :
		print 'the iteration number was ',i
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
		print output
		prob_match = probability_re.search(output)
		if prob_match == None :
			print'we should have found a probability'
		else :
			print 'the probability is %s'%prob_match.group(1)
		temp_corpus_probability = float(prob_match.group(1)[2:len(prob_match.group(1))])
		total_corpus_probability = 0.693147181 * temp_corpus_probability
		print 'reading language fractional counts'
	
		emMethods.readCarmelFractionalCounts('cipher.fst.trained',fractional_counts_channel,'channel')
		print 'read the fst'
			
		print' the probability of the corpus was %f' %total_corpus_probability
		print 'we are now checking the accuracies'
		noe_command = 'cat tagging.fsa | sed \'s/*e*//g\' > tagging.fsa.noe'
		(status,output) = commands.getstatusoutput(noe_command)
		print 'we wrote the noe fsa'
		viterbi_command = r'cat cipher.data.quotes | ./carmel -srbk 1 -QEWI cipher.wfsa.noe cipher.fst > decipherment_output'
		(status,output) = commands.getstatusoutput(viterbi_command)
		#tagged_sequence = emMethods.readTaggingOutput('tagging_output')	
		deciphered_sequence = emMethods.readCipherFile('decipherment_output')
		accuracy = emMethods.calcAccuracy(gold_cipher,deciphered_sequence)

		print 'The accuracy was %s and the objective function value was %s'%(str(accuracy),str(evaluateObjectiveFuncValue(total_corpus_probability,probabilities_channel,probabilities_language,alpha,beta)))
			

		#emMethods.reEstimateProbabilities(probabilities_channel,probabilities_language,fractional_counts_channel,fractional_counts_language)
		#fractional_counts_language = copy.deepcopy(free_parameters_language)
		#fractional_counts_channel = copy.deepcopy(free_parameters_channel)

		#first optimizing the tag bigrams and doing it per parameter
		for tag in initial_parameter_args.keys() :
			print 'we are starting a new tag'

			if len(initial_parameter_args[tag].keys()) == 1 :
				continue
			print 'we are currently optimizing for tag', tag
			#current_initial_parameter_args = initial_parameter_args[tag]
			current_optimization_tag = tag


			parameter_counter = len(initial_parameter_args[tag].keys())			
			current_fractional_counts = fractional_counts_channel
			constraint_tags_dict[1] = tag
			temp_language_probs = dict(probabilities_channel)
			#optimizing per constraint
			init_option = 'current_prob'
			current_optimization_params = 'channel'
			#creating the initial vector to pass to the projected gradient descent algorithm

			x=zeros(parameter_counter)
			expected_counts=zeros(parameter_counter)
			expected_counts_sum = 0.
			for next_tag in initial_parameter_args[current_optimization_tag].keys() :
				if initial_parameter_args[current_optimization_tag][next_tag] ==  'ONE' :
					continue
				else :
					parameter_number = initial_parameter_args[current_optimization_tag][next_tag]
					x[parameter_number] = probabilities_channel[current_optimization_tag][next_tag]	
					expected_counts[parameter_number] = current_fractional_counts[current_optimization_tag][next_tag]
					expected_counts_sum +=  current_fractional_counts[current_optimization_tag][next_tag]
					
			print 'expected counts sum was ',expected_counts_sum
			if expected_counts_sum <= 0.5 :
				#set all the probabilities to 0
				#for next_tag in initial_parameter_args[current_optimization_tag] :
				#	new_probabilities = zeros(len(initial_parameter_args[current_optimization_tag].keys()))
				#	assignProbs(new_probabilities,probabilities_channel)
				continue
			#raw_input()			
			#the optimization steps have to come here
			#print 'we are optimizing the tag bigrams now ' 

			#current_eta = eta_0 #/sqrt(num_iterations+1)
			print 'Doing projected gradient descent'
			new_probabilities = projectedGradientDescentWithArmijoRule(x = x,expected_counts = expected_counts,num_pgd_iterations = num_pgd_iterations,eta = eta,lower_bound = lower_bound,armijo_beta = armijo_beta,armijo_sigma = armijo_sigma)
			print 'finished projected gradient descent'
			#new_probabilities = newtonProjectedGradientDescentWithArmijoRule(x,num_pgd_iterations,eta,lower_bound,armijo_beta,armijo_sigma)
			#print new_probabilities
			#raw_input()
			if current_optimization_params == 'tag_bigrams' :
				#print 'we are replacing bigram probs'
				assignProbs(new_probabilities,probabilities_language)
		#	checkZeros(probabilities_language)
			else :
				print 'we are replacing channel probs'
				assignProbs(new_probabilities,probabilities_channel)
	
			
		#now writing the fsa back again
#		emMethods.writeFsa('tagging.fsa',probabilities_language)
#		emMethods.writeFst('tagging.fst',probabilities_channel)

		emMethods.writeFst('cipher.fst',probabilities_channel)
		#emMethods.writeFstLnFormat('tagging.fst',probabilities_channel)


		print 'checking the initial zeros in channel model'
		checkZeros(probabilities_channel)
	
		#fractional_counts_language.clear()
		#fractional_counts_channel.clear()

		fractional_counts_channel = copy.deepcopy(free_parameters_channel)
		final_probabilities_channel = copy.deepcopy(probabilities_channel)
		print 'at the end of the iteration'
	#raw_input()	

	elapsed_time = time.clock() - start_time
	print 'the elapsed time was ',elapsed_time

	#print 'the value of the optimization function was ',evaluateOptimizationFunction(initial_parameter_args,final_probabilities_language,final_fractional_counts_language,alpha,beta) 
	


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



def	evaluateObjectiveFuncValue(corpus_prob,probabilities_language,probabilities_channel,alpha,sigma) :
	#corpus_prob = 
	l0_norm_terms  =float (0)
	for tag in probabilities_language.keys() :
		for next_tag in probabilities_language[tag].keys() :
			l0_norm_terms +=  alpha * math.exp(-probabilities_language[tag][next_tag]/sigma)

	#for tag in probabilities_channel.keys() :
	#	for next_tag in probabilities_channel[tag].keys() :
	#		l0_norm_terms +=  alpha * math.exp(-probabilities_channel[tag][next_tag]/sigma)
	
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

#this goes through the tag bigram sequence and checks the number of tag bigram types
def	getTagBigramTypes(tagged_sequence) :
	tag_bigram_types = {}
	for i in range(1,len(tagged_sequence)) :
		tag_bigram = '%s %s'%(tagged_sequence[i],tagged_sequence[i-1])
		if not tag_bigram_types.has_key(tag_bigram) :
			tag_bigram_types[tag_bigram] = 1
	print 'The number of tag bigram types in the viterbi tag sequence were %d'%len(tag_bigram_types.keys())

def approximateProjectOntoSimplex(v) :
	array_size = v.size
	projection = array(v)
	sum =0.
	for i in range(array_size) :
		if projection[i] < 0 :
			projection[i] = 0.
		else :
			sum += projection[i]
	if sum > 1. : #the project
		return(projectOntoSimplex(projection))
	else :
		return(projection)

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

def projectOntoSimplexWithLowerBound(v,lower_bound) :
	global parameter_counter
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
		#print 'temp rho',temp_rho
		#raw_input()
		if temp_rho > 0 :
			max_index = i
			max_index_sum = sum
		
	theta = 1/float(max_index+1) * (max_index_sum -1)

	#print 'theta is ',theta			
	final_vector = array([getMax(value-theta,0) for value in v])
	num_zero_elements = 0
	for value in final_vector :
		if value == 0 :
			num_zero_elements += 1
	total_zero_mass = lower_bound * num_zero_elements
	mass_to_remove = total_zero_mass/(final_vector.size - num_zero_elements)
	for index in range(0,parameter_counter) :
		if final_vector[index] == 0:
			final_vector[index] = lower_bound
		else :
			final_vector[index] -= mass_to_remove
	#print 'the projected point is' 
	#print final_vector 
	return(final_vector)


def getMax(value,theta) :
	return(max(value-theta,0))

'''
			try :
				func_value += current_fractional_counts[current_optimization_tag][next_tag]  * math.log(x[parameter_number]) + alpha * math.exp(-x[parameter_number]/sigma)
			except ValueError :
				print 'it was a value error in the func value function the value of the current parameter was ',x[parameter_number]
'''

#we're doing projected gradient descent with the armijo rule along feasible direction. Bertsekas: NonLinear programming. page 230
def projectedGradientDescentWithArmijoRule(x,expected_counts,num_pgd_iterations,eta,lower_bound,armijo_beta,armijo_sigma) :
	
	global parameter_counter,slack_option
	#current_point = array([0.0]*parameter_counter)
	current_point = array(x)
	#print 'the function value at the current point  is ',current_function_value
	for time in range(1,num_pgd_iterations +1) :
		current_function_value= evalFunction(current_point=current_point,expected_counts=expected_counts)
		#print 'the function value at the current point  is ',current_function_value

		grad = evalGradient(current_point=current_point,expected_counts=expected_counts) #getting the gradient at the current point
		#print grad
		new_point = array([0.0]*x.size)
		'''
		for index in range (0,parameter_counter) :
			new_point[index]= current_point[index]-eta*grad[index]
		'''
		new_point = current_point - eta*grad
		#print new_point
		new_feasible_point = array([])
		#current_function_value = evalFunction()
		if slack_option == True :
			new_feasible_point = projectOntoSimplexWithLowerBound(new_point,lower_bound)
		else :
			new_feasible_point = projectOntoSimplex(new_point)
		#print current_point
		armijo_bound = 0.0
		'''
		for index in range (0,parameter_counter) :
			bound_term = armijo_sigma * armijo_beta * grad[index] * ( new_feasible_point[index] - current_point[index])
			#print 'the bound term is ',bound_term
			armijo_bound -= bound_term
			print 'temp armijo bound ',armijo_bound	
		'''
		armijo_bound = - armijo_sigma * armijo_beta * dot(grad,(new_feasible_point-current_point))
		#raw_input()
		terminate_line_srch = False
		num_steps = 1
		#current_beta = 1.0 #armijo_beta
		current_beta = armijo_beta
		final_beta = 0
		current_armijo_bound = armijo_bound
		no_update = True 
		best_func_value = current_function_value
		#print 'the current function value is %.16f and the current tag is %s'%(current_function_value,current_optimization_tag)
		while(terminate_line_srch != True) :
			#print 'num steps is ',num_steps
			temp_point = array([0.0]*x.size)
			#for index in range (0,parameter_counter) :
			#	temp_point[index]= (1.0 - current_beta) * current_point[index] + current_beta * new_feasible_point[index]
			temp_point = current_point * (1.0 - current_beta) + current_beta * new_feasible_point
			#print 'the temp point is '
			#print temp_point
			func_value_at_temp_point = evalFunction(current_point= temp_point,expected_counts = expected_counts)

			#print 'the function value at the current point  is %.16f'%current_function_value
			#raw_input()
			if func_value_at_temp_point < best_func_value :
				best_func_value = func_value_at_temp_point
				final_beta = current_beta
				no_update = False
				#print 'we arrived at a better function value'
				#raw_input()
			#	print 'we just updated thef final beta to ',final_beta
			if (current_function_value - func_value_at_temp_point >= current_armijo_bound) :
				terminate_line_srch = True
			#elif current_function_value - func_value_at_temp_point < 0:
			#	terminate_line_srch = True
			if num_steps >= 20 :
				terminate_line_srch = True
			
			current_beta = armijo_beta * current_beta
			current_armijo_bound =  armijo_bound * current_beta
			num_steps += 1
			
		#next_point = array([0.0]*x.size)
		#else :
		#print 'the number of steps was ',num_steps
		if no_update == False :
			for index in range (0,parameter_counter) :
			#	next_point[index] = (1.0 - final_beta) * current_point[index] + final_beta * new_feasible_point[index]
				temp_coordinate = (1.0 - final_beta) * current_point[index] + final_beta * new_feasible_point[index]
				#if next_point[index] == 0 :
				if temp_coordinate == 0 :
					print 'next point was 0 and final beta was ',final_beta, ' and the current point was ',current_point[index], ' and the feasible point was' ,new_feasible_point[index]
				current_point[index] = temp_coordinate
		
		if no_update == True :#x.all() == current_point.all() :
			print 'not update was true'
			break;
	return(current_point)

def evalFunction(expected_counts,current_point) :
	global alpha,beta
	func_value = float(0.0)
	for i in range(len(expected_counts)) :
		if current_point[i] == 0.0 and expected_counts[i] != 0 :
			print 'the probability in position ',i,' was 0 and the fractional count was not 0'
			exit(1)

		if current_point[i] != 0.0 :
			func_value += expected_counts[i] * math.log(current_point[i]) + alpha * math.exp(-current_point[i]/beta)
		else :
			func_value += alpha * math.exp(-current_point[i]/beta);
	
	return(-func_value)

def evalGradient(expected_counts,current_point) :
	global alpha,beta
	#print alpha
	#print beta
	#print expected_counts
	#print current_point
	gradient = zeros(len(current_point))
	for i in range(len(expected_counts)) :
		if current_point[i] == 0.0 and expected_counts[i] != 0 :
			print 'the probability in position ',i,' was 0 and the fractional count was not 0'
			exit(1)
		if current_point[i] != 0.0 :
			gradient[i] = -1 *(expected_counts[i]/current_point[i] - alpha * exp(-current_point[i]/beta)/beta)
			#print 'value of reg term is ', alpha * exp(-current_point[i]/beta)/beta
		else :
			gradient[i] = -1 * (-alpha * exp(-current_point[i]/beta)/beta) 	

	return(gradient)


def assignProbs(x,probabilities) :
	global initial_parameter_args,current_optimization_tag
	for next_tag in initial_parameter_args[current_optimization_tag].keys() :
		if initial_parameter_args[current_optimization_tag][next_tag] == 'ONE' :
			break
		else :
			parameter_number = 	initial_parameter_args[current_optimization_tag][next_tag]
			probability = x[parameter_number]
			probabilities[current_optimization_tag][next_tag] = probability

def createRandomPoint(probabilities_channel) :
	for tag in probabilities_channel.keys() :
		if len(probabilities_channel[tag].keys()) == 1 :
			for word in probabilities_channel[tag].keys(): 
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


if __name__ == "__main__" :
	main()
