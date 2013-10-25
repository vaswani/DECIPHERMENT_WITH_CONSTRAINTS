import math

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
	#for tag in probabilities_language.keys() :
	#	for next_tag in probabilities_language[tag].keys() :
	#		l0_norm_terms +=  alpha * math.exp(-probabilities_language[tag][next_tag]/sigma)

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

#this goes through the tag bigram sequence and checks the number of tag bigram types
def	getTagBigramTypes(tagged_sequence) :
	tag_bigram_types = {}
	for i in range(1,len(tagged_sequence)) :
		tag_bigram = '%s %s'%(tagged_sequence[i],tagged_sequence[i-1])
		if not tag_bigram_types.has_key(tag_bigram) :
			tag_bigram_types[tag_bigram] = 1
	print 'The number of tag bigram types in the viterbi tag sequence were %d'%len(tag_bigram_types.keys())

