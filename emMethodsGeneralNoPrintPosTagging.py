#these em methods will apply for a general lattice and do forward backward on a general lattice
# the bigram separation charachter is >
import latticeNode
import random
import math
from math import e
import copy
import re

negative_infinity = float(-999999)
#(0 (0 "NNP" "Institut" 0.999999999999943))
#(0 (: *e* ":" 4.00000000000004))
carmel_wfst_re = re.compile(r'\((.+) \((.+) "(.+)" "(.+)" (.+)\)\)')
carmel_wfsa_re = re.compile(r'\((.+) \((.+) \*e\* (.+) (.+)\)\)')

#writing the fsa for the next iteration

def writeFsa(fsa_file,probabilities_language) :
	fsa = open(fsa_file,'w')
	fsa.write('F\n');
	#first writing the transitions from the start state
	for next_tag in probabilities_language['&start&'].keys() :
		line  = "(0 (%s *e* \"%s\" %s))\n"%(next_tag,next_tag,str(probabilities_language['&start&'][next_tag]))
		fsa.write(line)
		#print line
	for tag in probabilities_language.keys() :
		if tag == '&start&' :
			continue
		for next_tag in probabilities_language[tag].keys() :
			transition_label = "\"%s\""%next_tag
			next_state  = next_tag
			if next_tag == r'&end&' :
				transition_label = '*e*'
				next_state = 'F'
			line = "(%s (%s *e* %s %s))\n"%(tag,next_state,transition_label,str(probabilities_language[tag][next_tag]))
			fsa.write(line)
			#print line
			#raw_input()
	fsa.close()

#it will write the FSA in ln format 
def writeFsaLnFormat(fsa_file,probabilities_language) :
	fsa = open(fsa_file,'w')
	fsa.write('F\n');
	#first writing the transitions from the start state
	for next_tag in probabilities_language['&start&'].keys() :
		try :
			log_base_e = math.log(probabilities_language['&start&'][next_tag])
		except :
			print 'It was a value error in writing the fsa and the probability was ', probabilities_language['&start&'][next_tag],' tag1 was start and tag2 was ',next_tag

		line  = "(0 (%s *e* \"%s\" e^%s))\n"%(next_tag,next_tag,str(log_base_e))
		fsa.write(line)
		#print line
	for tag in probabilities_language.keys() :
		if tag == '&start&' :
			continue
		for next_tag in probabilities_language[tag].keys() :
			transition_label = "\"%s\""%next_tag
			next_state  = next_tag
			if next_tag == r'&end&' :
				transition_label = '*e*'
				next_state = 'F'
			try :
				log_base_e = math.log(probabilities_language[tag][next_tag])
			except ValueError :
				print 'It was a value error in writing the fsa and the probability was ', probabilities_language[tag][next_tag],' tag1 was ',tag,' and tag2 was ',next_tag,' and log base e was ',log_base_e
			line = "(%s (%s *e* %s e^%s))\n"%(tag,next_state,transition_label,str(log_base_e))
			fsa.write(line)
			#print line
			#raw_input()
	fsa.close()

			
def writeFst(fst_file,probabilities_channel) :
	fst = open(fst_file,'w')
	fst.write('0\n');
	#first writing the transitions from the start state
	for tag in probabilities_channel.keys() :
		for word in probabilities_channel[tag].keys() :
			line = "(0 (0 \"%s\" \"%s\" %s))\n"%(tag,word,str(probabilities_channel[tag][word]))
			fst.write(line)
			#print line
			#raw_input()
	fst.close()

			
def writeFstLnFormat(fst_file,probabilities_channel) :
	fst = open(fst_file,'w')
	fst.write('0\n');
	#first writing the transitions from the start state
	for tag in probabilities_channel.keys() :
		for word in probabilities_channel[tag].keys() :
			try :
				log_base_e = math.log(probabilities_channel[tag][word])
			except ValueError:
				print 'it was a value error in the fst and the probability was ',probabilities_channel[tag][word], ' and the tag was ',tag,' and the word was ',word,' and log base e was ',log_base_e
				log_base_e = -800 #this is more or less 0 and I'm keeping it fixed at taht for now
			line = "(0 (0 \"%s\" \"%s\" e^%s))\n"%(tag,word,str(log_base_e))
			fst.write(line)
			#print line
			#raw_input()
	fst.close()


#this reads the carmel file that carries the fractional counts and then stores them in the dictionary. right now this is only for bigrams
def readCarmelFractionalCounts(carmel_file,dictionary,parameters) :
	counter=0
	file = open(carmel_file)
	dump_line = file.readline()
	for line in file :
		#print line
		match = ''
		if parameters == 'bigram' :
			match = carmel_wfsa_re.match(line)
		else :
			match = carmel_wfst_re.match(line)
		#print "the match is %s"%match
		if match == None :
			print 'the line did not match while reading the fractional counts'		
			print line
			exit(0)
		else :
			if parameters == 'bigram' :
				tag = match.group(1)
				if tag == '0' :
					tag = '&start&'
				
				next_tag = match.group(2)
				#print 'the next tag is %s'%next_tag			
				#if next_tag == r'*e*' :
				#	next_tag = r'&end&'
				#else :
				#	next_tag = next_tag[1:len(next_tag)-1]
				#print 'temp tag is %s'%temp_tag
				#raw_input()
				if next_tag == 'F' :
					next_tag = '&end&'
				fractional_count = float(0)
				if match.group(4)[0] == 'e' :
					fractional_count = eval(match.group(4).replace('^','**'))
				else :
					fractional_count = float(match.group(4))
				#print 'the fractional count was ',fractional_count
				if not dictionary.has_key(tag) :
					print 'the language dictionary did not have the tag %s'%tag
					exit(0)
				elif not dictionary[tag].has_key(next_tag) :
					print 'the language dictionary did not have the next tag %s'%next_tag
					exit(0)
				dictionary[tag][next_tag] = fractional_count
			else :
				tag = match.group(3)
				word = match.group(4)
#				fractional_count = float(match.group(5))
				fractional_count = float(0)
				#print "\"%s\""%match.group(5)
				if match.group(5)[0] == 'e' :
					fractional_count = eval(match.group(5).replace('^','**'))
				else :
					fractional_count = float(match.group(5))

				#print 'the fractional count was ',fractional_count
				if not dictionary.has_key(tag) :
					print 'the channel dictionary did not have the tag %s'%tag
				elif not dictionary[tag].has_key(word) :
					print 'the channel dictionary did not have the word %s'%word
				dictionary[tag][word] = fractional_count
			#print 'the line matched'
			#print match.group(1)
			#print match.group(2)
			#print match.group(3)
			#print match.group(4)
			#if parameters == 'channel' :
			#	print match.group(5)
		#raw_input()
	return(counter)

#this reads the carmel file that carries the fractional counts and then stores them in the dictionary. right now this is only for bigrams
def readCarmelFile(carmel_file,dictionary,parameters) :
	counter=0
	file = open(carmel_file)
	dump_line = file.readline()
	for line in file :
		#print line
		match = ''
		if parameters == 'bigram' :
			match = carmel_wfsa_re.match(line)
		else :
			match = carmel_wfst_re.match(line)
		#print "the match is %s"%match
		if match == None :
			print 'the line did not match'		
			exit(0)
		else :
			if parameters == 'bigram' :
				tag = match.group(1)
				if tag == '0' :
					tag = '&start&'
				
				next_tag = match.group(2)
				#print 'the next tag is %s'%next_tag			
				#if next_tag == r'*e*' :
				#	next_tag = r'&end&'
				#else :
				#	next_tag = next_tag[1:len(next_tag)-1]
				#print 'temp tag is %s'%temp_tag
				#raw_input()
				if next_tag == 'F' :
					next_tag = '&end&'
				fractional_count = float(match.group(4))
				if not dictionary.has_key(tag) :
					temp = {}
					temp[next_tag] = fractional_count
					dictionary[tag] = temp
					counter += 1
				elif not dictionary[tag].has_key(next_tag) :
					dictionary[tag][next_tag] = fractional_count
					counter += 1
			else :
				tag = match.group(3)
				word = match.group(4)
				fractional_count = float(match.group(5))
				if not dictionary.has_key(tag) :
					temp = {}
					temp[word] = fractional_count
					dictionary[tag] = temp
					counter += 1
				elif not dictionary[tag].has_key(word) :
					dictionary[tag][word] = fractional_count
					counter += 1
			#print 'the line matched'
			#print match.group(1)
			#print match.group(2)
			#print match.group(3)
			#print match.group(4)
			#if parameters == 'channel' :
			#	print match.group(5)
		#raw_input()
	return(counter)

def createDictionary(dictionary_file) :
	dictionary = {}
	dict_file_handle = open(dictionary_file,'r')
	for line in dict_file_handle :
		line = line.strip()
		list = re.split(r'\t',line)
		word = list[0]
		tag_list = []
		for i in range(1,len(list)) :
			tag_list.append(list[i])
	
		dictionary[word] = tag_list	
	return (dictionary)
	#advertising NN VBG 

def readWordList (word_file) :
	word_list = []
	word_file_handle = open(word_file,'r')
	for line in word_file_handle :
		line = line.strip()
		#words = line.split('\t')
		word_list.append(line)
	return (word_list)

def readWordLines (word_file) :
	word_list = []
	word_file_handle = open(word_file,'r')
	for line in word_file_handle :
		line = line.strip()
		words = line.split(' ')
		word_list.append(words)
	return (word_list)


def getNumTaggings(word_list,dictionary) :
	lattice_size = int(0)
	num_taggings = int(1)
	for word in word_list :
		if dictionary.has_key(word) :
			num_taggings *= len(dictionary[word])
			lattice_size += len(dictionary[word])

		else :
			print word
	print 'the lattice size is %d' %lattice_size
	return(num_taggings)		

def makeLatticeBigram(start_node,end_node,word_list,dictionary,language_model_parameters,channel_model_parameters) :
	previous_column = [] #this will include all the nodes that have to be connected to the 
	previous_column.append(start_node)
	prev_word_tag_list = ['&start&']
	for word in word_list :
	#	print 'the word is %s'%word
	#	print 'the number of tags this word has is %d'%len(dictionary[word])
		next_iteration_previous_column = []
	#	print 'the length of the previous word tag list is %d'%len(prev_word_tag_list)
		for tag in dictionary[word] :
			#i have to create a node for every tag 
			edge_list = [] #this is the edge list that will connect to the node that is created for this tag for this word
			for previous_node in previous_column : #then create the edge from every node in the previous column to this node. each edge corresponds to p(a|b,c)
			#	if  len(previous_column) == 1 :
			#		fixed = True
				#now create the language parameter that is to be this edge
				parameter_id = "%s|%s"%(tag,previous_node.parameter_id) #creating the n gram probability parameter id
				if language_model_parameters.has_key(previous_node.parameter_id) :
					if not language_model_parameters[previous_node.parameter_id].has_key(tag) :
						language_model_parameters[previous_node.parameter_id][tag] =float(0)
				else :
					temp_dict = {}
					temp_dict[tag] =  float(0)
					language_model_parameters[previous_node.parameter_id] = temp_dict
				language_edge = latticeNode.edge(False,'language',parameter_id)
				language_edge.head = previous_node
				previous_node.forward_pointers.append(language_edge)
				edge_list.append(language_edge)
			#now to create the node. this node will sort of represent the trigram in the generative story
			bigram_prob_node = latticeNode.node('dummy',word)
			#now, i have to connect the edges in the edge list to the temp node
			#print 'the size of the edge list is %d'%len(edge_list)
			for language_edge in edge_list :
				language_edge.tail = bigram_prob_node
				bigram_prob_node.backward_pointers.append(language_edge)
				#print 'the number of backward pointers of this node is %d'%len(trigram_prob_node.backward_pointers)
			#raw_input()
			#now, from the node, i have to create the edge that will carry the channel probability. only one edge
			channel_edge_parameter_id = "%s|%s"%(word,tag)
			if channel_model_parameters.has_key(tag) :
				if not channel_model_parameters[tag].has_key(word) :
					channel_model_parameters[tag][word] = float(0)
			else :
				temp_dict = {}
				temp_dict[word] =  float(0)
				channel_model_parameters[tag] = temp_dict

			channel_edge = latticeNode.edge(False,'channel',channel_edge_parameter_id)
			channel_edge.head = bigram_prob_node
			bigram_prob_node.forward_pointers.append(channel_edge)
			#now i have to create the dummy node that will lead to the bigram nodes. this node sort of represents that the word is being spit out in the generative story
			unigram_node = latticeNode.node('language',tag)
			unigram_node.backward_pointers.append(channel_edge)
			channel_edge.tail = unigram_node
			next_iteration_previous_column.append(unigram_node)
		previous_column = list(next_iteration_previous_column)
	#now i have to connect everything to the end node
	for previous_node in previous_column :
		tag = r'&end&'
		bigram_parameter_id = "&end&|%s"%previous_node.parameter_id #creating the n gram probability parameter id
		if language_model_parameters.has_key(previous_node.parameter_id) :
			if not language_model_parameters[previous_node.parameter_id].has_key(tag) :
				language_model_parameters[previous_node.parameter_id][tag] = float(0)
		else :
			temp_dict = {}
			temp_dict[tag] =  float(0)
			language_model_parameters[previous_node.parameter_id] = temp_dict

		bigram_edge = latticeNode.edge(False,'language',bigram_parameter_id)
		bigram_edge.head = previous_node
		bigram_edge.tail = end_node
		previous_node.forward_pointers.append(bigram_edge)
		end_node.backward_pointers.append(bigram_edge)

def makeLatticeTrigram(start_node,end_node,word_list,dictionary,language_model_parameters,channel_model_parameters) :
	previous_column = [] #this will include all the nodes that have to be connected to the 
	previous_column.append(start_node)
	prev_word_tag_list = ['&start&']
	for word in word_list :
		next_iteration_previous_column = []
		for tag in dictionary[word] :
			for prev_tag in prev_word_tag_list :
				bigram = "%s>%s"%(tag,prev_tag)
				#now making the trigram node
				trigram_node = latticeNode.node('language',bigram)
				for prev_node in previous_column :
					trigram = "%s|%s"%(tag,prev_node.parameter_id)
					first_tag,second_tag = prev_node.parameter_id.split(r'>')
					temp_bigram = "%s>%s"%(tag,first_tag)
					if not bigram == temp_bigram : #then there needs to be no edge from teh previous node to this node created
						continue
					#getting the language model parameter		
					if language_model_parameters.has_key(prev_node.parameter_id) :
						if not language_model_parameters[prev_node.parameter_id].has_key(tag) :
							language_model_parameters[prev_node.parameter_id][tag] =float(0)
					else :
						temp_dict = {}
						temp_dict[tag] =  float(0)
						language_model_parameters[prev_node.parameter_id] = temp_dict
				
					trigram_edge = latticeNode.edge(False,'language',trigram)
					prev_node.forward_pointers.append(trigram_edge)
					trigram_edge.head = prev_node
					trigram_edge.tail = trigram_node
					trigram_node.backward_pointers.append(trigram_edge)
				#now do the channel stuff
				#creating the channel model parrameter
				if channel_model_parameters.has_key(tag) :
					if not channel_model_parameters[tag].has_key(word) :
						channel_model_parameters[tag][word] = float(0)
				else :
					temp_dict = {}
					temp_dict[word] =  float(0)
					channel_model_parameters[tag] = temp_dict

				channel_parameter = "%s|%s"%(word,tag)
				channel_edge = latticeNode.edge(False,'channel',channel_parameter)
				channel_edge.head = trigram_node
				trigram_node.forward_pointers.append(channel_edge)
				channel_node = latticeNode.node('channel',bigram)
				channel_edge.tail = channel_node
				channel_node.backward_pointers.append(channel_edge)
				next_iteration_previous_column.append(channel_node) 
		previous_column = list(next_iteration_previous_column)
		prev_word_tag_list = list(dictionary[word])	
	#now, i have to create links to the end node
	fixed_flag = False
	#if (len(previous_column)) == 1 :
	#	fixed_flag = True
	#print 'the size of the previous column before creating the end nodes was %d'%len(previous_column)
	for previous_node in previous_column :
		tag = r'&end&'
		trigram_parameter_id = "&end&|%s"%previous_node.parameter_id #creating the n gram probability parameter id
		if language_model_parameters.has_key(previous_node.parameter_id) :
			if not language_model_parameters[previous_node.parameter_id].has_key(tag) :
				language_model_parameters[previous_node.parameter_id][tag] = float(0)
		else :
			temp_dict = {}
			temp_dict[tag] =  float(0)
			language_model_parameters[previous_node.parameter_id] = temp_dict

		trigram_edge = latticeNode.edge(fixed_flag,'language',trigram_parameter_id)
		trigram_edge.head = previous_node
		trigram_edge.tail = end_node
		previous_node.forward_pointers.append(trigram_edge)
		end_node.backward_pointers.append(trigram_edge)

'''
def makeLatticeTrigram(start_node,end_node,word_list,dictionary,language_model_parameters,channel_model_parameters) :
	previous_column = [] #this will include all the nodes that have to be connected to the 
	previous_column.append(start_node)
	prev_word_tag_list = ['&start&']
	for word in word_list :
	#	print 'the word is %s'%word
	#	print 'the number of tags this word has is %d'%len(dictionary[word])
		next_iteration_previous_column = []
	#	print 'the length of the previous word tag list is %d'%len(prev_word_tag_list)
		for tag in dictionary[word] :
			#i have to create a node for every tag 
			edge_list = [] #this is the edge list that will connect to the node that is created for this tag for this word
			for previous_node in previous_column : #then create the edge from every node in the previous column to this node. each edge corresponds to p(a|b,c)
				fixed = False
			#	if  len(previous_column) == 1 :
			#		fixed = True
				#now create the language parameter that is to be this edge
				parameter_id = "%s|%s"%(tag,previous_node.parameter_id) #creating the n gram probability parameter id
v				if language_model_parameters.has_key(previous_node.parameter_id) :
					if not language_model_parameters[previous_node.parameter_id].has_key(tag) :
						language_model_parameters[previous_node.parameter_id][tag] =float(0)
				else :
					temp_dict = {}
					temp_dict[tag] =  float(0)
					language_model_parameters[previous_node.parameter_id] = temp_dict
				language_edge = latticeNode.edge(fixed,'language',parameter_id)
				language_edge.head = previous_node
				previous_node.forward_pointers.append(language_edge)
				edge_list.append(language_edge)
			#now to create the node. this node will sort of represent the trigram in the generative story
			trigram_prob_node = latticeNode.node('dummy',word)
			#now, i have to connect the edges in the edge list to the temp node
			#print 'the size of the edge list is %d'%len(edge_list)
			for language_edge in edge_list :
				language_edge.tail = trigram_prob_node
				trigram_prob_node.backward_pointers.append(language_edge)
				#print 'the number of backward pointers of this node is %d'%len(trigram_prob_node.backward_pointers)
			#raw_input()
			#now, from the node, i have to create the edge that will carry the channel probability. only one edge
			channel_edge_parameter_id = "%s|%s"%(word,tag)
			if channel_model_parameters.has_key(tag) :
				if not channel_model_parameters[tag].has_key(word) :
					channel_model_parameters[tag][word] = float(0)
			else :
				temp_dict = {}
				temp_dict[word] =  float(0)
				channel_model_parameters[tag] = temp_dict

			channel_edge = latticeNode.edge(False,'channel',channel_edge_parameter_id)
			channel_edge.head = trigram_prob_node
			trigram_prob_node.forward_pointers.append(channel_edge)
			#now i have to create the dummy node that will lead to the bigram nodes. this node sort of represents that the word is being spit out in the generative story
			channel_prob_node = latticeNode.node('dummy',word)
			channel_prob_node.backward_pointers.append(channel_edge)
			channel_edge.tail = channel_prob_node
			#now i have to create the bigram nodes . one bigram node for each bigram. also the edges, which will be fixed
	#		print 'the size of the previous word tag list is %d'%len(prev_word_tag_list)
			for prev_word_tag in prev_word_tag_list :
				bigram_node_parameter_id = "%s>%s"%(tag,prev_word_tag)
	#			print 'we just created the bigram node %s'%bigram_node_parameter_id
				fixed_edge = latticeNode.edge(True,'fixed','')
				fixed_edge.head = channel_prob_node
				channel_prob_node.forward_pointers.append(fixed_edge)
				bigram_node = latticeNode.node('language',bigram_node_parameter_id)
				bigram_node.backward_pointers.append(fixed_edge)
				fixed_edge.tail = bigram_node
				#i also have to add this bigram node to the previous nodes that will be considered for the next word in the word list	
				next_iteration_previous_column.append(bigram_node)
		#previous_column = [] 
		previous_column = list(next_iteration_previous_column)
		#print 'the word was %s'%word
	#	print 'the number of nodes in the previous column now is %d'%len(previous_column)
		#raw_input()
		#prev_word_tag_list = []
		prev_word_tag_list = list(dictionary[word])
	
	#now, i have to create links to the end node
	fixed_flag = False
	#if (len(previous_column)) == 1 :
	#	fixed_flag = True
	#print 'the size of the previous column before creating the end nodes was %d'%len(previous_column)
	for previous_node in previous_column :
		tag = r'&end&'
		trigram_parameter_id = "&end&|%s"%previous_node.parameter_id #creating the n gram probability parameter id
		if language_model_parameters.has_key(previous_node.parameter_id) :
			if not language_model_parameters[previous_node.parameter_id].has_key(tag) :
				language_model_parameters[previous_node.parameter_id][tag] = float(0)
		else :
			temp_dict = {}
			temp_dict[tag] =  float(0)
			language_model_parameters[previous_node.parameter_id] = temp_dict

		trigram_edge = latticeNode.edge(fixed_flag,'language',trigram_parameter_id)
		trigram_edge.head = previous_node
		trigram_edge.tail = end_node
		previous_node.forward_pointers.append(trigram_edge)
		end_node.backward_pointers.append(trigram_edge)
'''

def printLatticeNew (start_node,end_node) :
	reached_end_node_flag = False
	column_node_dictionary = {start_node:1}
	while len(column_node_dictionary.keys()) != 0 :
		print 'the length of the column is %d'%len(column_node_dictionary.keys())
		print 'we are in print lattice new'
		column_nodes = column_node_dictionary.keys()
		column_node_dictionary.clear()
		for node in column_nodes :
			if node == end_node :
				reached_end_node_flag = True
			print 'the node parameter type is %s'%node.parameter_type
			print 'the node id is %s'%node.parameter_id
			print 'the number of edges this node has is %d'%len(node.forward_pointers)
			for edge in node.forward_pointers :
				print 'the edge type is %s'%edge.parameter_type
				print 'the edge id is %s'%edge.parameter_id
				column_node_dictionary[edge.tail] = 1  #addint the next node to the dictionary
			print 'the number of edges connecting to this node are %d'%len(node.backward_pointers)
			for edge in node.backward_pointers:
				print 'the edge type is %s'%edge.parameter_type
				print 'the edge id is %s'%edge.parameter_id


	if reached_end_node_flag == False :
		print 'We didnt reach the end node in the lattice printing !'
		exit(0)


#printing the lattice to see if it has been created properly		
def printLattice(node) :
	print 'we are in print lattice'
	print 'The parameter type is %s'%node.parameter_type
	print 'The parameter id is %s'%node.parameter_id
	print 'the number of edges it has is %s'%(len(node.forward_pointers)) 	
	print 'now printing the edges'
	for edge in node.forward_pointers :
		print 'the edge type is %s'%edge.parameter_type
		print 'the edge id is %s'%edge.parameter_id
	for edge in node.forward_pointers :
		new_node = edge.tail
		print 'we are going to the node with the edge type %s and parameter id %s'%(edge.parameter_type,edge.parameter_id)
		printLattice(new_node)

def initUniformProbs (probabilities_channel,probabilities_language) :
	for tag in probabilities_channel :
		num_keys =  len(probabilities_channel[tag].keys())
		probability = float(1)/float(num_keys)
		for word in probabilities_channel[tag] :
			probabilities_channel[tag][word] = probability

	for tag in probabilities_language :
		num_keys = len(probabilities_language[tag].keys())
#		print num_keys
		probability = float(1)/float(num_keys)
		for next_tag in probabilities_language[tag] :
			probabilities_language[tag][next_tag] = probability

def forwardPass(start_node,end_node,probabilities_channel,probabilities_language) :
	global negative_infinity
	reached_end_node_flag =  False
	#print 'we are in the forward pass'
	start_node.alpha = 0
	column_node_dictionary = {start_node:1}
	while len(column_node_dictionary.keys()) != 0 :
		column_nodes = column_node_dictionary.keys()
		column_node_dictionary.clear()
		for node in column_nodes :
			if node == end_node :
				reached_end_node_flag = True	
			node_alpha_cost = node.alpha
			if node_alpha_cost <= negative_infinity :
				print 'the alpha cost of the node that you were at was negative infinity'
				exit(0)
			#print 'the alpha cost of the node with id %s was %f'%(node.parameter_id,node_alpha_cost)
			for edge in node.forward_pointers :
				#current_path_cost = 0
				if edge.parameter_type == 'language' :
					tag,bigram = edge.parameter_id.split(r'|')
	#				print 'the trigram is %s'%edge.parameter_id
					if not probabilities_language[bigram][tag] == float(0) :
						current_path_cost = node_alpha_cost + math.log10(probabilities_language[bigram][tag])
						next_node = edge.tail
						next_node.alpha = logSum(next_node.alpha,current_path_cost) #you're adding the current path cost to the node
						column_node_dictionary[next_node] = 1
				elif edge.parameter_type == 'channel' :
					word,tag = edge.parameter_id.split(r'|')
					if not probabilities_channel[tag][word] == float(0) :
						current_path_cost = node_alpha_cost + math.log10(probabilities_channel[tag][word])
						next_node = edge.tail
						next_node.alpha = logSum(next_node.alpha,current_path_cost) #you're adding the current path cost to the node
						column_node_dictionary[next_node] = 1
				elif edge.parameter_type == 'fixed' :
						current_path_cost = node_alpha_cost 
						next_node = edge.tail
						next_node.alpha = logSum(next_node.alpha,current_path_cost) #you're adding the current path cost to the node
						column_node_dictionary[next_node] = 1
				else :
					print 'the edge type was not recognized !'
					exit(0)
	#print 'we are done with the forward pass'	
	if reached_end_node_flag == False :
		print 'We didnt reach the end node in the forward pass !'
		exit(0)

'''
def forwardPass(start_node,end_node,probabilities_channel,probabilities_language) :
	global negative_infinity
	reached_end_node_flag =  False
	print 'we are in the forward pass'
	start_node.alpha = float(1)
	column_node_dictionary = {start_node:1}
	while len(column_node_dictionary.keys()) != 0 :
		column_nodes = column_node_dictionary.keys()
		column_node_dictionary.clear()
		for node in column_nodes :
			if node == end_node :
				reached_end_node_flag = True	
			node_alpha_cost = node.alpha
			print 'the node alpha cost was %f'%node_alpha_cost
			if node_alpha_cost <= negative_infinity :
				print 'the alpha cost of the node that you were at was negative infinity'
				exit(0)
	#		print 'the alpha cost of the node was %f'%node_alpha_cost
			for edge in node.forward_pointers :
				#current_path_cost = 0
				if edge.parameter_type == 'language' :
					tag,bigram = edge.parameter_id.split(r'|')
	#				print 'the trigram is %s'%edge.parameter_id
					if not probabilities_language[bigram][tag] == float(0) :
						current_path_cost = node_alpha_cost + probabilities_language[bigram][tag]
						next_node = edge.tail
						next_node.alpha = next_node.alpha + current_path_cost #logSum(next_node.alpha,current_path_cost) #you're adding the current path cost to the node
						column_node_dictionary[next_node] = 1
				elif edge.parameter_type == 'channel' :
					word,tag = edge.parameter_id.split(r'|')
					if not probabilities_channel[tag][word] == float(0) :
						current_path_cost = node_alpha_cost * probabilities_channel[tag][word]
						next_node = edge.tail
						next_node.alpha = next_node.alpha,current_path_cost #you're adding the current path cost to the node
						column_node_dictionary[next_node] = 1
				elif edge.parameter_type == 'fixed' :
						current_path_cost = node_alpha_cost 
						next_node = edge.tail
						next_node.alpha = next_node.alpha,current_path_cost #you're adding the current path cost to the node
						column_node_dictionary[next_node] = 1
				else :
					print 'the edge type was not recognized !'
					exit(0)
	print 'we are done with the forward pass'	
	if reached_end_node_flag == False :
		print 'We didnt reach the end node in the forward pass !'
		exit(0)
'''


def backwardPass(start_node,end_node,probabilities_channel,probabilities_language) :
	global negative_infinity
	reached_start_node_flag = False
	#print 'we are in the backward pass'
	end_node.beta = 0
	column_node_dictionary = {end_node:1}
	while len(column_node_dictionary.keys()) != 0 :
		column_nodes = column_node_dictionary.keys()
		column_node_dictionary.clear()
		for node in column_nodes :
			if node == start_node :
				reached_start_node_flag = True	
			node_beta_cost = node.beta
		#	print 'the node with the id %s had the beta %f'%(node.parameter_id,node_beta_cost)
			if node_beta_cost <= negative_infinity :
				print 'the beta cost of the node that you were at was negative infinity'
				exit(0)
			#print 'the beta cost of the node for the node id %s was %f'%(node.parameter_id,node_beta_cost)
			for edge in node.backward_pointers :
				#current_path_cost = 0
				if edge.parameter_type == 'language' :
					tag,bigram = edge.parameter_id.split(r'|')
					if not probabilities_language[bigram][tag] == float(0) :
						current_path_cost = node_beta_cost + math.log10(probabilities_language[bigram][tag])
						next_node = edge.head
						next_node.beta = logSum(next_node.beta,current_path_cost) #you're adding the current path cost to the node
						column_node_dictionary[next_node] = 1

				elif edge.parameter_type == 'channel' :
					word,tag = edge.parameter_id.split(r'|')
					if not probabilities_channel[tag][word] == float(0) :
						current_path_cost = node_beta_cost + math.log10(probabilities_channel[tag][word])
						next_node = edge.head
						next_node.beta = logSum(next_node.beta,current_path_cost) #you're adding the current path cost to the node
						column_node_dictionary[next_node] = 1
				elif edge.parameter_type == 'fixed' :
						current_path_cost = node_beta_cost 
						next_node = edge.head
						next_node.beta = logSum(next_node.beta,current_path_cost) #you're adding the current path cost to the node
						column_node_dictionary[next_node] = 1
				else :
					print 'the edge type was not recognized !'
					exit(0)
	#print 'we are done with the backward pass'	
	if reached_start_node_flag == False :
		print 'We didnt reach the start node in the backward pass!'
		exit(0)


def getFractionalCounts(start_node,end_node,probabilities_channel,probabilities_language,fractional_counts_channel,fractional_counts_language) :
	global negative_infinity
	reached_end_node_flag = False
	end_node.beta = 0
	column_node_dictionary = {start_node:1}
	while len(column_node_dictionary.keys()) != 0 :
		column_nodes = column_node_dictionary.keys()
		column_node_dictionary.clear()
		for node in column_nodes :
			if node == end_node :
				reached_end_node_flag = True
			for edge in node.forward_pointers :
				if edge.parameter_type == 'language' :
					(tag,bigram)  = edge.parameter_id.split('|')
					if not probabilities_language[bigram][tag] == 0.0 :
						log_fractional_count = 	edge.head.alpha + math.log10(probabilities_language[bigram][tag]) + edge.tail.beta - end_node.alpha
						fractional_count = math.pow(10,log_fractional_count)
						edge.fractional_count = fractional_count
						fractional_counts_language[bigram][tag] += fractional_count
						#print 'the fractional count of the parameter %s is %f'%(edge.parameter_id,edge.fractional_count)
						#print 'the fractional count was '
						column_node_dictionary[edge.tail] = 1  #addint the next node to the dictionary
				elif edge.parameter_type == 'channel' :
					#print 'the channel parameter is %s'%edge.parameter_id
					(word,tag)  = edge.parameter_id.split('|')
					if not probabilities_channel[tag][word] == 0.0 :
						log_fractional_count = 	edge.head.alpha + math.log10(probabilities_channel[tag][word]) + edge.tail.beta - end_node.alpha
						#print 'the log fractional count was %f'%log_fractional_count
						fractional_count = math.pow(10,log_fractional_count)
						edge.fractional_count = fractional_count
						fractional_counts_channel[tag][word] += fractional_count
						#print 'the fractional count of the channel parameter %s was %f'%(edge.parameter_id,edge.fractional_count)
						column_node_dictionary[edge.tail] = 1  #addint the next node to the dictionary
				elif edge.parameter_type == 'fixed' :
						next_node = edge.tail
						column_node_dictionary[next_node] = 1
				else :
					print 'the edge type was not recognized !'
					print 'the edge parameter type was %s'%edge.parameter_type
					exit(0)
	if reached_end_node_flag == False :
		print 'We didnt reach the end node in the get fractional counts method!'
		exit(0)


#given fractional counts, we are getting the probabilities again			
def reEstimateProbabilities(probabilities_channel,probabilities_language,fractional_counts_channel,fractional_counts_language) :
	for tag in probabilities_channel :
		sum = float(0)
		for word in probabilities_channel[tag] :
			sum += fractional_counts_channel[tag][word]
		if sum == float(0) :
			print 'the sum of the fractional counts was zero'
			break
		for word in probabilities_channel[tag] :
			probabilities_channel[tag][word] = fractional_counts_channel[tag][word] / sum
	'''	
	for tag in probabilities_language :
		sum = float (0)
#		print num_keys
		for next_tag in probabilities_language[tag] :
			sum += fractional_counts_language[tag][next_tag]
		if sum == float(0) :
			print 'the sum of the fractional counts was zero'
			break
		for next_tag in probabilities_language[tag] :
			probabilities_language[tag][next_tag] = fractional_counts_language[tag][next_tag] / sum
	'''
'''
def viterbiDecoding(start_node,end_node,probabilities_channel,probabilities_language) :
	print 'we are starting the viterbi decoding'
	global negative_infinity
	reached_end_node_flag = False
	start_node.viterbi_cost = 0
	start_node.viterbi_tag_sequence = []
	column_node_dictionary = {start_node:1}
	while len(column_node_dictionary.keys()) != 0 :
		column_nodes = column_node_dictionary.keys()
		column_node_dictionary.clear()
		for node in column_nodes :
			if node == end_node :
				reached_end_node_flag = True
			#I have to pick the viterbi cost and the tag sequencei			
			for edge in node.forward_pointers :
				if edge.parameter_type == 'language' :
					(tag,bigram)  = edge.parameter_id.split('|')
					if not probabilities_language[bigram][tag] == float(0) :
					#	continue
					#else :
						viterbi_cost = edge.head.viterbi_cost + math.log10(probabilities_language[bigram][tag]) 
						edge.viterbi_tag_sequence = list(edge.head.viterbi_tag_sequence)#edge.viterbi_tag_sequence.extend(edge.head.viterbi_tag_sequence)
						edge.viterbi_cost = viterbi_cost
						#now to update the viterbi tag sequence of the tail node
						if edge.tail.viterbi_cost > edge.viterbi_cost : #then we need to update
							edge.tail.viterbi_cost = edge.viterbi_cost
							#print edge
							edge.tail.viterbi_tag_sequence = list(edge.viterbi_tag_sequence)
						column_node_dictionary[edge.tail] = 1  #addint the next node to the dictionary
				elif edge.parameter_type == 'channel' :
					#print 'the channel parameter is %s'%edge.parameter_id
					(word,tag)  = edge.parameter_id.split('|')
					if not probabilities_channel[tag][word] == 0 :
					#	continue
					#else :
						viterbi_cost = edge.head.viterbi_cost + math.log10(probabilities_channel[tag][word]) 
						edge.viterbi_tag_sequence = list(edge.head.viterbi_tag_sequence)#edge.viterbi_tag_sequence.extend(edge.head.viterbi_tag_sequence)
						edge.viterbi_tag_sequence.append(tag) #= edge.viterbi_tag_sequence.append(tag) #appending the tag
						edge.viterbi_cost = viterbi_cost
						#print edge.tail
						#print edge.viterbi_tag_sequence
						#print edge.head.viterbi_tag_sequence

						#now to update the viterbi tag sequence of the tail node
						if edge.tail.viterbi_cost > edge.viterbi_cost : #then we need to update
							edge.tail.viterbi_cost = edge.viterbi_cost
							#print edge.tail
							#print edge.viterbi_tag_sequence
							#print edge.tail.viterbi_tag_sequence
							edge.tail.viterbi_tag_sequence = list(edge.viterbi_tag_sequence)

						column_node_dictionary[edge.tail] = 1  #addint the next node to the dictionary
				elif edge.parameter_type == 'fixed' :
						edge.viterbi_cost = edge.head.viterbi_cost
						edge.viterbi_tag_sequence = (edge.head.viterbi_tag_sequence)
						if edge.tail.viterbi_cost > edge.viterbi_cost : #then we need to update
							edge.tail.viterbi_cost = edge.viterbi_cost
							edge.tail.viterbi_tag_sequence = list(edge.viterbi_tag_sequence)
						
						#next_node = edge.tail
						column_node_dictionary[edge.tail] = 1
				else :
					print 'the edge type was not recognized !'
					print 'the edge parameter type was %s'%edge.parameter_type
					exit(0)
	#print 'we are done with the viterbi decoding'
	if reached_end_node_flag == False :
		print 'We didnt reach the end node in the viterbi decoding !'
		exit(0)
'''	
def viterbiDecoding(start_node,end_node,probabilities_channel,probabilities_language) :
	#print 'we are starting the viterbi decoding'
	global negative_infinity
	reached_end_node_flag = False
	start_node.viterbi_cost = 0
	start_node.viterbi_tag_sequence = []
	column_node_dictionary = {}
	for edge in start_node.forward_pointers :
		column_node_dictionary[edge.tail] = 1
	while len(column_node_dictionary.keys()) != 0 :
		column_nodes = column_node_dictionary.keys()
		column_node_dictionary.clear()
		for node in column_nodes :
			if node == end_node :
				reached_end_node_flag = True
			#I have to pick the viterbi cost and the tag sequencei			
			for edge in node.backward_pointers :
				if edge.parameter_type == 'language' :
					(tag,bigram)  = edge.parameter_id.split('|')
					if not probabilities_language[bigram][tag] == float(0) :
					#	continue
					#else :
						viterbi_cost = edge.head.viterbi_cost + math.log10(probabilities_language[bigram][tag]) 
						edge.viterbi_cost = viterbi_cost
						#now to update the viterbi tag sequence of the tail node
						if node.viterbi_cost < edge.viterbi_cost : #then we need to update
							node.viterbi_cost = edge.viterbi_cost
							node.viterbi_edge = edge
							#print edge
							#edge.tail.viterbi_tag_sequence = list(edge.viterbi_tag_sequence)
				elif edge.parameter_type == 'channel' :
					#print 'the channel parameter is %s'%edge.parameter_id
					(word,tag)  = edge.parameter_id.split('|')
					if not probabilities_channel[tag][word] == 0 :
					#	continue
					#else :
						viterbi_cost = edge.head.viterbi_cost + math.log10(probabilities_channel[tag][word]) 
						#edge.viterbi_tag_sequence = list(edge.head.viterbi_tag_sequence)#edge.viterbi_tag_sequence.extend(edge.head.viterbi_tag_sequence)
						#edge.viterbi_tag_sequence.append(tag) #= edge.viterbi_tag_sequence.append(tag) #appending the tag
						edge.viterbi_cost = viterbi_cost
						#node.viterbi_cost = edge.viterbi_cost
						#node.viterbi_edge = edge
						#print edge.tail
						#print edge.viterbi_tag_sequence
						#print edge.head.viterbi_tag_sequence
#probabilities_language[tag][next_tag]
						#now to update the viterbi tag sequence of the tail node
						if node.viterbi_cost < edge.viterbi_cost : #then we need to update
							node.viterbi_cost = edge.viterbi_cost
							node.viterbi_edge = edge
							#print edge.tail
							#print edge.viterbi_tag_sequence
							#print edge.tail.viterbi_tag_sequence

				elif edge.parameter_type == 'fixed' :
						edge.viterbi_cost = edge.head.viterbi_cost
						if node.viterbi_cost < edge.viterbi_cost : #then we need to update
							node.viterbi_cost = edge.viterbi_cost
							node.viterbi_edge = edge
						
						#next_node = edge.tail
				else :
					print 'the edge type was not recognized !'
					print 'the edge parameter type was %s'%edge.parameter_type
					exit(0)
			for edge in node.forward_pointers :
				column_node_dictionary[edge.tail] = 1
	
	#print 'we are done with the viterbi decoding'
	if reached_end_node_flag == False :
		print 'We didnt reach the end node in the viterbi decoding !'
		exit(0)

#after doing the viterbi decoding, we have to get the viterbi path
def getViterbiPath(start_node,node,viterbi_path) : 
	current_node = node
	while (1) :
		if current_node == start_node :
			break
		else :
			if current_node.viterbi_edge == '':
				print 'the viterbi edge of the node with id %s was empty!!'%current_node.parameter_id
				exit(0)
			else :
				viterbi_edge = current_node.viterbi_edge
				if viterbi_edge.parameter_type == 'channel' :
					word,tag = viterbi_edge.parameter_id.split('|')
					viterbi_path.append(tag)
				current_node = viterbi_edge.head
				#getViterbiPath(start_node,viterbi_edge.head,viterbi_path)

	

def initializeParameters(start_node,end_node) :
	reached_end_node_flag = False
	column_node_dictionary = {start_node:1}
	while len(column_node_dictionary.keys()) != 0 :
		column_nodes = column_node_dictionary.keys()
		column_node_dictionary.clear()
		for node in column_nodes :
			if node == end_node :
				reached_end_node_flag = True
			node.alpha = negative_infinity
			node.beta = negative_infinity
			node.viterbi_cost=negative_infinity
			node.viterbi_tag_sequence=[]
			node.viterbi_edge = ''
			for edge in node.forward_pointers :
				column_node_dictionary[edge.tail] = 1  #addint the next node to the dictionary
				edge.fractional_count=0.0
				edge.viterbi_cost=negative_infinity
				edge.viterbi_tag_sequence = []
	if reached_end_node_flag == False :
		print 'We didnt reach the end node in the initialization of parameters!'
		exit(0)


def logSum(x,y)	:
#	print 'x is %f' %x
#	print 'y is %f' %y
	min = float(0)
	if x < y :
		min = x
	else :
		min = y

	if y - x > 30 :
#		print' the first if was true'
#		print' we are returning %f' %y
		return(y)
	elif x -y > 30 :
		return(x)
	else :
#		print 'we are in the last else and min is %f' %min
		sum = min + math.log10(math.pow(10,(x-min)) + math.pow(10,(y-min)))
#		print 'sum is %f' %sum
		return sum
#probabilities_language[tag][next_tag]
def calcAccuracy (gold,best_tag) :
	if not len(gold) == len(best_tag) :
		print 'the gold tag sequence length was %d and the viterbi tagged sequence length was %d'%(len(gold),len(best_tag))
#	print len(gold)
#	print len(best_tag)
	match = int(0)	
	for i in range (0,len(best_tag)) :
		if gold[i] == best_tag[i] :
			match += 1
	accuracy = float(match)/float(len(best_tag))
	return(accuracy)
	print 'accuracy = %f' %accuracy

	
def getFreeParametersBigram(word_list,tag_dictionary,free_parameters_language,free_parameters_channel) :
	num_language_parameters = 0
	num_channel_parameters =  0
	prev_word_tag_list = [r'&start&']
	for word in word_list :
		#print 'the word is %s'%word
		for tag in prev_word_tag_list :
			for next_tag in tag_dictionary[word] :
				if not free_parameters_language.has_key(tag) :
					temp_dict = {next_tag:0.0}
					free_parameters_language[tag] = temp_dict
					num_language_parameters += 1
					#print 'bigram created %s %s'%(tag,next_tag)
				elif not free_parameters_language[tag].has_key(next_tag) :
					free_parameters_language[tag][next_tag] = 0.0
					num_language_parameters += 1
					#print 'bigram created %s %s'%(tag,next_tag)
		#now making the end tags
		#now creating the channel parameters
		#prev_word_tag_list = []
		prev_word_tag_list = list(tag_dictionary[word])
		for tag in tag_dictionary[word] :
			if not free_parameters_channel.has_key(tag) :
				temp_dict = {word:0.0}
				free_parameters_channel[tag] = temp_dict
				num_channel_parameters += 1
			elif not free_parameters_channel[tag].has_key(word) :
				free_parameters_channel[tag][word] = 0.0
				num_channel_parameters += 1
	last_word = word_list[len(word_list) -1]
	#print 'the last word is %s'%last_word
	for tag in tag_dictionary[last_word] :
		#print 'the tag is %s' %tag
		if not free_parameters_language.has_key(tag) :
			temp_dict = {r'&end&':0.0}
			free_parameters_language[tag] = temp_dict
			num_language_parameters += 1
#			print 'bigram created &end& %s'%(tag)
		elif not free_parameters_language[tag].has_key(r'&end&') :
			free_parameters_language[tag][r'&end&'] = 0.0
			num_language_parameters += 1
#			print 'bigram created &end& %s'%(tag)

	#print num_language_parameters
	#print num_channel_parameters
	return ((num_language_parameters,num_channel_parameters))
	
def readTaggingOutput(tagged_file) :
	file = open(tagged_file,'r')
	tag_seq = []
	for line in file :
		line = line.strip()
		temp = line.split(' ')
		tag_seq.extend(temp)
	return(tag_seq)

