#!/usr/bin/python

negative_infinity = float(-999999)

'''
class node :
#a node has a label which is is the id of the node 
        def __init__ (self,word,tag) : #,backward_link,forward_link,fixed_flag) :
		self.word = word
		self.tag = tag
		self.p_word_given_tag = float(0)
		self.alpha = float(0)
		self.beta = float(0)
		self.alpha_language = float(0) #this basically will hold the alphas and the betas for the free language parameter models. we can obtain alpha langugage by muliptlying the alpha with the channel probability
		self.beta_language = float(0)  
		self.fixed_flag = 0 #this will indicate if the parameter is fixed or not
		self.count_tag_word = float(0)
		self.forward_pointers = [] #list of arcs from node in n col to n+1 col  nodes lablelled wiht p_tag_given_tag
		self.backward_pointers = []  #arcs from node in n col to n-1 col
		self.best_tag_sequence = [] #this stores the best tagging up to now
		self.best_tag_score = float(-99999999) #this stores the score of the best tagging sequence
		self.in_greedy_path = False
		self.to_delete = False
'''

class edge :
	global negative_infinity
	def __init__(self,fixed_flag,parameter_type,parameter_id) :
		self.fixed_flag = False
		self.parameter_type = parameter_type # 'language' or 'channel' parameter or dummy
		self.fractional_count=0.0
		self.parameter_id = parameter_id #what bigram/trigram/channel etc. parameter it is
		self.head = ''
		self.tail = ''
		self.viterbi_cost=negative_infinity
		self.viterbi_tag_sequence = []

class node :
	global negative_infinity
	def __init__(self,parameter_type,parameter_id) :
		self.alpha = negative_infinity
		self.parameter_type = parameter_type #language or channel or dummy
		self.parameter_id = parameter_id #now, the node can hold the language model parameter. Nodes will hold bigrams or unigrams, which will then later help in creating the edges that will carry the bigram or trigram probabilities
		self.beta = negative_infinity
		self.viterbi_cost=negative_infinity
		self.viterbi_tag_sequence=[]
		self.viterbi_edge = '' #this will hold the id of the best viterbi edge
		self.backward_pointers = []
		self.forward_pointers = []
