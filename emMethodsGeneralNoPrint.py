#these em methods will apply for a general lattice and do forward backward on a general lattice
# the bigram separation charachter is >
import latticeNode
import random
import math
from math import e
import copy
import re
import sys
from collections import defaultdict

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


def getFreeCipherParametersChannel(plaintext,free_parameters_channel,cipher_letter_dict=None) :
  if cipher_letter_dict != None:
    for plaintext_letter in plaintext :
      for cipher_letter in cipher_letter_dict.keys() :
        if cipher_letter == "_" :
          continue 
        if not free_parameters_channel.has_key(plaintext_letter) :
          temp = {}
          temp[cipher_letter] = 0.0
          free_parameters_channel[plaintext_letter] =  temp
        else :
          free_parameters_channel[plaintext_letter][cipher_letter] = 0.0
  else :
    for k in range(65, 91):
      plain_letter = chr(k)
      for j in range(65, 91):
        cipher_letter = chr(j).lower()
        free_parameters_channel[plain_letter][cipher_letter] = 0.0

def writePosteriorDecodingFST(ciphertext,probabilities_channel,posterior_decoding_file_name):
  posterior_decoding_file = open(posterior_decoding_file_name,'w')
  posterior_decoding_file.write("F\n")
  current_state = 0
  next_state = 1
  for i,cipher_letter in enumerate(ciphertext):
    #print 'cipher letter is',cipher_letter
    if i == len(ciphertext)-1:
      string_next_state = 'F'
    else :
      string_next_state = str(next_state)
    if cipher_letter != "_":
      for j,k in enumerate(range(65, 91)):
        plain_letter = chr(k)
        posterior_decoding_file.write("(%d (%s \"%s\" \"%s\" %s))\n"%(current_state,string_next_state,plain_letter,cipher_letter,repr(probabilities_channel[plain_letter][cipher_letter])))
    else:
      posterior_decoding_file.write("(%d (%s \"_\" \"_\" 1))\n"%(current_state,string_next_state))
    current_state += 1
    next_state += 1
  posterior_decoding_file.close()


#this reads the carmel file that carries the fractional counts and then stores them in the dictionary. right now this is only for bigrams
def getPosteriorDecode(posterior_decode_file):
  file = open(posterior_decode_file)
  dump_line = file.readline()
  current_start_state = "0"
  current_end_state = "1"
  argmax = "A"
  decode = []
  max_marginal = -99999
  for line in file :
    #print line
    match = ''
    match = carmel_wfst_re.match(line)
    #print "the match is %s"%match
    if match == None :
      print 'the line did not match while reading the fractional counts'    
      print line
      exit(0)
    else :
      start_state = match.group(1)
      end_state = match.group(2)
      plain_letter = match.group(3)
      cipher_letter = match.group(4)
#        fractional_count = float(match.group(5))
      fractional_count = float(0)
      #print "\"%s\""%match.group(5)
      if match.group(5)[0] == 'e' :
        fractional_count = eval(match.group(5).replace('^','**'))
      else :
        fractional_count = float(match.group(5))
      '''
      if not dictionary.has_key(tag) :
        print 'the channel dictionary did not have the tag %s'%tag
      elif not dictionary[tag].has_key(word) :
        print 'the channel dictionary did not have the word %s'%word
      '''
      if fractional_count <= 10E-10:
        fractional_count = 10E-10
      #print 'plain letter is ',plain_letter
      if start_state == current_start_state and end_state == current_end_state:
        if max_marginal < fractional_count:
          max_marginal = fractional_count
          argmax = plain_letter
      else:
        current_start_state = start_state
        current_end_state = end_state
        decode.append(argmax)
        max_marginal = -99999
        argmax = plain_letter
    #print 'the line matched'
    #print match.group(1)
    #print match.group(2)
    #print match.group(3)
    #print match.group(4)
    #if parameters == 'channel' :
    #  print match.group(5)
  #raw_input()
  return(decode)

def readCipherFile(file_name) :
  file = open(file_name)

  data = []

  for line in file :
    line = line.strip()
    items = line.split(' ')
    for item in items :
      if item != "_" :
        data.append(item)

  file.close()
  return data

#this reads the carmel file that carries the fractional counts and then stores them in the dictionary. right now this is only for bigrams
def initFromWfst(carmel_file,dictionary,parameters) :
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
        #  next_tag = r'&end&'
        #else :
        #  next_tag = next_tag[1:len(next_tag)-1]
        #print 'temp tag is %s'%temp_tag
        #raw_input()
        if next_tag == 'F' :
          next_tag = '&end&'
        fractional_count = float(0)
        if match.group(4)[0] == 'e' :
          fractional_count = eval(match.group(4).replace('^','**'))
        else :
          fractional_count = float(match.group(4))
        if not dictionary.has_key(tag) :
          print 'the language dictionary did not have the tag %s'%tag
          exit(0)
        elif not dictionary[tag].has_key(next_tag) :
          print 'the language dictionary did not have the next tag %s'%next_tag
          exit(0)

        if fractional_count <= 10E-10:
          fractional_count = 10E-10
        dictionary[tag][next_tag] = fractional_count
      else :
        tag = match.group(3)
        word = match.group(4)
#        fractional_count = float(match.group(5))
        fractional_count = float(0)
        #print "\"%s\""%match.group(5)
        if match.group(5)[0] == 'e' :
          fractional_count = eval(match.group(5).replace('^','**'))
        else :
          fractional_count = float(match.group(5))

        if not dictionary.has_key(tag) :
          print 'the channel dictionary did not have the tag %s'%tag
        elif not dictionary[tag].has_key(word) :
          print 'the channel dictionary did not have the word %s'%word
        if fractional_count <= 10E-10:
          fractional_count = 10E-10
        dictionary[tag][word] = fractional_count
      #print 'the line matched'
      #print match.group(1)
      #print match.group(2)
      #print match.group(3)
      #print match.group(4)
      #if parameters == 'channel' :
      #  print match.group(5)
    #raw_input()
  return(counter)

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
        #  next_tag = r'&end&'
        #else :
        #  next_tag = next_tag[1:len(next_tag)-1]
        #print 'temp tag is %s'%temp_tag
        #raw_input()
        if next_tag == 'F' :
          next_tag = '&end&'
        fractional_count = float(0)
        if match.group(4)[0] == 'e' :
          fractional_count = eval(match.group(4).replace('^','**'))
        else :
          fractional_count = float(match.group(4))
        if not dictionary.has_key(tag) :
          print 'the language dictionary did not have the tag %s'%tag
          exit(0)
        elif not dictionary[tag].has_key(next_tag) :
          print 'the language dictionary did not have the next tag %s'%next_tag
          exit(0)

        if fractional_count <= 10E-10:
          fractional_count = 10E-10
        dictionary[tag][next_tag] = fractional_count
      else :
        tag = match.group(3)
        word = match.group(4)
#        fractional_count = float(match.group(5))
        fractional_count = float(0)
        #print "\"%s\""%match.group(5)
        if match.group(5)[0] == 'e' :
          fractional_count = eval(match.group(5).replace('^','**'))
        else :
          fractional_count = float(match.group(5))

        if not dictionary.has_key(tag) :
          print 'the channel dictionary did not have the tag %s'%tag
        elif not dictionary[tag].has_key(word) :
          print 'the channel dictionary did not have the word %s'%word
        if fractional_count <= 10E-10:
          fractional_count = 10E-10
        dictionary[tag][word] = fractional_count
      #print 'the line matched'
      #print match.group(1)
      #print match.group(2)
      #print match.group(3)
      #print match.group(4)
      #if parameters == 'channel' :
      #  print match.group(5)
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
        #  next_tag = r'&end&'
        #else :
        #  next_tag = next_tag[1:len(next_tag)-1]
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
      #  print match.group(5)
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

def readCipherLines (cipher_file) :
  cipher_letter_dict = {}
  cipher_list = []
  cipher_file_handle = open(cipher_file,'r')
  for line in cipher_file_handle :
    line = line.strip()
    if line == '' :
      continue
    #words = line.split('\t')
    cipher_letters = line.split(' ')
    for letter in cipher_letters :
      cipher_letter_dict[letter] = 1
    cipher_list.append(line)
  
  return (cipher_list,cipher_letter_dict)


def getUniqCipherLetters(cipher_file) :
  cipher_letter_dict = {}
  cipher_file_handle = open(cipher_file,'r')
  for line in cipher_file_handle :
    line = line.strip()
    if line == '' :
      continue
    #words = line.split('\t')
    cipher_letters = line.split(' ')
    for letter in cipher_letters :
      cipher_letter_dict[letter] = 1
  
  return (cipher_letter_dict)


def getCipherLetterProbs(probs,cipher_file):
  for line in open(cipher_file):
    line = line.strip()
    if line == '':
      continue
    sum = 0.
    for letter in line.split():
      if letter == '_': 
        continue
      probs[letter] += 1
      sum += 1
    for letter in probs:
      probs[letter] /= sum

def initUniformProbs (probabilities_channel) :
  for tag in probabilities_channel :
    num_keys =  len(probabilities_channel[tag].keys())
    probability = float(1)/float(num_keys)
    for word in probabilities_channel[tag] :
      probabilities_channel[tag][word] = probability

def readGaussians(gaussians,gaussian_params_file):
  for line in open(gaussian_params_file) :
    line = line.strip()
    if line == '':
      continue
    items = line.split()
    gaussians[items[0]].append(float(items[1]))
    gaussians[items[0]].append(float(items[2]))


def initFromGaussians(probabilities_channel,gaussians,prob_ciphertext,std_mult) :
  W = defaultdict(lambda:defaultdict(float))
  sqrt_two_pi = 1/math.sqrt(2*math.pi)
  for plain_letter in probabilities_channel:
    if len(probabilities_channel[plain_letter].keys()) > 1:
      #constant = (1./gaussians[plain_letter][1])*sqrt_two_pi
      sum = 0.
      adjusted_std = gaussians[plain_letter][1]*std_mult
      for cipher_letter in probabilities_channel[plain_letter] :
        number = math.exp(-(gaussians[plain_letter][0]-prob_ciphertext[cipher_letter])**2/(2*adjusted_std**2))
        if number <= 10E-10:
          number = 10E-10
        W[plain_letter][cipher_letter] = number
        #print 'plain_prob for letter',plain_letter,':',prob_plaintext[plain_letter],' cipher_prob for letter',cipher_letter,':',prob_ciphertext[cipher_letter]
        #print 'plain:',plain_letter,' cipher:',cipher_letter,' pdf:',W[plain_letter][cipher_letter] 
        sum += W[plain_letter][cipher_letter]
      for cipher_letter in probabilities_channel[plain_letter]:
        probabilities_channel[plain_letter][cipher_letter] =  W[plain_letter][cipher_letter]/sum
    else :
      probabilities_channel[plain_letter][probabilities_channel[plain_letter].keys()[0]] = 1.

def dictionaryToArray(array_to_populate,dictionary_to_populate_from,parameter_to_index):
  for i,k in enumerate(range(65, 91)):
    plain_letter = chr(k)
    #print plain_letter
    #print 'number of cipher letters is ',len(probabilities_channel[plain_letter].keys())
    for cipher_letter in dictionary_to_populate_from[plain_letter] :
      parameter_number = parameter_to_index[plain_letter][cipher_letter]
      #print 'parameter number is ',parameter_number
      #print 'i is ',i
      array_to_populate[i][parameter_number] = dictionary_to_populate_from[plain_letter][cipher_letter]  

def initFromGaussiansSingleStd(probabilities_channel,prob_plaintext,prob_ciphertext,std_mult) :
  W = defaultdict(lambda:defaultdict(float))
  sqrt_two_pi = 1/math.sqrt(2*math.pi)
  for plain_letter in probabilities_channel:
    if len(probabilities_channel[plain_letter].keys()) > 1:
      #constant = (1./gaussians[plain_letter][1])*sqrt_two_pi
      sum = 0.
      adjusted_std = 0.01*std_mult
      for cipher_letter in probabilities_channel[plain_letter] :
        number = math.exp(-(prob_plaintext[plain_letter]-prob_ciphertext[cipher_letter])**2/(2*adjusted_std**2))
        if number <= 10E-10:
          number = 10E-10
        W[plain_letter][cipher_letter] = number
        #print 'plain_prob for letter',plain_letter,':',prob_plaintext[plain_letter],' cipher_prob for letter',cipher_letter,':',prob_ciphertext[cipher_letter]
        #print 'plain:',plain_letter,' cipher:',cipher_letter,' pdf:',W[plain_letter][cipher_letter] 
        sum += W[plain_letter][cipher_letter]
      for cipher_letter in probabilities_channel[plain_letter]:
        probabilities_channel[plain_letter][cipher_letter] =  W[plain_letter][cipher_letter]/sum
    else :
      probabilities_channel[plain_letter][probabilities_channel[plain_letter].keys()[0]] = 1.

def initIdentity(probabilities_channel) :
  for plain_letter in probabilities_channel:
    if len(probabilities_channel[plain_letter].keys()) > 1:
      sum = 0.
      for cipher_letter in probabilities_channel[plain_letter] :
        if cipher_letter.upper() == plain_letter :
          probabilities_channel[plain_letter][cipher_letter] = 1.
        else:
          probabilities_channel[plain_letter][cipher_letter] = 10E-5
        sum += probabilities_channel[plain_letter][cipher_letter]
      for cipher_letter in probabilities_channel[plain_letter]:
        probabilities_channel[plain_letter][cipher_letter] /= sum
    else :
      probabilities_channel[plain_letter][probabilities_channel[plain_letter].keys()[0]] = 1.

def generateNoiseCipher(ciphertext):
  #first get positions of spaces
  space_dict = {}
  for i,item in enumerate(ciphertext):
    if item == '_':
      space_dict[i] = 1
  #randomly permute four pairs of cipher letters
  current_positions = {}
  for i in range(10) :
    position1 = random.randint(0,len(ciphertext)-1)
    while position1 in space_dict or position1 in current_positions:
      position1 = random.randint(0,len(ciphertext)-1)
    current_positions[position1] = 1
    position2 = random.randint(0,len(ciphertext)-1)
    while position2 in current_positions or position2 in space_dict:
      position2 = random.randint(0,len(ciphertext)-1)
    current_positions[position2] = 1.
    ciphertext[position1],ciphertext[position2] = ciphertext[position2],ciphertext[position1]

  print 'after swapping ciphertext is ',' '.join(ciphertext)
  '''  
  position3 = random.randint(0,len(ciphertext)-1)
  while position3 in current_positions or position3 in space_dict:
    position3 = random.randint(0,len(ciphertext)-1)
  current_positions[position3] = 1
  position4 = random.randint(0,len(ciphertext)-1)
  while position4 in current_positions or position4 in space_dict:
    position4 = random.randint(0,len(ciphertext)-1)

  ciphertext[position3],ciphertext[position4] = ciphertext[position4],ciphertext[position3]
  '''



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
#    print num_keys
    for next_tag in probabilities_language[tag] :
      sum += fractional_counts_language[tag][next_tag]
    if sum == float(0) :
      print 'the sum of the fractional counts was zero'
      break
    for next_tag in probabilities_language[tag] :
      probabilities_language[tag][next_tag] = fractional_counts_language[tag][next_tag] / sum
  '''

def normalizeTable(probabilities) :
  for given in probabilities:
    sum = 0.
    for event in probabilities[given] :
      sum += probabilities[given][event]
    if sum == 1. :
      continue
    else :
      print 'the sum of probabilities was not 1 and was ',sum
    for event in probabilities[given] :
      probabilities[given][event] = probabilities[given][event]/sum 
  

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


def logSum(x,y)  :
#  print 'x is %f' %x
#  print 'y is %f' %y
  min = float(0)
  if x < y :
    min = x
  else :
    min = y

  if y - x > 30 :
#    print' the first if was true'
#    print' we are returning %f' %y
    return(y)
  elif x -y > 30 :
    return(x)
  else :
#    print 'we are in the last else and min is %f' %min
    sum = min + math.log10(math.pow(10,(x-min)) + math.pow(10,(y-min)))
#    print 'sum is %f' %sum
    return sum
'''
def calcAccuracy (gold,best_tag) :
  if not len(gold) == len(best_tag) :
    print 'the gold tag sequence length was %d and the viterbi tagged sequence length was %d'%(len(gold),len(best_tag))
#  print len(gold)
#  print len(best_tag)
  match = int(0)  
  for i in range (0,len(best_tag)) :
    if gold[i] == best_tag[i] :
      match += 1
  accuracy = float(match)/float(len(best_tag))
  return(accuracy)
  print 'accuracy = %f' %accuracy
'''
def calcAccuracy(gold_data,output_data) :
  space_cleaned_gold_data = [item for item in gold_data if item != "_"]
  space_cleaned_output_data = [item for item in output_data if item != "_" ]
  if not len(space_cleaned_output_data) == len(space_cleaned_gold_data) :
    print 'The gold file is of length %d and the output file is of length %d'%(len(space_cleaned_gold_data),len(space_cleaned_output_data))
    sys.exit()

  total_items = 0
  matches   = 0

  for i in range(0,len(space_cleaned_gold_data)) :
    total_items += len(space_cleaned_gold_data[i])
    for j in range(0, len(space_cleaned_gold_data[i])) :
      if space_cleaned_gold_data[i][j] == space_cleaned_output_data[i][j]:
        matches += 1
  
  accuracy = float(matches)/total_items
  #print 'the acccuracy was %f'%accuracy
  return(accuracy)

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
#      print 'bigram created &end& %s'%(tag)
    elif not free_parameters_language[tag].has_key(r'&end&') :
      free_parameters_language[tag][r'&end&'] = 0.0
      num_language_parameters += 1
#      print 'bigram created &end& %s'%(tag)

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

