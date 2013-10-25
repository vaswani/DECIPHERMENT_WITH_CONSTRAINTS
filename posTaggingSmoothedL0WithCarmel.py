#!/usr/bin/env python

#the change in this one is going to be that we are going to include both tag bigrams and tags
#I have to change the optimization here. I will only optimze every set of conditional probabilities on its own. So, there will be 44 optimizations in each iteration
#i dump theta
#3000 pgd iterations gives 77.7 % acc
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
from emUtil import *
from pgdMethods import *


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



#initializing the parameter vectors. Here, each parameters is given a number which is the index in the parameter vector that algencan optimizes. The constraints are also given a number . We need these constraint numbers for algencan
#the new technique here does the optimization for each constraint separately. Therefore,it makes sense to prepare these in advance and then send them off to algencan when they are needed
#here, we will create the indexing of the parameters needed by algencan so that we can scale it. Instead of throwing the entire problem to algencan, we will #also, I will call this just once because the only code that is affected by being called in different places is the one that checks if the probability is zero. Since we will not make any probabilities zero in algencan, we just need to call this once
def createParametersForScaling(probabilities,parameter_to_index) :
  #first make a list of all the parametersa
  #constraint_parameters.clear()
  #so, I have to take the current probabilities and then create the constraint parameters and the initial parameters
  #you should also calculate the current number of zero parameters and pass it to the argmax function
  for tag in probabilities :
    rowwise_parameter_to_index = {}
    paramter_counter = 0
    for next_tag in probabilities[tag]:
      rowwise_parameter_to_index[next_tag] = paramter_counter
      paramter_counter += 1
    parameter_to_index[tag] = rowwise_parameter_to_index
  #print 'after indexing, the paramter to index is '
  #print parameter_to_index

def main() :
  # setting up options
  parser = OptionParser()
  parser.add_option("--num_iter", action="store", type="int", dest="num_iterations",default=100,help="number of iterations you would like to run em+smoothed l0. Default is 50")
  parser.add_option("--alpha", action="store", type="float", dest="alpha",default = 0.0,help="alpha,the weight of the smoothed l0 penalty. Defaut is 0 ")
  parser.add_option("--beta", action="store", type="float", dest="beta",default = 0.5,help="beta, the smoothness of the l0 prior, smaller the sigma, closer the approximation to true L0. Default is 0.5")
  parser.add_option("--slack", action="store_true", dest="slack_option",default = False,help="if you want to project on the simplex with slack.")
  parser.add_option("--noslack", action="store_false", dest="slack_option",default = False, help="if you want to project on the simplex with no slack. This is the regular projection approach")
  parser.add_option("--num_pgd_iterations", action="store", type="int", dest="num_pgd_iterations",default = 10,help="Number of Projected Gradient Descent to run. Default is 100")
  parser.add_option("--eta", action="store", type="float", dest="eta",default = 0.1,help="Eta, the constant step size for PGD using armijo line search. Default is 0.1")
  parser.add_option("--armijo_beta", action="store", type="float", dest="armijo_beta",default = 0.5,help="Set value for Armijo beta, the beta used in in armijo line search. Default value is 0.2")
  parser.add_option("--armijo_sigma", action="store", type="float", dest="armijo_sigma",default = 0.5,help="Set value for Armijo sigma, the sigma used in in armijo line search. Lower bound is 0.0001")
  parser.add_option("--lower_bound", action="store", type="float", dest="lower_bound",default = 0.000001,help="Set value for the lower bound on the probability.Default is 10E-6")
  (options, args) = parser.parse_args()

  print options
  print args
  
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
  print 'alpha is ',options.alpha
  print 'beta is ',options.beta
  print 'eta is ',options.eta

  #setting up paramters
  fractional_counts_language = {}
  fractional_counts_channel = {}
  probabilities_channel = {}
  probabilities_language = {}
  current_probabilities = {}
  current_fractional_counts = {}
  #constraint_parameters = []
  initial_parameter_args = {}
  current_initial_parameter_args = {}
  parameter_to_index = {}
  parameter_counter = 0
  num_constraints = 1
  constraint_tags_dict = {} #this will hold the tag for the constraint. the key is the constraint id and the value is the tag corresponding to this constraint
  #beta = float(0)
  #alpha = float(0)
  global_num_parameters = 0
  init_option = ''
  current_optimization_tag = '' #this will hold the of the constraint for which we are doing the optimization

  #adding parser options
  '''
  print 'beta is ',beta
  print 'alpha is ',alpha
  print 'eta is ',eta
  raw_input()
  '''
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
  #now, we will build all the lattices, and create a special start node and end node for every sentence
  start_node_end_node_list = []

  fractional_counts_channel = copy.deepcopy(free_parameters_channel)
  probabilities_channel = copy.deepcopy(free_parameters_channel)

  #createRandomPoint(probabilities_channel)
  emMethods.initUniformProbs(probabilities_channel)
#  emMethods.initUniformProbs(probabilities_language,probabilities_channel)
  emMethods.writeFst('cipher.fst',probabilities_channel)
  #sys.exit()
  final_probabilities_channel = copy.deepcopy(free_parameters_channel)

  run_training = r'./carmel  --train-cascade -u -M 0 -m -HJ cipher.data cipher.wfsa cipher.fst'

  #running the EM iterations
  #we are creating the indexes for algencan . Notice that here, the probabilities language is already uniform and therefore none of them will be zero
  createParametersForScaling(probabilities_channel,parameter_to_index)
  start_time = time.clock()
  print 'start time was ',start_time
  #fractional_counts_dump_file = open('fractional.params','w')
  #probabilities_dump_file = open('probs.params','w')
  #optimized_probabilities_dump_file = open('probs.optimized.params','w')
  for i in range (0,num_iterations) :
    print 'the iteration number was ',i
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
    viterbi_command = r'cat cipher.data.quotes | ./carmel  -u  -srbk 1 -QEWI cipher.wfsa.noe cipher.fst > decipherment_output'
    (status,output) = commands.getstatusoutput(viterbi_command)
    #tagged_sequence = emMethods.readTaggingOutput('tagging_output')  
    deciphered_sequence = emMethods.readCipherFile('decipherment_output')
    accuracy = emMethods.calcAccuracy(gold_cipher,deciphered_sequence)

    print 'The accuracy was %s and the objective function value was %s'%(str(accuracy),str(evaluateObjectiveFuncValue(total_corpus_probability,probabilities_language,probabilities_channel,alpha,beta)))
      

    #first optimizing the tag bigrams and doing it per parameter
    for tag in probabilities_channel :
      print 'we are starting a new tag'
      if len(probabilities_channel[tag]) == 1 :
        continue
      print 'we are currently optimizing for tag', tag
      #current_initial_parameter_args = initial_parameter_args[tag]
      current_optimization_tag = tag
      parameter_counter = len(probabilities_channel[tag].keys())      
      current_fractional_counts = fractional_counts_channel
      #optimizing per constraint

      x=zeros(parameter_counter)
      expected_counts=zeros(parameter_counter)
      expected_counts_sum = 0.
      for next_tag in probabilities_channel[current_optimization_tag].keys() :
        parameter_number = parameter_to_index[current_optimization_tag][next_tag]
        x[parameter_number] = probabilities_channel[current_optimization_tag][next_tag]  
        expected_counts[parameter_number] = current_fractional_counts[current_optimization_tag][next_tag]
        expected_counts_sum +=  current_fractional_counts[current_optimization_tag][next_tag]
        
      print 'expected counts sum was ',expected_counts_sum
      #print parameter_to_index

      #current_eta = eta_0 #/sqrt(num_iterations+1)
      print 'Doing projected gradient descent'
      new_probabilities = projectedGradientDescentWithArmijoRule(x = x,expected_counts = expected_counts,num_pgd_iterations = num_pgd_iterations,eta = eta,lower_bound = lower_bound,armijo_beta = armijo_beta,armijo_sigma = armijo_sigma,alpha = alpha,beta = beta,slack_option = slack_option, parameter_counter = parameter_counter)
      print 'finished projected gradient descent'
      #new_probabilities = newtonProjectedGradientDescentWithArmijoRule(x,num_pgd_iterations,eta,lower_bound,armijo_beta,armijo_sigma)
      #print new_probabilities
      #raw_input()
      print 'we are replacing channel probs'
      assignProbs(new_probabilities,probabilities_channel[current_optimization_tag],parameter_to_index[current_optimization_tag])
  
      
    #now writing the fsa back again
    emMethods.writeFst('cipher.fst',probabilities_channel)


    print 'checking the initial zeros in channel model'
    checkZeros(probabilities_channel)
  
    fractional_counts_channel = copy.deepcopy(free_parameters_channel)
    final_probabilities_channel = copy.deepcopy(probabilities_channel)
    print 'at the end of the iteration'

  elapsed_time = time.clock() - start_time
  print 'the elapsed time was ',elapsed_time




def assignProbs(x,probabilities,parameter_to_index) :
  for next_tag in parameter_to_index :
      parameter_number =   parameter_to_index[next_tag]
      probability = x[parameter_number]
      probabilities[next_tag] = probability

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
