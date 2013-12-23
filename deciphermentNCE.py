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
from collections import defaultdict


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
  #for tag in probabilities :
  for k in range(65, 91):
    plain_letter = chr(k)
    rowwise_parameter_to_index = {}
    paramteter_counter = 0
    for cipher_letter in probabilities[plain_letter]:
      rowwise_parameter_to_index[cipher_letter] = paramteter_counter
      paramteter_counter += 1
    parameter_to_index[plain_letter] = rowwise_parameter_to_index
  #print 'after indexing, the paramter to index is '
  #print parameter_to_index

def getProb(output):
  prob_match = probability_re.search(output)
  if prob_match == None :
    print'we should have found a probability'
  else :
    print 'the probability is %s'%prob_match.group(1)
  temp_corpus_probability = float(prob_match.group(1)[2:len(prob_match.group(1))])
  total_corpus_probability = 0.693147181 * temp_corpus_probability
  return(total_corpus_probability)

def computePositiveWeight(log_true_prob,log_noise_prob,num_noise_samples):
    print 'true prob',log_true_prob
    print 'noise prob',log_noise_prob
    log_noise_prob += math.log(num_noise_samples) 
    log_positive_weight  =  log_noise_prob - logaddexp(log_noise_prob,log_true_prob)
    return(math.exp(log_positive_weight))

def computeNegativeWeight(log_true_prob,log_noise_prob,num_noise_samples):
    print 'true prob',log_true_prob
    print 'noise prob',log_noise_prob
    log_noise_prob += math.log(num_noise_samples) 
    log_negative_weight =  log_true_prob - logaddexp(log_noise_prob,log_true_prob)
    return(math.exp(log_negative_weight))

def runViterbiAndGetAccuracy(gold_cipher,viterbi_command):
  print 'we are now checking the accuracies'
  (status,output) = commands.getstatusoutput(viterbi_command)
  print 'status',status
  print 'output',output
  #tagged_sequence = emMethods.readTaggingOutput('tagging_output')  
  deciphered_sequence = emMethods.readCipherFile('decipherment_output')
  print 'length of deciphered sequence was ',len(deciphered_sequence)
  accuracy = emMethods.calcAccuracy(gold_cipher,deciphered_sequence)
  print 'The accuracy was %s'%str(accuracy)
  #print 'The accuracy was %s and the objective function value was %s'%(str(accuracy),str(evaluateObjectiveFuncValue(total_corpus_probability,probabilities_language,probabilities_channel,current_alpha,beta)))

def runCarmel(run_training):
  print 'the run training command is', run_training
  #this will create the parameters
  #getting the probability of the true data under the true model
  total_true_corpus_probability = 0.0
  (status,output) = commands.getstatusoutput(run_training)  
  print 'output:',output
  print 'status:',status
  corpus_probability = getProb(output)
  print' the probability of the corpus was ', repr(corpus_probability)

  return(corpus_probability)

def runCarmelAndGetFracCounts(run_true_training,run_noise_training,fractional_counts_channel,trained_cipher_file):
  total_true_corpus_probability = runCarmel(run_true_training)
  total_noise_corpus_probability = runCarmel(run_noise_training)

  print 'reading language fractional counts from',trained_cipher_file
  emMethods.readCarmelFractionalCounts(trained_cipher_file,fractional_counts_channel,'channel')
  print 'read the fst'

  return(total_true_corpus_probability,total_noise_corpus_probability)
  
def main() :
  # setting up options
  parser = OptionParser()
  parser.add_option("--num_iter", action="store", type="int", dest="num_iterations",default=100,help="number of iterations you would like to run em+smoothed l0. Default is 50")
  parser.add_option("--num_noise_samples", action="store", type="int", dest="num_noise_samples",default=10,help="number of noise samples. Default is 10")
  parser.add_option("--initial_alpha", action="store", type="float", dest="initial_alpha",default = 0.0,help="initial_alpha,the weight of the smoothed l0 penalty. Defaut is 0 ")
  parser.add_option("--final_alpha", action="store", type="float", dest="final_alpha",default = 0.0,help="final_alpha,the weight of the smoothed l0 penalty. Defaut is 0 ")
  parser.add_option("--beta", action="store", type="float", dest="beta",default = 0.5,help="beta, the smoothness of the l0 prior, smaller the sigma, closer the approximation to true L0. Default is 0.5")
  parser.add_option("--slack", action="store_true", dest="slack_option",default = False,help="if you want to project on the simplex with slack.")
  parser.add_option("--noslack", action="store_false", dest="slack_option",default = False, help="if you want to project on the simplex with no slack. This is the regular projection approach")
  parser.add_option("--num_pgd_iterations", action="store", type="int", dest="num_pgd_iterations",default = 10,help="Number of Projected Gradient Descent to run. Default is 100")
  parser.add_option("--eta", action="store", type="float", dest="eta",default = 0.1,help="Eta, the constant step size for PGD using armijo line search. Default is 0.1")
  parser.add_option("--armijo_beta", action="store", type="float", dest="armijo_beta",default = 0.5,help="Set value for Armijo beta, the beta used in in armijo line search. Default value is 0.2")
  parser.add_option("--armijo_sigma", action="store", type="float", dest="armijo_sigma",default = 0.5,help="Set value for Armijo sigma, the sigma used in in armijo line search. Lower bound is 0.0001")
  parser.add_option("--lower_bound", action="store", type="float", dest="lower_bound",default = 0.000001,help="Set value for the lower bound on the probability.Default is 10E-6")
  parser.add_option("--cipher_data_file", action="store", type="string", dest="cipher_data_file",default ='cipher.data',help="Cipher data file for training")
  parser.add_option("--cipher_noq_file", action="store", type="string", dest="cipher_noq_file",default ='cipher.noq',help="Cipher data without quotes")
  parser.add_option("--cipher_decode_file", action="store", type="string", dest="cipher_decode_file",default ='cipher.decode',help="Cipher file for decoding")
  parser.add_option("--cipher_gold_file", action="store", type="string", dest="cipher_gold_file",default ='cipher.gold',help="The correct decipherment")
  parser.add_option("--lm", action="store", type="string", dest="lm",default ='lm.carmel',help="The lm file")
  parser.add_option("--noe_lm", action="store", type="string", dest="noe_lm",default ='lm.carmel',help="The noe lm file")
  parser.add_option("--gaussian_params_file", action="store", type="string", dest="gaussian_params_file",default =None,help="The means and variances for each gaussian letter")

  parser.add_option("--unigram_probs_file", action="store", type="string", dest="unigram_probs_file",default =None,help="The means and variances for each gaussian letter")
  parser.add_option("--std_mult", action="store", type="float", dest="std_mult",default =1.,help="The multiplier for the std ")
  parser.add_option("--u", action="store_true", dest="uniform_init",default=False)
  parser.add_option("--g", action="store_true", dest="gaussian_init",default=False)
  parser.add_option("--i", action="store_true", dest="identity_init",default=False)
  parser.add_option("--hpc", action="store_true", dest="hpc",default=False)
  parser.add_option("--full_fst", action="store_true", dest="full_fst",default=False)
  parser.add_option("--fst_init", action="store_true", dest="fst_init",default=False)

  parser.add_option("--noise_lm", action="store", type="string", dest="noise_lm",default ='noise.lm.carmel',help="The noise lm file")
  parser.add_option("--noise_channel_fst", action="store", type="string", dest="noise_channel_fst",default ='noise.fst.carmel',help="The channel fst")
  parser.add_option("--noise_probs_file", action="store", type="string", dest="noise_probs_file",default ='noise.probs',help="The noise probs file")
  parser.add_option("--noise_samples_file", action="store", type="string", dest="noise_samples_file",default ='noise.samples',help="The noise samples file")
  
  (options, args) = parser.parse_args()

  print options
  print args
  
  #getting the values from optparse
  num_iterations = options.num_iterations
  initial_alpha = options.initial_alpha
  final_alpha = options.final_alpha
  beta = options.beta
  slack_option = options.slack_option
  num_pgd_iterations = options.num_pgd_iterations
  eta = options.eta
  armijo_beta = options.armijo_beta
  armijo_sigma = options.armijo_sigma
  lower_bound = options.lower_bound
  cipher_data_file = options.cipher_data_file
  cipher_decode_file = options.cipher_decode_file
  cipher_noq_file = options.cipher_noq_file
  cipher_gold_file = options.cipher_gold_file
  lm = options.lm
  noe_lm = options.noe_lm
  gaussian_params_file = options.gaussian_params_file
  unigram_probs_file = options.unigram_probs_file
  uniform_init = options.uniform_init
  gaussian_init = options.gaussian_init
  identity_init = options.identity_init
  hpc = options.hpc
  std_mult = options.std_mult
  full_fst = options.full_fst
  num_noise_samples = options.num_noise_samples
  noise_lm = options.noise_lm
  noise_channel_fst = options.noise_channel_fst
  fst_init = options.fst_init
  noise_probs_file = options.noise_probs_file
  noise_samples_file = options.noise_samples_file
  noise_probs_file = options.noise_probs_file

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
  gold_cipher = emMethods.readCipherFile(cipher_gold_file)
  print gold_cipher
  #dictionary = emMethods.createDictionary('complete.dict.new-formatted')#.small')
  #word_lines= emMethods.readWordLines('test.words.new-formatted')
  cipher_letter_dict = emMethods.getUniqCipherLetters(cipher_noq_file)
  cipher_letters = open(cipher_noq_file).readline().split()


  #sys.exit()

  cipher_probs = defaultdict(float)
  emMethods.getCipherLetterProbs(cipher_probs,cipher_noq_file)
  print 'cipher probs are '
  #gaussians = defaultdict(list)
  #emMethods.readGaussians(gaussians,gaussian_params_file)
  plain_unigram_probs = dict((line.strip().split()[0],float(line.strip().split()[1])) for line in open(unigram_probs_file))

  #get cipher letter probs
  del cipher_letter_dict['_']
  #word_list_five = emMethods.readWordList('TEXT.3.linear')
  #plaintext = map(chr, range(97, 123))
  plaintext = []
  for k in range(65, 91):
    plaintext.append(chr(k))
  print plaintext
  print 'the number of unique cipher letter is %d'%len(cipher_letter_dict.keys())
  print cipher_letter_dict
  num_cipher_letters = len(cipher_letter_dict.keys()) 
  num_plain_letters = 26
  #gold_tag_sequence = emMethods.readWordList('test.tags.new-formatted.linear')
     
  free_parameters_channel = defaultdict(lambda:defaultdict(float))
  free_parameters_language = defaultdict(lambda:defaultdict(float))
  print 'starting to create parameters'
  total_language_parameters = 0
  total_channel_parameters = 0
  #for line in cipher_lines :
    #print 'created parameters for a line'
    #(language_parameters,channel_parameters) = emMethods.getFreeParametersBigram(line,dictionary,free_parameters_language,free_parameters_channel)
  if full_fst == True:
    emMethods.getFreeCipherParametersChannel(plaintext,free_parameters_channel)
    num_cipher_letters = 26
  else :
    emMethods.getFreeCipherParametersChannel(plaintext,free_parameters_channel,cipher_letter_dict = cipher_letter_dict)
  temp = {'_':0.0}
  free_parameters_channel['_'] = temp
  #now, we will build all the lattices, and create a special start node and end node for every sentence
  start_node_end_node_list = []

  fractional_counts_channel = copy.deepcopy(free_parameters_channel)
  probabilities_channel = copy.deepcopy(free_parameters_channel)
  #print 'gaussians'
  #print gaussians
  #createRandomPoint(probabilities_channel)
  if (uniform_init == True) :
    print 'uniform initialization'
    emMethods.initUniformProbs(probabilities_channel)
#  emMethods.initUniformProbs(probabilities_language,probabilities_channel)
  if (gaussian_init == True) :
    print 'gaussian initialization'
    #emMethods.initFromGaussians(probabilities_channel,gaussians,cipher_probs,std_mult)

    emMethods.initFromGaussiansSingleStd(probabilities_channel,plain_unigram_probs,cipher_probs,std_mult)
  if (identity_init == True) :
    emMethods.initIdentity(probabilities_channel)

  if fst_init == True:
    emMethods.initFromWfst('init.fst',probabilities_channel,'channel') 
  #print 'channel probabilities after weighting are '
  #print probabilities_channel
  #raw_input()

  emMethods.writeFst('cipher.fst',probabilities_channel)
  #sys.exit()
  final_probabilities_channel = copy.deepcopy(free_parameters_channel)


  #running the EM iterations
  #we are creating the indexes for algencan . Notice that here, the probabilities language is already uniform and therefore none of them will be zero
  createParametersForScaling(probabilities_channel,parameter_to_index)
  start_time = time.clock()
  print 'start time was ',start_time
  #fractional_counts_dump_file = open('fractional.params','w')
  #probabilities_dump_file = open('probs.params','w')
  #optimized_probabilities_dump_file = open('probs.optimized.params','w')
  alpha_delta = (final_alpha-initial_alpha)/(num_iterations-1)
  current_alpha = initial_alpha
   
  ###########GENERATING THE NOISE FILE LIST##########################################
  noise_file_list = []
  for k in range(num_noise_samples):
    noise_file_list.append("%s.noise%d"%(cipher_noq_file,k))
  
  #reading the noise probsa
  noise_probs = []
  for line in open(noise_probs_file):
    if line.strip() == '':
      continue
    else:
      noise_probs.append(float(line.strip()))
  print 'number of noise probs is ',len(noise_probs)
  #this will be needed for computations
  log_num_noise_samples = math.log(num_noise_samples)

  #reading noise samples
  total_noise_samples = 0
  noise_samples = []
  for line in open(noise_samples_file):
    if line.strip() == '':
      continue
    else:
      noise_samples.append(line.strip())
    total_noise_samples += 1
  print 'total number of noise samples is',total_noise_samples
  ###############GENERATE CARMEL TRAINING AND DECODING COMMANDS##############################
  print 'Generating the training and decoding commands...'
  run_true_training = ''
  run_noise_training = ''
  if hpc == True:
    run_true_training = "/home/nlg-05/vaswani/graehl/carmel/bin/linux64/carmel --train-cascade -u -M 0 -m -HJ %s %s cipher.fst"%(cipher_data_file,lm)
    run_noise_training = "/home/nlg-05/vaswani/graehl/carmel/bin/linux64/carmel --train-cascade -u -M 0 -m -HJ %s %s %s"%(cipher_data_file,noise_lm,noise_channel_fst)
  else :
    run_true_training = "/Users/avaswani/graehl/carmel/bin/macosx/carmel --train-cascade -u -M 0 -m -HJ %s %s cipher.fst"%(cipher_data_file,lm)
    run_noise_training = "/Users/avaswani/graehl/carmel/bin/macosx/carmel --train-cascade -u -M 0 -m -HJ %s %s %s"%(cipher_data_file,noise_lm,noise_channel_fst)
  
  noise_sample_training_commands_for_model = []
  #noise_sample_training_commands_for_noise = []
  for k in range(num_noise_samples):
    run_noise_sample_true_training = ''
    #run_noise_sample_noise_training = ''
    if hpc == True:
      run_noise_sample_true_training = "/home/nlg-05/vaswani/graehl/carmel/bin/linux64/carmel --train-cascade -u -M 0 -m -HJ %s %s cipher.fst"%(noise_file_list[k],lm)

      #run_noise_sample_noise_training = "/home/nlg-05/vaswani/graehl/carmel/bin/linux64/carmel --train-cascade -u -M 0 -m -HJ %s %s %s"%(noise_file_list[k],noise_lm,noise_channel_fst)
    else :
      run_noise_sample_true_training = "/Users/avaswani/graehl/carmel/bin/macosx/carmel --train-cascade -u -M 0 -m -HJ %s %s cipher.fst"%(noise_file_list[k],lm)

      #run_noise_sample_noise_training = "/Users/avaswani/graehl/carmel/bin/macosx/carmel --train-cascade -u -M 0 -m -HJ %s %s %s"%(noise_file_list[k],noise_lm,noise_channel_fst)
    noise_sample_training_commands_for_model.append(run_noise_sample_true_training)
    #noise_sample_training_commands_for_noise.append(run_noise_sample_noise_training)

  viterbi_command = '' 
  if hpc == True:
    viterbi_command = "cat %s | /home/nlg-05/vaswani/graehl/carmel/bin/linux64/carmel  -u  -srbk 1 -QEWI %s  cipher.fst > decipherment_output"%(cipher_decode_file,noe_lm)
  else :
    viterbi_command = "cat %s | /Users/avaswani/graehl/carmel/bin/macosx/carmel  -u  -srbk 1 -QEWI %s  cipher.fst > decipherment_output"%(cipher_decode_file,noe_lm)
  
  print 'Running Training ...'
  print 'command for running noise training is ',run_noise_training
  log_prob_under_noise = runCarmel(run_noise_training)
  '''
  log_noise_sample_probs_under_noise = []
  for k in range(num_noise_samples):
    #log_noise_sample_prob_under_noise = runCarmel(noise_sample_training_commands_for_noise[k])
    log_noise_sample_probs_under_noise.append(log_noise_sample_prob_under_noise)
  print 'the log noise sample probs under noise are ',log_noise_sample_probs_under_noise 
  '''
  ##########################RUN NCE TRAINING###############################
  for i in range (0,num_iterations) :

    ###############GENERATING K NOISE SAMPLES#########################
    print 'Generating noise samples...'
    noise_ciphertext = list(cipher_letters)
    for k in range(num_noise_samples):
      noise_file = open(noise_file_list[k],'w')
      #emMethods.generateNoiseCipher(noise_ciphertext)
      noise_ciphertext = noise_samples[i*num_noise_samples+k]
      #print 'noise ciphertext is ',' '.join(noise_ciphertext)
      print 'noise ciphertext is ',noise_ciphertext
      #noise_file.write("\n%s\n"%' '.join(["\"%s\""%item for item in noise_ciphertext]))
      noise_file.write("\n%s"%noise_ciphertext)
      noise_file.close()

    print 'Getting Accuracy'
    ################RUNNING VITERBI AND GETTING ACCURCY###################
    runViterbiAndGetAccuracy(gold_cipher,viterbi_command)


    #dampening the learning rate
    eta = eta/(1+i)
    current_function_value = 0.
    print 'the iteration number was ',i
    print 'current alpha is ',current_alpha
  
    print 'Computing expected counts...'
    ###############COMPUTING EXPECTED COUNTS ###############################
    log_prob_under_model = runCarmel(run_true_training)
    print 'reading language fractional counts from cipher.fst.trained'
    emMethods.readCarmelFractionalCounts('cipher.fst.trained',fractional_counts_channel,'channel')

    print 'Running the noise model'
    log_prob_under_noise = runCarmel(run_noise_training)


    positive_weight = computePositiveWeight(log_prob_under_model,log_prob_under_noise,num_noise_samples)

    #compute the 
    p_d_equals_one = log_prob_under_model - (logaddexp(log_num_noise_samples+log_prob_under_noise,log_prob_under_model))
    current_function_value += p_d_equals_one
    print 'p of d equals 1 at the current point is ', p_d_equals_one
    print 'The positive weight was',positive_weight
    if positive_weight <= 10E-10:
      positive_weight = 10E-10
    expected_counts = zeros(shape=(num_plain_letters,num_plain_letters))
    emMethods.dictionaryToArray(expected_counts,fractional_counts_channel,parameter_to_index)
    expected_counts *= positive_weight


    total_noise_expected_counts = zeros(shape=(num_plain_letters,num_plain_letters))

    log_noise_sample_probs_under_noise = []
    #log_noise_sample_probs_under_noise = []
    for k in range(num_noise_samples):
      temp_noise_fractional_counts = copy.deepcopy(free_parameters_channel)
      log_noise_sample_prob_under_model = runCarmel(noise_sample_training_commands_for_model[k])
      print 'reading noise language fractional counts from cipher.fst.trained'
      emMethods.readCarmelFractionalCounts('cipher.fst.trained',temp_noise_fractional_counts,'channel')
      #log_noise_sample_prob_under_noise = runCarmel(noise_sample_training_commands_for_noise[k])
      log_noise_sample_prob_under_noise = noise_probs[i*num_noise_samples+k]
      log_noise_sample_probs_under_noise.append(log_noise_sample_prob_under_noise)
      negative_weight = computeNegativeWeight(log_noise_sample_prob_under_model,log_noise_sample_probs_under_noise[k],num_noise_samples)
      if negative_weight <= 10E-10 :
        negative_weight = 10E-10
      print 'The negative weight was',negative_weight
      p_d_equals_zero = log_noise_sample_probs_under_noise[k]+ log_num_noise_samples -\
          (logaddexp(log_num_noise_samples+log_noise_sample_probs_under_noise[k],log_noise_sample_prob_under_model))
      current_function_value += p_d_equals_zero
      print 'p of d equals 0 at the current point',k,' is ', p_d_equals_zero

      temp_noise_expected_counts = zeros(shape=(num_plain_letters,num_plain_letters))
      emMethods.dictionaryToArray(temp_noise_expected_counts,temp_noise_fractional_counts,parameter_to_index)
      #print 'temp noise expected counts are ',temp_noise_expected_counts
      temp_noise_expected_counts *= negative_weight 
      #print 'temp noise exptected counts are ',temp_noise_expected_counts
      total_noise_expected_counts += temp_noise_expected_counts

    print 'the log noise sample probs under noise are ',log_noise_sample_probs_under_noise 
    #print 'total noise expected counts are ',total_noise_expected_counts
    #getting the final exptected counts
    expected_counts -= total_noise_expected_counts
  
    print 'Current function value is ',current_function_value
    ###############PROJECTING THE POINT ONTO THE DOUBLY STOCHASTIC MATRIX###########################
    new_feasible_point,current_point,grad = nceUpdate(eta,expected_counts,probabilities_channel,parameter_to_index,num_plain_letters)
 
    print 'Running line search....'
    #############LINE SEARCH#########################
    armijo_bound = 0.0
    #print 'new feasible point is ',new_feasible_point
    #raw_input()
    #print 'gradient is ',grad
    #raw_input()
    armijo_bound = -armijo_sigma * armijo_beta * (grad*(new_feasible_point-current_point)).sum()
    print 'Armijo bound is ',armijo_bound
    terminate_line_srch = False
    num_steps = 1
    #current_beta = 1.0 #armijo_beta
    current_beta = armijo_beta
    final_beta = 0
    current_armijo_bound = armijo_bound
    no_update = True 
    best_func_value = current_function_value

    while(terminate_line_srch != True) :
      #print 'num steps is ',num_steps
      temp_point = zeros(shape=current_point.shape)
      temp_point = current_point * (1.0 - current_beta) + current_beta * new_feasible_point
      #print 'the temp point is '
      #print temp_point
      #evaluating the function at the current point
      #first writing the fst
      #assigning the probabilities and writing the fsa
      temp_probabilities_channel = copy.deepcopy(free_parameters_channel)
      assignProbsMatrix(temp_point,temp_probabilities_channel,parameter_to_index)
      temp_probabilities_channel['_']['_'] = 1.0
      emMethods.writeFst('cipher.fst',temp_probabilities_channel)
      #print 'temp point is ',temp_point
      #sys.exit()
      func_value_at_temp_point = 0.
      temp_log_prob_under_model = runCarmel(run_true_training)
      temp_p_d_equals_one = temp_log_prob_under_model - (logaddexp(log_num_noise_samples+log_prob_under_noise,temp_log_prob_under_model))
      print 'p of d equals 1 at the temp point is ', temp_p_d_equals_one
      func_value_at_temp_point += temp_p_d_equals_one

      #then go over the noise samples 
      for k in range(num_noise_samples):
        temp_log_noise_sample_prob_under_model = runCarmel(noise_sample_training_commands_for_model[k])
        temp_p_d_equals_zero = log_noise_sample_probs_under_noise[k] + log_num_noise_samples -\
            (logaddexp(log_num_noise_samples+log_noise_sample_probs_under_noise[k],temp_log_noise_sample_prob_under_model))

        func_value_at_temp_point += temp_p_d_equals_zero
        print 'p of d equals 0 at the temp point',k,' is ', temp_p_d_equals_zero
      
      print 'the function value at the temp point  is %.16f'%func_value_at_temp_point
      #raw_input()
      if func_value_at_temp_point > best_func_value :
        best_func_value = func_value_at_temp_point
        final_beta = current_beta
        no_update = False
        #print 'we arrived at a better function value'
        #raw_input()
      #  print 'we just updated thef final beta to ',final_beta
      #if (func_value_at_temp_point - current_function_value >= current_armijo_bound) :
      #  terminate_line_srch = True
      #elif current_function_value - func_value_at_temp_point < 0:
      #  terminate_line_srch = True
      if num_steps >= 5 :
        terminate_line_srch = True
      
      current_beta = armijo_beta * current_beta
      current_armijo_bound =  armijo_bound * current_beta
      num_steps += 1
 
    if no_update == False :
      current_point = (1.0-final_beta)*current_point + final_beta*new_feasible_point
    if no_update == True :#x.all() == current_point.all() :
      print 'not update was true'
      break;
    #print 'the current point is ',current_point
    #raw_input()
    #assigning the probabilities and writing the fsa
    assignProbsMatrix(current_point,probabilities_channel,parameter_to_index)
    emMethods.writeFst('cipher.fst',probabilities_channel)


    #print 'checking the initial zeros in channel model'
    #checkZeros(probabilities_channel)
  
    fractional_counts_channel = copy.deepcopy(free_parameters_channel)
    final_probabilities_channel = copy.deepcopy(probabilities_channel)
    print 'at the end of the iteration'
    current_alpha += alpha_delta

  elapsed_time = time.clock() - start_time
  print 'the elapsed time was ',elapsed_time


if __name__ == "__main__" :
  #asd
  main()
