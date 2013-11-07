from numpy import *
import numpy
from emUtil import *
import math
from dstoc import *

def approximateProjectOntoSimplex(v) :
  array_size = v.size
  projection = array(v)
  sum =0.
  for i in range(array_size) :
    if projection[i] < 0 :
      projection[i] = 0.
    else :
      sum += projection[i]
  print 'sum was ',sum
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

#we're doing projected gradient descent with the armijo rule along feasible direction. Bertsekas: NonLinear programming. page 230
def projectedGradientDescentWithArmijoRule(x,expected_counts,num_pgd_iterations,eta,lower_bound,armijo_beta,armijo_sigma,alpha,beta,slack_option,parameter_counter) :
  
  #global parameter_counter,slack_option
  #current_point = array([0.0]*parameter_counter)
  current_point = array(x)
  #print 'the function value at the current point  is ',current_function_value
  for time in range(1,num_pgd_iterations +1) :
    current_function_value= evalFunction(current_point=current_point,expected_counts=expected_counts,alpha=alpha,beta=beta)
    #print 'the function value at the current point  is ',current_function_value

    grad = evalGradient(current_point=current_point,expected_counts=expected_counts,alpha = alpha,beta = beta) #getting the gradient at the current point
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
      #  temp_point[index]= (1.0 - current_beta) * current_point[index] + current_beta * new_feasible_point[index]
      temp_point = current_point * (1.0 - current_beta) + current_beta * new_feasible_point
      #print 'the temp point is '
      #print temp_point
      func_value_at_temp_point = evalFunction(current_point= temp_point,expected_counts = expected_counts,alpha = alpha,beta = beta)

      #print 'the function value at the current point  is %.16f'%current_function_value
      #raw_input()
      if func_value_at_temp_point < best_func_value :
        best_func_value = func_value_at_temp_point
        final_beta = current_beta
        no_update = False
        #print 'we arrived at a better function value'
        #raw_input()
      #  print 'we just updated thef final beta to ',final_beta
      if (current_function_value - func_value_at_temp_point >= current_armijo_bound) :
        terminate_line_srch = True
      #elif current_function_value - func_value_at_temp_point < 0:
      #  terminate_line_srch = True
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
      #  next_point[index] = (1.0 - final_beta) * current_point[index] + final_beta * new_feasible_point[index]
        temp_coordinate = (1.0 - final_beta) * current_point[index] + final_beta * new_feasible_point[index]
        #if next_point[index] == 0 :
        if temp_coordinate == 0 :
          print 'next point was 0 and final beta was ',final_beta, ' and the current point was ',current_point[index], ' and the feasible point was' ,new_feasible_point[index]
        current_point[index] = temp_coordinate
    
    if no_update == True :#x.all() == current_point.all() :
      print 'not update was true'
      break;
  return(current_point)

def evalFunction(expected_counts,current_point,alpha,beta) :
  #global alpha,beta
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

def evalGradient(expected_counts,current_point,alpha,beta) :
  #global alpha,beta
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

def pdgOnlyRowConstraints(probabilities_channel,fractional_counts_channel,parameter_to_index,num_pgd_iterations,alpha,beta,eta,lower_bound,armijo_beta,armijo_sigma,slack_option):
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

def pdgRowAndColumnConstraints(probabilities_channel,fractional_counts_channel,parameter_to_index,num_pgd_iterations,alpha,beta,eta,lower_bound,armijo_beta,armijo_sigma,num_plain_letters,num_cipher_letters):
  
  p = zeros(shape=(num_plain_letters,num_plain_letters))
  expected_counts = zeros(shape=(num_plain_letters,num_plain_letters))
  current_fractional_counts = fractional_counts_channel
  #first populating the expected counts and probabilities matrix
  #for i,plain_letter in enumerate(probabilities_channel.keys()) :
  for i,k in enumerate(range(65, 91)):
    plain_letter = chr(k)
    expected_counts_sum = 0.
    #print plain_letter
    #print 'number of cipher letters is ',len(probabilities_channel[plain_letter].keys())
    for cipher_letter in probabilities_channel[plain_letter].keys() :
      parameter_number = parameter_to_index[plain_letter][cipher_letter]
      #print 'parameter number is ',parameter_number
      #print 'i is ',i
      p[i][parameter_number] = probabilities_channel[plain_letter][cipher_letter]  
      expected_counts[i][parameter_number] = current_fractional_counts[plain_letter][cipher_letter]
      expected_counts_sum +=  current_fractional_counts[plain_letter][cipher_letter]
    #print 'expected counts sum was ',expected_counts_sum
    #print parameter_to_index

  #current_eta = eta_0 #/sqrt(num_iterations+1)
  print 'Doing projected gradient descent'
  new_probabilities = projectedGradientDescentWithArmijoRuleMatrix(p = p,expected_counts = expected_counts,num_pgd_iterations = num_pgd_iterations,eta = eta,lower_bound = lower_bound,armijo_beta = armijo_beta,armijo_sigma = armijo_sigma,alpha = alpha,beta = beta,num_cipher_letters = num_cipher_letters)
  print 'finished projected gradient descent'
  #new_probabilities = newtonProjectedGradientDescentWithArmijoRule(x,num_pgd_iterations,eta,lower_bound,armijo_beta,armijo_sigma)
  #print new_probabilities
  #raw_input()
  print 'we are replacing channel probs'
  assignProbsMatrix(new_probabilities,probabilities_channel,parameter_to_index)
  '''
  print 'the probabilities after em step are '
  for i in range(0,p.shape[0]):
    #print ' '.join([item[:3] for item in map(str,p[i])])
    numpy.set_printoptions(precision=3)
    #max = numpy.max(p,)
    print p[i]
    print 'max for letter ',chr(i+65),' is ',numpy.max(p[i])
  raw_input()
  for i in range(0,p.shape[0]):
    print 'max for letter ',chr(i+65),' is ',numpy.max(p[i])
  raw_input()
  '''
def projectedGradientDescentWithArmijoRuleMatrix(p,expected_counts,num_pgd_iterations,eta,lower_bound,armijo_beta,armijo_sigma,alpha,beta,num_cipher_letters) :
  current_point = array(p)
  #print 'on entering, the current point is ',
  #print current_point
  #print current_point.shape
  #raw_input()
  #print 'the expected counts are ',
  #print expected_counts
  #raw_input()
  #print 'the function value at the current point  is ',current_function_value
  for time in range(1,num_pgd_iterations +1) :
    current_function_value=evalFunctionMatrix(current_point=current_point,expected_counts=expected_counts,alpha=alpha,beta=beta,num_cipher_letters=num_cipher_letters)
    #print 'the function value at the current point  is ',current_function_value
    #raw_input()
    grad = evalGradientMatrix(current_point=current_point,expected_counts=expected_counts,alpha = alpha,beta = beta,num_cipher_letters = num_cipher_letters) #getting the gradient at the current point
    #print grad
    #print 'that was gradient'
    #raw_input()
    new_point = zeros(shape=current_point.shape)
    '''
    for index in range (0,parameter_counter) :
      new_point[index]= current_point[index]-eta*grad[index]
    '''
    new_point = current_point - eta*grad
    #print new_point
    #current_function_value = evalFunction()
    new_feasible_point,blah = dykstra(new_point,10E-6)
    new_feasible_point = array(new_feasible_point)
    #print 'new feasible point is '
    #print new_feasible_point
    #raw_input()
    #print current_point
    armijo_bound = 0.0
    armijo_bound = -armijo_sigma * armijo_beta * (grad*(new_feasible_point-current_point)).sum()
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
      temp_point = zeros(shape=current_point.shape)
      temp_point = current_point * (1.0 - current_beta) + current_beta * new_feasible_point
      #print 'the temp point is '
      #print temp_point
      func_value_at_temp_point = evalFunctionMatrix(current_point= temp_point,expected_counts = expected_counts,alpha= alpha,beta=beta,num_cipher_letters=num_cipher_letters)

      #print 'the function value at the current point  is %.16f'%current_function_value
      #raw_input()
      if func_value_at_temp_point < best_func_value :
        best_func_value = func_value_at_temp_point
        final_beta = current_beta
        no_update = False
        #print 'we arrived at a better function value'
        #raw_input()
      #  print 'we just updated thef final beta to ',final_beta
      if (current_function_value - func_value_at_temp_point >= current_armijo_bound) :
        terminate_line_srch = True
      #elif current_function_value - func_value_at_temp_point < 0:
      #  terminate_line_srch = True
      if num_steps >= 20 :
        terminate_line_srch = True
      
      current_beta = armijo_beta * current_beta
      current_armijo_bound =  armijo_bound * current_beta
      num_steps += 1
      
    #next_point = array([0.0]*x.size)
    #else :
    #print 'the number of steps was ',num_steps
    if no_update == False :
      current_point = (1.0-final_beta)*current_point + final_beta*new_feasible_point
      '''
      for index in range (0,parameter_counter) :
      #  next_point[index] = (1.0 - final_beta) * current_point[index] + final_beta * new_feasible_point[index]
        temp_coordinate = (1.0 - final_beta) * current_point[index] + final_beta * new_feasible_point[index]
        #if next_point[index] == 0 :
        if temp_coordinate == 0 :
          print 'next point was 0 and final beta was ',final_beta, ' and the current point was ',current_point[index], ' and the feasible point was' ,new_feasible_point[index]
        current_point[index] = temp_coordinate
      '''
    if no_update == True :#x.all() == current_point.all() :
      print 'not update was true'
      break;
  return(current_point)

def evalFunctionMatrix(current_point,expected_counts,alpha,beta,num_cipher_letters):
  #print expected_counts[:,0:num_cipher_letters].shape
  #print current_point[:,0:num_cipher_letters].shape
  #print expected_counts[:,0:num_cipher_letters]
  #print current_point[:,0:num_cipher_letters]
  #raw_input()
  #print 'sum of first term is ',(expected_counts[:,0:num_cipher_letters]*numpy.log(current_point[:,0:num_cipher_letters])).sum()
  #raw_input()
  func_val = (expected_counts[:,0:num_cipher_letters]*numpy.log(current_point[:,0:num_cipher_letters])).sum()  + (alpha*numpy.exp(-current_point[:,0:num_cipher_letters]/beta)).sum()
  return(-func_val)

def evalGradientMatrix(current_point,expected_counts,alpha,beta,num_cipher_letters) :

  gradient = zeros(current_point.shape)
  '''
  print 'num cipher letters are ',num_cipher_letters
  print 'first part of gradient is ',
  print expected_counts[:,0:num_cipher_letters]/current_point[:,0:num_cipher_letters]
  raw_input()
  '''
  gradient[:,0:num_cipher_letters] = expected_counts[:,0:num_cipher_letters]/current_point[:,0:num_cipher_letters] - alpha*numpy.exp(-current_point[:,0:num_cipher_letters]/beta)/beta
  return(-gradient)

