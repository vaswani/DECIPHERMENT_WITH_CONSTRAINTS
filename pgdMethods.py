from numpy import *
import math

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
			#	temp_point[index]= (1.0 - current_beta) * current_point[index] + current_beta * new_feasible_point[index]
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

