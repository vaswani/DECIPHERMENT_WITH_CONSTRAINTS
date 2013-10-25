#!/usr/bin/env python

import sys
data = []
dict = {}
file  = open('accuracies_file')
for line in file :
	line = line.strip()
	data.append(line)

for i in range(0,len(data)) :
	if i%2 == 1 :
		continue
	else :
		(alpha,sigma) = data[i].split(' ')
		print alpha
		print sigma
		accuracy = data[i+1]
		if dict.has_key(alpha) :
			if dict[alpha].has_key(sigma) :
				print 'could not have had teh same configuration repeat'
			else :
				#print 'we are adding' 
				#temp[sigma] = float(accuracy)
				dict[alpha][sigma] = float(accuracy)
		else :
			temp = {}
			temp[sigma] = float(accuracy)
			dict[alpha] = temp

#print dict
#raw_input()
sigma_list = ['0.5','0.25','0.075','0.05','0.025','0.0075','0.005']
alpha_list = ['10','20','30','40','50','60','70','80','90','100']

temp_string = ' & '.join(["%s"%sigma for sigma in sigma_list])
print temp_string

for alpha in alpha_list :
		temp_string = ' & '.join(["%.1f"%(dict[alpha][sigma]*100) for sigma in sigma_list])
		final_string = "%s &%s \\\\ "%(alpha,temp_string)
		print final_string
			#print dict[alpha][sigma]	
