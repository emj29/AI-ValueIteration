import copy
import numpy as np
import sys
import getopt
import math as mt
from random import randint
import time as tm

def usage():
  print "\nUsage function\n"

  print 'python ValueIteration.py [-h] [-d num] [-l num] [-m num] [-s num] [-o num] [-c num]\
   [-b num] [-N num] [-M num] [-W num] [-P num]'

  print '-h, --help \n\t invoke the usage information'
  print '-d num, --dist=num\n\t 1 for poisson, 0 for gaussian'
  print '-l num, --lambda=num\n\t Mean for the poisson distribution'
  print '-m num, --mu=num\n\t Mean for the Gaussian distribution'
  print '-s num, --sigma=num\n\t Standard Deviation for the Gaussian distribution'
  print '-o num, --E1=num\n\t Cost of opening a room'
  print '-c num, --E2=num\n\t Cost of closing a room'
  print '-b num, --E3=num\n\t Cost of operating a room'
  print '-N num\n\t Number of seats in the waiting room'
  print '-M num\n\t Maximum number of rooms that could be operating'
  print '-W num\n\t Cost of waiting per hour per patient'
  print '-P num\n\t Cost of rejecting a patient'

def state_reward(open_rooms, people):
	cost=open_rooms*E3
	waiting=(people>open_rooms)*min(people-open_rooms,N)
	rejected=people-open_rooms-waiting
	rejected=rejected*(rejected>0)
	cost=cost+(waiting*W)+rejected*P
	#print 'W is ',(waiting/float(N)*W)
	return cost

def action_reward(open_rooms,free_rooms):
	return abs(open_rooms)*((open_rooms>0)*(E1+E3)+(open_rooms<0)*(E2))\
	-min(free_rooms,abs(open_rooms))*(open_rooms<0)*E3

def get_pmf(poisson,lamb,mu,sigma,epi):
	value=1
	pmf=[]
	k=0
	sum_n=0
	before_mean=1;
	while (value>=0.001 or before_mean):
		if (poisson):
			value=lamb**k*mt.exp(-lamb)/mt.factorial(k)
		else:
			value=1/(mt.sqrt(2*mt.pi*sigma**2))*mt.exp((-(k-mu)**2)/(2*sigma**2))
		#print value
		#if (value>epi):
		if (value>=0.001):
			pmf.append([k,value])
		if (value>0.001): before_mean=0;
		k=k+1
		sum_n=0.0
		for x in pmf:
			sum_n=sum_n+x[1]
	sum_n=0.0
	for x in pmf:
		sum_n=sum_n+x[1]
	for i,item in enumerate(pmf):
		#print i,item
		pmf[i][1]=item[1]/sum_n
	return pmf

def get_poisson_pmf(lamb,epi):
	value=1
	poisson_pmf=[]
	k=0
	sum_n=0
	before_mean=1;
	while (value>=0.001 or before_mean):
		value=lamb**k*mt.exp(-lamb)/mt.factorial(k)
		print value
		#if (value>epi):
		if (value>=0.001):
			poisson_pmf.append([k,value])
		if (value>0.001): before_mean=0;
		k=k+1
		sum_n=0.0
		for x in poisson_pmf:
			sum_n=sum_n+x[1]
	sum_n=0.0
	for x in poisson_pmf:
		sum_n=sum_n+x[1]
	for i,item in enumerate(poisson_pmf):
		print i,item
		poisson_pmf[i][1]=item[1]/sum_n
	#poisson_pmf[:]=[x/sum(poisson_pmf) for x in poisson_pmf]

	return poisson_pmf

def get_normal_pmf(mu,sigma,epi):
	normal_pmf=[]
	value=1
	k=0
	before_mean=1;
	while(value>=0.001 or before_mean):
		value=1/(mt.sqrt(2*mt.pi*sigma**2))*mt.exp((-(k-mu)**2)/(2*sigma**2))
		print 'bottom', 1/(mt.sqrt(2*mt.pi*sigma**2))
		print 'top',(-((k-mu)**2)),'bottom',(2*sigma**2)
		print value
		if (value>0.001): before_mean=0
		if (value>=0.001):
			normal_pmf.append([k,value])
		k=k+1
	sum_n=0.0
	for x in normal_pmf:
		sum_n=sum_n+x[1]
	for i,item in enumerate(normal_pmf):
		print i,item
		normal_pmf[i][1]=item[1]/sum_n
	#normal_pmf[:,1]=[x[1]/1 for x in normal_pmf]

	return normal_pmf

def main(argv):

	try:
		opts,args=getopt.getopt(argv,"hd:l:m:s:o:c:b:N:M:W:P:",["help","dist=","lambda=","mu=","sigma=","E1=","E2=","E3="])
	except getopt.GetoptError:
		print 'bad input'
		sys.exit(2)

	for opt,arg in opts:
		if opt in ("-h","--help"):
			usage()
			sys.exit()
		elif opt in("-d","--dist"):
			global use_poisson
			use_poisson=int(arg)
		elif opt in("-l","--lambda"):
			global lambda_
			lambda_=float(arg)
		elif opt in("-m","--mu"):
			global mu
			mu=float(arg)
		elif opt in("-s","--sigma"):
			global sigma
			sigma=float(arg)
		elif opt in ("-N"):
			global N
			N=int(arg)
		elif opt in ("-M"):
			global M
			M=int(arg)
		elif opt in ("-o","--E1"):
			print 'we got here'
			global E1
			E1=float(arg)
		elif opt in ("-c","--E2"):
			global E2
			E2=float(arg)
		elif opt in ("-b","--E3"):
			global E3
			E3=float(arg)
		elif opt in ("-W"):
			global W
			W=float(arg)
		elif opt in ("-P"):
			global P
			P=float(arg)


	global gamma
	gamma=0.997


	#create the poisson probability mass functions
	#look at the true mean of the Gaussian being created
	# for i in range(1,10):
	# 	for j in range(1,10):
	# 		x=get_pmf(0,0,i/2.0,j/2.0,0.001)
	# 		t_m=0
	# 		for s in x:
	# 			t_m=t_m+s[0]*s[1]
	# 		print 'mean ',i/2.0, 'std ',j/2.0, 'true mean ',t_m
	# sys.exit(2)
	#print poisson.pmf(range(0,poisson.ppf(0.999,mu).astype(int)+1),mu)
	#print get_normal_pmf(float(11),float(1),0.001)
	#sys.exit(2)
	if (use_poisson):
		prob_dist=get_pmf(1,lambda_,3,2,0.001)
		mean=lambda_
	else:
		prob_dist=get_pmf(0,0,mu,sigma,0.001)
		mean=mu
	print prob_dist
	print prob_dist[-1][0]
	#sys.exit(2)
	#prob_dist=poisson_pmf
	#print prob_dist
	#sys.exit(2)
	#print range(0,poisson.ppf(0.999,mu).astype(int))

	Nm=M+1
	#print Nm
	Np=N+prob_dist[-1][0]+1
	print Np
	#print 'length of dist',prob_dist, len(prob_dist)
	#sys.exit(2)
	#create the state space 

	current_utility=np.zeros([3,Np,Nm])
	#print current_utility
	#initialie with current 

	##code for the value iteration
	##iterate through each state and perform the update
	#max_error=0.01*(1-gamma)/gamma #max error we tolerate for value iteration.
	max_error=mean*min(E3,P)/2 
	print mean,E3,P,max_error
	#sys.exit(2)
	#max_error=100
	delta=1e10
	iteration=0
	while (delta>max_error):
		print 'iteration: ',iteration,' delta: ',delta
		#print np.mean(current_utility[0,:,:]/(iteration+1))
		print 'utility mean %r and STD %r'% (np.mean(current_utility[0,:,:]/(iteration+1)),np.std(current_utility[0,:,:]/(iteration+1)))
		iteration=iteration+1
		for i in xrange(Np):
			for j in xrange(Nm):
				#print 'now on state: ',i,j
				delta=0
				min_val=1e100
				wait_room=(i>j)*min(i-j,N)
				#max_close=max(j-i,0)
				max_close=j
				free_rooms=max(0,j-i)
				#for a in xrange(-j,M-j+1,1):
				#print 'people ',i,' and rooms ', j, ' actions are ',xrange(-max_close,M-j+1,1)
				for a in xrange(-max_close,M-j+1,1):
					cost=0.0
					#for state in xrange(len(poisson_pmf)):
					for state in prob_dist:
						#print state
						prob=state[1]
						#print state[0]
						#print 'wait_room', wait_room
						#print current_utility.shape
						cost=cost+prob*current_utility[0,wait_room+state[0],j+a]
					cost=cost*gamma
					cost=cost+state_reward(j,i)+action_reward(a,free_rooms)

					if (cost<min_val):
						#print cost
						min_val=cost
						min_action=a
						current_utility[2,i,j]=a
				current_utility[1,i,j]=min_val
				difference=(abs(current_utility[1,i,j]-current_utility[0,i,j]))
				#print 'difference is: ', difference
				if (difference>delta):
					delta=difference

			

		current_utility[0,:,:]=current_utility[1,:,:]

	for i in xrange(Np):
		for j in xrange(Nm):
			print('For rooms %r and people %r, action is %1.0f'%(j,i,current_utility[2,i,j]))
	#print the average utility and STD
	print current_utility[0,:,:]
	print 'Average utility is %r with and STD of %r' % (np.mean(current_utility[0,:,:]*(1-gamma)),np.std(current_utility[0,:,:]*(1-gamma)))
	value_iteration_policy=current_utility[2,:,:]
	#write code to query


	# while (True):
	# 	Nm_input=int(raw_input("Enter the number of open rooms (0,%r):" % M))
	# 	Np_input=int(raw_input("Enter the number of total people (0,%r): "%(Np-1)))
	# 	#print('For rooms %r and people %r, action is %2.2f'%(Nm_input,Np_input,current_utility[2,Np_input,Nm_input]))
	# 	action_r=int(current_utility[2,Np_input,Nm_input])
	# 	#print action_r
	# 	if (action_r>0):
	# 		action_s="to open %r rooms" % action_r
	# 	elif(action_r<0):
	# 		action_s="to close %r rooms" % abs(action_r)
	# 	else:
	# 		action_s="to stay put"
	# 	print 'Best action is ',action_s


#create function to perform policy iteration

	##Initialize all utilities and actions
	current_utility=np.zeros([3,Np,Nm])
	##five a random policy
	# for i in xrange(Np):
	# 	for j in xrange(Nm):
	# 		current_utility[2,i,j]=randint(-j,M-j)
	iteration=0;
	converged=False
	k=4
	while (not converged):
		##for each state update the utilities
		print 'iteration', iteration
		print 'converged', converged
		iteration=iteration+1
		delta=1e10
		#for i in xrange(k):
		iter_p=0
		while (delta>3*max_error/float(iteration)):
			iter_p=iter_p+1
			delta=0
			max_diff=0
			for j in xrange(Np):
				for q in xrange(Nm):
					
					wait_room=(j>q)*min(j-q,N)
					cost=0.0
					#for state in xrange(len(poisson_pmf)):
					a=current_utility[2,j,q]
					for state in prob_dist:
						prob=state[1]

						cost=cost+prob*current_utility[0,wait_room+state[0],q+a]
					free_rooms=max(0,q-j)
					cost=cost*gamma
					cost=cost+state_reward(q,j)+action_reward(a,free_rooms)
					current_utility[1,j,q]=cost
					diff=abs(current_utility[1,j,q]-current_utility[0,j,q])
					if (diff>delta):
						#print diff
						delta=diff
			print iter_p,delta
			current_utility[0,:,:]=current_utility[1,:,:]
		##now we check for the new updated best policy 
		converged=True ##the moment we have an update, set it to false
		print 'find new best policy'
		for i in xrange(Np):
			for j in xrange(Nm):
				min_val=1e100
				wait_room=(i>j)*min(i-j,N)
				#max_close=max(j-i,0)
				max_close=j
				free_rooms=max(0,j-i)
				#for a in xrange(-j,M-j+1,1):
				#print 'people ',i,' and rooms ', j, ' actions are ',xrange(-max_close,M-j+1,1)
				for a in xrange(-max_close,M-j+1,1):
					cost=0.0
					#for state in xrange(len(poisson_pmf)):
					for state in prob_dist:
						prob=state[1]

						cost=cost+prob*current_utility[0,wait_room+state[0],j+a]
					cost=cost+state_reward(j,i)+action_reward(a,free_rooms)

					if (cost<min_val):
						#print cost
						min_val=cost
						min_action=a
						#current_utility[2,i,j]=a
				if (min_action!=current_utility[2,i,j]):
					print min_action, current_utility[2,i,j]
					current_utility[2,i,j]=min_action
					converged=False
	for i in xrange(Np):
		for j in xrange(Nm):
			print('For rooms %r and people %r, action is %1.0f'%(j,i,current_utility[2,i,j]))
	policy_iteration_policy=current_utility[2,:,:]
	##check for equality
	print 'Value and Policy iteration converged to the same optimal policy: ',np.array_equal(value_iteration_policy,policy_iteration_policy)


if __name__=="__main__":
	main(sys.argv[1:])

