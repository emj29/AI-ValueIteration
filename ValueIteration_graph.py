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
			#print 'we got here'
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
	gamma=0.999
	M_max=11
	N_max=31
	value_time_ar=np.zeros([M_max,N_max])
	policy_time_ar=np.zeros([M_max,N_max])
	for M in xrange(1,M_max):
		for N in xrange(1,N_max):
			print 'Running for M %r and N %r' % (M,N)
			if (use_poisson):
				prob_dist=get_pmf(1,lambda_,3,2,0.001)
				mean=lambda_
				#print 'Using Poisson distribution with mean ',lambda_
				#print prob_dist
				#sys.exit(2)
			else:
				prob_dist=get_pmf(0,0,mu,sigma,0.001)
				t_m=0
				for s in prob_dist:
					t_m=t_m+s[0]*s[1]
				print 'Using Gaussian distribution with original mean',mu,'. \nAfter discretizing it, the new mean is %1.2f' % t_m
				mean=t_m


			Nm=M+1 #number of possible rooms that can be open (0,1,...,M)

			Np=N+prob_dist[-1][0]+1 #number of people in the waiting room + new arrivals (0,1,...,N+max(prob_dist))


			current_utility=np.zeros([3,Np,Nm]) #store the values of the utilities and actions
			#current_utility[0,:,:] --> utilities at time step i
			#current_utility[1,:,:] --> utilities at time step i+1
			#current_utility[2,:,:] --> best possible action so far for each state


			##VALUE ITERATION CODE

			max_error=mean*min(E3,P)/2 #error value used for convergence of value iteration
			#print max_error
			delta=1e10
			iteration=0
			value_start=tm.time()
			initial_delta=1e10
			print 'Performing Value Iteration...'
			while (delta>max_error):
				print 'iteration: ',iteration,' Algorithm Completion: ', (iteration>40)*(initial_delta-delta)/(initial_delta-max_error),'%'
				sys.stdout.write("\033[F") 
				#print np.mean(current_utility[0,:,:]/(iteration+1))
				#print 'utility mean %r and STD %r'% (np.mean(current_utility[0,:,:]/(iteration+1)),np.std(current_utility[0,:,:]/(iteration+1)))
				iteration=iteration+1
				if (iteration==40): initial_delta=delta
				for i in xrange(Np):
					for j in xrange(Nm):
						#print 'now on state: ',i,j
						delta=0
						min_val=1e100
						wait_room=(i>j)*min(i-j,N)
						#max_close=max(j-i,0)
						max_close=j
						free_rooms=max(0,j-i)
						for a in xrange(-max_close,M-j+1,1):
							cost=0.0
							for state in prob_dist:
								#print state
								prob=state[1]
								cost=cost+prob*current_utility[0,wait_room+state[0],j+a]
							cost=cost*gamma

							cost=cost+state_reward(j,i)+action_reward(a,free_rooms)
							#if (i==4 and j==4): print 'action ',a, 'value',cost
							if (cost<min_val):

								min_val=cost
								min_action=a
								current_utility[2,i,j]=a
							
						current_utility[1,i,j]=min_val ##here it is bitches?

						difference=(abs(current_utility[1,i,j]-current_utility[0,i,j]))
						#print 'difference is: ', difference
						if (difference>delta):
							delta=difference

				current_utility[0,:,:]=current_utility[1,:,:]
			print 'iteration: ',iteration,' Value Iteration Algorithm Complete!: '

			value_time=tm.time()-value_start
			value_time_ar[M,N]=value_time
			value_iteration_policy=current_utility[2,:,:]
			#uncomment below to print out the final optimal policy (actions for each state)
			# for i in xrange(Np):
			# 	for j in xrange(Nm):
			# 		print('For rooms %r and people %r, action is %1.0f'%(j,i,current_utility[2,i,j]))
			# #print the average utility and STD
			# print current_utility[0,:,:]
			# print 'Average utility is %r with and STD of %r' % (np.mean(current_utility[0,:,:]*(1-gamma)),np.std(current_utility[0,:,:]*(1-gamma)))




			#POLICY ITERATION CODE

			##Initialize all utilities and actions
			current_utility=np.zeros([3,Np,Nm])

			##Uncomment below if you want to start with a random policy
			# for i in xrange(Np):
			# 	for j in xrange(Nm):
			# 		current_utility[2,i,j]=M-j
			iteration=0;
			converged=False #converged will tell you if our policy has converged or not
			k=50
			policy_start=tm.time()
			print 'Performing Policy Iteration....'
			while (not converged):
				if (iteration>30): break
				iteration=iteration+1
				delta=1e10
				iter_p=0
				#First step is to perform a policy evaluation given our current policy
				for i in xrange(k):
				#while (delta>4*max_error/float(iteration)):
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
			
								delta=diff

					current_utility[0,:,:]=current_utility[1,:,:]
				##now we check for the new updated best policy 
				converged=True ##the moment we have an update, set it to false
				print 'Policy Update # ',iteration
				sys.stdout.write("\033[F") 
				for i in xrange(Np):
					for j in xrange(Nm):
						min_val=1e100
						wait_room=(i>j)*min(i-j,N)
						max_close=j
						free_rooms=max(0,j-i)

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
							current_utility[2,i,j]=min_action
							converged=False
			
			#uncoment below if you wish to print the optimal policy
			# for i in xrange(Np):
			# 	for j in xrange(Nm):
			# 		print('For rooms %r and people %r, action is %1.0f'%(j,i,current_utility[2,i,j]))
			policy_iteration_policy=current_utility[2,:,:]
			policy_time=tm.time()-policy_start
			policy_time_ar[M,N]=policy_time
			##check for equality
			print '                               '
			print 'Policy Iteration Complete!!!!'
			print 'Elapsed time for value iteration: %3.2f' % value_time,'seconds'
			print 'Elapsed time for policy iteration: %3.2f' % policy_time,'seconds'
			print 'Value and Policy iteration converged to the same optimal policy: ',np.array_equal(value_iteration_policy,policy_iteration_policy)
			#tm.sleep(2)

			#Code to allow for querying a particular state and returning the action
			# while (True):
			# 	Nm_input=int(raw_input("Enter the number of open rooms (0,%r): ('-1' to quit)" % M))
			# 	if (Nm_input==-1): break
			# 	Np_input=int(raw_input("Enter the number of total people (0,%r): "%(Np-1)))
			# 	#print('For rooms %r and people %r, action is %2.2f'%(Nm_input,Np_input,current_utility[2,Np_input,Nm_input]))
			# 	action_r=int(policy_iteration_policy[Np_input,Nm_input])

			# 	#print action_r
			# 	if (action_r>0):
			# 		action_s="to open %r rooms" % action_r
			# 	elif(action_r<0):
			# 		action_s="to close %r rooms" % abs(action_r)
			# 	else:
			# 		action_s="to stay put"
			# 	print 'Best action is ',action_s

			for i in xrange(Np):
				for j in xrange(Nm):
					if (value_iteration_policy[i,j]!=policy_iteration_policy[i,j]):
						print 'for %r and %r, vlaue gives %r and policy gives %r' % (i,j,value_iteration_policy[i,j],policy_iteration_policy[i,j])
	print value_time_ar,policy_time_ar
	np.save('val_time',value_time_ar)
	np.save('pol_time',policy_time_ar)


if __name__=="__main__":
	main(sys.argv[1:])

