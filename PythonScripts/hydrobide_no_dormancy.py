# This model does not include dormancy.
# Hydrobide: An individual-based community simulation of birth, immigration, death,
# emigration, resource limitation, and dormancy in a chemostat environment.", 
    

# Import some modules"
import sys               
import  matplotlib.pyplot as plt
import numpy as np                                     
#import pylab as plt
from scipy import stats                                  
import random                                            
from random import choice, randrange 

# Main abiotic constraints of interest: Volume (V), inflow rate (r)    
V = 1000.0      # volume 
r = 10.0        # inflow rate

# Other empirical constraints, i.e. not tunable
R_cons = 0.37  # inflowing resource concentration 
I_cons = 0.36  # propagule concentration of inflowing medium

# Resulting variables
I = int(round(I_cons * r)) # propagules immigrating per unit time 
R = R_cons # initial resource concentration in the chemostat equals the inflowing
TR = R*V   # total resources in the local environment

# Tunable biotic constraints ...too many of these and we can draw an elephant
Qmin = 0.0 # min cell quota
K = 1.0     # half saturation constant, i.e. R when v = 0.5
m = 0.01   # resources needed for maintenance

# Declare the community; it will be a list of lists. Declare lists to track community features over time. There can also be lists to track diversity, evenness, etc.
COM = [] 
Nlist = []   # track N 
Rlist = []   # track resource concentration 
TRlist = []  # track total resources

# Declare the number of time steps and the initial number of steps to discard (for burn-in)
time = 1000  # length of the experiment 
burnin = int(0.5*time)  # number of initial time steps to discard

# BEGIN THE SIMULATION
t = 0
tlist = []

while t <= time: 
    """ inflow of resources, Immigration of propagules """
    lgp = 0.70 # arbitrarily set log-series parameter
    props = np.random.logseries(lgp, I)  # An initial set of propagules; a list of log-series distributed integers. Assume the source community is infinite. Because large communities are approximately log-series distributed (most things are rare and relatively few things are abundant) immigration of propagules will occur in a log-series distributed fashion.
    
    for prop in props: 
        Q = float(random.randrange(100))/100 # Q of propagules is a uniform random variable between 0.0 and 1.0
        i = [prop, Q] # an individual and its cell quota
        COM.append(i)
    N = len(COM) # size of the community
    TR += r * R_cons # Total Resources increase due to inflow
    R = TR/V # Resource concentration increases due to inflow
    
    """ Reproduction/Growth/Death """
    random.shuffle(COM) # randomize the community
    for i, p in enumerate(COM):  # loop through the randomized community
        v = R/(K+R)  # uptake (nutrients per cell per day)
        TR -= v      # total resources decreases
        R = TR/V     # new resource concentration  
        Q = p[1] + v  # increase cell quota by uptake
        Q -= m  # decrease cell quota by resources needed to maintain activity
        #print Q
        COM[i][1] = Q # change the individual's cell quota
        
        if Q >= 1.0: # active, reproducing, i.e. Q reaches or exceeds normalized Qmax
            COM.append([p[0],Q/2.0]) # individuals produce sister cells with half Q
            COM[i][1] = Q/2.0          
        elif Q <= 0: # dead
            COM.pop(i) # remove the individual

    """ outflow of resources, Emigration """
    num_out = int(round((N/V) * r)) # no. individuals lost per unit time
    if num_out > 0:
        ct = 0
        while ct < num_out:
            i = randrange(0,len(COM)) # randomly pick an individual
            COM.pop(i) # the individuals dies or emigrates
            ct+=1
    N = len(COM) # community size is reduced
    
    """ recording community info from time-steps """
    if t >= burnin and t%20 == 0: # allow a burn-in and only print/record info every so many steps
        print 'COM size: ',N,' total resources: ',round(TR,3),' time:',t
        """ the above print statement is useful for checking how the community changes as constraint values are changed """
        tlist.append(t)
        Nlist.append(N) 
        Rlist.append(R)
        TRlist.append(TR)
       
    t += 1
    
    
fs = 12 # fontsize
fig = pl.figure(figsize=(11.0,11.0))
ax = fig.add_subplot(2,2,1)
plt.plot(tlist,Nlist,c='k')
plt.ylim(0,max(Nlist)+0.1*max(Nlist))
plt.xlabel("time",fontsize=fs)
plt.ylabel("N",fontsize=fs)

ax = fig.add_subplot(2,2,2)
plt.plot(tlist,Rlist,c='m')
plt.xlabel("time",fontsize=fs)
plt.ylabel("Resource density",fontsize=fs)

ax = fig.add_subplot(2,2,3)
plt.plot(tlist,TRlist,c='g')
plt.xlabel("time",fontsize=fs)
plt.ylabel("Total resources",fontsize=fs)
plt.savefig('test.png')
print 'done'
