# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <markdowncell>

# # This is the simplest hydrobide model. It does not include dormancy.
# ## It includes resource limited growth, but does not include a parameter for cell maintenance. As a result, cell quota does not decrease; it can only increase. 
# 
# # Hydrobide: An individual-based community simulation of birth, immigration, death, emigration, resource limitation in a chemostat environment.
# 
# There is 1 source, 1 local community, and 1 outlet. All BIDE processes occur via random sampling, i.e. individuals are picked at random to reproduce, immigrate, etc.
# The model is spatially implicit and individual-based. Information on every individual's Q and species identity is tracked. Resource concentration in the environment is also tracked. 
# 
# Some things are dynamic and will vary or arise naturally: Community size (N), composition, demography, cell quota (Q), resource concentration (R),
# nutrient uptake rate, biomass, biomass turnover. Others can remain constant across time: inflow rate (r), Volume (V). R influence resource uptake (v). v influences cell quota Q,
# and Q influences the per capita probability of reproduction. The rates of these processes are typiclally conceived to rely on the following:
#      
#     GROWTH/REPRODUCTION: mu = mu_max*(1 - Qmin/Q)/(1 - Qmin/Qmax)
#         mu is growth rate
#         mu_max is maximal (nutrient saturated) growth rate
#         Q is the cell quota (nutrient per cell, cannot equal 0)
#         Qmin is the minimal cell quota needed for life
#              
#     UPTAKE: v/vmax = R/(K+R)  (basically monod)
#         v is the rate of uptake (what we want to know)
#         vmax = maximal rate of uptake
#         R is external resource concentration (which our model tracks)
#         K is the half saturation constant, i.e. value of R when v/vmax = 0.5

# <markdowncell>

# ### Import some modules

# <codecell>

import sys                                            
import  matplotlib.pyplot as plt
import numpy as np                                    
from scipy import stats                                 
import random                                           
from random import choice, randrange

# <markdowncell>

# ### Main abiotic constraints of interest (empirical, not tunable): Volume (V), inflow rate (r)    

# <codecell>

Vs = [100.0, 120.0, 130.0]      # volume                                                                                           

# <markdowncell>

# ### Other empirical constraints, i.e. not tunable

# <codecell>

R_cons = 0.4  # inflowing resource concentration
I_cons = 0.4  # propagule concentration of inflowing medium

# <markdowncell>

# ### Tunable biotic constraints ...too many of these and we can draw an elephant

# <codecell>

Qmin = 0.0 # min cell quota

# <markdowncell>

# ### Declare the number of time steps and the initial number of steps to discard (for burn-in)

# <codecell>

time = 1000  # length of the experiment
burnin = int(0.5*time)  # number of initial time steps to discard

Nlists = []
Rlists = []
TRlists = []
TOlists = []
rlists = []
REStimes = []
Srichness = []

# <markdowncell>

# # BEGIN THE SIMULATION

# <codecell>

for V in Vs:
    rs = range(int(V/20),int(V),int(V/20))
    rs.reverse()
    AVGN = []
    AVGR = []
    AVGTR = []
    AVGTO = []
    AVGS = []
    REStime = []
    for r in rs:
        Q = 0
        I = int(round(I_cons * r)) # propagules immigrating per unit time  
        R = R_cons # initial resource concentration in the chemostat equals the inflowing concentration
        TR = R*V   # total resources in the local environment
        t = 0
        tlist = []
        COM = [] # the community will be a list of lists
        Nlist = []   # track N
        Rlist = []   # track resource concentration
        TRlist = []  # track total resources
        TOlist = []
        Slist = []
        while t <= time: 
            _in = I
            Srich = []
            """ inflow of resources, Immigration of propagules """
            lgp = 0.70 # arbitrarily set log-series parameter
            props = np.random.logseries(lgp, I)  # An initial set of propagules; a list of log-series distributed integers.
            # Assume the source community is infinite. Because large communities are
            # approximately log-series distributed (most things are rare and relatively
            # few things are abundant) immigration of propagules will occur in a log-series
            # distributed fashion.
            for prop in props: 
                Q = 0#float(random.randrange(100))/100 # Q of propagules is a uniform random variable between 0.0 and 1.0
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
                COM[i][1] = Q # change the individual's cell quota
                    
                if Q >= 1.0: # active, reproducing, i.e. Q reaches or exceeds normalized Qmax
                    COM.append([p[0],Q/2.0]) # individuals produce sister cells with half Q
                    COM[i][1] = Q/2.0
                    _in += 1
                   
            """ outflow of resources, Emigration """
            _out = int(round((N/V) * r)) # no. individuals lost per unit time
            if _out > 0:
                ct = 0
                while ct < _out:
                    i = randrange(0,len(COM)) # randomly pick an individual
                    COM.pop(i) # the individuals dies or emigrates
                    ct+=1
            N = len(COM) # community size is reduced
            
            for i, p in enumerate(COM): # Find species richness
                Srich.append(p[0])
            Srich = set(Srich)
            S = len(Srich)
            
            """ recording community info from time-steps """
            if t >= burnin: # allow a burn-in and only print/record info every so many steps
                #print 'Volume',V,'inflow',r,'COM size: ',N,' total resources: ',round(TR,3),' time:',t
                """ the above print statement is useful for checking how the community changes as constraint values are changed """
                tlist.append(t)
                Nlist.append(N) 
                Rlist.append(R)
                TRlist.append(TR)
                Slist.append(S)
                TOlist.append(abs(_in - _out))
                #print I,_in,' ',_out,' ',N
            t+=1
            
        avgN = np.mean(Nlist)
        AVGN.append(np.log(avgN))
        print 'V:',V,'r:',r,' avgN:',int(avgN),
        
        avgR = np.mean(Rlist)
        AVGR.append(avgR)
        #print ' avgR:',round(avgR,3),
        
        avgTR = np.mean(TRlist)
        AVGTR.append(avgTR)
        #print ' avgTR:',round(avgTR,3),
        
        REStime.append(np.log(V/r))
        print 'RT:',round(V/r,2),
        
        avgTO = np.mean(TOlist)
        AVGTO.append(avgTO)
        print ' avgTO:',round(avgTO,2),
        
        avgS = np.mean(Slist)
        AVGS.append(avgS)
        print ' avgs:',int(avgS)
        
        
    Nlists.append(AVGN)
    Rlists.append(AVGR)
    TRlists.append(AVGTR)
    TOlists.append(AVGTO)
    rlists.append(rs)
    REStimes.append(REStime)
    Srichness.append(AVGS)
    print len(rs),len(AVGN),len(AVGR),len(AVGTO),len(REStime),len(AVGS) # these lists should have equal lengths

# <markdowncell>

# ### generate figures of change over time

# <codecell>

print 'generating figures'
fig = plt.figure(figsize=(11.0,11.0))
fs = 12 # fontsize
colors = ['b','r','g']

ax = fig.add_subplot(2,2,1)
for i, Nlist in enumerate(Nlists):
    plt.plot(REStimes[i],Nlist,c=colors[i],label=str(int(Vs[i])))
    plt.ylim(0,max(Nlist)+0.1*max(Nlist))
plt.xlabel("log Residence time",fontsize=fs)
plt.ylabel("log(Total abundance)",fontsize=fs)
leg = plt.legend(loc=1, prop={'size':12})
leg.draw_frame(False)

ax = fig.add_subplot(2,2,2)
for i, Rlist in enumerate(Rlists):
    plt.plot(REStimes[i],Rlist,c=colors[i],label=str(int(Vs[i])))
plt.xlabel("log Residence time",fontsize=fs)
plt.ylabel("Resource density",fontsize=fs)
leg = plt.legend(loc=1, prop={'size':12})
leg.draw_frame(False)

ax = fig.add_subplot(2,2,3)
for i, Slist in enumerate(Srichness):
    plt.plot(REStimes[i],Slist,c=colors[i],label=str(int(Vs[i])))
plt.xlabel("log Residence time",fontsize=fs)
plt.ylabel("Species richness",fontsize=fs)
leg = plt.legend(loc=1, prop={'size':12})
leg.draw_frame(False)

ax = fig.add_subplot(2,2,4)
for i, TOlist in enumerate(TOlists):
    plt.plot(REStimes[i],TOlist,c=colors[i],label=str(int(Vs[i])))
plt.xlabel("log Residence time",fontsize=fs)
plt.ylabel("Biomass turnover",fontsize=fs)
leg = plt.legend(loc=1, prop={'size':12})
leg.draw_frame(False)

plt.show

# <codecell>


