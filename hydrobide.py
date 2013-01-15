#!/usr/local/bin/python                                                                   

import sys                                            
import os                                                  
import  matplotlib.pyplot as plt     
sys.path.append("/home/kenlocey/hydrobide_modules")
import hydrobide_modules as hm
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from pylab import *                                       
import numpy as np                                    
from scipy import stats                                 
import random                                           
from random import choice                          
import re                                             
from decimal import *                                 
import math                                    
from random import randrange                


"""  This script runs nicely, and appears to operate bug free. yippie!
     It will need to be checked for artifacts and mis-encoded processes. """

"""  The following is a simulation-based ecological neutral model for a                                                       
     chemostat/organ/bioreactor (COB) scenario. There is 1 source, 1 community, and
     1 outlet. This model incorporates dormancy and per capita resource limitation.
     The eventual goal will be to include species, taxa, and trait differences.
     The model is spatially implicit
     
     Volume (V) and influent rate (r) are modeled as influencing multiple resource
     and community-related rates and dynamics. The incoming concentration of
     resources and density of propagules are also accounted for. Influent rate (r)
     is the constraining rate on all passively in-flowing stuff. 

     This model is not zero-sum, i.e. community size (N) is not fixed and can
     fluctuate wildly. There is no mutation rate and no speciation, but there could be.
     Individuals are modeled such that they can grow to a reproductively viable state.
     The rate of per capita growth is then influenced by resource availability (R)
     and R is influenced by V, r, and the concentration of resources in the inflowing medium. 
     In this way, the model is similar to the Monod growth model. """     

    
V = 1000.0      # volume                                                                                           
r = 10.0        # influent rate                                                             
res_dens = 1.0 # growth limiting resource concentration of inflowing medium
prop_dens = 1.0 # propagule density, (cells per unit volume of inflowing medium)                     

dorm_lim = 0.01  # dormancy threshold; dormancy is undertaken if per capita resource availability
bp = 0.5 # arbitrarily set binomial probability that a propagule is active

rates = hm.get_in_rates(r,prop_dens,res_dens) # get immigration rate and resource delivery rate
im_rate = rates[0]  # immigration rate, (cells/unit time)                                     
res_rate = rates[1] # resource delivery rate, (unit resource/unit time) 

R = V * res_dens # Amount of resources in the environment before inoculation.

""" The community (COBcom) will be a list of lists: """
COBcom = []
comlist = hm.immigration(COBcom, im_rate, bp,lgp = 0.7) # inoculate the community with propagules
COBcom = comlist[0] # The community
N = len(COBcom)     # size of the community
num_A = comlist[1]  # number of active individuals.                                      

""" Things that will change as the community changes """         
c = 1.0 # constant of proportionality; will eventually be used to create taxa differences
R += res_rate # R increases due to inflow
R -= c*num_A  # Starting assumption: total resources (R) decreases
              # in direct proportion to the number of active individuals

a = 1.0 # constant of proportionality; will eventually be used to create taxa differences    
ind_res = R/num_A # per capita resource availability
ind_grow = a*ind_res # proportion growth towards reproductive viability achieved per unit time (not constant)
                     # This makes per capita growth rate directly proportional to per capita resource availability.
                     # So, in general, individuals grow according to their share of resources
                
""" Having set up the community, it's time to turn it loose. """
    
time = 400   # length of the experiment
burnin = int(0.75*time)  # number of initial time steps to discard

N_COBcom = []   # list to track N over time
A_COBcom = []   # list to track number of active individuals over time
pcr_COBcom = [] # list to tack per capita resources over time
R_COBcom = []   # list to track total resources over time
D_COBcom = []   # list to track dormancy

RAD_Ahigh = []  # lists to hold the community at various stages
RAD_Alow = [] 
RAD_Amedium = [] 

t = 0
while t <= time: 
    #print 'time',t,'immigrants',im_rate,' ','size=',N,'per capita resources=',ind_res,'active',num_A#,ct
    
    """ inflow of individuals, i.e., immigration """
    comlist = hm.immigration(COBcom, im_rate, bp, lgp=0.7) # add some propagules to the community
    COBcom = comlist[0]
    
    """ recalculate parameter values """
    N = len(COBcom)     
    num_A += comlist[1] 
    R += res_rate        
    ind_res = R/num_A  # per capita resource availability may change
                       # according to a change in the active portion
                       # of the COB community
    ind_grow = a*ind_res 
    
    """" The community responds to the inflow & changes """
    # Simulation should reflect that the flow of individuals and resources
    # into the COB occurs independently of the community dynamics inside.
       
    for i, v in enumerate(COBcom):
        
        if v[1] == 1:  # if the individual is active
            if ind_res <= dorm_lim: # if per capita resource availability <= dormancy threshold, 
                                    # then the individual can go dormant or starve and die
                
                x = choice([1,2]) # assume 50/50 chance of going dormant 
                                  # this could be made to vary among taxa
                if x == 2:
                    COBcom[i][1] = 2 # go dormant
                    num_A -= 1
                elif x == 1: 
                    COBcom.pop(i) # starve and die
                    N -= 1
                    num_A -= 1
                    
            else: # if there are enough resources to grow or reproduce
                if v[2] >= 100.0:
                    COBcom.append([v[0],1,0.0]) # reproduce if mature, offspring are active
                    COBcom[i][2] = 0.0          # one individual produces two sister cells at growth level 0  
                    num_A += 1
                    
                else:
                    COBcom[i][2] += ind_grow # grow if not mature
                    R -= ind_grow
                    if num_A != 0: ind_res = R/num_A
                    ind_grow = a*ind_res
                    
        elif v[1] == 2: # if the individual is dormant
            if ind_res > dorm_lim: # if per capita resource availability > the dormancy threshold
                COBcom[i][1] = 1 # go active
                num_A += 1
    
    """ outflow of individuals, i.e., death/emigration """
    N = len(COBcom)
    num_out = int(round(N/V * r)) # no. individuals lost per unit time
    num_a = 0
    if num_out >= 1:
        comlist = hm.death_emigration(COBcom, num_out)
        COBcom = comlist[0] 
        num_a = comlist[1]  
    
    """ recalculate parameter values """
    N = len(COBcom)  
    num_A -= num_a 
    R -= (R/V) * r 
    ind_res = R/num_A  
    ind_grow = a*ind_res 
    
    if t >= burnin:# and t%10 == 0: # allow a burn-in
        N_COBcom.append(np.log(N)) # using natural logs when values can be enormous
        A_COBcom.append(np.log(num_A))
        D_COBcom.append(np.log(N-num_A))
        pcr_COBcom.append(ind_res)
        R_COBcom.append(np.log(R))
        
        percent_A = num_A/float(N)
        if percent_A >= 0.66:
            RAD = hm.get_rad(COBcom)
            RAD_Ahigh.append(RAD)
        elif percent_A < 0.33:
           RAD = hm.get_rad(COBcom)
           RAD_Alow.append(RAD)
        else:
           RAD = hm.get_rad(COBcom)
           RAD_Amedium.append(RAD)

    t += 1
    random.shuffle(COBcom) # randomize the community, prevent artifacts from arising due to list order
    """ Here, we have completed one time interval of inflow/outflow """

nRADs = 10
if len(RAD_Ahigh) > nRADs: RAD_Ahigh = random.sample(RAD_Ahigh,nRADs) 
if len(RAD_Amedium) > nRADs: RAD_Amedium = random.sample(RAD_Amedium,nRADs) 
if len(RAD_Alow) > nRADs: RAD_Alow = random.sample(RAD_Alow,nRADs)

sets_of_RADS = [RAD_Ahigh,RAD_Amedium,RAD_Alow] # A list to hold lists that have captured the community
                                                # at times of high, low, and medium activity
    
""" Here, we end the experiment and the result is a huge list of small lists. Each small list
    looks something like [1, 2, 50.0] with the first index representing the taxa label, the
    second index representing active/dormant, & the third index representing % growth to reproductive
    viability. """

# plot fig1
hm.fig1(N_COBcom,A_COBcom,D_COBcom,R_COBcom,pcr_COBcom,time,burnin,sets_of_RADS,V,r,res_dens,prop_dens,dorm_lim,nRADs)
# write the community to a file
#hm.comm_to_file(COBcom)
