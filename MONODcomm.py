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
                                          
                                  

"""  UNDER HEAVY CONSTRUCTION. FUTURE SITE OF BIOGEOGRAPHIC MONOD MODEL  """

"""  The following is a simulation-based combination of ecological neutral theory and
     the Monod chemostat model of microbial growth. There is 1 source, 1 community,
     and 1 outlet. Like the Monod model, growth is related to the concentration of 
     a limiting substrate in the medium.
     
     This model is biogeographic in that it incorporate immigration and emigration
     in an ecologically neutral way, pulling from a log-series distributed source
     community.
     
     Volume (V), influent rate (r), and hence residence time (V/r) are modeled
     as influencing multiple resource and community-related rates and dynamics.
     Dilution rate also needs to be accounted for. Influent rate (r) is the
     constraining rate on all passively in-flowing stuff (i.e. resources, propagules). 

     There is no mutation rate and no speciation, but there could be.
     This model is not zero-sum, i.e. community size (N) is not fixed and can fluctuate. 
     Community size = total cells
          
     Unlike ecological neutral models, individuals do not have to be 'born' into a
     state of reproductive maturity. Individuals are modeled such that they can
     grow to a reproductively viable state. Like the Monod model per capita growth
     is influenced by total available resources (R) and R is influenced by V, r, and
     the concentration of inflowing resources. """

     
V = 1000.0      # volume of the COB                                                                                                  
r = 10.0        # influent rate (unit volume/unit time)                                                                  

res_dens = 1.0 # growth limiting resource concentration of inflowing medium,
                # e.g. (grams cellulose + grams x + grams y)/ (liter of inflowing medium) 

prop_dens = 1.0 # propagule density, (cells per unit volume of inflowing medium)                  
                # Assume initially that resource concentration of the influent equals
                # the resource concentration of the COB. This makes sense if we're 
                # starting with a community of zero individuals.

dorm_lim = 0.01  # dormancy threshold; dormancy is undertaken if per capita resource availability
                # is below some threshhold (low resources -> low metabolism -> slow growth = go dormant)
                # This could be made to vary among species

bp = 0.5 # arbitrarily set binomial probability
         # the probability that a propagule is active

im_rate = int(round(prop_dens * r)) # immigration rate (cells/unit time)                                     
res_rate = res_dens * r # resource delivery rate, (unit resource/unit time) 

R = V * res_dens # Amount of resources in the COB before inoculation. If inflow = outflow and the COB is empty,
                 # then the initial resource concentration would equal that of the inflow and outflow. Eventually,
                 # total resources (R) will change as a result of a changing community.

""" The community (COBcom) will be a list of lists: """
COBcom = []
comlist = hm.immigration(COBcom, im_rate, bp,lgp = 0.7) # inoculate the community with propagules using the above function
COBcom = comlist[0] # The community
N = len(COBcom)     # size of the community
num_A = comlist[1]  # number of active individuals in COBcom 
                    # Each list will contain the individual species label and reveal whether the          
                    # individual is dormant or active, and how close the individuals is to          
                    # being reproductively viable.                                                  

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
                     
                
""" Having set up the community, it's time to turn it loose. From here on, community size, biomass,
    dormancy, compositional and noncompositional community structure, population structure,
    replacement, & turnover will all ride on random drift. """
    
time = 400   # length of the experiment
burnin = 300  # number of initial time steps to discard

N_COBcom = []   # list to track N over time
A_COBcom = []   # list to track number of active individuals over time
pcr_COBcom = [] # list to tack per capita resources over time
R_COBcom = []   # list to track total resources over time
D_COBcom = []   # list to track dormancy

RAD_Ahigh = []  # lists to hold the community at various stages
RAD_Alow = [] 
RAD_Amedium = [] 

t = 0
while t <= time: # looping one time unit at a time
    #ct = 0
    #for i in COBcom:
    #    if i[1] == 1: ct += 1
    print 'time',t,'immigrants',im_rate,' ','size=',N,'per capita resources=',ind_res,'active',num_A#,ct
    
    """ inflow of individuals, i.e., immigration """
    comlist = hm.immigration(COBcom, im_rate, bp, lgp=0.7) # add some propagules to the community
    COBcom = comlist[0]
    
    """ recalculate parameter values """
    N = len(COBcom)     # Total community size will increase
    num_A += comlist[1]  # total number of active individuals may increase
    
    R += res_rate        # total resources increases
    ind_res = R/num_A  # per capita resource availability may change
                       # according to a change in the active portion
                       # of the COB community
    ind_grow = a*ind_res # per capita growth rate will change according
                         # to a change in per capita resources
    
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
                    # individual growth decreases R, which decreases
                    # per capita resource availability and
                    # per capita growth rate
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
        COBcom = comlist[0] # the newly decreased community
        num_a = comlist[1]  # account for the number of lost active individuals
    
    """ recalculate parameter values """
    N = len(COBcom)  # total community size will have decreased
    num_A -= num_a # number of active individuals may decrease
    R -= (R/V) * r # total resources decreases according to the
                   # amount of resources lost to outflow 
    ind_res = R/num_A  # per capita resource availability may change
    ind_grow = a*ind_res # per capita growth rate may change
    
    
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
    random.shuffle(COBcom) # randomize the community, prevent artifacts
                           # from arising in the next iteration of the 
                           # for loop
    """ Here, we have completed one time interval of inflow/outflow """
#print len(RAD_Ahigh),len(RAD_Amedium),len(RAD_Alow)
nRADs = 10
if len(RAD_Ahigh) > nRADs: RAD_Ahigh = random.sample(RAD_Ahigh,nRADs) 
if len(RAD_Amedium) > nRADs: RAD_Amedium = random.sample(RAD_Amedium,nRADs) 
if len(RAD_Alow) > nRADs: RAD_Alow = random.sample(RAD_Alow,nRADs)

sets_of_RADS = [RAD_Ahigh,RAD_Amedium,RAD_Alow] # A list to hold lists that have captured the community
                                                # at times of high, low, and medium activity
    
""" Here, we end the experiment and the result is
    a huge list of small lists. Each small list
    looks something like [1, 2, 50.0] with the 
    first index representing the taxa label, the
    second index representing active/dormant, & the
    third index representing % growth to reproductive
    viability. """

# plot fig1
hm.fig1(N_COBcom,A_COBcom,D_COBcom,R_COBcom,pcr_COBcom,time,burnin,sets_of_RADS,V,r,res_dens,prop_dens,dorm_lim,nRADs)
# write the community to a file
#hm.comm_to_file(COBcom)
