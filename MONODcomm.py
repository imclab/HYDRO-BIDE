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


"""  UNDER HEAVY CONSTRUCTION. BIOGEOGRAPHIC MONOD MODEL  """

"""  The following is a simulation-based combination of ecological neutral theory and
     the Monod chemostat model of microbial growth. There is 1 source, 1 community,
     and 1 outlet. The model is spatially implicit. Like the Monod model, growth is
     influenced by the concentration of a limiting substrate in the medium.
     
     MONOD EQUATION: mu = mu_max * ( S / (Ks + S))
     
     S = concentration of a limiting substrate or resource
     mu = specific growth rate, increase in cell mass per unit time, grams of cells per gram of cell per hour
     mu_max = max. specific growth rate, nutrient-saturated reproductive rate
     Ks = half-saturation coefficient; value of S when mu/mu_max = 0.5., the concentration at which µ is one-half
     of its maximum
     
     note: mu_max and Ks are empirical coefficients that can differ among species and ambient environmental conditions 
     
     Here, growth as in many neutral and community ecological models is per capita. population and community dynamics
     are NOT necessarily neutral but they ARE stochastic. This maintains a non-deterministic form of competition for
     limiting resources. This means that individuals are still picked at random to immigrate, die, reproduce, emigrate,
     etc. but according to probabilities that may or may not differ among species (i.e. neutral or non-neutral).  
     
     The model will simulate a community under varying resource concentration. The inflow of resources will remain
     constant, and so will the volume. But, due to changing community size, composition, and the number of dormant
     and active individuals, the concentration of limiting resource or substrate will change from one time period
     to the next. 
     
     Since the community will have active and dormant portions that vary in size and composition through time, and
     since the community will not be at a steady state, specific growth rate and dilution rate (volume/flow rate)
     will not be equivalent.
     
     The specific growth rate (mu) will actually be a probability of reproducing. This probability will be influenced
     by the concentration of resources in the environment (S). S will then be driven by the rate of inflow (r),
     volume of the environment (V), and concentration of resources in the inflowing medium, as well as by the size of
     the active portion of the community in the previous time step.  """

    
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


""" more variables """       
c = 1.0 # constant of proportionality; will eventually be used to create taxa differences
R += res_rate # R increases due to inflow
R -= c*num_A  # Starting assumption: total resources (R) decreases
              # in direct proportion to the number of active individuals

a = 1.0 # constant of proportionality; will eventually be used to create taxa differences    
ind_res = R/num_A # per capita resource availability
ind_grow = a*ind_res # proportion growth towards reproductive viability achieved per unit time (not constant)
                     # This makes per capita growth rate directly proportional to per capita resource availability.
                     # So, in general, individuals grow according to their share of resources
                
time = 400   # length of the experiment
burnin = int(0.75*time)  # number of initial time steps to discard

""" lists to track changes over time """

N_COBcom = []   # list to track N over time
A_COBcom = []   # list to track number of active individuals over time
pcr_COBcom = [] # list to tack per capita resources over time
R_COBcom = []   # list to track total resources over time
D_COBcom = []   # list to track dormancy

RAD_Ahigh = []  # lists to hold the community at various stages
RAD_Alow = [] 
RAD_Amedium = [] 


""" Having set up the community, it's time to turn it loose. """
t = 0
while t <= time: 
    #print 'time',t,'immigrants',im_rate,' ','size=',N,'per capita resources=',ind_res,'active',num_A#,ct
    
    """ inflow of individuals and resources """
    COBcom, num_a = hm.immigration(COBcom, im_rate, bp, lgp=0.7) # add some propagules to the community
    N = len(COBcom)
    # recalculate parameter values
    in_or_out = 'in'
    num_A, R, ind_res, ind_grow = hm.params_neutral(V, r, a, num_A, num_a, R, ind_res, res_rate, in_or_out)
    
    """" Community responds to inflow """
    # Simulation should reflect that the flow of individuals and resources
    # into the COB occurs independently of the community dynamics inside.
    COBcom, ind_res, num_A, N, ind_grow, R = hm.loop_thru_neutral_comm(COBcom, ind_res, dorm_lim, num_A, N, ind_grow, R, a)
    N = len(COBcom)
    
    """ outflow of individuals and resources """
    COBcom, num_a = hm.death_emigration(COBcom, N, V, r)
    N = len(COBcom)
    #recalculate parameter values """
    in_or_out = 'out'
    num_A, R, ind_res, ind_grow = hm.params_neutral(V, r, a, num_A, num_a, R, ind_res, res_rate, in_or_out)
    
    """ recording community info from time-steps """
    if t >= burnin:# and t%10 == 0: # allow a burn-in
        # call function to add information to lists, will be used for plotting changes, e.g. across time
        N_COBcom, A_COBcom, D_COBcom, pcr_COBcom, R_COBcom, RAD_Ahigh, RAD_Alow, RAD_Amedium = hm.add_to_lists_neutral(COBcom, N, num_A, ind_res, R, N_COBcom, A_COBcom, D_COBcom, pcr_COBcom, R_COBcom, RAD_Ahigh, RAD_Alow, RAD_Amedium)
        
    t += 1
    random.shuffle(COBcom) # randomize the community, prevent artifacts from arising due to list order
    """ Here, we have completed one time interval of inflow/outflow """


""" Here, we end the experiment and the result is a huge list of small lists. Each small list
    looks something like [1, 2, 50.0] with the first index representing the taxa label, the
    second index representing active/dormant, & the third index representing % growth to reproductive
    viability. """

nRADs = 10
if len(RAD_Ahigh) > nRADs: RAD_Ahigh = random.sample(RAD_Ahigh,nRADs) 
if len(RAD_Amedium) > nRADs: RAD_Amedium = random.sample(RAD_Amedium,nRADs) 
if len(RAD_Alow) > nRADs: RAD_Alow = random.sample(RAD_Alow,nRADs)

sets_of_RADS = [RAD_Ahigh,RAD_Amedium,RAD_Alow] # A list to hold lists that have captured the community
                                                # at times of high, low, and medium activity
    

# plot fig1
hm.fig1(N_COBcom,A_COBcom,D_COBcom,R_COBcom,pcr_COBcom,time,burnin,sets_of_RADS,V,r,res_dens,prop_dens,dorm_lim,nRADs)
# write the community to a file
#hm.comm_to_file(COBcom)
