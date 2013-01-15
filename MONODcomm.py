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
     Ks = half-saturation coefficient; S when mu/mu_max = 0.5., concentration at which µ is one-half its maximum
     
     Here, community dynamics are NOT necessarily neutral but they DO have a stochastic component that maintains a
     non-deterministic competition for limiting resources among many distinct species. This means that
     individuals are still picked at random to immigrate, reproduce, emigrate, etc. but according to
     probabilities that may or may not differ among species (i.e. neutral or non-neutral).  
     
     The inflow of a limiting resource and propagules will remain constant through time and so will the volume
     of the environment. Volume, inflow rate, proportion of active individuals, and concentration of inflowing
     resources all influence resource concentration. Resource concentration influences specific growth rate (mu).
     Mu can be translated to a per capita probability of reproducing in a given time step (during which all
     individuals could theoretically reproduce if S was not limiting).   
     
     In a given time period when individuals can either reproduce or not, and the lack of reproduction
     implies dormancy and not death (outflow captures death + emigration):
         
         1. Let Pr be the per capita probability of reproducing
                  
         2. Let Smax be the smallest resource concentration required for mu_max
         
         3. At a concentration S at or above Smax, mu = mu_max & Pr = 1.0
               i.e. everybody is expected to reproduce during a period when everybody
                    has the opportunity, Malthusian growth
         
         4. Let there be a threshold on resource concentration (Smin) at which, growth/reproduction
              does not occur. In this case, individuals persist in the environment without growing/reproducing
              but do not die (i.e. effectively dormant).
              
         In a given time period when individuals can either reproduce or not, and the lack of reproduction
         implies dormancy and not death, this means:
                 
             The probability of being dormant (Pd) is one minus the probability of reproducing, Pd = 1-Pr 
                 
             In an infinitely large population (or neutral community), the probability of reproducing would 
             equal the portion of active individuals. Likewise, the probability of being dormant would equal
             one minus the probability of reproducing, i.e.
                 Pd = (D/N) = 1 - (A/N) = 1 - Pr   
         
             
             Knowing that:
                 at Smax, mu = mu_max & Pr = 1
                 at Smin, mu = 0 & Pr = 0
                 at Smin < S < Smax, 0 < mu < mu_max
                 
                 Proposition:
                 We can use S, Smax, Smin to derive Pr.
                 We can use Pr to 1.) sample the community, inducing dormancy and activity
                                  2.) derive expected size of dormant and active portions 
                                  3.) induce immigration (see below).
                 
                 Q. How do we find Pr?
                 A. Relate S, Smax, and Smin to Pr
                 
                 Account for the floor (Pr = 0 when S = Smin)
                     subtract Smin from each term:
                         S' = S - Smin
                         Smax' = Smax - Smin
                 
                 Account for the ceiling (Pr = 1 when S = Smax)
                     if S' > Smax':
                         Pr = (S' - (S' - Smax')) / Smax'
                     if Smin < S <= Smax: 
                         Pr = S'/Smax'
                         
                 Finally:
                     Pr = [(S' + Smax') - max(S', Smax')]/Smax' 
                     mu = mu_max * Pr
                     
     
     Dilution rate will determine rates of emigration and death, which will vary with the community size.
     When mu drops below a certain probability, things go dormant (or die)
     When mu is above a certain probability, dormant things can become active
     
     This model will capture aspects of:
     
     Tilman (2004). Niche tradeoffs, neutrality, and community structure: A stochastic theory of resource
     competition, invasion, and community assembly. PNAS, 101:10854-10861.
     
         Tilman (2004) recognizes the importance of niche differences and stochastic dynamics
     
     Pueyo S, Fangliang H & Zillio T (2007) The maximum entropy formalism and the idiosyncratic theory
     of biodiversity. Ecology Letters, 10:1017-1028.
     
         Whereas neutral theory ignores all species differences, the idiosyncratic theory of Pueyo et al. (2007) 
         ignores all species similarities. These represent the two ends of a continuum of possible models that
         produce realistic log-series like community structure. 
     
       """

    
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
