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


"""  UNDER HEAVY CONSTRUCTION. BIOGEOGRAPHIC MONOD MODEL

     None of the Python code below actually accounts (yet) for any of the following comments.  """

"""  The following is (or will be) a simulation-based combination of ecological neutral theory and
     the Monod chemostat model of microbial growth. There is 1 source, 1 community, and 1 outlet.
     The model is spatially implicit. Like the Monod model, growth is influenced by the concentration
     of a limiting substrate in the medium. It would be nice if this model wasn't so constrained as to only
     apply to microbes and aquatic organisms.
     
     MONOD EQUATION: mu = mu_max * S/(Ks + S)
     
     S = concentration of a limiting substrate or resource
     mu = specific growth rate, increase in cell mass per unit time, grams of cells per gram of cell per hour
     mu_max = max. specific growth rate, nutrient-saturated reproductive rate
     Ks = half-saturation coefficient; S when mu/mu_max = 0.5., concentration at which µ is one-half its maximum
     
     Here, community dynamics are NOT necessarily neutral but they DO have a stochastic component that maintains a
     non-deterministic competition for limiting resources among many potentially distinct species. This means that
     individuals are still picked at random to immigrate, reproduce, emigrate, etc. but according to
     probabilities that may or may not differ among species (i.e. neutral or non-neutral).  
     
     The inflow of a limiting resource and propagules will remain constant through time and so will the volume
     of the environment. Volume, inflow rate, and concentration of incoming resources will influence S.
     S influences specific growth rate (mu).
     
     In a time period when metabolically active growing individuals would reproduce, during which, the lack of
     reproduction implies dormancy and not death (i.e. outflow will capture death & emigration):
         
         1. Let Pr be the per capita probability of reproducing. Then the probability of being dormant (Pd)
            during the time period is one minus the probability of reproducing, Pd = 1-Pr
                  
         2. Let Smax be the smallest resource concentration required for mu_max
                Smax = 2*Ks
             
         3. At a concentration S at or above Smax, mu = mu_max & Pr = 1.0
               i.e. everybody is expected to grow and reproduce during a period when everybody has
                    the opportunity, i.e. Malthusian growth
         
         4. Let there be a threshold on resource concentration (Smin) at which, growth/reproduction cannot occur.
            Individuals persist without growing/reproducing (i.e. effectively dormant).
              
         Now, in an infinitely large population, the frequency of active individuals could equal the 
         probability of being active individual. Likewise, the frequency of dormant individuals 
         would equal one minus the active frequency, i.e. Pd = (D/N) = 1 - (A/N) = 1 - Pr   
         
         if S = Smax, then mu = mu_max & Pr = 1
         if S = Smin, then mu = 0 & Pr = 0
         if Smin < S < Smax, then  0 < mu < mu_max  &  0 < Pr < 1
         So, let Pr = mu/mu_max 
                 
         Proposition:
         We can use the concentrations of S, Smax, & Smin to derive Pr.
         We can use Pr to 1.) sample the community to induce dormancy and reproduction
                          2.) derive expected size of dormant and active portions 
                          3.) link constraints like V and r to 1 & 2
                          4.) derive mu
                          5.) induce immigration from the source
         
         Q. How do we find Pr?
         A. Through S, Smax, & Smin
                 
         Account for the floor (Pr = 0 when S = Smin):
             subtract Smin from each term:
                 S' = S - Smin
                 Smax' = Smax - Smin
                 
         Account for the ceiling (Pr = 1 when S = Smax):
             if S' > Smax':
                 Pr = (S' - (S' - Smax')) / Smax'
             if Smin < S <= Smax: 
                 Pr = S'/Smax'
                         
         Finally:
             Pr = [(S' + Smax') - max(S', Smax')]/Smax' 
             mu = mu_max * Pr
                     
     The presumed benefit to dormancy, in this Monod-type biogeographic/community model,
     is that not having enough resources to grow/reproduce does not lead to death.
     Actually, in this model, the only thing that leads to death is getting flushed out.
     
     Dilution rate (r/V) will determine the rate of emigration+death (ED). That is,
     all individuals dormant and active have the same probability of leaving the 
     community. The number of individuals flowing out will vary with community size.
     Highly active populations are more buffered against dilution than highly dormant ones.
     
     Consequently, this model (as is) does not allow an active/dormant tradeoff:
         i.e.  active => maybe reproduce but maybe die in a harsh environment 
              dormant => don't reproduce but don't die in a harsh environment
     
     Perhaps, a harsh environment should be built-in because, right now, there are no drawback to being active. 
     Perhaps, we should allow something in the environment that, at a certain level, kills active things.
     
     Regardless of examining the dormancy/activity tradeoff, this model could capture aspects of:
     
         Tilman (2004). Niche tradeoffs, neutrality, and community structure: A stochastic theory of resource
         competition, invasion, and community assembly. PNAS, 101:10854-10861.
     
             Tilman (2004) recognizes the importance of niche differences and stochastic dynamics. This Monod-type
             model allows niche differences and random drift.
     
         Pueyo S, Fangliang H & Zillio T (2007) The maximum entropy formalism and the idiosyncratic theory
         of biodiversity. Ecology Letters, 10:1017-1028.
     
             Whereas neutral theory ignores all species differences, the idiosyncratic theory of Pueyo et al. (2007) 
             ignores all species similarities. These represent the two ends of a continuum of possible models that
             produce realistic log-series like community structure. This Monod-type model can vary between totally
             neutral and totally idiosyncratic.
     
       None of the code below actually accounts (yet) for any of the above comments. """

    
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
