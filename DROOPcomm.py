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


"""  FUTURE SITE OF A BIOGEOGRAPHIC DROOP/CELL QUOTA MODEL. CHECK BACK LATER  """

"""  The following is (or will be) a simulation-based combination of ecological neutral theory and
     the Droop/Cell Quota chemostat model of microbial growth. There is 1 source, 1 community, and 1 outlet.
     The model is spatially implicit. Like the Droop model, growth is influenced by the intracellular
     concentration of a limiting substrate. It would be nice if this model wasn't so constrained as to only
     apply to microbes and aquatic organisms.
     
     DROOP EQUATION: mu = mu_max * ( 1 - kq/Q)
     
     mu = specific growth rate, increase in cell mass per unit time, grams of cells per gram of cell per hour
     mu_max = max. specific growth rate, growth rate at infinite cell quota, theoretical maximal reproductive rate
     Q = cell quota, intracellular concentration of a growth limiting resource
     kq or Qmin = the minimum quota necessary for life (or activity)
     
     A normalized form of the Droop equation makes the parameter mu_max meaningful, such that when Q = Qmax then
     mu = mu_max (max growth achieved at max cell quota):
     
     RELATIVIZED DROOP EQUATION: mu_rel = mu/mu_max = (1 - Qmin/Q) / (1 - Qmin/Qmax)
     
     Question: for purely theoretical work, should it matter whether C quota or cell quota is used?
     
     Other: Nr = specific nutrient uptake rate, relates concentration in the environment to cell quota
     
     
     The importance of the Droop/Cell Quota model is not only in that growth is primarily influenced by
     intracellular concentrations, but that growth can continue after resource concentrations outside the
     cell have been greatly depleted.
     
     
     
     Here, community dynamics are NOT necessarily neutral but they DO have a stochastic component that maintains a
     non-deterministic competition for limiting resources among many potentially distinct species. This means that
     individuals are still picked at random to immigrate, reproduce, emigrate, etc. but according to
     probabilities that may or may not differ among species (i.e. neutral or non-neutral).  
     
     The inflow of a limiting resource and propagules will remain constant through time and so will the volume
     of the environment. Volume, inflow rate, and concentration of incoming resources will influence the
     concentration of resource in the evironment. 
     S influences specific growth rate (mu).
     
     In a time period when metabolically active growing individuals would reproduce, during which, the lack of
     reproduction implies dormancy and not death (i.e. outflow will capture death & emigration):

"""



""" Things that will remain constant through time, the values of which, must be reasonable
    or the COB will crash, explode, or worse 

    Some values for cow rumen obtainable here: http://microbewiki.kenyon.edu/index.php/Bovine_Rumen """
    
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
comlist = immigration(COBcom, im_rate, bp) # inoculate the community with propagules using the above function
COBcom = comlist[0] # The community
N = len(COBcom)     # size of the community
num_A = comlist[1]  # number of active individuals in COBcom 
                    # Each list will contain the individual species label and reveal whether the          
                    # individual is dormant or active, and how close the individuals is to          
                    # being reproductively viable.                                                  

#print COBcom
#sys.exit()
            
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
    
time = 800   # length of the experiment
burnin = 700  # number of initial time steps to discard

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
    comlist = immigration(COBcom, im_rate, bp) # add some propagules to the community
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
        comlist = death_emigration(COBcom, num_out)
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
            RAD = get_rad(COBcom)
            RAD_Ahigh.append(RAD)
        elif percent_A < 0.33:
           RAD = get_rad(COBcom)
           RAD_Alow.append(RAD)
        else:
           RAD = get_rad(COBcom)
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



""" NOTE: Everything below pertains to plotting and graphing.
    This code will eventually be moved to a module KL will
    create for this project """

fig = plt.figure(figsize=(10.0,8.0))


ax = plt.subplot2grid((2,2), (0,0), rowspan=1) # plotting N, dormancy, & activity through time
plt.plot(N_COBcom,'0.5',label='Total')
plt.plot(A_COBcom,'r',label='Active')
plt.plot(D_COBcom,'b',label='Dormant')
ymax = max(N_COBcom)+0.1
ymin = min(min(A_COBcom),min(D_COBcom))-0.1
plt.ylim(ymin,ymax) 
plt.xlabel("Time",fontsize=12)
plt.ylabel("ln(abundance)",fontsize=12)
plt.title("EQ: N ~constant, D & A change in sync",fontsize=12) 
# Add legend
leg = plt.legend(loc=10,prop={'size':12})
leg.draw_frame(False)

ax = plt.subplot2grid((2,2), (0,1), rowspan=1) # plotting N, R, & per capita resources through time
plt.ylim(-0.1,max(R_COBcom)+0.5)
plt.plot(R_COBcom, 'b', label='ln(total resources)')
plt.plot(pcr_COBcom, 'r', label='per capita resources')
plt.xlabel("Time",fontsize=12)
plt.ylabel("Value",fontsize=12)
# add some vertical gridlines
t_range = range(0,(time - burnin),10)
for t in t_range:
    plt.axvline(x=t,color='0.80',ls='--',lw=1) # plot a vertical line at the mode
# Add legend
leg = plt.legend(loc=10,prop={'size':12})
leg.draw_frame(False)

ax = plt.subplot2grid((2,2), (1,0), rowspan=1) # plotting activity vs. dormancy
plt.scatter(A_COBcom, D_COBcom, c='0.4', lw=0.5)#, label='active vs. dormant')
plt.title("Each point represents a time step",fontsize=12) 
plt.xlabel("ln(active abundance)",fontsize=12)
plt.ylabel("ln(dormant abundance)",fontsize=12)

ax = plt.subplot2grid((2,2), (1,1), rowspan=1) # plotting activity vs. dormancy
ct = 0
colors = ['r','0.4','b']
series = ['high activity', 'medium activity', 'high dormancy'] 
for RADs in sets_of_RADS:
    plt.plot([0],[0],color=colors[ct],label=series[ct],lw=3)
    for RAD in RADs:
        rank = range(1,len(RAD)+1)
        plt.plot(rank, RAD, c=colors[ct], lw=0.5)
    ct += 1
# add labels
plt.xlabel("Rank",fontsize=12)
plt.ylabel("ln(abundance)",fontsize=12)
plt.title(str(nRADs)+' randomly chosen RADs per activity level',fontsize=12) 
# Add legend
leg = plt.legend(loc=1,prop={'size':12})
leg.draw_frame(False)


plt.subplots_adjust(wspace=0.2, hspace=0.3)
plt.savefig('results/COBcom V='+str(V)+' r'+str(r)+' resdens='+str(res_dens)+' propdens='+str(prop_dens)+' dormlim='+str(dorm_lim)+'.png', dpi=400, bbox_inches = 'tight', pad_inches=0.1) 

sys.exit()    
# write the list to a file
OUT = open('/results/COBcom V='+str(V)+' r'+str(r)+' resdens='+str(res_dens)+' propdens='+str(prop_dens)+' dormlim='+str(dorm_lim),'w+')
for _list in COBcom:
    print>>OUT, _list[0],_list[1],_list[2]    
OUT.close() 
