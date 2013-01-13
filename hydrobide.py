#!/usr/local/bin/python                                                                    _____ 
#                                                                                        //  _  \
#                                                                                        || | |_|
import sys#                                              Esophagus                       || |  _  CHEMOSTAT
import numpy as np#                                         |                            || |_| |__ 
import  matplotlib.pyplot as plt#      (inflow of resources & immigration )              \\____/   \
import re#                                                  |                               || | | |
from pylab import *#                                      | V |                             || | | | ORGAN
from random import choice#                                |   |_____                        || |_| |__           
from random import randrange#                             |   |     \                       \\____/   \
import math#                                              |         |                          || |_| |
from scipy import stats#                                  / birth   | -biomass turnover?       ||  _ <  BIOREACTOR
#                                                        /    &     |                          || |_| |
#                                               ________/   death  /  -community structure?    ||____/
#                                              /  ___             /
#                                             /  /  \            /    -dormancy dynamic?
#                                              |     \__________/ 
#                                              |       
#                           (death, emigration, & loss of resources)      
#                                              |
#                                              V
#                                    To intestines & beyond  

"""  This script runs nicely, but needs to be checked for bugs.  """

"""  The following is a simulation-based ecological neutral model for a                                                       
     chemostat/organ/bioreactor (COB) scenario. There is 1 source, 1 community, and
     1 outlet. This model incorporates dormancy and per capita resource limitation.
     The eventual goal will be to include species, taxa, and trait differences.
     Accounting for temperature and pH would be fun, too.
     
     The model is spatially implicit, which makes sense because COBs are often
     well-mixed and we have no interest in tracking individual bacteria
     
     Volume (V), influent rate (r), and hence residence time (V/r) are modeled
     as influencing multiple resource and community-related rates and dynamics.
     The incoming concentration of resources and density of propagules are also
     accounted for. Influent rate (r) is the constraining rate on all passively
     in-flowing stuff. 

     Unlike some neutral models, this model is not zero-sum, i.e. community size (N)     
     is not fixed and can fluctuate wildly. As in real COBs, the community should
     go extinct if the inflow rate (r) is to high. Likewise, N should grow to
     capacity if r is to low and if the resource concentration is to high. 
     Like any good model, this model should yield realistic behavior.
     
     There is no mutation rate and no speciation, but there could be.
     Unlike any other ecological neutral model, individuals do not have to be
     'born' into a state of reproductive maturity. Individuals are modeled such
     that they can grow to a reproductively viable state. The rate of per capita
     growth is then influenced by resource availability (R) and R is influenced
     by V, r, and the concentration of resources in the inflowing medium. """

""" A function to find the RAD of the community """
def get_rad(CODcom):
    rad = []
    tx_labels = []
    for i in CODcom:
        tx_labels.append(i[0])
    taxa = set(tx_labels)
    for t in taxa:
        ab = tx_labels.count(t)
        rad.append(np.log(ab))
    rad.sort()
    rad.reverse()
    
    return rad

""" Some functions to simulate immigration and death/emigration """                      
def immigration(COBcom,im_rate):
    p = 0.7 # arbitrarily set log-series parameter
    props = np.random.logseries(p,im_rate) # An initial set of propagules; a list of log-series distributed integers.
                                              # Assume the source community is infinite. Because large communities are
                                              # approximately log-series distributed (most things are rare and relatively
                                              # few things are abundant) immigration of propagules will occur in a log-series
                                              # distributed fashion (i.e. not all species have the same chance of contributing 
                                              # propagules. Still neutral in the per capita sense.
    num_A = 0 # number of active individuals
    for p in props: 
        state = choice([1,2]) #   starting assumptions: 1. propagules are as likely to be active (1) as dormant (2)
        if state == 1: num_A += 1
        growth = float(np.random.randint(0,101))#       2. propagules have equal chances of being 0 to 100% reproductively viable  
        i = [p,state,growth] 
        COBcom.append(i) # adding the propagule's taxa label, activity state, and growth state to the community 
    
    return [COBcom,num_A]

    
def death_emigration(COBcom,num_out):
    ct = 0
    num_A = 0
    while ct < num_out:
        random_i = randrange(0,len(COBcom)) # randomly pick an individual
        if COBcom[random_i][1] == 1: num_A += 1 # count the number of active individuals lost
        COBcom.pop(random_i) # die/emigrate
        ct+=1
    return [COBcom,num_A]



""" Things that will remain constant through time, the values of which, must be reasonable
    or the COB will crash, explode, or worse 

    Some values for cow rumen obtainable here: http://microbewiki.kenyon.edu/index.php/Bovine_Rumen """
    
V = 1000.0       # volume of the COB                                                                                                  
r = 100.0         # influent rate (unit volume/unit time)                                                                  
prop_dens = 10.0 # propagule density, (cells or biomass per unit volume of inflowing medium)

res_dens = 0.1  # growth limiting resource concentration of inflowing medium,
                # e.g. (grams cellulose + grams x + grams y) / (grams of medium flowing in) 
                
                # Assume initially that resource concentration of the influent equals
                # the resource concentration of the COB. This makes sense if we're 
                # starting with a community of zero individuals.

dorm_lim = 0.05 # dormancy threshold; dormancy is undertaken if per capita resource availability
                # is below some threshhold (low resources -> low metabolism -> slow growth = go dormant)
                # This could be made to vary among species

im_rate = int(round(prop_dens * r)) # immigration rate (cells/unit time)                                     
res_rate = res_dens * r # resource delivery rate, (unit resource/unit time) 

R = V * res_dens # Amount of resources in the COB before inoculation. If inflow = outflow and the COB is empty,
                 # then the initial resource concentration would equal that of the inflow and outflow. Eventually,
                 # total resources (R) will change as a result of a changing community.

""" The community (COBcom) will be a list of lists: """
COBcom = []
comlist = immigration(COBcom, im_rate) # inoculate the community with propagules using the above function
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
    
time = 100   # length of the experiment

N_COBcom = []   # list to track N over time
A_COBcom = []   # list to track number of active individuals over time
pcr_COBcom = [] # list to tack per capita resources over time
R_COBcom = []   # list to track total resources over time
D_COBcom = []   # list to track dormancy

t = 0
while t <= time: # looping one time unit at a time
    #ct = 0
    #for i in COBcom:
    #    if i[1] == 1: ct += 1
    print 'time',t,' ','size =',N,'active',num_A#,ct
    
    """ inflow of individuals, i.e., immigration """
    comlist = immigration(COBcom, im_rate) # add some propagules to the community
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
    
    burnin = 20
    if t >= burnin: #and t%10 == 0: # allow a burn-in
        N_COBcom.append(np.log(N)) # using natural logs when values can be enormous
        A_COBcom.append(np.log(num_A))
        D_COBcom.append(np.log(N-num_A))
        pcr_COBcom.append(ind_res)
        R_COBcom.append(np.log(R))
    
    t += 1
    
    """ Here, we have completed one time interval of inflow/outflow """
    
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
plt.plot(A_COBcom,'b',label='Active')
plt.plot(D_COBcom,'r',label='Dormant')
plt.ylim(3.0,9.0) 
plt.xlabel("Time",fontsize=12)
plt.ylabel("ln(abundance)",fontsize=12)
plt.text(40,10.0,'Volume = '+str(V)+', inflow rate = '+str(r),fontsize=14,color='k',weight='heavy')
plt.text(-10.0,9.4,'Resource density = '+str(res_dens)+', Propagule density = '+str(prop_dens)+', Dormancy threshold = '+str(dorm_lim),fontsize=14,color='k',weight='heavy')
# Add legend
leg = plt.legend(loc=8,prop={'size':12})
leg.draw_frame(False)

ax = plt.subplot2grid((2,2), (0,1), rowspan=1) # plotting N, R, & per capita resources through time
plt.plot(N_COBcom, '0.5', label='ln(total abundance)')
plt.plot(R_COBcom, 'b', label='ln(total resources)')
plt.plot(pcr_COBcom, 'r', label='per capita resources')

t_range = range(0,(time - burnin),10)
for t in t_range:
    plt.axvline(x=t,color='0.80',ls='--',lw=1) # plot a vertical line at the mode

plt.xlabel("Time",fontsize=12)
plt.ylabel("Value",fontsize=12)
# Add legend
leg = plt.legend(loc=8,prop={'size':12})
leg.draw_frame(False)

ax = plt.subplot2grid((2,2), (1,0), rowspan=1) # plotting activity vs. dormancy
plt.scatter(A_COBcom, D_COBcom, c='0.4', lw=0.5)#, label='active vs. dormant')
plt.xlabel("ln(active abundance)",fontsize=12)
plt.ylabel("ln(dormant abundance)",fontsize=12)
# Add legend
#leg = plt.legend(loc=3,prop={'size':12})
#leg.draw_frame(False)

ax = plt.subplot2grid((2,2), (1,1), rowspan=1) # plotting activity vs. dormancy
RAD = get_rad(COBcom)
rank = range(1,len(RAD)+1)
plt.scatter(rank, RAD, c='0.4', lw=0.5)#, label='RAD')
plt.xlabel("Rank",fontsize=12)
plt.ylabel("ln(abundance)",fontsize=12)
# Add legend
#leg = plt.legend(loc=1,prop={'size':12})
#leg.draw_frame(False)

plt.subplots_adjust(wspace=0.2, hspace=0.2)
plt.savefig('COBcom.png', dpi=400, bbox_inches = 'tight', pad_inches=0.1) 


sys.exit()    
# write the list to a file
OUT = open('/home/kenlocey/COBcom.txt','w+')
for _list in COBcom:
    print>>OUT, _list[0],_list[1],_list[2]    
OUT.close() 
