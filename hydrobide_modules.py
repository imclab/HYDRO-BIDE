#!/usr/local/bin/python                                                                   
                                                                       
import sys                                          
import os                                                                  
import  matplotlib.pyplot as plt   
from pylab import *                                        
import numpy as np                                     
from scipy import stats                                      
import random                                            
from random import choice                               
import re                                               
from decimal import *                               
import math                                 
from random import randrange                
             
             
"""  UNDER HEAVY CONSTRUCTION. FUTURE SITE OF MANY COOL FUNCTIONS

     This file contains functions for plotting, statistical analysis,
     examinations of community structure, and other fun stuff """

  
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

""" A function to get rates of immigration and resource delivery """
def get_in_rates(r,prop_dens,res_dens):
    
    im_rate = int(round(prop_dens * r)) # immigration rate (cells/unit time)                                     
    res_rate = res_dens * r # resource delivery rate, (unit resource/unit time) 

    return [im_rate,res_rate]

""" Some functions to simulate immigration and death/emigration """                      
def immigration(COBcom,im_rate,bp,lgp):
    p = 0.70 # arbitrarily set log-series parameter
    props = np.random.logseries(p,im_rate) # An initial set of propagules; a list of log-series distributed integers.
                                              # Assume the source community is infinite. Because large communities are
                                              # approximately log-series distributed (most things are rare and relatively
                                              # few things are abundant) immigration of propagules will occur in a log-series
                                              # distributed fashion (i.e. not all species have the same chance of contributing 
                                              # propagules. Still neutral in the per capita sense.
    num_A = 0 # number of active individuals
    
    for prop in props: 
        state = np.random.binomial(1,bp,1) #   
        if state == 1: num_A += 1
        else: state == 2
        growth = float(np.random.randint(0,101))#       2. propagules have equal chances of being 0 to 100% reproductively viable  
        i = [prop,state,growth] 
        COBcom.append(i) # adding the propagule's taxa label, activity state, and growth state to the community 
    
    return [COBcom,num_A]

    
def death_emigration(COBcom, N, V, r):
    num_out = int(round(N/V * r)) # no. individuals lost per unit time
    
    if num_out <= 0: 
        num_a = 0
        return [COBcom,num_a]
    
    else:
        ct = 0
        num_a = 0
        while ct < num_out:
            random_i = randrange(0,len(COBcom)) # randomly pick an individual
            if COBcom[random_i][1] == 1: num_a += 1 # count the number of active individuals lost
            COBcom.pop(random_i) # die/emigrate
            ct+=1
    return [COBcom,num_a]


""" A function to recalculate parameter values for the neutral hydrobide model """
def params_neutral(V, r, a, num_A, num_a, R, ind_res, res_rate, in_or_out):

    if in_or_out == 'in':
        num_A += num_a 
        R += res_rate        
        ind_res = R/num_A  
        ind_grow = a*ind_res 
    
    elif in_or_out == 'out':
        num_A -= num_a
        R -= (R/V) * r 
        ind_res = R/num_A  
        ind_grow = a*ind_res 
    
    return [num_A, R, ind_res, ind_grow]


""" A function to add info about dormancy, activity, community size, per capita resources, and total resources to lists  """
def add_to_lists_neutral(COBcom,N,num_A,ind_res,R,N_COBcom,A_COBcom,D_COBcom,pcr_COBcom,R_COBcom,RAD_Ahigh,RAD_Alow,RAD_Amedium):

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
    
    return [N_COBcom, A_COBcom, D_COBcom, pcr_COBcom, R_COBcom, RAD_Ahigh, RAD_Alow, RAD_Amedium]

""" A function that loops through a list containing a neutral-style community """
def loop_thru_neutral_comm(COBcom, ind_res, dorm_lim, num_A, N, ind_grow, R, a):
    
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

    return(COBcom, ind_res, num_A, N, ind_grow, R)



    
    
""" Some plotting function and figures """

def fig1(N_COBcom,A_COBcom,D_COBcom,R_COBcom,pcr_COBcom,time,burnin,sets_of_RADS,V,r,res_dens,prop_dens,dorm_lim,nRADs):

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
    
    return
    

def comm_to_file(COBcom,V,r,res_dens,prop_dens,dorm_lim):
    # write the list to a file
    OUT = open('/results/COBcom V='+str(V)+' r'+str(r)+' resdens='+str(res_dens)+' propdens='+str(prop_dens)+' dormlim='+str(dorm_lim),'w+') 
    for _list in COBcom:
        print>>OUT, _list[0],_list[1],_list[2]    
    OUT.close() 
    return
