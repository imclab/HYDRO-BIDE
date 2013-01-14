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

    
def death_emigration(COBcom,num_out):
    ct = 0
    num_A = 0
    while ct < num_out:
        random_i = randrange(0,len(COBcom)) # randomly pick an individual
        if COBcom[random_i][1] == 1: num_A += 1 # count the number of active individuals lost
        COBcom.pop(random_i) # die/emigrate
        ct+=1
    return [COBcom,num_A]


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
