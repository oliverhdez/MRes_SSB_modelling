from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
import random
  

# To do:
# add transcription and translation
    # mRNA degradation rate in yeast?
# run through several generations
# think of a way of presenting the data



def find_index_from_time(t_obs,time,start_index=0):  
    i=start_index
    while i+1<len(t_obs):  
        if t_obs[i+1]>time:
            break
        i=i+1
    return i 
      
def resample_observations(t_obs_in, s_obs_in, t_obs_out):
    s_obs_out=[] 
    pos=0 
    for time in t_obs_out:  
        i=find_index_from_time(t_obs_in,time, start_index=pos)
        si = s_obs_in[i]  
        s_obs_out.append(si) 
        pos = i
    return s_obs_out


def gen_next_event_time(rate):
    t=random.expovariate(rate)
    return t

def calc_dillution(time_of_phase, initial_volume, final_volume):
    kdil=(math.log(final_volume/initial_volume))/time_of_phase
    return kdil


def gillespie(s0,t_obs_out,params):

    #--0--# Unpack parameters and species variables

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kexpression, percentage_sd, avogadro = params
    
    AUXIN, PIN2, VOLUME  = s0

    print(kexpression)

    #--0--#

    # create arrays for output
    s_obs=[]
    t_obs=[]

    # read in start time and end time
    t_init=t_obs_out[0]
    t_final=t_obs_out[-1]

    t=t_init
    t_obs.append(t)
    s_obs.append(s0)

    while t < t_final:

        types=['influx','efflux','expression','growth']
        rate_influx = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
        if rate_influx < 0:
            rate_influx = 0
        rate_efflux = permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN
        if rate_efflux < 0:
            rate_efflux = 0
        rate_expression = kexpression
        rate_growth = VOLUME * kdil
        rates=[rate_influx,rate_efflux,rate_expression,rate_growth]
        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        next_event=gen_next_event_time(rate_all_events)
        print(next_event)
        t=t+next_event

        AUXIN+= np.random.poisson(rate_influx*next_event)
        AUXIN-= np.random.poisson(rate_efflux*next_event)
        PIN2+= (np.random.poisson(rate_expression*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = np.random.normal((rate_growth*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution



        s = (AUXIN,PIN2,VOLUME)

        t_obs.append(t)
        s_obs.append(s)

    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)

#setting seed to make results reproducible

np.random.seed=1
random.seed=1

AUXIN0 = 0
PIN20 = 0
VOLUME0 = 26*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 0.5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081 # um*s^-1 per umol/L of PIN2
Nu = 5.02
kexpression = 1.2*(10**(-4)) # based on the average expression rate of proteins in yeast
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)


# Mother 0 G1

G1_length = 60*20
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=26*(10**(-15))
G1_final_volume=38*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)
print(kdil)

s0 = (AUXIN0,PIN20,VOLUME0)
params = (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kexpression,percentage_sd,avogadro)

s_obs=gillespie(s0,t_G1,params)

results = s_obs[:]

#n_runs=1
#result_runs=[]

#for i in range(n_runs):
#    print('Simulating {} of {} runs...'.format(i+1,n_runs))
#    s_obs=gillespie(s0,t_G1,params)
#    results = s_obs[:]
#    result_runs.append(results)

#result_runs = np.array(result_runs)

fig=plt.figure(figsize=(12,8))

ax1=fig.add_subplot(3,1,1)
ax2=fig.add_subplot(3,1,2)
ax3=fig.add_subplot(3,1,3)


# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_G1, results[:,0],'r-',label='AUXIN')
ax2.plot(t_G1, results[:,1],'g-',label='PIN2')
ax3.plot(t_G1, results[:,2],'b-',label='VOLUME')


#ax1.plot(t_obs, result_runs[1],'g-',label='run 2')
#ax1.plot(t_obs, result_runs[2],'y-',label='run 3')
#ax1.plot(t_obs, result_runs[3],'r-',label='run 4')
#ax1.plot(t_obs, result_runs[4],'k-',label='run 5')


#ax1.plot(t_obs, RNA_mean,'k',label='Average concentration of mRNA')
ax1.set_xlabel('Time, in seconds')
ax1.set_ylabel('Concentration in arbitrary units')
ax1.legend()

ax2.set_xlabel('Time, in seconds')
ax2.set_ylabel('Concentration in arbitrary units')
ax2.legend()

ax3.set_xlabel('Time, in seconds')
ax3.set_ylabel('Concentration in arbitrary units')
ax3.legend()

plt.show()
