from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
import random
  
# code for stochastic runs

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


def random_choice_from_pdf(pdf):
    cdf=[]
    cumulative_p=0
    for p in pdf:
        cumulative_p+=p
        cdf.append(cumulative_p)
    rand=random.random()

    for i in range(len(cdf)):
        if rand<cdf[i]:
            return i
    # last cdf should be 1.0 so the following should never happen!
    print("Error generating choice, check PDF")
    return None

def gillespie(s0,t_obs_out,params):

    #--0--# Unpack parameters and species variables

    k0,k1 = params
    
    RNA = s0

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

        types=["birth","death"]      
        rate_birth = k0
        rate_death= k1*RNA
        rates=[rate_birth,rate_death]
        rate_all_events=sum(rates)
        if rate_all_events==0:
            break
        next_event=gen_next_event_time(rate_all_events)
        pdf=[]
        for event_rate in rates:
            p_event = event_rate/sum(rates)
            pdf.append(p_event)

        rand_i =  random_choice_from_pdf(pdf)
        event_type=types[rand_i]
        t=t+next_event

        if event_type=="birth":
            RNA+=1
        elif event_type=="death":
            RNA-=1
        else:
            print("error unknown event type!!")

        s = (RNA)

        t_obs.append(t)
        s_obs.append(s)

    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)

#setting seed to make results reproducible

np.random.seed=1
random.seed=1

RNA0 = 0

k0 = 0.2
k1 = 0.01

s0 = (RNA0)

params = (k0,k1)

t_max = 1000
t_obs = np.linspace(0,t_max,t_max+1)

s_obs=gillespie(s0,t_obs,params)

RNA_obs = s_obs[:]

n_runs=5
RNA_runs=[]

for i in range(n_runs):
    print("Simulating {} of {} runs...".format(i+1,n_runs))
    s_obs=gillespie(s0,t_obs,params)
    RNA_obs = s_obs[:]
    RNA_runs.append(RNA_obs)
RNA_runs = np.array(RNA_runs)

fig=plt.figure(figsize=(12,8))

ax1=fig.add_subplot(1,1,1)


# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_obs, RNA_runs[0],'b-',label="run 1")
ax1.plot(t_obs, RNA_runs[1],'g-',label="run 2")
ax1.plot(t_obs, RNA_runs[2],'y-',label="run 3")
ax1.plot(t_obs, RNA_runs[3],'r-',label="run 4")
ax1.plot(t_obs, RNA_runs[4],'k-',label="run 5")

np.savetxt("RNA_runs.txt",RNA_runs)

#ax1.plot(t_obs, RNA_mean,'k',label="Average concentration of mRNA")
ax1.set_xlabel("Time, in seconds")
ax1.set_ylabel("Concentration in arbitrary units")
ax1.legend()

plt.show()
