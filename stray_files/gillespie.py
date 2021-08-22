from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
import random
  

# To do:
# what about ace2 dynamics
# run through several generations
# think of a way of presenting the data. What about an overlay of all daughter cells on one graph and mother cells on another graph over 10/20 generations


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


def gen_next_event_time(rate):
    t=random.expovariate(rate)
    return t

def calc_dillution(time_of_phase, initial_volume, final_volume):
    kdil=(math.log(final_volume/initial_volume))/time_of_phase
    return kdil


def gillespie(s0,t_obs_out,params):

    #--0--# Unpack parameters and species variables

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay = params
    
    AUXIN, mRNA_PIN2, PIN2, VOLUME  = s0

    AUXIN *= VOLUME*avogadro
    mRNA_PIN2 *= VOLUME*avogadro
    PIN2 *= VOLUME*avogadro

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

    percent=10

    while t < t_final:

        types=['influx','efflux','expression','RNA_decay','synthesis']
        rate_influx = (permeability_IAAH/pm_thickness)*(fIAAH_w*conc_out*VOLUME*avogadro - fIAAH_c*AUXIN*VOLUME*avogadro)
        rate_efflux = permeability_IAA*PIN2*VOLUME*avogadro*Nu*(1-fIAAH_c)*AUXIN*VOLUME*avogadro
        rate_expression = k_rnaexpression*VOLUME*avogadro*VOLUME*avogadro
        rate_decay = k_rnadecay*mRNA_PIN2*VOLUME*avogadro
        rate_synthesis = k_PIN2expression*mRNA_PIN2*VOLUME*avogadro
        rates=[rate_influx,rate_efflux,rate_expression,rate_decay,rate_synthesis]
        rate_growth = VOLUME * kdil
        position=0
        for i in rates:
            if i < 0:
                rates[position]=0
            position += 1

        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        next_event=gen_next_event_time(rate_all_events)
        pdf=[]
        for event_rate in rates:
            p_event = event_rate/sum(rates)
            pdf.append(p_event)

        rand_i =  random_choice_from_pdf(pdf)
        event_type=types[rand_i]
        t+=next_event*100000

        if event_type=="influx":
            AUXIN+=100000
        elif event_type=="efflux":
            AUXIN-=100000
        elif event_type=="expression":
            PIN2_mRNA+=100000
        elif event_type=="rna_decay":
            PIN2_mRNA-=100000
        elif event_type=="synthesis":
            PIN2+=100000
        else:
            print("error unknown event type!!")

        #AUXIN+= (np.random.poisson(rates[0]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        #AUXIN-= (np.random.poisson(rates[1]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        #mRNA_PIN2+= (np.random.poisson(rates[2]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        #mRNA_PIN2-= (np.random.poisson(rates[3]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        #PIN2+= (np.random.poisson(rates[4]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = np.random.normal((rate_growth*next_event*100000),(rate_growth*next_event*100000*percentage_sd))
        VOLUME+= event_dilution
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution

        if t/t_final*100 >= percent:
            print("{}% of generation".format(percent))
            percent+=10

        s = (AUXIN/(VOLUME*avogadro), mRNA_PIN2/(VOLUME*avogadro), PIN2/(VOLUME*avogadro), VOLUME)

        t_obs.append(t)
        s_obs.append(s)

    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)

#setting seed to make results reproducible

np.random.seed=1
random.seed=1

AUXIN0 = 0
mRNA_PIN20 = 9.09364077747e-10
PIN20 = 0.15
VOLUME0 = 30*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081*(10**(-1)) # um*s^-1 per umol/L of PIN2
Nu = 5.02
k_rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_PIN2expression = 52769.9671394 # um per liter per mRNA per second 
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**17) # as units are in uM
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536


# Mother 0 G1

G1_length = 91*60
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=30*(10**(-15))
G1_final_volume=60*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)


s0 = (AUXIN0,mRNA_PIN20,PIN20,VOLUME0)
params_daughter = (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)
params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, 0,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)

#s_obs=gillespie(s0,t_G1,params)

#results = s_obs[:]

print("stuck")
gen0=gillespie(s0,t_G1,params_daughter)
print("not stuck")

ori_results=[]
daughter_results=[]

ori_results.append(gen0)
s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.4)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))
s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.6)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))


n_generations=2

for i in range(n_generations):
    print("Simulating {} of {} generations...".format(i+1,n_generations))
    d2=gillespie(s0_daughter,t_G1,params_daughter)
    daughter_results.append(d2)
    d1=gillespie(s0_ori,t_G1,params)
    ori_results.append(d1)
    s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.4)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))
    s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.6)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))


ori_results=np.array(ori_results)
daughter_results=np.array(daughter_results)
#print(daughter_results[0])
#print(daughter_results[0][...,0])

fig=plt.figure(figsize=(12,8))
fig2=plt.figure(figsize=(12,8))

ax1=fig.add_subplot(4,1,1)
ax2=fig.add_subplot(4,1,2)
ax3=fig.add_subplot(4,1,3)
ax4=fig.add_subplot(4,1,4)

# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_G1, daughter_results[:][...,0].T,'g:')
ax2.plot(t_G1, daughter_results[:][...,1].T,'b:')
ax3.plot(t_G1, daughter_results[:][...,2].T,'y:')
ax4.plot(t_G1, daughter_results[:][...,3].T,'r:')

auxin_mean = np.average(daughter_results[:][...,0], axis=0,)
mrna_mean = np.average(daughter_results[:][...,1], axis=0)
pin2_mean = np.average(daughter_results[:][...,2], axis=0)
volume_mean = np.average(daughter_results[:][...,3], axis=0)

ax1.plot(t_G1, auxin_mean,'k',label="Auxin")
ax2.plot(t_G1, mrna_mean,'k',label="PIN2 mRNA")
ax3.plot(t_G1, pin2_mean,'k',label="PIN2")
ax4.plot(t_G1, volume_mean,'k',label="Volume")

ax1.set_xlabel("Time, in minutes")
ax1.set_ylabel("Concentration in umol per litre")
ax1.legend()

ax2.set_xlabel("Time, in minutes")
ax2.set_ylabel("Concentration in umol per litre")
ax2.legend()

ax3.set_xlabel("Time, in minutes")
ax3.set_ylabel("Concentration in umol per litre")
ax3.legend()

ax4.set_xlabel("Time, in minutes")
ax4.set_ylabel("Volume in litres")
ax4.legend()



ax1=fig2.add_subplot(4,1,1)
ax2=fig2.add_subplot(4,1,2)
ax3=fig2.add_subplot(4,1,3)
ax4=fig2.add_subplot(4,1,4)

t_all = np.linspace(0,G1_length*n_generations+G1_length,G1_length*n_generations+G1_length)

ax1.plot(t_all, ori_results[:][...,0].flatten().tolist(),'g',label="Auxin")
ax2.plot(t_all, ori_results[:][...,1].flatten().tolist(),'b',label="PIN2_mRNA")
ax3.plot(t_all, ori_results[:][...,2].flatten().tolist(),'y',label="PIN2")
ax4.plot(t_all, ori_results[:][...,3].flatten().tolist(),'r',label="Volume")

ax1.set_xlabel("Time, in minutes")
ax1.set_ylabel("Concentration in umol per litre")
ax1.legend()

ax2.set_xlabel("Time, in minutes")
ax2.set_ylabel("Concentration in umol per litre")
ax2.legend()

ax3.set_xlabel("Time, in minutes")
ax3.set_ylabel("Concentration in umol per litre")
ax3.legend()

ax4.set_xlabel("Time, in minutes")
ax4.set_ylabel("Volume in Litres")
ax4.legend()





plt.show()

