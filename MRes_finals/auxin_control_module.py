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


def gen_next_event_time(rate):
    t=random.expovariate(rate)
    return t

def calc_dillution(time_of_phase, initial_volume, final_volume):
    kdil=(math.log(final_volume/initial_volume))/time_of_phase
    return kdil


def gillespie(s0,t_obs_out,params):

    #--0--# Unpack parameters and species variables

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_PIN2rnaexpression, k_ACE2rnaexpression, k_PIN2expression,percentage_sd,avogadro,k_rnadecay = params
    
    AUXIN, mRNA_PIN2, PIN2, mRNA_ACE2, ACE2, VOLUME  = s0

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

        types=['influx','efflux','expression','RNA_decay','synthesis',"ace2_expression","ace2_decay","ace2_synthesis",'growth']
        rate_influx = (permeability_IAAH/pm_thickness)*(fIAAH_w*conc_out - fIAAH_c*AUXIN)#*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
        rate_efflux = permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN
        rate_expression = (k_PIN2rnaexpression)*ACE2
        rate_decay = k_rnadecay*mRNA_PIN2
        rate_synthesis = k_PIN2expression*mRNA_PIN2
        rate_ace2_expression =  k_ACE2rnaexpression
        rate_ace2_decay = k_rnadecay*mRNA_ACE2
        rate_ace2_synthesis = k_PIN2expression*mRNA_ACE2
        rate_growth = VOLUME * kdil
        rates=[rate_influx,rate_efflux,rate_expression,rate_decay,rate_synthesis,rate_ace2_expression,rate_ace2_decay,rate_ace2_synthesis,rate_growth]
        position=0
        for i in rates:
            if i < 0:
                rates[position]=0
            position += 1

        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        next_event=gen_next_event_time(rate_all_events)
        next_event=1
        t=t+next_event

        AUXIN+= (np.random.poisson(rates[0]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        AUXIN-= (np.random.poisson(rates[1]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2+= (np.random.poisson(rates[2]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2-= (np.random.poisson(rates[3]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        PIN2+= (np.random.poisson(rates[4]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_ACE2+= (np.random.poisson(rates[5]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_ACE2-= (np.random.poisson(rates[6]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        ACE2+= (np.random.poisson(rates[7]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = np.random.normal((rates[-1]*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution



        s = (AUXIN, mRNA_PIN2, PIN2, mRNA_ACE2, ACE2, VOLUME)

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
mRNA_ACE20 = 9.09364077747e-10*0.8
ACE20 = .389
VOLUME0 = 30*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081*(10**(-1)) # um*s^-1 per umol/L of PIN2
Nu = 5.02
k_PIN2rnaexpression = 2.18862203686e-12/.2 # based on average expression rate of RNAs in yeast
k_ACE2rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_PIN2expression = 52769.9671394 # um per liter per mRNA per second 
percentage_sd = 0.0636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536


# Mother 0 G1

G1_length = 91*60
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=30*(10**(-15))
G1_final_volume=60*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)


s0 = (AUXIN0,mRNA_PIN20,PIN20,mRNA_ACE20,ACE20,VOLUME0)
params_daughter = (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_PIN2rnaexpression,k_ACE2rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)
params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_PIN2rnaexpression,k_ACE2rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)

#s_obs=gillespie(s0,t_G1,params)

#results = s_obs[:]


gen0=gillespie(s0,t_G1,params_daughter)

ori_results=[]
daughter_results=[]

ori_results.append(gen0)
s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.3)/(avogadro*ori_results[-1][-1][-1]/2)),ori_results[-1][-1][3],((ori_results[-1][-1][4]*ori_results[-1][-1][-1]*avogadro*0)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))
s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.7)/(avogadro*ori_results[-1][-1][-1]/2)),ori_results[-1][-1][3],((ori_results[-1][-1][4]*ori_results[-1][-1][-1]*avogadro*2)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))


n_generations=30

for i in range(n_generations):
    print("Simulating {} of {} generations...".format(i+1,n_generations))
    d2=gillespie(s0_daughter,t_G1,params_daughter)
    daughter_results.append(d2)
    d1=gillespie(s0_ori,t_G1,params)
    ori_results.append(d1)
    s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.3)/(avogadro*ori_results[-1][-1][-1]/2)),ori_results[-1][-1][3],((ori_results[-1][-1][4]*ori_results[-1][-1][-1]*avogadro*0)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))
    s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.7)/(avogadro*ori_results[-1][-1][-1]/2)),ori_results[-1][-1][3],((ori_results[-1][-1][4]*ori_results[-1][-1][-1]*avogadro*1)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))



ori_results=np.array(ori_results)
daughter_results=np.array(daughter_results)
#print(daughter_results[0])
#print(daughter_results[0][...,0])

fig=plt.figure(figsize=(12,8))
fig.suptitle("Figure 4: Resulting daughter cells overlayed")
fig2=plt.figure(figsize=(12,8))
fig2.suptitle("Figure 3: Auxin dynamics in a mother cell over 30 cell divisions, using the endocytic recycling signal peptide")


ax1=fig.add_subplot(3,2,1)
ax2=fig.add_subplot(3,2,2)
ax3=fig.add_subplot(3,2,3)
ax4=fig.add_subplot(3,2,4)
ax5=fig.add_subplot(3,2,5)
ax6=fig.add_subplot(3,2,6)

t_G1/=60


# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_G1, daughter_results[:][...,0].T,'g:')
ax2.plot(t_G1, daughter_results[:][...,1].T,'b:')
ax3.plot(t_G1, daughter_results[:][...,2].T,'y:')
ax4.plot(t_G1, daughter_results[:][...,3].T,'m:')
ax5.plot(t_G1, daughter_results[:][...,4].T,'c:')
ax6.plot(t_G1, daughter_results[:][...,-1].T,'r:')


auxin_mean = np.average(daughter_results[:][...,0], axis=0,)
mrna_mean = np.average(daughter_results[:][...,1], axis=0)
pin2_mean = np.average(daughter_results[:][...,2], axis=0)
ace2_mrna_mean = np.average(daughter_results[:][...,3], axis=0)
ace2_mean = np.average(daughter_results[:][...,4], axis=0)
volume_mean = np.average(daughter_results[:][...,-1], axis=0)

ax1.plot(t_G1, auxin_mean,'k',label="Auxin")
ax2.plot(t_G1, mrna_mean,'k',label="PIN2_mRNA")
ax3.plot(t_G1, pin2_mean,'k',label="PIN2")
ax4.plot(t_G1, ace2_mrna_mean,'k',label="ACE2_mRNA")
ax5.plot(t_G1, ace2_mean,'k',label="ACE2")
ax6.plot(t_G1, volume_mean,'k',label="Volume")

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
ax4.set_ylabel("Concentration in umol per litre")
ax4.legend()

ax5.set_xlabel("Time, in minutes")
ax5.set_ylabel("Concentration in umol per litre")
ax5.legend()

ax6.set_xlabel("Time, in minutes")
ax6.set_ylabel("Volume in litres")
ax6.legend()

ax7=fig2.add_subplot(3,2,1)
ax8=fig2.add_subplot(3,2,2)
ax9=fig2.add_subplot(3,2,3)
ax10=fig2.add_subplot(3,2,4)
ax11=fig2.add_subplot(3,2,5)
ax12=fig2.add_subplot(3,2,6)

t_all = np.linspace(0,G1_length*n_generations+G1_length,G1_length*n_generations+G1_length)
coefficients=[]
degrees = [6,5,5,1,5,1]
for i in range(len(degrees)):
    c = np.polyfit(t_all,ori_results[:][...,i].flatten().tolist(),degrees[i])
    coefficients.append(c)
auxin_fit = []
pin2_mrna_fit = []
pin2_fit = []
ace2_mrna_fit = []
ace2_fit = []
volume_fit = []
y=[auxin_fit,pin2_mrna_fit,pin2_fit,ace2_mrna_fit,ace2_fit,volume_fit]
for i in range(6):
    poly=np.poly1d(coefficients[i])
    y[i]=poly(t_all)

t_all /= 91*60

ax7.plot(t_all, ori_results[:][...,0].flatten().tolist(),'g:',label="Auxin")
ax8.plot(t_all, ori_results[:][...,1].flatten().tolist(),'b:',label="PIN2_mRNA")
ax9.plot(t_all, ori_results[:][...,2].flatten().tolist(),'y:',label="PIN2")
ax10.plot(t_all, ori_results[:][...,3].flatten().tolist(),'m:',label="ACE2_mRNA")
ax11.plot(t_all, ori_results[:][...,4].flatten().tolist(),'c:',label="ACE2")
ax12.plot(t_all, ori_results[:][...,-1].flatten().tolist(),'r:',label="Volume")

ax7.plot(t_all, y[0],'k')
ax8.plot(t_all, y[1],'k')
ax9.plot(t_all, y[2],'k')
ax10.plot(t_all, y[3],'k')
ax11.plot(t_all, y[4],'k')
ax12.plot(t_all, y[-1],'k')



ax7.set_xlabel("Number of divisions")
ax7.set_ylabel("Concentration in umol per litre")
ax7.legend()

ax8.set_xlabel("Number of divisions")
ax8.set_ylabel("Concentration in umol per litre")
ax8.legend()

ax9.set_xlabel("Number of divisions")
ax9.set_ylabel("Concentration in umol per litre")
ax9.legend()

ax10.set_xlabel("Number of divisions")
ax10.set_ylabel("Concentration in umol per litre")
ax10.legend()

ax11.set_xlabel("Number of divisions")
ax11.set_ylabel("Concentration in umol per litre")
ax11.legend()

ax12.set_xlabel("Number of divisions")
ax12.set_ylabel("Volume in Litres")
ax12.legend()





plt.show()

