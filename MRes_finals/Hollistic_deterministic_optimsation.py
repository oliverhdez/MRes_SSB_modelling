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

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_PIN2rnaexpression,k_ACE2rnaexpression,k_PIN2expression,k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad,percentage_sd,avogadro = params
    
    AUXIN, mRNA_PIN2, PIN2, mRNA_ACE2, ACE2, mRNA_TIR1, TIR1, mRNA_RFP, RFP, mRNA_BFP, BFP, VOLUME  = s0
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

        types=['influx','efflux','expression','RNA_decay','synthesis',"ace2_expression","ace2_decay","ace2_synthesis",'TIR1_transcription','TIR1_rna_decay','TIR1_translation','RFP_transcription','RFP_rna_decay','RFP_translation',"leaky_degradation",'RFP_degradation','BFP_transcription','BFP_rna_decay','BFP_translation','growth']
        rate_influx = (permeability_IAAH/pm_thickness)*(fIAAH_w*conc_out - fIAAH_c*AUXIN)#*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
        rate_efflux = permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN
        rate_PIN2mRNA_expression = (k_PIN2rnaexpression)*ACE2
        rate_PIN2mRNA_decay = k_rnadecay*mRNA_PIN2
        rate_PIN2_synthesis = k_PIN2expression*mRNA_PIN2
        rate_ace2_expression =  k_ACE2rnaexpression
        rate_ace2_decay = k_rnadecay*mRNA_ACE2
        rate_ace2_synthesis = k_PIN2expression*mRNA_ACE2
        rate_TIR1_transcription = k_TIR1rnaexpression
        rate_TIR1_rna_decay = k_rnadecay*mRNA_TIR1
        rate_TIR1_translation = k_translation*mRNA_TIR1
        rate_RFP_transcription = k_RFPrnaexpression
        rate_RFP_rna_decay =k_rnadecay*mRNA_RFP
        rate_RFP_translation = k_translation*mRNA_RFP
        rate_leaky_degradation = k_leakydegrad*RFP*TIR1
        rate_RFP_degradation = k_degrad*RFP*AUXIN*TIR1
        rate_BFP_transcription = k_BFPrnaexpression
        rate_BFP_rna_decay =k_rnadecay*mRNA_BFP
        rate_BFP_translation =k_translation*mRNA_BFP
        rate_growth = VOLUME * kdil
        rates=[rate_influx,rate_efflux,rate_PIN2mRNA_expression,rate_PIN2mRNA_decay,rate_PIN2_synthesis,rate_ace2_expression,rate_ace2_decay,rate_ace2_synthesis,rate_TIR1_transcription,rate_TIR1_rna_decay,rate_TIR1_translation,rate_RFP_transcription,rate_RFP_rna_decay,rate_RFP_translation,rate_leaky_degradation,rate_RFP_degradation,rate_BFP_transcription,rate_BFP_rna_decay,rate_BFP_translation,rate_growth]
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

        AUXIN+= (rates[0]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        AUXIN-= (rates[1]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_PIN2+= (rates[2]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_PIN2-= (rates[3]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        PIN2+= (rates[4]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_ACE2+= (rates[5]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_ACE2-= (rates[6]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        ACE2+= (rates[7]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_TIR1+= (rates[8]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_TIR1-= (rates[9]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        TIR1+= (rates[10]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_RFP+= (rates[11]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_RFP-= (rates[12]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        RFP+= (rates[13]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        RFP-= (rates[14]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        RFP-= (rates[15]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_BFP+= (rates[16]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_BFP-= (rates[17]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        BFP+= (rates[18]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)

        event_dilution = (rates[-1]*next_event)
        mRNA_PIN2-= mRNA_PIN2-((mRNA_PIN2*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        mRNA_ACE2 -= mRNA_ACE2-((mRNA_ACE2*VOLUME)/(VOLUME+event_dilution))
        ACE2 -= ACE2-((ACE2*VOLUME)/(VOLUME+event_dilution))
        mRNA_TIR1-= mRNA_TIR1-((mRNA_TIR1*VOLUME)/(VOLUME+event_dilution))
        TIR1-= TIR1-((TIR1*VOLUME)/(VOLUME+event_dilution))
        mRNA_RFP-= mRNA_RFP-((mRNA_RFP*VOLUME)/(VOLUME+event_dilution))
        RFP-= RFP-((RFP*VOLUME)/(VOLUME+event_dilution))
        mRNA_BFP-= mRNA_BFP-((mRNA_BFP*VOLUME)/(VOLUME+event_dilution))
        BFP-= BFP-((BFP*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution

        s = (AUXIN, mRNA_PIN2, PIN2, mRNA_ACE2, ACE2, mRNA_TIR1, TIR1, mRNA_RFP, RFP, mRNA_BFP, BFP, VOLUME)

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
ACE20 = 0
mRNA_TIR10 = 9.09364077747e-10
TIR10 = 0.389
mRNA_RFP0 = 1.8*9.09364077747e-10
RFP0 = 0.47
mRNA_BFP0 = 1.8*9.09364077747e-10
BFP0 = 0.378*1.8
VOLUME0 = 30*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 750  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081*(10**(-1)) # um*s^-1 per umol/L of PIN2
Nu = 5.02
k_PIN2rnaexpression = 2.18862203686e-12/.2 # based on average expression rate of RNAs in yeast
k_ACE2rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_PIN2expression = 52769.9671394 # um per liter per mRNA per second 
k_TIR1rnaexpression = 2.43e-12 # based on average expression rate of RNAs in yeast
k_RFPrnaexpression = 4e-12 # based on average expression rate of RNAs in yeast
k_BFPrnaexpression = 4e-12 # based on average expression rate of RNAs in yeast
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536
k_translation = 52769.9671394 # um per liter per mRNA per second 
k_leakydegrad = 0.000150354734602 # umolar per second
k_degrad = 4.1718216004e-09 # umolar per second
percentage_sd = 0.0636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536


# Mother 0 G1

G1_length = 91*60
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=30*(10**(-15))
G1_final_volume=60*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)


s0 = (AUXIN0,mRNA_PIN20,PIN20,mRNA_ACE20,ACE20, mRNA_TIR10, TIR10, mRNA_RFP0, RFP0, mRNA_BFP0, BFP0, VOLUME0)
params_daughter = (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_PIN2rnaexpression,k_ACE2rnaexpression,k_PIN2expression,k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad,percentage_sd,avogadro)
params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_PIN2rnaexpression,k_ACE2rnaexpression,k_PIN2expression,k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad,percentage_sd,avogadro)

#s_obs=gillespie(s0,t_G1,params)

#results = s_obs[:]


gen0=gillespie(s0,t_G1,params_daughter)

conc_out=0
params_daughter = (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_PIN2rnaexpression,k_ACE2rnaexpression,k_PIN2expression,k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad,percentage_sd,avogadro)
params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_PIN2rnaexpression,k_ACE2rnaexpression,k_PIN2expression,k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad,percentage_sd,avogadro)

# Make the phase plot of auxin induction vs mother daughter assymmetry
mother_rfp=[]
daughter_rfp=[]
conc_outs=[]
step=0

while conc_out <= 100:
    ori_results=[]
    daughter_results=[]
    print("step:",step)
    step+=1
    params_daughter = (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_PIN2rnaexpression,k_ACE2rnaexpression,k_PIN2expression,k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad,percentage_sd,avogadro)
    params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_PIN2rnaexpression,k_ACE2rnaexpression,k_PIN2expression,k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad,percentage_sd,avogadro)


    ori_results.append(gen0)
    AUXIN0=ori_results[-1][-1][0]
    s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.3)/(avogadro*ori_results[-1][-1][-1]/2)),ori_results[-1][-1][3],((ori_results[-1][-1][4]*ori_results[-1][-1][-1]*avogadro*0)/(avogadro*ori_results[-1][-1][-1]/2)),ori_results[-1][-1][5],ori_results[-1][-1][6],ori_results[-1][-1][7],ori_results[-1][-1][8],ori_results[-1][-1][9],ori_results[-1][-1][10],(ori_results[-1][-1][-1]/2))
    s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.7)/(avogadro*ori_results[-1][-1][-1]/2)),ori_results[-1][-1][3],((ori_results[-1][-1][4]*ori_results[-1][-1][-1]*avogadro*1)/(avogadro*ori_results[-1][-1][-1]/2)),ori_results[-1][-1][5],ori_results[-1][-1][6],ori_results[-1][-1][7],ori_results[-1][-1][8],ori_results[-1][-1][9],ori_results[-1][-1][10],(ori_results[-1][-1][-1]/2))

    ori_results=[]


    d2=gillespie(s0_daughter,t_G1,params_daughter)
    daughter_results.append(d2)
    d1=gillespie(s0_ori,t_G1,params)
    ori_results.append(d1)

    ori_results=np.array(ori_results)
    daughter_results=np.array(daughter_results)

    mother_rfp.append(ori_results[-1][-1][8])
    print(ori_results[-1][-1][8])
    daughter_rfp.append(daughter_results[-1][-1][8])
    print(daughter_results[-1][-1][8])
    conc_outs.append(conc_out)

    conc_out+=1

rfp_differences = []
zip_rfps=zip(mother_rfp,daughter_rfp)
for i, j in zip_rfps:
    rfp_differences.append(j-i)
optimum_expression=conc_outs[rfp_differences.index(max(rfp_differences))]
print(max(rfp_differences))
print(optimum_expression)

fig=plt.figure(figsize=(12,8))

ax1=fig.add_subplot(1,1,1)

ax1.plot(conc_outs, rfp_differences,'r')

ax1.set_xlabel("Expression level of TIR1 mRNA, in umol per litre")
ax1.set_ylabel("Differences between mother and daughter concentration of mRUBY2, in umol per litre")
ax1.legend()

plt.show()


