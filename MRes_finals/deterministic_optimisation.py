# Use this to find the theoretically perfect level of TIR1 expressions

from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
import random
  

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


def assymetrical_cell_division(s0,t_obs_out,params):

    #--0--# Unpack parameters and species variables

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad, percentage_sd,avogadro, kdil, induction_level = params
    
    AUXIN, mRNA_TIR1, TIR1, mRNA_RFP, RFP, mRNA_BFP, BFP, VOLUME  = s0

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
    induction= "false"

    while t < t_final:

        types=['influx','TIR1_transcription','TIR1_rna_decay','TIR1_translation','RFP_transcription','RFP_rna_decay','RFP_translation',"leaky_degradation",'RFP_degradation','BFP_transcription','BFP_rna_decay','BFP_translation','growth']
        rate_influx = (permeability_IAAH/pm_thickness)*((fIAAH_w*conc_out - fIAAH_c*AUXIN))#*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
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
        rates=[rate_influx,rate_TIR1_transcription,rate_TIR1_rna_decay,rate_TIR1_translation,rate_RFP_transcription,rate_RFP_rna_decay,rate_RFP_translation,rate_leaky_degradation,rate_RFP_degradation,rate_BFP_transcription,rate_BFP_rna_decay,rate_BFP_translation,rate_growth]
        position=0
        for i in rates:
            if i < 0:
                rates[position]=0
            position += 1

        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        #next_event=gen_next_event_time(rate_all_events)
        next_event=1

        t+=1



        AUXIN+= (rates[0]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_TIR1+=(rates[1]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_TIR1-= (rates[2]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        TIR1+= (rates[3]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_RFP+= (rates[4]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_RFP-= (rates[5]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        RFP+= (rates[6]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        RFP-= (rates[7]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        RFP-= (rates[8]*VOLUME*avogadro*next_event)/(VOLUME*avogadro)
        mRNA_BFP+= ((rates[9]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_BFP-= ((rates[10]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        BFP+= ((rates[11]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = (rates[12]*next_event)
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        mRNA_TIR1-= mRNA_TIR1-((mRNA_TIR1*VOLUME)/(VOLUME+event_dilution))
        TIR1-= TIR1-((TIR1*VOLUME)/(VOLUME+event_dilution))
        mRNA_RFP-= mRNA_RFP-((mRNA_RFP*VOLUME)/(VOLUME+event_dilution))
        RFP-= RFP-((RFP*VOLUME)/(VOLUME+event_dilution))
        mRNA_BFP-= mRNA_BFP-((mRNA_BFP*VOLUME)/(VOLUME+event_dilution))
        BFP-= BFP-((BFP*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution

        if t >= 60*20 and induction == "false":
            rfp_start=RFP
            conc_out=induction_level
            induction="true"

        if t == 60*65:
            rfp_end=RFP

        s = (AUXIN, mRNA_TIR1, TIR1, mRNA_RFP, RFP, mRNA_BFP, BFP, VOLUME)

        t_obs.append(t)
        s_obs.append(s)

    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out), rfp_start, rfp_end

#setting seed to make results reproducible

np.random.seed=1
random.seed=1

AUXIN0 = 0
mRNA_TIR10 = 0
TIR10 = 0
mRNA_RFP0 = 1.8*9.09364077747e-10
RFP0 = 0.378*1.8
mRNA_BFP0 = 1.8*9.09364077747e-10
BFP0 = 0.378*1.8
VOLUME0 = 30*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 0  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
k_TIR1rnaexpression = 2.37e-12 # based on average expression rate of RNAs in yeast
k_RFPrnaexpression = 4e-12 # based on average expression rate of RNAs in yeast
k_BFPrnaexpression = 4e-12 # based on average expression rate of RNAs in yeast
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536
k_translation = 52769.9671394 # um per liter per mRNA per second 
k_leakydegrad = 0.000150354734602 # umolar per second
k_degrad = 4.1718216004e-09 # umolar per second
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)


# Mother 0 G1

G1_length = 91*60
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=30*(10**(-15))
G1_final_volume=60*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)

s0 = (AUXIN0, mRNA_TIR10, TIR10, mRNA_RFP0, RFP0, mRNA_BFP0, BFP0, VOLUME0)
params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad, percentage_sd,avogadro, kdil)
#params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, 0,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)

#s_obs=assymetrical_cell_division(s0,t_G1,params)

#results = s_obs[:]

k_TIR1rnaexpression=1e-13
s0 = (AUXIN0, mRNA_TIR10, TIR10, mRNA_RFP0, RFP0, mRNA_BFP0, BFP0, VOLUME0)
rfp_starts=[]
rfp_ends=[]
expression_levels=[]
step=0
#for i in range(10):
while k_TIR1rnaexpression <= 1e-11:
    induction_level=0
    params=(permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad, percentage_sd,avogadro, kdil,induction_level)
    calibrate=[]
    round_results=[]
    cal,dum,dummy = assymetrical_cell_division(s0,t_G1,params)
    calibrate.append(cal)
    calibrate = np.array(calibrate)
    induction_level=750
    params=(permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad, percentage_sd,avogadro, kdil,induction_level)
    mRNA_TIR10=np.average(calibrate[:][...,1], axis=0)[-1]
    TIR10=np.average(calibrate[:][...,2], axis=0)[-1]
    mRNA_RFP0=np.average(calibrate[:][...,3], axis=0)[-1]
    RFP0=np.average(calibrate[:][...,4], axis=0)[-1]
    rfp_start_list=[]
    rfp_end_list=[]
    s0 = (AUXIN0, mRNA_TIR10, TIR10, mRNA_RFP0, RFP0, mRNA_BFP0, BFP0, VOLUME0)
    gen,a,b=assymetrical_cell_division(s0,t_G1,params)
    rfp_start_list.append(a)
    rfp_end_list.append(b)
    rfp_starts.append(sum(rfp_start_list)/len(rfp_start_list))
    rfp_ends.append(sum(rfp_end_list)/len(rfp_end_list))
    expression_levels.append(k_TIR1rnaexpression)
    k_TIR1rnaexpression+=1e-14
    step+=1
    print("step: ",step)


rfp_differences = []
zip_rfps=zip(rfp_starts,rfp_ends)
for i, j in zip_rfps:
    rfp_differences.append(i-j)
optimum_expression=expression_levels[rfp_differences.index(max(rfp_differences))]
print(max(rfp_differences))
print(optimum_expression)

fig=plt.figure(figsize=(12,8))

ax1=fig.add_subplot(1,1,1)

ax1.plot(expression_levels, rfp_differences,'k',label="Drop in RFP concentration")

ax1.set_xlabel("Expression level of TIR1 mRNA, in umol per litre")
ax1.set_ylabel("Drop in RFP concentration over 45 mins, in umol per litre")
ax1.legend()

plt.show()
