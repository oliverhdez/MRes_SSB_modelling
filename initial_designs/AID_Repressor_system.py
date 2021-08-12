from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
import random

def sdot_mother(s,t,params):
	IAAH, IAM, REPRESSOR, TIR1, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinhibition, permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil # 1.38 comes from the average surface area to volume ratio in haploid yeast: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3910951/
	dREPRESSOR = kconstitutive - kdegrad*REPRESSOR*AUXIN*TIR1 - REPRESSOR*kdil
	dTIR1= kconstitutive - TIR1*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = kconstitutive - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dIAAH, dIAM, dREPRESSOR, dTIR1, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_budding(s,t,params):
	IAAH, IAM, REPRESSOR, TIR1, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinhibition,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil
	dREPRESSOR = kconstitutive - kdegrad*REPRESSOR*AUXIN*TIR1 - REPRESSOR*kdil
	dTIR1= kconstitutive - TIR1*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = kconstitutive - PIN2*kdil
	dVOLUME = VOLUME*kdil


	ds = [dIAAH, dIAM, dREPRESSOR, dTIR1, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_daughter1(s,t,params):
	IAAH, IAM, REPRESSOR, TIR1, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinhibition,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil
	dREPRESSOR = kconstitutive - kdegrad*REPRESSOR*AUXIN*TIR1 - REPRESSOR*kdil
	dTIR1= kconstitutive - TIR1*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = kconstitutive - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dIAAH, dIAM, dREPRESSOR, dTIR1, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_daughter2(s,t,params):
	IAAH, IAM, REPRESSOR, TIR1, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinhibition,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil
	dREPRESSOR = kconstitutive - kdegrad*REPRESSOR*AUXIN*TIR1 - REPRESSOR*kdil
	dTIR1= kconstitutive - TIR1*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dIAAH, dIAM, dREPRESSOR, dTIR1, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def calc_dillution(time_of_phase, initial_volume, final_volume):
	kdil=(math.log(final_volume/initial_volume))/time_of_phase
	return kdil


def assymetrical_cell_division():

# Some useful yeast bionumbers:
#	- Number of genes in yeast: 6275
#	- Total number of mRNAs per cell: 15000
#	- average number of mRNAs per transcript: 4000
#	- Therefore, average number of proteins per gene: 9500
#	- That means on average 0.378 umol/L of each protein per gene

	IAAH_0 = 0.378 # start at steady state
	IAM_0 = 0
	REPRESSOR_0 = 0.378 # start at steady state
	TIR1_0 = 0.378
	AUXIN_0 = 0
	VIOLACEIN_0 = 0 # start at steady state
	PIN2_0 = 0.378
	VOLUME_0 = 26*(10**(-15))

	G1_length = 20*60
	t_G1 = np.linspace(0,G1_length,G1_length)
	G1_initial_volume=26*(10**(-15))
	G1_final_volume=38*(10**(-15))

	kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)
	kconstitutive = 9*(10**(-5)) # based on the average expression rate of proteins in yeast
	kinflux = 0.00356/((4.4*(9**(-14)))*1.38) # Based on the target steady state of 10uM of IAM inside cells. Denominator is the average surface area of a haploid yeast cell.
	km = 1.2 #uM https://watermark.silverchair.com/35-3-329.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAs0wggLJBgkqhkiG9w0BBwagggK6MIICtgIBADCCAq8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMVAMlhJRLBt7XRuUuAgEQgIICgOY-zSswl1IkJ0msTGpSQDV0gEXkTVDmHnSdHfqEtEPv84lI5vXRN9uAfIC2xGEHKzx_7T9_7z8arhnKuV_EKYPOika6o27cl6Q9QBs1-pmhq4Hz_6LU2S2XGU4DmPKbgPzmkQUi2WOQsKntvAkLlNQWY4Ur-98HNhl4s0C2ZBtr2l5ay79PESYKcnBvQIR3ihYWzV-AGbW6zrTreDkjCn5NY5vWh-jCh4Fru8ZHgTRDLaqUs2lGicB6_2H4A4mhgOnF0usWH1kgQ4BnZsJDOX9x52b0TLFqxy0g-QadE4MTzdQ1iPUszyD6wmwIAg7l6JxYfmnv0HOdwnUYctBLC1niMNh6bgDkgPG-n5ympzK47KJyOwIog32c2pdnvk4P936-fkle0R5IU-wu6DtKNvvbc2UekYptB2E8a0r5sKabGFwac6tO75b84EE2b0ek6A7p5ru741dJJGntE46Afdm6hL99oWifR0P23EBZEfFW3EkFSo_lwB4agOdqkNsWCllYrZ-855O0uJB_cjAa5F9Jhcl5K1lFxZL4flGB1AV_VhIVz9tO-NlMv9uNrgSPrn6iKhcp8wGuCEygUdtQa7manU3ki_mAfa5EQ5xdCWPTx9ixzSWx9uFvsvQlXSTldmI07FVaxODdSpaRDhG4m3PTyxhBwJz8m8A2W8SVWrpHiLmSXvcbCFcnCtwG6qXrF2wna7s8IcIVfUZD4OYVKTmIC1G5CN2dVR2fiVfxZ-49m5bJuIoi_Bc4IkPtnQ0gRBQ-FiJuUUwaWspidbNlOM6iI7kQuUObFs5b75iavm7ILRm2MzvXsWkBeEGY6uxg-yYMreGVlOfxDTleYAGS20s
	kcat = 0.2202/60
	kdegrad = (5.78*(10**(-2)))# LAST ONE TO FIGURE OUT, READ https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3323547/#:~:text=The%20protein%20degradation%20kinetics%20affects,state%20is%20maintained%20the%20same.
	kinhibition = 0.003818 # Made up "realistic value", if repressor at steady state, promoter is 99% repressed
	permeability_IAAH = 3.89*(10**(-7)) # m*s^-1 
	permeability_IAA =  2.25*(10**(-5)) # m*s^-1 per umol/L of PIN2
	fIAAH_w = 0.25
	fIAAH_c = 0.0004
	Nu = 5.02



	s0_G1 = IAAH_0, IAM_0, REPRESSOR_0, TIR1_0, AUXIN_0, VIOLACEIN_0, PIN2_0, VOLUME_0
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinhibition,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_G1 = odeint(sdot_mother, s0_G1, t_G1, args=(params,))

	IAAH_G1 = s_G1[:,0]
	IAM_G1 = s_G1[:,1]
	REPRESSOR_G1 = s_G1[:,2]
	TIR1_G1 = s_G1[:,3]
	AUXIN_G1 = s_G1[:,4]
	VIOLACEIN_G1 = s_G1[:,5]
	PIN2_G1 = s_G1[:,6]

	bud_length = 71*60
	t_bud = np.linspace(0,bud_length,bud_length)
	bud_initial_volume=38*(10**(-15))
	bud_final_volume=60*(10**(-15))

	kdil = calc_dillution(bud_length,bud_initial_volume,bud_final_volume)

	s0_bud = IAAH_G1[-1], IAM_G1[-1], REPRESSOR_G1[-1], TIR1_G1[-1], AUXIN_G1[-1], VIOLACEIN_G1[-1], PIN2_G1[-1], bud_initial_volume
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinhibition,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_bud = odeint(sdot_budding, s0_bud, t_bud, args=(params,))

	IAAH_bud = s_bud[:,0]
	IAM_bud = s_bud[:,1]
	REPRESSOR_bud = s_bud[:,2]
	TIR1_bud = s_bud[:,3]
	AUXIN_bud = s_bud[:,4]
	VIOLACEIN_bud = s_bud[:,5]
	PIN2_bud = s_bud[:,6]

	daughter1_length = 6.6*60
	t_daughter1 = np.linspace(0,daughter1_length,daughter1_length)
	daughter1_initial_volume = 34*(10**(-15)) 
	daughter1_final_volume = 38*(10**(-15)) 

	kdil = calc_dillution(daughter1_length,daughter1_initial_volume,daughter1_final_volume)

	daughter1_split=((bud_final_volume)*(daughter1_initial_volume/bud_final_volume))/daughter1_initial_volume

	s0_daughter1 = IAAH_bud[-1]*daughter1_split, IAM_bud[-1]*daughter1_split, REPRESSOR_bud[-1]*daughter1_split, TIR1_bud[-1]*daughter1_split, AUXIN_bud[-1]*daughter1_split, VIOLACEIN_bud[-1]*daughter1_split, ((PIN2_bud[-1]*bud_final_volume)/daughter1_initial_volume), daughter1_initial_volume
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinhibition,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_daughter1 = odeint(sdot_daughter1, s0_daughter1, t_daughter1, args=(params,))

	IAAH_daughter1 = s_daughter1[:,0]
	IAM_daughter1 = s_daughter1[:,1]
	REPRESSOR_daughter1 = s_daughter1[:,2]
	TIR1_daughter1 = s_daughter1[:,3]
	AUXIN_daughter1 = s_daughter1[:,4]
	VIOLACEIN_daughter1 = s_daughter1[:,5]
	PIN2_daughter1 = s_daughter1[:,6]

	daughter2_length = 20*60
	t_daughter2 = np.linspace(0,daughter2_length,daughter2_length)
	daughter2_initial_volume = 26*(10**(-15))
	daughter2_final_volume = 38*(10**(-15))

	kdil = calc_dillution(daughter2_length,daughter2_initial_volume,daughter2_final_volume)

	daughter2_split=((bud_final_volume)*(daughter2_initial_volume/bud_final_volume))/daughter2_initial_volume


	s0_daughter2 = IAAH_bud[-1]*daughter2_split, IAM_bud[-1]*daughter2_split, REPRESSOR_bud[-1]*daughter2_split, TIR1_bud[-1]*daughter2_split, AUXIN_bud[-1]*daughter2_split, VIOLACEIN_bud[-1]*daughter2_split, 0, daughter2_initial_volume
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinhibition,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_daughter2 = odeint(sdot_daughter2, s0_daughter2, t_daughter2, args=(params,))

	IAAH_daughter2 = s_daughter2[:,0]
	IAM_daughter2 = s_daughter2[:,1]
	REPRESSOR_daughter2 = s_daughter2[:,2]
	TIR1_daughter2 = s_daughter2[:,3]
	AUXIN_daughter2 = s_daughter2[:,4]
	VIOLACEIN_daughter2 = s_daughter2[:,5]
	PIN2_daughter2 = s_daughter2[:,6]

	fig1 = plt.figure(figsize=(5,6))
	ax1 = fig1.add_subplot(2,1,1)
	ax2 = fig1.add_subplot(2,1,2)
	fig2 = plt.figure(figsize=(5,6))
	ax3 = fig2.add_subplot(2,1,1)
	ax4 = fig2.add_subplot(2,1,2)


	#ax1.plot(t_G1, IAAH_G1,'b-',label="IAAH")
	#ax1.plot(t_G1, IAM_G1,'g-',label="IAM")
	ax1.plot(t_G1, REPRESSOR_G1,'c-',label="REPRESSOR")
	ax1.plot(t_G1, TIR1_G1,'y-',label="TIR1")
	ax1.plot(t_G1, AUXIN_G1,'r.',label="AUXIN")
	ax1.plot(t_G1, VIOLACEIN_G1,'m-',label="VIOLACEIN")
	ax1.plot(t_G1, PIN2_G1,'g--',label="PIN2")

	#ax2.plot(t_bud, IAAH_bud,'b-',label="IAAH")
	#ax2.plot(t_bud, IAM_bud,'g-',label="IAM")
	ax2.plot(t_bud, REPRESSOR_bud,'c-',label="REPRESSOR")
	ax2.plot(t_bud, TIR1_bud,'y-',label="TIR1")
	ax2.plot(t_bud, AUXIN_bud,'r.',label="AUXIN")
	ax2.plot(t_bud, VIOLACEIN_bud,'m-',label="VIOLACEIN")
	ax2.plot(t_bud, PIN2_bud,'g--',label="PIN2")

	#ax3.plot(t_daughter1, IAAH_daughter1,'b-',label="IAAH")
	#ax3.plot(t_daughter1, IAM_daughter1,'g-',label="IAM")
	ax3.plot(t_daughter1, REPRESSOR_daughter1,'c-',label="REPRESSOR")
	ax3.plot(t_daughter1, TIR1_daughter1,'y-',label="TIR1")
	ax3.plot(t_daughter1, AUXIN_daughter1,'r.',label="AUXIN")
	ax3.plot(t_daughter1, VIOLACEIN_daughter1,'m-',label="VIOLACEIN")
	ax3.plot(t_daughter1, PIN2_daughter1,'g--',label="PIN2")

	#ax4.plot(t_daughter2, IAAH_daughter2,'b-',label="IAAH")
	#ax4.plot(t_daughter2, IAM_daughter2,'g-',label="IAM")
	ax4.plot(t_daughter2, REPRESSOR_daughter2,'c-',label="REPRESSOR")
	ax4.plot(t_daughter2, TIR1_daughter2,'y-',label="TIR1")
	ax4.plot(t_daughter2, AUXIN_daughter2,'r.',label="AUXIN")
	ax4.plot(t_daughter2, VIOLACEIN_daughter2,'m-',label="VIOLACEIN")
	ax4.plot(t_daughter2, PIN2_daughter2,'g--',label="PIN2")

	ax1.legend()
	ax2.legend()
	ax3.legend()
	ax4.legend()

	ax1.set_xlabel("time (s)")
	ax2.set_xlabel("time (s)")
	ax3.set_xlabel("time (s)")
	ax4.set_xlabel("time (s)")

	ax1.set_ylabel("Concentration, in M")
	ax2.set_ylabel("Concentration, in M")
	ax3.set_ylabel("Concentration, in M")
	ax4.set_ylabel("Concentration, in M")

	ax1.set_title("Mother cell in G1")
	ax2.set_title("Budding cell")
	ax3.set_title("Daughter 1 in G1 (Used to be mother cell)")
	ax4.set_title("Daughter 2 in G1 (Used to be budding cell)")

	plt.show()

assymetrical_cell_division()