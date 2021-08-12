# When I transitioned to github, some older versions of files ended up without a parent folder, so I am pasting them here:

#ADCTF model
from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
import random

def sdot_mother(s,t,params):
	IAAH, IAM, ADCTF, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced, permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil # 1.38 comes from the average surface area to volume ratio in haploid yeast: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3910951/
	dADCTF = kconstitutive - kdegrad*ADCTF*AUXIN - ADCTF*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - 10000*permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = kinduced*ADCTF - VIOLACEIN*kdil
	dPIN2 = kconstitutive - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dIAAH, dIAM, dADCTF, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_budding(s,t,params):
	IAAH, IAM, ADCTF, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil
	dADCTF = kconstitutive* - kdegrad*ADCTF*AUXIN - ADCTF*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - 10000*permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = kinduced*ADCTF - VIOLACEIN*kdil
	dPIN2 = kconstitutive - PIN2*kdil
	dVOLUME = VOLUME*kdil


	ds = [dIAAH, dIAM, dADCTF, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_daughter1(s,t,params):
	IAAH, IAM, ADCTF, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil
	dADCTF = kconstitutive - kdegrad*ADCTF*AUXIN - ADCTF*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - 10000*permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = kinduced*ADCTF - VIOLACEIN*kdil
	dPIN2 = kconstitutive - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dIAAH, dIAM, dADCTF, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_daughter2(s,t,params):
	IAAH, IAM, ADCTF, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil
	dADCTF = kconstitutive - kdegrad*ADCTF*AUXIN - ADCTF*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = kinduced*ADCTF - VIOLACEIN*kdil
	dPIN2 = - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dIAAH, dIAM, dADCTF, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

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
	ADCTF_0 = 0.378 # start at steady state
	AUXIN_0 = 0
	VIOLACEIN_0 = 1 # start at steady state
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
	kdegrad = (5.78*(10**(-4)))# LAST ONE TO FIGURE OUT, READ https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3323547/#:~:text=The%20protein%20degradation%20kinetics%20affects,state%20is%20maintained%20the%20same.
	kinduced = 0.000743 #based on a violacein target steady state of 1uM
	permeability_IAAH = 3.89*(10**(-7)) # m*s^-1 
	permeability_IAA =  2.25*(10**(-5)) # m*s^-1 per umol/L of PIN2
	fIAAH_w = 0.25
	fIAAH_c = 0.0004
	Nu = 5.02



	s0_G1 = IAAH_0, IAM_0, ADCTF_0, AUXIN_0, VIOLACEIN_0, PIN2_0, VOLUME_0
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_G1 = odeint(sdot_mother, s0_G1, t_G1, args=(params,))

	IAAH_G1 = s_G1[:,0]
	IAM_G1 = s_G1[:,1]
	ADCTF_G1 = s_G1[:,2]
	AUXIN_G1 = s_G1[:,3]
	VIOLACEIN_G1 = s_G1[:,4]
	PIN2_G1 = s_G1[:,5]

	bud_length = 71*60
	t_bud = np.linspace(0,bud_length,bud_length)
	bud_initial_volume=38*(10**(-15))
	bud_final_volume=60*(10**(-15))

	kdil = calc_dillution(bud_length,bud_initial_volume,bud_final_volume)

	s0_bud = IAAH_G1[-1], IAM_G1[-1], ADCTF_G1[-1], AUXIN_G1[-1], VIOLACEIN_G1[-1], PIN2_G1[-1], bud_initial_volume
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_bud = odeint(sdot_budding, s0_bud, t_bud, args=(params,))

	IAAH_bud = s_bud[:,0]
	IAM_bud = s_bud[:,1]
	ADCTF_bud = s_bud[:,2]
	AUXIN_bud = s_bud[:,3]
	VIOLACEIN_bud = s_bud[:,4]
	PIN2_bud = s_bud[:,5]

	daughter1_length = 6.6*60
	t_daughter1 = np.linspace(0,daughter1_length,daughter1_length)
	daughter1_initial_volume = 34*(10**(-15)) 
	daughter1_final_volume = 38*(10**(-15)) 

	kdil = calc_dillution(daughter1_length,daughter1_initial_volume,daughter1_final_volume)

	daughter1_split=((bud_final_volume)*(daughter1_initial_volume/bud_final_volume))/daughter1_initial_volume

	s0_daughter1 = IAAH_bud[-1]*daughter1_split, IAM_bud[-1]*daughter1_split, ADCTF_bud[-1]*daughter1_split, AUXIN_bud[-1]*daughter1_split, VIOLACEIN_bud[-1]*daughter1_split, PIN2_bud[-1]*daughter1_split, daughter1_initial_volume
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_daughter1 = odeint(sdot_daughter1, s0_daughter1, t_daughter1, args=(params,))

	IAAH_daughter1 = s_daughter1[:,0]
	IAM_daughter1 = s_daughter1[:,1]
	ADCTF_daughter1 = s_daughter1[:,2]
	AUXIN_daughter1 = s_daughter1[:,3]
	VIOLACEIN_daughter1 = s_daughter1[:,4]
	PIN2_daughter1 = s_daughter1[:,5]

	daughter2_length = 20*60
	t_daughter2 = np.linspace(0,daughter2_length,daughter2_length)
	daughter2_initial_volume = 26*(10**(-15))
	daughter2_final_volume = 38*(10**(-15))

	kdil = calc_dillution(daughter2_length,daughter2_initial_volume,daughter2_final_volume)

	daughter2_split=((bud_final_volume)*(daughter2_initial_volume/bud_final_volume))/daughter2_initial_volume


	s0_daughter2 = IAAH_bud[-1]*daughter2_split, IAM_bud[-1]*daughter2_split, ADCTF_bud[-1]*daughter2_split, AUXIN_bud[-1]*daughter2_split, VIOLACEIN_bud[-1]*daughter2_split, PIN2_bud[-1]*daughter2_split, daughter2_initial_volume
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_daughter2 = odeint(sdot_daughter2, s0_daughter2, t_daughter2, args=(params,))

	IAAH_daughter2 = s_daughter2[:,0]
	IAM_daughter2 = s_daughter2[:,1]
	ADCTF_daughter2 = s_daughter2[:,2]
	AUXIN_daughter2 = s_daughter2[:,3]
	VIOLACEIN_daughter2 = s_daughter2[:,4]
	PIN2_daughter2 = s_daughter2[:,5]

	fig1 = plt.figure(figsize=(5,6))
	ax1 = fig1.add_subplot(2,1,1)
	ax2 = fig1.add_subplot(2,1,2)
	fig2 = plt.figure(figsize=(5,6))
	ax3 = fig2.add_subplot(2,1,1)
	ax4 = fig2.add_subplot(2,1,2)


	#ax1.plot(t_G1, IAAH_G1,'b-',label="IAAH")
	ax1.plot(t_G1, IAM_G1,'g-',label="IAM")
	ax1.plot(t_G1, ADCTF_G1,'c-',label="ADCTF")
	ax1.plot(t_G1, AUXIN_G1,'r.',label="AUXIN")
	ax1.plot(t_G1, VIOLACEIN_G1,'m-',label="VIOLACEIN")
	ax1.plot(t_G1, PIN2_G1,'g--',label="PIN2")

	#ax2.plot(t_bud, IAAH_bud,'b-',label="IAAH")
	ax2.plot(t_bud, IAM_bud,'g-',label="IAM")
	ax2.plot(t_bud, ADCTF_bud,'c-',label="ADCTF")
	ax2.plot(t_bud, AUXIN_bud,'r.',label="AUXIN")
	ax2.plot(t_bud, VIOLACEIN_bud,'m-',label="VIOLACEIN")
	ax2.plot(t_bud, PIN2_bud,'g--',label="PIN2")

	#ax3.plot(t_daughter1, IAAH_daughter1,'b-',label="IAAH")
	ax3.plot(t_daughter1, IAM_daughter1,'g-',label="IAM")
	ax3.plot(t_daughter1, ADCTF_daughter1,'c-',label="ADCTF")
	ax3.plot(t_daughter1, AUXIN_daughter1,'r.',label="AUXIN")
	ax3.plot(t_daughter1, VIOLACEIN_daughter1,'m-',label="VIOLACEIN")
	ax3.plot(t_daughter1, PIN2_daughter1,'g--',label="PIN2")

	#ax4.plot(t_daughter2, IAAH_daughter2,'b-',label="IAAH")
	ax4.plot(t_daughter2, IAM_daughter2,'g-',label="IAM")
	ax4.plot(t_daughter2, ADCTF_daughter2,'c-',label="ADCTF")
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

#design 1

from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
import random

# units: umol/l for concentration, litre for volume , um for distance (diffusions)

def sdot_mother(s,t,params):
	AUXIN, REPRESSOR, VIOLACEIN, PIN2, VOLUME  = s
	permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition = params

	dAUXIN = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)/VOLUME - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN - AUXIN*kdil
	dREPRESSOR = kconstitutive - (kdegrad*AUXIN*REPRESSOR) - REPRESSOR*kdil #/(AUXIN+Ksaturation when you figureout hill coefficients)
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dAUXIN, dREPRESSOR, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_budding(s,t,params):
	AUXIN, REPRESSOR, VIOLACEIN, PIN2, VOLUME  = s
	permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition = params

	dAUXIN = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)/VOLUME - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN - AUXIN*kdil # Permeability should increase with surface area, not decrease
	dREPRESSOR = kconstitutive - kdegrad*AUXIN*REPRESSOR - REPRESSOR*kdil
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = - PIN2*kdil
	dVOLUME = VOLUME*kdil


	ds = [dAUXIN, dREPRESSOR, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_daughter1a(s,t,params):
	AUXIN, REPRESSOR, VIOLACEIN, PIN2, VOLUME  = s
	permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition = params

	dAUXIN = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)/VOLUME - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN - AUXIN*kdil
	dREPRESSOR = kconstitutive - kdegrad*AUXIN*REPRESSOR - REPRESSOR*kdil
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dAUXIN, dREPRESSOR, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_daughter1b(s,t,params):
	AUXIN, REPRESSOR, VIOLACEIN, PIN2, VOLUME  = s
	permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition = params

	dAUXIN = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)/VOLUME - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN - AUXIN*kdil
	dREPRESSOR = kconstitutive - kdegrad*AUXIN*REPRESSOR - REPRESSOR*kdil
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = kconstitutive - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dAUXIN, dREPRESSOR, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def calc_dillution(time_of_phase, initial_volume, final_volume):
	kdil=(math.log(final_volume/initial_volume))/time_of_phase
	print(kdil)
	return kdil


def assymetrical_cell_division():

# Some useful yeast bionumbers:
#	- Number of genes in yeast: 6275
#	- Total number of mRNAs per cell: 15000
#	- average number of proteins per mrna: 4000
#	- Therefore, average number of proteins per gene: 9500
#	- That means on average 0.378 umol/L of each protein per gene

	permeability_IAAH = .389 # um*s^-1  
	fIAAH_w = 0.25
	conc_out = 0.5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
	fIAAH_c = 0.0004
	pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
	permeability_IAA = 1.081*(10**3) # um*s^-1 per umol/L of PIN2
	Nu = 5.02
	kconstitutive = 9*(10**(-5)) # based on the average expression rate of proteins in yeast
	kdegrad = 5.28 * (10**(-6)) # umolar per second
	kinhibition = 0.003818 # Made up "realistic value", if repressor at steady state, promoter is 99% repressed

	AUXIN_0 = 0
	REPRESSOR_0 = 0.378 # start at steady state
	VIOLACEIN_0 = 0 
	PIN2_0 = 0
	VOLUME_0 = 26*(10**(-15))

	# Mother 0 G1

	G1_length = 20*60
	t_G1 = np.linspace(0,G1_length,G1_length)
	G1_initial_volume=26*(10**(-15))
	G1_final_volume=38*(10**(-15))

	kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)


	s0_G1 = AUXIN_0, REPRESSOR_0, VIOLACEIN_0, PIN2_0, VOLUME_0
	params = permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition

	s_G1 = odeint(sdot_mother, s0_G1, t_G1, args=(params,))

	AUXIN_G1 = s_G1[:,0]
	REPRESSOR_G1 = s_G1[:,1]
	VIOLACEIN_G1 = s_G1[:,2]
	PIN2_G1 = s_G1[:,3]

	# Mother 0 Bud length

	bud_length = 71*60
	t_bud = np.linspace(0,bud_length,bud_length)
	bud_initial_volume=38*(10**(-15))
	bud_final_volume=60*(10**(-15))

	kdil = calc_dillution(bud_length,bud_initial_volume,bud_final_volume)

	s0_bud = AUXIN_G1[-1], REPRESSOR_G1[-1], VIOLACEIN_G1[-1], PIN2_G1[-1], bud_initial_volume
	params = permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition

	s_bud = odeint(sdot_budding, s0_bud, t_bud, args=(params,))

	AUXIN_bud = s_bud[:,0]
	REPRESSOR_bud = s_bud[:,1]
	VIOLACEIN_bud = s_bud[:,2]
	PIN2_bud = s_bud[:,3]

	# Daughter1a up to citokynesis, used to be mother cell

	daughter1a_length = (6.6+71)*60
	t_daughter1a = np.linspace(0,daughter1a_length,daughter1a_length)
	daughter1a_initial_volume = 34*(10**(-15)) 
	daughter1a_final_volume = 60*(10**(-15)) 

	kdil = calc_dillution(daughter1a_length,daughter1a_initial_volume,daughter1a_final_volume)

	daughter1a_split=((bud_final_volume)*(daughter1a_initial_volume/bud_final_volume))/daughter1a_initial_volume

	s0_daughter1a = AUXIN_bud[-1]*daughter1a_split, REPRESSOR_bud[-1]*daughter1a_split, VIOLACEIN_bud[-1]*daughter1a_split, ((PIN2_bud[-1]*bud_final_volume)/daughter1a_initial_volume), daughter1a_initial_volume
	params = permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition

	s_daughter1a = odeint(sdot_daughter1a, s0_daughter1a, t_daughter1a, args=(params,))

	AUXIN_daughter1a = s_daughter1a[:,0]
	REPRESSOR_daughter1a = s_daughter1a[:,1]
	VIOLACEIN_daughter1a = s_daughter1a[:,2]
	PIN2_daughter1a = s_daughter1a[:,3]

	# Daughter1b up to citokinesis, used to be bud

	daughter1b_length = 91*60
	t_daughter1b = np.linspace(0,daughter1b_length,daughter1b_length)
	daughter1b_initial_volume = 26*(10**(-15))
	daughter1b_final_volume = 60*(10**(-15))

	kdil = calc_dillution(daughter1b_length,daughter1b_initial_volume,daughter1b_final_volume)

	daughter1b_split=((bud_final_volume)*(daughter1b_initial_volume/bud_final_volume))/daughter1b_initial_volume


	s0_daughter1b = AUXIN_bud[-1]*daughter1b_split, REPRESSOR_bud[-1]*daughter1b_split, VIOLACEIN_bud[-1]*daughter1b_split, 0, daughter1b_initial_volume
	params = permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition

	s_daughter1b = odeint(sdot_daughter1b, s0_daughter1b, t_daughter1b, args=(params,))

	# Daughter2ba, used to be mother, keeps PIN2 on its membrane until it dilutes out

	# Daugther2bb, used to be bud, keeps 


	AUXIN_daughter1b = s_daughter1b[:,0]
	REPRESSOR_daughter1b = s_daughter1b[:,1]
	VIOLACEIN_daughter1b = s_daughter1b[:,2]
	PIN2_daughter1b = s_daughter1b[:,3]

	fig1 = plt.figure(figsize=(5,6))
	ax1 = fig1.add_subplot(2,1,1)
	ax2 = fig1.add_subplot(2,1,2)
	fig2 = plt.figure(figsize=(5,6))
	ax3 = fig2.add_subplot(2,1,1)
	ax4 = fig2.add_subplot(2,1,2)
	fig3 = plt.figure(figsize=(5,6))
	ax5 = fig3.add_subplot(1,1,1)


	ax1.plot(t_G1, REPRESSOR_G1,'c-',label="REPRESSOR")
	#ax1.plot(t_G1, AUXIN_G1,'r-',label="AUXIN")
	ax1.plot(t_G1, VIOLACEIN_G1,'m-',label="VIOLACEIN")
	ax1.plot(t_G1, PIN2_G1,'g-',label="PIN2")

	ax2.plot(t_bud, REPRESSOR_bud,'c-',label="REPRESSOR")
	#ax2.plot(t_bud, AUXIN_bud,'r-',label="AUXIN")
	ax2.plot(t_bud, VIOLACEIN_bud,'m-',label="VIOLACEIN")
	ax2.plot(t_bud, PIN2_bud,'g-',label="PIN2")

	ax3.plot(t_daughter1a, REPRESSOR_daughter1a,'c-',label="REPRESSOR")
	#ax3.plot(t_daughter1a, AUXIN_daughter1a,'r-',label="AUXIN")
	ax3.plot(t_daughter1a, VIOLACEIN_daughter1a,'m-',label="VIOLACEIN")
	ax3.plot(t_daughter1a, PIN2_daughter1a,'g-',label="PIN2")

	ax4.plot(t_daughter1b, REPRESSOR_daughter1b,'c-',label="REPRESSOR")
	#ax4.plot(t_daughter1b, AUXIN_daughter1b,'r-',label="AUXIN")
	ax4.plot(t_daughter1b, VIOLACEIN_daughter1b,'m-',label="VIOLACEIN")
	ax4.plot(t_daughter1b, PIN2_daughter1b,'g-',label="PIN2")

	ax5.plot(t_G1, AUXIN_G1,'c-',label="AUXIN of mother")
	ax5.plot(t_bud, AUXIN_bud,'r-',label="AUXIN of bud")
	ax5.plot(t_daughter1a, AUXIN_daughter1a,'m-',label="AUXIN of daughter 1")
	ax5.plot(t_daughter1b, AUXIN_daughter1b,'g-',label="AUXIN of daughter 2")

	ax1.legend()
	ax2.legend()
	ax3.legend()
	ax4.legend()
	ax5.legend()

	ax1.set_xlabel("time (s)")
	ax2.set_xlabel("time (s)")
	ax3.set_xlabel("time (s)")
	ax4.set_xlabel("time (s)")
	ax5.set_xlabel("time (s)")


	ax1.set_ylabel("Concentration, in M")
	ax2.set_ylabel("Concentration, in M")
	ax3.set_ylabel("Concentration, in M")
	ax4.set_ylabel("Concentration, in M")
	ax5.set_ylabel("Concentration, in M")


	ax1.set_title("Mother cell in G1")
	ax2.set_title("Budding cell")
	ax3.set_title("Daughter 1a in G1 (Used to be mother cell)")
	ax4.set_title("Daughter 1b in G1 (Used to be budding cell)")

	plt.show()

assymetrical_cell_division()

# Stochastic trial 3

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

# auxin reporter module:
	# Use this to find the theoretically perfect level of TIR1 expressions

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

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad, percentage_sd,avogadro, kdil = params
    
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
        rate_influx = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
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

        t=t+next_event

        AUXIN+= (np.random.poisson(rates[0]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_TIR1+= (np.random.poisson(rates[1]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_TIR1-= (np.random.poisson(rates[2]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	TIR1+= (np.random.poisson(rates[3]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	mRNA_RFP+= (np.random.poisson(rates[4]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_RFP-= (np.random.poisson(rates[5]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	RFP+= (np.random.poisson(rates[6]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	RFP-= (np.random.poisson(rates[7]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	RFP-= (np.random.poisson(rates[8]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	mRNA_BFP+= (np.random.poisson(rates[9]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_BFP-= (np.random.poisson(rates[10]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	BFP+= (np.random.poisson(rates[11]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	event_dilution = np.random.normal((rates[12]*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        mRNA_TIR1-= mRNA_TIR1-((mRNA_TIR1*VOLUME)/(VOLUME+event_dilution))
        TIR1-= TIR1-((TIR1*VOLUME)/(VOLUME+event_dilution))
        mRNA_RFP-= mRNA_RFP-((mRNA_RFP*VOLUME)/(VOLUME+event_dilution))
        RFP-= RFP-((RFP*VOLUME)/(VOLUME+event_dilution))
        mRNA_BFP-= mRNA_BFP-((mRNA_BFP*VOLUME)/(VOLUME+event_dilution))
        BFP-= BFP-((BFP*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution

        s = (AUXIN, mRNA_TIR1, TIR1, mRNA_RFP, RFP, mRNA_BFP, BFP, VOLUME)

        t_obs.append(t)
        s_obs.append(s)

        if t >= 1200*2.5 and induction == "false":
        	conc_out=0.5
        	induction="true"


    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)

#setting seed to make results reproducible

np.random.seed=1
random.seed=1

AUXIN0 = 0
mRNA_TIR10 = 9.09364077747e-10
TIR10 = 0.378
mRNA_RFP0 = 9.09364077747e-10
RFP0 = 0.378
mRNA_BFP0 = 9.09364077747e-10
BFP0 = 0.378
VOLUME0 = 30*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 0  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
k_TIR1rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_RFPrnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_BFPrnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536
k_translation = 52769.9671394 # um per liter per mRNA per second 
k_leakydegrad = 3.056 * (10**(-4)) # umolar per second
k_degrad = 5.28 * (10**(-6)) # umolar per second
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)


# Mother 0 G1

G1_length = 91*60
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=30*(10**(-15))
G1_final_volume=60*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)
print(kdil)

s0 = (AUXIN0, mRNA_TIR10, TIR10, mRNA_RFP0, RFP0, mRNA_BFP0, BFP0, VOLUME0)
params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad, percentage_sd,avogadro, kdil)
#params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, 0,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)

#s_obs=gillespie(s0,t_G1,params)

#results = s_obs[:]


results=[]
n_runs=20

for i in range(n_runs):
    print("Simulating {} of {} cells...".format(i+1,n_runs))
    gen=gillespie(s0,t_G1,params)
    results.append(gen)

results=np.array(results)
#print(results[0])
#print(results[0][...,0])

fig=plt.figure(figsize=(12,8))

ax1=fig.add_subplot(5,1,1)
ax2=fig.add_subplot(5,1,2)
ax3=fig.add_subplot(5,1,3)
ax4=fig.add_subplot(5,1,4)
ax5=fig.add_subplot(5,1,5)

# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_G1, results[:][...,0].T,'k:')
ax2.plot(t_G1, results[:][...,2].T,'y:')
ax3.plot(t_G1, results[:][...,4].T,'r:')
ax4.plot(t_G1, results[:][...,6].T,'b:')
ax5.plot(t_G1, results[:][...,7].T,'g:')

auxin_mean = np.average(results[:][...,0], axis=0,)
TIR1_mean = np.average(results[:][...,2], axis=0)
RFP_mean = np.average(results[:][...,4], axis=0)
BFP_mean = np.average(results[:][...,6], axis=0)
volume_mean = np.average(results[:][...,7], axis=0)

ax1.plot(t_G1, auxin_mean,'k',label="Auxin")
ax2.plot(t_G1, TIR1_mean,'k',label="TIR1")
ax3.plot(t_G1, RFP_mean,'k',label="RFP")
ax4.plot(t_G1, BFP_mean,'k',label="BFP")
ax5.plot(t_G1, volume_mean,'k',label="Volume")

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
ax5.set_ylabel("Volume in litres")
ax5.legend()



#ax1=fig2.add_subplot(5,1,1)
#ax2=fig2.add_subplot(5,1,2)
#ax3=fig2.add_subplot(5,1,3)
#ax4=fig2.add_subplot(5,1,4)
#ax5=fig2.add_subplot(5,1,5)


#t_all = np.linspace(0,G1_length*n_generations+G1_length,G1_length*n_generations+G1_length)

#ax1.plot(t_all, ori_results[:][...,0].flatten().tolist(),'g',label="Auxin")
#ax2.plot(t_all, ori_results[:][...,2].flatten().tolist(),'b',label="TIR1")
#ax3.plot(t_all, ori_results[:][...,4].flatten().tolist(),'y',label="RFP")
#ax4.plot(t_all, ori_results[:][...,6].flatten().tolist(),'r',label="BFP")
#ax5.plot(t_all, ori_results[:][...,7].flatten().tolist(),'r',label="Volume")

#ax1.set_xlabel("Time, in minutes")
#ax1.set_ylabel("Concentration in umol per litre")
#ax1.legend()

#ax2.set_xlabel("Time, in minutes")
#ax2.set_ylabel("Concentration in umol per litre")
#ax2.legend()

#ax3.set_xlabel("Time, in minutes")
#ax3.set_ylabel("Concentration in umol per litre")
#ax3.legend()

#ax4.set_xlabel("Time, in minutes")
#ax4.set_ylabel("Concentration in umol per litre")
#ax4.legend()

#ax5.set_xlabel("Time, in minutes")
#ax5.set_ylabel("Volume in litres")
#ax5.legend()





plt.show()


# ADCTF system # 2

from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
import random

def sdot_mother(s,t,params):
	IAAH, IAM, ADCTF, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced, permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil # 1.38 comes from the average surface area to volume ratio in haploid yeast: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3910951/
	dADCTF = kconstitutive - kdegrad*ADCTF*AUXIN - ADCTF*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - 10000*permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = kinduced*ADCTF - VIOLACEIN*kdil
	dPIN2 = kconstitutive - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dIAAH, dIAM, dADCTF, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_budding(s,t,params):
	IAAH, IAM, ADCTF, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil
	dADCTF = kconstitutive* - kdegrad*ADCTF*AUXIN - ADCTF*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - 10000*permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = kinduced*ADCTF - VIOLACEIN*kdil
	dPIN2 = kconstitutive - PIN2*kdil
	dVOLUME = VOLUME*kdil


	ds = [dIAAH, dIAM, dADCTF, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_daughter1(s,t,params):
	IAAH, IAM, ADCTF, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil
	dADCTF = kconstitutive - kdegrad*ADCTF*AUXIN - ADCTF*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - 10000*permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = kinduced*ADCTF - VIOLACEIN*kdil
	dPIN2 = kconstitutive - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dIAAH, dIAM, dADCTF, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_daughter2(s,t,params):
	IAAH, IAM, ADCTF, AUXIN, VIOLACEIN, PIN2, VOLUME  = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu  = params

	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil
	dADCTF = kconstitutive - kdegrad*ADCTF*AUXIN - ADCTF*kdil
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - AUXIN*kdil - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN # Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVIOLACEIN = kinduced*ADCTF - VIOLACEIN*kdil
	dPIN2 = - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dIAAH, dIAM, dADCTF, dAUXIN, dVIOLACEIN, dPIN2, dVOLUME]

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
	ADCTF_0 = 0.378 # start at steady state
	AUXIN_0 = 0
	VIOLACEIN_0 = 1 # start at steady state
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
	kdegrad = (5.78*(10**(-4)))# LAST ONE TO FIGURE OUT, READ https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3323547/#:~:text=The%20protein%20degradation%20kinetics%20affects,state%20is%20maintained%20the%20same.
	kinduced = 0.000743 #based on a violacein target steady state of 1uM
	permeability_IAAH = 3.89*(10**(-7)) # m*s^-1 
	permeability_IAA =  2.25*(10**(-5)) # m*s^-1 per umol/L of PIN2
	fIAAH_w = 0.25
	fIAAH_c = 0.0004
	Nu = 5.02



	s0_G1 = IAAH_0, IAM_0, ADCTF_0, AUXIN_0, VIOLACEIN_0, PIN2_0, VOLUME_0
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_G1 = odeint(sdot_mother, s0_G1, t_G1, args=(params,))

	IAAH_G1 = s_G1[:,0]
	IAM_G1 = s_G1[:,1]
	ADCTF_G1 = s_G1[:,2]
	AUXIN_G1 = s_G1[:,3]
	VIOLACEIN_G1 = s_G1[:,4]
	PIN2_G1 = s_G1[:,5]

	bud_length = 71*60
	t_bud = np.linspace(0,bud_length,bud_length)
	bud_initial_volume=38*(10**(-15))
	bud_final_volume=60*(10**(-15))

	kdil = calc_dillution(bud_length,bud_initial_volume,bud_final_volume)

	s0_bud = IAAH_G1[-1], IAM_G1[-1], ADCTF_G1[-1], AUXIN_G1[-1], VIOLACEIN_G1[-1], PIN2_G1[-1], bud_initial_volume
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_bud = odeint(sdot_budding, s0_bud, t_bud, args=(params,))

	IAAH_bud = s_bud[:,0]
	IAM_bud = s_bud[:,1]
	ADCTF_bud = s_bud[:,2]
	AUXIN_bud = s_bud[:,3]
	VIOLACEIN_bud = s_bud[:,4]
	PIN2_bud = s_bud[:,5]

	daughter1_length = 6.6*60
	t_daughter1 = np.linspace(0,daughter1_length,daughter1_length)
	daughter1_initial_volume = 34*(10**(-15)) 
	daughter1_final_volume = 38*(10**(-15)) 

	kdil = calc_dillution(daughter1_length,daughter1_initial_volume,daughter1_final_volume)

	daughter1_split=((bud_final_volume)*(daughter1_initial_volume/bud_final_volume))/daughter1_initial_volume

	s0_daughter1 = IAAH_bud[-1]*daughter1_split, IAM_bud[-1]*daughter1_split, ADCTF_bud[-1]*daughter1_split, AUXIN_bud[-1]*daughter1_split, VIOLACEIN_bud[-1]*daughter1_split, PIN2_bud[-1]*daughter1_split, daughter1_initial_volume
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_daughter1 = odeint(sdot_daughter1, s0_daughter1, t_daughter1, args=(params,))

	IAAH_daughter1 = s_daughter1[:,0]
	IAM_daughter1 = s_daughter1[:,1]
	ADCTF_daughter1 = s_daughter1[:,2]
	AUXIN_daughter1 = s_daughter1[:,3]
	VIOLACEIN_daughter1 = s_daughter1[:,4]
	PIN2_daughter1 = s_daughter1[:,5]

	daughter2_length = 20*60
	t_daughter2 = np.linspace(0,daughter2_length,daughter2_length)
	daughter2_initial_volume = 26*(10**(-15))
	daughter2_final_volume = 38*(10**(-15))

	kdil = calc_dillution(daughter2_length,daughter2_initial_volume,daughter2_final_volume)

	daughter2_split=((bud_final_volume)*(daughter2_initial_volume/bud_final_volume))/daughter2_initial_volume


	s0_daughter2 = IAAH_bud[-1]*daughter2_split, IAM_bud[-1]*daughter2_split, ADCTF_bud[-1]*daughter2_split, AUXIN_bud[-1]*daughter2_split, VIOLACEIN_bud[-1]*daughter2_split, PIN2_bud[-1]*daughter2_split, daughter2_initial_volume
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, kinduced,  permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nu

	s_daughter2 = odeint(sdot_daughter2, s0_daughter2, t_daughter2, args=(params,))

	IAAH_daughter2 = s_daughter2[:,0]
	IAM_daughter2 = s_daughter2[:,1]
	ADCTF_daughter2 = s_daughter2[:,2]
	AUXIN_daughter2 = s_daughter2[:,3]
	VIOLACEIN_daughter2 = s_daughter2[:,4]
	PIN2_daughter2 = s_daughter2[:,5]

	fig1 = plt.figure(figsize=(5,6))
	ax1 = fig1.add_subplot(2,1,1)
	ax2 = fig1.add_subplot(2,1,2)
	fig2 = plt.figure(figsize=(5,6))
	ax3 = fig2.add_subplot(2,1,1)
	ax4 = fig2.add_subplot(2,1,2)


	#ax1.plot(t_G1, IAAH_G1,'b-',label="IAAH")
	ax1.plot(t_G1, IAM_G1,'g-',label="IAM")
	ax1.plot(t_G1, ADCTF_G1,'c-',label="ADCTF")
	ax1.plot(t_G1, AUXIN_G1,'r.',label="AUXIN")
	ax1.plot(t_G1, VIOLACEIN_G1,'m-',label="VIOLACEIN")
	ax1.plot(t_G1, PIN2_G1,'g--',label="PIN2")

	#ax2.plot(t_bud, IAAH_bud,'b-',label="IAAH")
	ax2.plot(t_bud, IAM_bud,'g-',label="IAM")
	ax2.plot(t_bud, ADCTF_bud,'c-',label="ADCTF")
	ax2.plot(t_bud, AUXIN_bud,'r.',label="AUXIN")
	ax2.plot(t_bud, VIOLACEIN_bud,'m-',label="VIOLACEIN")
	ax2.plot(t_bud, PIN2_bud,'g--',label="PIN2")

	#ax3.plot(t_daughter1, IAAH_daughter1,'b-',label="IAAH")
	ax3.plot(t_daughter1, IAM_daughter1,'g-',label="IAM")
	ax3.plot(t_daughter1, ADCTF_daughter1,'c-',label="ADCTF")
	ax3.plot(t_daughter1, AUXIN_daughter1,'r.',label="AUXIN")
	ax3.plot(t_daughter1, VIOLACEIN_daughter1,'m-',label="VIOLACEIN")
	ax3.plot(t_daughter1, PIN2_daughter1,'g--',label="PIN2")

	#ax4.plot(t_daughter2, IAAH_daughter2,'b-',label="IAAH")
	ax4.plot(t_daughter2, IAM_daughter2,'g-',label="IAM")
	ax4.plot(t_daughter2, ADCTF_daughter2,'c-',label="ADCTF")
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

# Auxin control module:

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

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay = params
    
    AUXIN, mRNA_PIN2, PIN2, VOLUME  = s0

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

        types=['influx','efflux','expression','RNA_decay','synthesis','growth']
        rate_influx = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
        rate_efflux = permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN
        rate_expression = k_rnaexpression
        rate_decay = k_rnadecay*mRNA_PIN2
        rate_synthesis = k_PIN2expression*mRNA_PIN2
        rate_growth = VOLUME * kdil
        rates=[rate_influx,rate_efflux,rate_expression,rate_decay,rate_synthesis,rate_growth]
        position=0
        for i in rates:
            if i < 0:
                rates[position]=0
            position += 1

        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        next_event=gen_next_event_time(rate_all_events)

        t=t+next_event

        AUXIN+= (np.random.poisson(rates[0]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        AUXIN-= (np.random.poisson(rates[1]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2+= (np.random.poisson(rates[2]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2-= (np.random.poisson(rates[3]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        PIN2+= (np.random.poisson(rates[4]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = np.random.normal((rates[5]*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution



        s = (AUXIN, mRNA_PIN2, PIN2, VOLUME)

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
conc_out = 0.5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081*(10**(-1)) # um*s^-1 per umol/L of PIN2
Nu = 5.02
k_rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_PIN2expression = 52769.9671394 # um per liter per mRNA per second 
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)
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


gen0=gillespie(s0,t_G1,params_daughter)

ori_results=[]
daughter_results=[]

ori_results.append(gen0)
s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.9)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))
s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.1)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))


n_generations=10

for i in range(n_generations):
    print("Simulating {} of {} generations...".format(i+1,n_generations))
    d2=gillespie(s0_daughter,t_G1,params_daughter)
    daughter_results.append(d2)
    d1=gillespie(s0_ori,t_G1,params)
    ori_results.append(d1)
    s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.9)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))
    s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.1)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))


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
ax2.plot(t_G1, mrna_mean,'k',label="mRNA")
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
ax4.set_ylabel("Concentration in umol per litre")
ax4.legend()



ax1=fig2.add_subplot(4,1,1)
ax2=fig2.add_subplot(4,1,2)
ax3=fig2.add_subplot(4,1,3)
ax4=fig2.add_subplot(4,1,4)

t_all = np.linspace(0,G1_length*n_generations+G1_length,G1_length*n_generations+G1_length)

ax1.plot(t_all, ori_results[:][...,0].flatten().tolist(),'g',label="Auxin")
ax2.plot(t_all, ori_results[:][...,1].flatten().tolist(),'b',label="mRNA")
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
ax4.set_ylabel("Concentration in umol per litre")
ax4.legend()





plt.show()


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

        types=['influx','efflux','expression','RNA_decay','synthesis','growth']
        rate_influx = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
        if rate_influx < 0:
            rate_influx = 0
        rate_efflux = permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN
        if rate_efflux < 0:
            rate_efflux = 0
        rate_expression = k_rnaexpression
        rate_decay = k_rnadecay*mRNA_PIN2
        rate_synthesis = k_PIN2expression*mRNA_PIN2
        rate_growth = VOLUME * kdil
        rates=[rate_influx,rate_efflux,rate_expression,rate_decay,rate_synthesis,rate_growth]
        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        next_event=gen_next_event_time(rate_all_events)

        t=t+next_event

        AUXIN+= (np.random.poisson(rate_influx*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        AUXIN-= (np.random.poisson(rate_efflux*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2+= (np.random.poisson(rate_expression*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2-= (np.random.poisson(rate_decay*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        PIN2+= (np.random.poisson(rate_synthesis*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = np.random.normal((rate_growth*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution



        s = (AUXIN, mRNA_PIN2, PIN2, VOLUME)

        t_obs.append(t)
        s_obs.append(s)

    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)

#setting seed to make results reproducible

np.random.seed=1
random.seed=1

AUXIN0 = 0
mRNA_PIN20 = 9.09364077747e-10
PIN20 = .378
VOLUME0 = 26*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 0.5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081*(10**(-1)) # um*s^-1 per umol/L of PIN2
Nu = 5.02
k_rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_PIN2expression = 131453.654001 # um per liter per mRNA per second 
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536


# Mother 0 G1

G1_length = 60*90
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=26*(10**(-15))
G1_final_volume=38*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)


s0 = (AUXIN0,mRNA_PIN20,PIN20,VOLUME0)
params = (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)

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

ax1=fig.add_subplot(4,1,1)
ax2=fig.add_subplot(4,1,2)
ax3=fig.add_subplot(4,1,3)
ax4=fig.add_subplot(4,1,4)


# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_G1, results[:,0],'r-',label='AUXIN')
ax2.plot(t_G1, results[:,1],'g-',label='PIN2_mRNA')
ax3.plot(t_G1, results[:,2],'g-',label='PIN2')
ax4.plot(t_G1, results[:,3],'b-',label='VOLUME')


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

# Auxin control module:

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

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay = params
    
    AUXIN, mRNA_PIN2, PIN2, VOLUME  = s0

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

        types=['influx','efflux','expression','RNA_decay','synthesis','growth']
        rate_influx = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
        rate_efflux = permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN
        rate_expression = k_rnaexpression
        rate_decay = k_rnadecay*mRNA_PIN2
        rate_synthesis = k_PIN2expression*mRNA_PIN2
        rate_growth = VOLUME * kdil
        rates=[rate_influx,rate_efflux,rate_expression,rate_decay,rate_synthesis,rate_growth]
        position=0
        for i in rates:
            if i < 0:
                rates[position]=0
            position += 1

        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        next_event=gen_next_event_time(rate_all_events)

        t=t+next_event

        AUXIN+= (np.random.poisson(rates[0]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        AUXIN-= (np.random.poisson(rates[1]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2+= (np.random.poisson(rates[2]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2-= (np.random.poisson(rates[3]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        PIN2+= (np.random.poisson(rates[4]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = np.random.normal((rates[5]*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution



        s = (AUXIN, mRNA_PIN2, PIN2, VOLUME)

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
conc_out = 0.5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081*(10**(-1)) # um*s^-1 per umol/L of PIN2
Nu = 5.02
k_rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_PIN2expression = 52769.9671394 # um per liter per mRNA per second 
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)
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


gen0=gillespie(s0,t_G1,params_daughter)

ori_results=[]
daughter_results=[]

ori_results.append(gen0)
s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.9)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))
s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.1)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))


n_generations=10

for i in range(n_generations):
    print("Simulating {} of {} generations...".format(i+1,n_generations))
    d2=gillespie(s0_daughter,t_G1,params_daughter)
    daughter_results.append(d2)
    d1=gillespie(s0_ori,t_G1,params)
    ori_results.append(d1)
    s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.9)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))
    s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.1)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))


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
ax2.plot(t_G1, mrna_mean,'k',label="mRNA")
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
ax4.set_ylabel("Concentration in umol per litre")
ax4.legend()



ax1=fig2.add_subplot(4,1,1)
ax2=fig2.add_subplot(4,1,2)
ax3=fig2.add_subplot(4,1,3)
ax4=fig2.add_subplot(4,1,4)

t_all = np.linspace(0,G1_length*n_generations+G1_length,G1_length*n_generations+G1_length)

ax1.plot(t_all, ori_results[:][...,0].flatten().tolist(),'g',label="Auxin")
ax2.plot(t_all, ori_results[:][...,1].flatten().tolist(),'b',label="mRNA")
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
ax4.set_ylabel("Concentration in umol per litre")
ax4.legend()





plt.show()


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

        types=['influx','efflux','expression','RNA_decay','synthesis','growth']
        rate_influx = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
        if rate_influx < 0:
            rate_influx = 0
        rate_efflux = permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN
        if rate_efflux < 0:
            rate_efflux = 0
        rate_expression = k_rnaexpression
        rate_decay = k_rnadecay*mRNA_PIN2
        rate_synthesis = k_PIN2expression*mRNA_PIN2
        rate_growth = VOLUME * kdil
        rates=[rate_influx,rate_efflux,rate_expression,rate_decay,rate_synthesis,rate_growth]
        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        next_event=gen_next_event_time(rate_all_events)

        t=t+next_event

        AUXIN+= (np.random.poisson(rate_influx*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        AUXIN-= (np.random.poisson(rate_efflux*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2+= (np.random.poisson(rate_expression*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2-= (np.random.poisson(rate_decay*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        PIN2+= (np.random.poisson(rate_synthesis*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = np.random.normal((rate_growth*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution



        s = (AUXIN, mRNA_PIN2, PIN2, VOLUME)

        t_obs.append(t)
        s_obs.append(s)

    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)

#setting seed to make results reproducible

np.random.seed=1
random.seed=1

AUXIN0 = 0
mRNA_PIN20 = 9.09364077747e-10
PIN20 = .378
VOLUME0 = 26*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 0.5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081*(10**(-1)) # um*s^-1 per umol/L of PIN2
Nu = 5.02
k_rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_PIN2expression = 131453.654001 # um per liter per mRNA per second 
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536


# Mother 0 G1

G1_length = 60*90
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=26*(10**(-15))
G1_final_volume=38*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)


s0 = (AUXIN0,mRNA_PIN20,PIN20,VOLUME0)
params = (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)

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

ax1=fig.add_subplot(4,1,1)
ax2=fig.add_subplot(4,1,2)
ax3=fig.add_subplot(4,1,3)
ax4=fig.add_subplot(4,1,4)


# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_G1, results[:,0],'r-',label='AUXIN')
ax2.plot(t_G1, results[:,1],'g-',label='PIN2_mRNA')
ax3.plot(t_G1, results[:,2],'g-',label='PIN2')
ax4.plot(t_G1, results[:,3],'b-',label='VOLUME')


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

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay = params
    
    AUXIN, mRNA_PIN2, PIN2, VOLUME  = s0

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

        types=['influx','efflux','expression','RNA_decay','synthesis','growth']
        rate_influx = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
        if rate_influx < 0:
            rate_influx = 0
        rate_efflux = permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN
        if rate_efflux < 0:
            rate_efflux = 0
        rate_expression = k_rnaexpression
        rate_decay = k_rnadecay*mRNA_PIN2
        rate_synthesis = k_PIN2expression*mRNA_PIN2
        rate_growth = VOLUME * kdil
        rates=[rate_influx,rate_efflux,rate_expression,rate_decay,rate_synthesis,rate_growth]
        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        next_event=gen_next_event_time(rate_all_events)

        t=t+next_event

        AUXIN+= (np.random.poisson(rate_influx*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        AUXIN-= (np.random.poisson(rate_efflux*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2+= (np.random.poisson(rate_expression*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2-= (np.random.poisson(rate_decay*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        PIN2+= (np.random.poisson(rate_synthesis*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = np.random.normal((rate_growth*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution



        s = (AUXIN, mRNA_PIN2, PIN2, VOLUME)

        t_obs.append(t)
        s_obs.append(s)

    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)

#setting seed to make results reproducible

np.random.seed=1
random.seed=1

AUXIN0 = 0
mRNA_PIN20 = 9.09364077747e-10
PIN20 = .378
VOLUME0 = 26*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 0.5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081*(10**(-1)) # um*s^-1 per umol/L of PIN2
Nu = 5.02
k_rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_PIN2expression = 131453.654001 # um per liter per mRNA per second 
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536


# Mother 0 G1

G1_length = 60*90
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=26*(10**(-15))
G1_final_volume=38*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)


s0 = (AUXIN0,mRNA_PIN20,PIN20,VOLUME0)
params = (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)

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

ax1=fig.add_subplot(4,1,1)
ax2=fig.add_subplot(4,1,2)
ax3=fig.add_subplot(4,1,3)
ax4=fig.add_subplot(4,1,4)


# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_G1, results[:,0],'r-',label='AUXIN')
ax2.plot(t_G1, results[:,1],'g-',label='PIN2_mRNA')
ax3.plot(t_G1, results[:,2],'g-',label='PIN2')
ax4.plot(t_G1, results[:,3],'b-',label='VOLUME')


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


def gillespie(s0,t_obs_out,params):

    #--0--# Unpack parameters and species variables

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad, percentage_sd,avogadro, kdil = params
    
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
        rate_influx = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
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

        t=t+next_event

        AUXIN+= (np.random.poisson(rates[0]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_TIR1+= (np.random.poisson(rates[1]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_TIR1-= (np.random.poisson(rates[2]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	TIR1+= (np.random.poisson(rates[3]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	mRNA_RFP+= (np.random.poisson(rates[4]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_RFP-= (np.random.poisson(rates[5]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	RFP+= (np.random.poisson(rates[6]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	RFP-= (np.random.poisson(rates[7]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	RFP-= (np.random.poisson(rates[8]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	mRNA_BFP+= (np.random.poisson(rates[9]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_BFP-= (np.random.poisson(rates[10]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	BFP+= (np.random.poisson(rates[11]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
       	event_dilution = np.random.normal((rates[12]*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        mRNA_TIR1-= mRNA_TIR1-((mRNA_TIR1*VOLUME)/(VOLUME+event_dilution))
        TIR1-= TIR1-((TIR1*VOLUME)/(VOLUME+event_dilution))
        mRNA_RFP-= mRNA_RFP-((mRNA_RFP*VOLUME)/(VOLUME+event_dilution))
        RFP-= RFP-((RFP*VOLUME)/(VOLUME+event_dilution))
        mRNA_BFP-= mRNA_BFP-((mRNA_BFP*VOLUME)/(VOLUME+event_dilution))
        BFP-= BFP-((BFP*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution

        s = (AUXIN, mRNA_TIR1, TIR1, mRNA_RFP, RFP, mRNA_BFP, BFP, VOLUME)

        t_obs.append(t)
        s_obs.append(s)

        if t >= 1200*2.5 and induction == "false":
        	conc_out=0.5
        	induction="true"


    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)

#setting seed to make results reproducible

np.random.seed=1
random.seed=1

AUXIN0 = 0
mRNA_TIR10 = 9.09364077747e-10
TIR10 = 0.378
mRNA_RFP0 = 9.09364077747e-10
RFP0 = 0.378
mRNA_BFP0 = 9.09364077747e-10
BFP0 = 0.378
VOLUME0 = 30*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 0  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
k_TIR1rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_RFPrnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_BFPrnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536
k_translation = 52769.9671394 # um per liter per mRNA per second 
k_leakydegrad = 3.056 * (10**(-4)) # umolar per second
k_degrad = 5.28 * (10**(-6)) # umolar per second
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)


# Mother 0 G1

G1_length = 91*60
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=30*(10**(-15))
G1_final_volume=60*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)
print(kdil)

s0 = (AUXIN0, mRNA_TIR10, TIR10, mRNA_RFP0, RFP0, mRNA_BFP0, BFP0, VOLUME0)
params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, k_TIR1rnaexpression, k_RFPrnaexpression, k_BFPrnaexpression, k_rnadecay, k_translation, k_leakydegrad, k_degrad, percentage_sd,avogadro, kdil)
#params= (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, 0,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)

#s_obs=gillespie(s0,t_G1,params)

#results = s_obs[:]


results=[]
n_runs=20

for i in range(n_runs):
    print("Simulating {} of {} cells...".format(i+1,n_runs))
    gen=gillespie(s0,t_G1,params)
    results.append(gen)

results=np.array(results)
#print(results[0])
#print(results[0][...,0])

fig=plt.figure(figsize=(12,8))

ax1=fig.add_subplot(5,1,1)
ax2=fig.add_subplot(5,1,2)
ax3=fig.add_subplot(5,1,3)
ax4=fig.add_subplot(5,1,4)
ax5=fig.add_subplot(5,1,5)

# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_G1, results[:][...,0].T,'k:')
ax2.plot(t_G1, results[:][...,2].T,'y:')
ax3.plot(t_G1, results[:][...,4].T,'r:')
ax4.plot(t_G1, results[:][...,6].T,'b:')
ax5.plot(t_G1, results[:][...,7].T,'g:')

auxin_mean = np.average(results[:][...,0], axis=0,)
TIR1_mean = np.average(results[:][...,2], axis=0)
RFP_mean = np.average(results[:][...,4], axis=0)
BFP_mean = np.average(results[:][...,6], axis=0)
volume_mean = np.average(results[:][...,7], axis=0)

ax1.plot(t_G1, auxin_mean,'k',label="Auxin")
ax2.plot(t_G1, TIR1_mean,'k',label="TIR1")
ax3.plot(t_G1, RFP_mean,'k',label="RFP")
ax4.plot(t_G1, BFP_mean,'k',label="BFP")
ax5.plot(t_G1, volume_mean,'k',label="Volume")

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
ax5.set_ylabel("Volume in litres")
ax5.legend()



#ax1=fig2.add_subplot(5,1,1)
#ax2=fig2.add_subplot(5,1,2)
#ax3=fig2.add_subplot(5,1,3)
#ax4=fig2.add_subplot(5,1,4)
#ax5=fig2.add_subplot(5,1,5)


#t_all = np.linspace(0,G1_length*n_generations+G1_length,G1_length*n_generations+G1_length)

#ax1.plot(t_all, ori_results[:][...,0].flatten().tolist(),'g',label="Auxin")
#ax2.plot(t_all, ori_results[:][...,2].flatten().tolist(),'b',label="TIR1")
#ax3.plot(t_all, ori_results[:][...,4].flatten().tolist(),'y',label="RFP")
#ax4.plot(t_all, ori_results[:][...,6].flatten().tolist(),'r',label="BFP")
#ax5.plot(t_all, ori_results[:][...,7].flatten().tolist(),'r',label="Volume")

#ax1.set_xlabel("Time, in minutes")
#ax1.set_ylabel("Concentration in umol per litre")
#ax1.legend()

#ax2.set_xlabel("Time, in minutes")
#ax2.set_ylabel("Concentration in umol per litre")
#ax2.legend()

#ax3.set_xlabel("Time, in minutes")
#ax3.set_ylabel("Concentration in umol per litre")
#ax3.legend()

#ax4.set_xlabel("Time, in minutes")
#ax4.set_ylabel("Concentration in umol per litre")
#ax4.legend()

#ax5.set_xlabel("Time, in minutes")
#ax5.set_ylabel("Volume in litres")
#ax5.legend()





plt.show()



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

# Auxin control module:

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

    permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay = params
    
    AUXIN, mRNA_PIN2, PIN2, VOLUME  = s0

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

        types=['influx','efflux','expression','RNA_decay','synthesis','growth']
        rate_influx = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
        rate_efflux = permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN
        rate_expression = k_rnaexpression
        rate_decay = k_rnadecay*mRNA_PIN2
        rate_synthesis = k_PIN2expression*mRNA_PIN2
        rate_growth = VOLUME * kdil
        rates=[rate_influx,rate_efflux,rate_expression,rate_decay,rate_synthesis,rate_growth]
        position=0
        for i in rates:
            if i < 0:
                rates[position]=0
            position += 1

        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        next_event=gen_next_event_time(rate_all_events)

        t=t+next_event

        AUXIN+= (np.random.poisson(rates[0]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        AUXIN-= (np.random.poisson(rates[1]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2+= (np.random.poisson(rates[2]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2-= (np.random.poisson(rates[3]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        PIN2+= (np.random.poisson(rates[4]*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = np.random.normal((rates[5]*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution



        s = (AUXIN, mRNA_PIN2, PIN2, VOLUME)

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
conc_out = 0.5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081*(10**(-1)) # um*s^-1 per umol/L of PIN2
Nu = 5.02
k_rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_PIN2expression = 52769.9671394 # um per liter per mRNA per second 
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)
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


gen0=gillespie(s0,t_G1,params_daughter)

ori_results=[]
daughter_results=[]

ori_results.append(gen0)
s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.9)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))
s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.1)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))


n_generations=10

for i in range(n_generations):
    print("Simulating {} of {} generations...".format(i+1,n_generations))
    d2=gillespie(s0_daughter,t_G1,params_daughter)
    daughter_results.append(d2)
    d1=gillespie(s0_ori,t_G1,params)
    ori_results.append(d1)
    s0_ori=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.9)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))
    s0_daughter=(ori_results[-1][-1][0],ori_results[-1][-1][1],((ori_results[-1][-1][2]*ori_results[-1][-1][-1]*avogadro*.1)/(avogadro*ori_results[-1][-1][-1]/2)),(ori_results[-1][-1][-1]/2))


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
ax2.plot(t_G1, mrna_mean,'k',label="mRNA")
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
ax4.set_ylabel("Concentration in umol per litre")
ax4.legend()



ax1=fig2.add_subplot(4,1,1)
ax2=fig2.add_subplot(4,1,2)
ax3=fig2.add_subplot(4,1,3)
ax4=fig2.add_subplot(4,1,4)

t_all = np.linspace(0,G1_length*n_generations+G1_length,G1_length*n_generations+G1_length)

ax1.plot(t_all, ori_results[:][...,0].flatten().tolist(),'g',label="Auxin")
ax2.plot(t_all, ori_results[:][...,1].flatten().tolist(),'b',label="mRNA")
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
ax4.set_ylabel("Concentration in umol per litre")
ax4.legend()





plt.show()


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

        types=['influx','efflux','expression','RNA_decay','synthesis','growth']
        rate_influx = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)#/VOLUME
        if rate_influx < 0:
            rate_influx = 0
        rate_efflux = permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN
        if rate_efflux < 0:
            rate_efflux = 0
        rate_expression = k_rnaexpression
        rate_decay = k_rnadecay*mRNA_PIN2
        rate_synthesis = k_PIN2expression*mRNA_PIN2
        rate_growth = VOLUME * kdil
        rates=[rate_influx,rate_efflux,rate_expression,rate_decay,rate_synthesis,rate_growth]
        rate_all_events=sum(rates)
        if rate_all_events <= 0:
            break
        next_event=gen_next_event_time(rate_all_events)

        t=t+next_event

        AUXIN+= (np.random.poisson(rate_influx*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        AUXIN-= (np.random.poisson(rate_efflux*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2+= (np.random.poisson(rate_expression*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        mRNA_PIN2-= (np.random.poisson(rate_decay*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        PIN2+= (np.random.poisson(rate_synthesis*VOLUME*avogadro*next_event))/(VOLUME*avogadro)
        event_dilution = np.random.normal((rate_growth*next_event),(rate_growth*next_event*percentage_sd))
        AUXIN-= AUXIN-((AUXIN*VOLUME)/(VOLUME+event_dilution))
        PIN2-= PIN2-((PIN2*VOLUME)/(VOLUME+event_dilution))
        VOLUME+= event_dilution



        s = (AUXIN, mRNA_PIN2, PIN2, VOLUME)

        t_obs.append(t)
        s_obs.append(s)

    s_obs_out=resample_observations(t_obs,s_obs,t_obs_out)
    return np.array(s_obs_out)

#setting seed to make results reproducible

np.random.seed=1
random.seed=1

AUXIN0 = 0
mRNA_PIN20 = 9.09364077747e-10
PIN20 = .378
VOLUME0 = 26*(10**(-15))

permeability_IAAH = .389 # um*s^-1  
fIAAH_w = 0.25
conc_out = 0.5  # Taken out of the AID2 paper Assuming an average intracellular concentration of auxin of 23.8 umol per litre, look at system 1 notes, 0.03808
fIAAH_c = 0.0004
pm_thickness = 0.0092 #um https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=0&id=108569#:~:text=%22The%20average%20thickness%20of%20the,al.%2C%201994).%22
permeability_IAA = 1.081*(10**(-1)) # um*s^-1 per umol/L of PIN2
Nu = 5.02
k_rnaexpression = 2.18862203686e-12 # based on average expression rate of RNAs in yeast
k_PIN2expression = 131453.654001 # um per liter per mRNA per second 
percentage_sd = 0.00636 # Standard deviation of Sc BY4741: https://www.sciencedirect.com/science/article/pii/S2468501120300201
avogadro = 6.022*(10**23)
k_rnadecay = 0.00240676104375 # per second, https://elifesciences.org/articles/32536


# Mother 0 G1

G1_length = 60*90
t_G1 = np.linspace(0,G1_length,G1_length)
G1_initial_volume=26*(10**(-15))
G1_final_volume=38*(10**(-15))

kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)


s0 = (AUXIN0,mRNA_PIN20,PIN20,VOLUME0)
params = (permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, k_rnaexpression,k_PIN2expression,percentage_sd,avogadro,k_rnadecay)

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

ax1=fig.add_subplot(4,1,1)
ax2=fig.add_subplot(4,1,2)
ax3=fig.add_subplot(4,1,3)
ax4=fig.add_subplot(4,1,4)


# plot all runs at once (the .T transposes the matrix)
ax1.plot(t_G1, results[:,0],'r-',label='AUXIN')
ax2.plot(t_G1, results[:,1],'g-',label='PIN2_mRNA')
ax3.plot(t_G1, results[:,2],'g-',label='PIN2')
ax4.plot(t_G1, results[:,3],'b-',label='VOLUME')


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

