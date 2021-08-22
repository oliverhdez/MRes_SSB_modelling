from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
import random

# units: umol/l for concentration, litre for volume , um for distance (diffusions)

def sdot_mother(s,t,params):
	AUXIN, REPRESSOR, VIOLACEIN, PIN2, VOLUME  = s
	permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition = params

	dAUXIN = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness) - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN - AUXIN*kdil
	dREPRESSOR = kconstitutive - (kdegrad*AUXIN*REPRESSOR) - REPRESSOR*kdil #/(AUXIN+Ksaturation when you figureout hill coefficients)
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dAUXIN, dREPRESSOR, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_budding(s,t,params):
	AUXIN, REPRESSOR, VIOLACEIN, PIN2, VOLUME  = s
	permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition = params

	dAUXIN = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness) - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN - AUXIN*kdil # Permeability should increase with surface area, not decrease
	dREPRESSOR = kconstitutive - kdegrad*AUXIN*REPRESSOR - REPRESSOR*kdil
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = - PIN2*kdil
	dVOLUME = VOLUME*kdil


	ds = [dAUXIN, dREPRESSOR, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_daughter1a(s,t,params):
	AUXIN, REPRESSOR, VIOLACEIN, PIN2, VOLUME  = s
	permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition = params

	dAUXIN = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness) - permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN - AUXIN*kdil
	dREPRESSOR = kconstitutive - kdegrad*AUXIN*REPRESSOR - REPRESSOR*kdil
	dVIOLACEIN = (kconstitutive/(1+(REPRESSOR/kinhibition))) - VIOLACEIN*kdil
	dPIN2 = - PIN2*kdil
	dVOLUME = VOLUME*kdil

	ds = [dAUXIN, dREPRESSOR, dVIOLACEIN, dPIN2, dVOLUME]

	return ds

def sdot_daughter1b(s,t,params):
	AUXIN, REPRESSOR, VIOLACEIN, PIN2, VOLUME  = s
	permeability_IAAH, fIAAH_w, conc_out, fIAAH_c, pm_thickness, permeability_IAA, Nu, kdil, kconstitutive, kdegrad, kinhibition = params

	dAUXIN = ((permeability_IAAH*((fIAAH_w*conc_out - fIAAH_c*AUXIN)*(4.83597586205*(VOLUME**(2/3)))))/pm_thickness)- permeability_IAA*PIN2*Nu*(1-fIAAH_c)*AUXIN - AUXIN*kdil
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