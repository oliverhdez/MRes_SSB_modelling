from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint
import random

def sdot_auxin(s,t,params):
	IAAH, IAM, ADCTF, AUXIN, VOLUME = s
	kdil, kconstitutive, kinflux, km, kcat, kdegrad, permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nµ = params
	
	dIAAH = kconstitutive - IAAH*kdil
	dIAM = kinflux*VOLUME*1.38 - ((IAAH*kcat*IAM)/(km + IAM)) - IAM*kdil # 1.38 comes from the average surface area to volume ratio in haploid yeast: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3910951/
	dADCTF = kconstitutive - ADCTF*kdil - kdegrad*ADCTF*AUXIN
	dAUXIN = ((IAAH*kcat*IAM)/(km + IAM)) - permeability_IAAH*(fIAAH_c*AUXIN-fIAAH_w*0) - permeability_IAA*0*Nµ*(1-fIAAH_c)*AUXIN - AUXIN*kdil# Assume auxin disperses too fast therefore no influx back into the cell. If looking at emergent properties, simply rearrange the efflux equation into influx
	dVOLUME = VOLUME*kdil

	ds = [dIAAH, dIAM, dADCTF, dAUXIN, dVOLUME]

	return ds

def calc_dillution(time_of_phase, initial_volume, final_volume):
	kdil=((math.log(final_volume/initial_volume))/time_of_phase)
	return kdil

def assymetrical_cell_division():

	IAAH_0 = 0.378 # start at steady state
	IAM_0 = 0
	ADCTF_0 = 0.378 # start at steady state
	AUXIN_0 = 0
	VOLUME_0 = 26*(10**(-15))

	G1_length = 20*60
	t_G1 = np.linspace(0,G1_length,G1_length)

	G1_initial_volume=26*(10**(-15))
	G1_final_volume=38*(10**(-15))

	kdil = calc_dillution(G1_length,G1_initial_volume,G1_final_volume)
	kconstitutive = 0.0001195 # based on the average expression rate of proteins in yeast
	kinflux = 0.00356/((4.4*(10**(-14)))*1.38) # Based on the target steady state of 10µM of IAM inside cells. Denominator is the average surface area of a haploid yeast cell.
	km = 1.2 #µM https://watermark.silverchair.com/35-3-329.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAs0wggLJBgkqhkiG9w0BBwagggK6MIICtgIBADCCAq8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMVAMlhJRLBt7XRuUuAgEQgIICgOY-zSswl1IkJ0msTGpSQDV0gEXkTVDmHnSdHfqEtEPv84lI5vXRN9uAfIC2xGEHKzx_7T9_7z8arhnKuV_EKYPOika6o27cl6Q9QBs1-pmhq4Hz_6LU2S2XGU4DmPKbgPzmkQUi2WOQsKntvAkLlNQWY4Ur-98HNhl4s0C2ZBtr2l5ay79PESYKcnBvQIR3ihYWzV-AGbW6zrTreDkjCn5NY5vWh-jCh4Fru8ZHgTRDLaqUs2lGicB6_2H4A4mhgOnF0usWH1kgQ4BnZsJDOX9x52b0TLFqxy0g-QadE4MTzdQ1iPUszyD6wmwIAg7l6JxYfmnv0HOdwnUYctBLC1niMNh6bgDkgPG-n5ympzK47KJyOwIog32c2pdnvk4P936-fkle0R5IU-wu6DtKNvvbc2UekYptB2E8a0r5sKabGFwac6tO75b84EE2b0ek6A7p5ru741dJJGntE46Afdm6hL99oWifR0P23EBZEfFW3EkFSo_lwB4agOdqkNsWCllYrZ-855O0uJB_cjAa5F9Jhcl5K1lFxZL4flGB1AV_VhIVz9tO-NlMv9uNrgSPrn6iKhcp8wGuCEygUdtQa7manU3ki_mAfa5EQ5xdCWPTx9ixzSWx9uFvsvQlXSTldmI07FVaxODdSpaRDhG4m3PTyxhBwJz8m8A2W8SVWrpHiLmSXvcbCFcnCtwG6qXrF2wna7s8IcIVfUZD4OYVKTmIC1G5CN2dVR2fiVfxZ-49m5bJuIoi_Bc4IkPtnQ0gRBQ-FiJuUUwaWspidbNlOM6iI7kQuUObFs5b75iavm7ILRm2MzvXsWkBeEGY6uxg-yYMreGVlOfxDTleYAGS20s
	kcat = 0.2202/60
	kdegrad = (5.78*(10**(-4)))# LAST ONE TO FIGURE OUT, READ https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3323547/#:~:text=The%20protein%20degradation%20kinetics%20affects,state%20is%20maintained%20the%20same.
	kinduced = 9*(10**(-5))
	kandgate = 9*(10**(-5))
	ki = 4.2*(10**(-7)) # Made up "realistic value" that we can tune later, inhibition constant
	permeability_IAAH = 3.89*(10**(-7)) # m*s^-1 
	permeability_IAA =  2.25*(10**(-5)) # m*s^-1 per µmol/L of PIN2
	fIAAH_w = 0.25
	fIAAH_c = 0.0004
	Nµ = 5.02
	print(kdil)



	s0_G1 = IAAH_0, IAM_0, ADCTF_0, AUXIN_0, VOLUME_0
	params = kdil, kconstitutive, kinflux, km, kcat, kdegrad, permeability_IAAH, permeability_IAA, fIAAH_w, fIAAH_c, Nµ

	s_G1 = odeint(sdot_auxin, s0_G1, t_G1, args=(params,))

	IAAH_G1 = s_G1[:,0]
	IAM_G1 = s_G1[:,1]
	ADCTF_G1 = s_G1[:,2]
	AUXIN_G1 = s_G1[:,3]

	fig1 = plt.figure(figsize=(5,6))
	ax1 = fig1.add_subplot(1,1,1)

	ax1.plot(t_G1, IAAH_G1,'b-',label="IAAH")
	#ax1.plot(t_G1, IAM_G1,'g-',label="IAM")
	#ax1.plot(t_G1, ADCTF_G1,'c-',label="ADCTF")
	#ax1.plot(t_G1, AUXIN_G1,'r.',label="AUXIN")

	ax1.legend()

	ax1.set_xlabel("time (s)")

	ax1.set_ylabel("Concentration, in µM")


	plt.show()


assymetrical_cell_division()

