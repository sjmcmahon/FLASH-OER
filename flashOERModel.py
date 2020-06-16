import numpy as np

# Choose oxygen depletion model
# True makes depletion rate independent of O2 concentraiton
# False makes depletion portional to O2 concentration
linearModel = False

# Set default parameters. O2 units are fraction of 100% O2. 
if linearModel:
	REFoxDep = 0.000518  # Depletion Per Gray
	REFoxRec = 1.0	     # Recovery Per Second
else:
	REFoxDep = 0.053
	REFoxRec = 1.0

# Set mid-point of OER curve to 1% O2
REFOERCenter = 0.010

# Determine oxygen concentration at a time t of exposure to doseRate, starting from initial 
# concentration of O2Mean. t can be single time or array-like list of values.
def oxygenCurve(t, doseRate, O2Mean, oxDepletion=REFoxDep, oxRecovery=REFoxRec):
	t=np.array(t)
	if linearModel:
		estO2 = O2Mean + doseRate*oxDepletion/oxRecovery*(np.exp(-oxRecovery*t)-1)
		return max(0,estO2)

	return ((O2Mean/(doseRate*oxDepletion+oxRecovery)) * 
    	    (doseRate*oxDepletion*np.exp(-(doseRate*oxDepletion+oxRecovery)*t)+oxRecovery) )

# Reference calculation for 'standard' OER curve
# A value of 1.0 is full sensitivity, fully de-oxyenated will return a value of 1/3.
def OER(oxygen, centralPoint):
	oerVal = (3*oxygen+centralPoint)/(oxygen+centralPoint)/3.0
	return oerVal

# Calculate a cumulative OER for an exposure to a dose at a given doserate. 
# A value of 1.0 is full sensitivity, fully de-oxyenated will return a value of 1/3.
def cumulativeOER(dose, doseRate, O2Mean, oxDepletion=REFoxDep, oxRecovery=REFoxRec, 
				  OERCenter=REFOERCenter):
	if dose==0:	return 0

	t = dose/doseRate
	# O2 depletion proportional to O2 concentration and dose rate
	if not linearModel:
		r = oxDepletion*doseRate
		
		# Calculate numerator and denominator for OER
		NumTermOne = t*(3*oxRecovery*O2Mean + OERCenter*(oxRecovery+r) )
		NumTermTwo = 2*OERCenter*np.log( (OERCenter+O2Mean)*(oxRecovery+r) )
		NumTermThree = -2*OERCenter*np.log( OERCenter*(oxRecovery+r)+
			   							    O2Mean*(oxRecovery+r*np.exp(-(oxRecovery+r)*t) ) )
		OERNum = NumTermOne + NumTermTwo + NumTermThree

		OERDen = t*(oxRecovery*O2Mean+OERCenter*(oxRecovery+r))

		return OERNum/OERDen/3

	# O2 depletion independent of O2 concentration
	# Calculate depletion threshold - if this is >1, we'll always have O2 concentration >0.
	# If not, we'll need to work out crossing time, as O2 concentration can't be negative.
	depletionThreshold = O2Mean * oxRecovery / (doseRate*oxDepletion)
	if depletionThreshold <1:
		crossTime = -np.log(1 - depletionThreshold)/oxRecovery
	else:
		crossTime = np.nan

	# Work out integrated OER*T product (E) to enable accomodation of crossover time.
	# If we'll not reach zero O2 during the exposure, just calculate directly.
	if depletionThreshold >= 1 or crossTime>t:
		NumTermOne = (3*doseRate*oxDepletion - OERCenter*oxRecovery - 3*oxRecovery*O2Mean)*t 
		NumTermTwo = 2*OERCenter*np.log( 1 + doseRate*oxDepletion*(np.exp(-oxRecovery*t)-1) / (oxRecovery*(OERCenter+O2Mean)) )
		Den = (3*doseRate*oxDepletion - 3*oxRecovery*(O2Mean+OERCenter) )
		E = (NumTermOne + NumTermTwo)/Den
	# Otherwise, split in two, calculate up to crosstime and then after.
	else:
		NumTermOne = (3*doseRate*oxDepletion - OERCenter*oxRecovery - 3*oxRecovery*O2Mean)*crossTime 
		NumTermTwo = 2*OERCenter*np.log( 1 + doseRate*oxDepletion*(np.exp(-oxRecovery*crossTime)-1) / (oxRecovery*(OERCenter+O2Mean)) )
		Den = (3*doseRate*oxDepletion - 3*oxRecovery*(O2Mean+OERCenter) )
		E = (NumTermOne + NumTermTwo)/Den

		E+=(t-crossTime)/3

	return E/t

# Calculate log survival for a given set of conditions
# Exposures is a list of [dose, doseRate, O2] sub-lists
def predictSurvival(exposures, alpha, beta, oxDep = REFoxDep, oxRec = REFoxDep, 
					OERCentre = REFOERCenter):
	retVals = []
	for dose, doseRate, o2 in exposures:
		OER = cumulativeOER(dose, doseRate, o2, oxDep, oxRec, OERCentre)
		scaledDose = dose*OER
		surv = np.exp(-scaledDose*alpha-beta*scaledDose*scaledDose)
		retVals.append(surv)
	return retVals
