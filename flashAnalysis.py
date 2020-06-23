# ############################################################################
# 
# This software is made freely available in accordance with the simplifed BSD
# license:
# 
# Copyright (c) <2020>, <Stephen McMahon>
# All rights reserved
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation 
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND ANY 
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contacts: Stephen McMahon,	stephen.mcmahon@qub.ac.uk
# 
# ############################################################################
# 
# Basic analytic function plotting different outputs of FLASH oxygen depletion
# model. 
#
# A Quantitative Analysis of the Role of Oxygen Tension in FLASH Radiation Therapy
# Petersson et al, IJROBP, 107(3), 539-547
# https://doi.org/10.1016/j.ijrobp.2020.02.63
# 
# Please cite this paper if you are building on this tool.
# 
# ############################################################################

from flashOERModel import *
import matplotlib.pyplot as plt 

# Plot of oxygen depletion
def oxygenDepletion(baseO2=0.05):
	# Build array of dose-rates, exposure times and O2 concentrations
	doseRates = [pow(10,x) for x in np.arange(-1,2.1,0.5)]
	times = [pow(10,x) for x in np.arange(-7,3,0.025)]
	oxCurves = []
	for dDot in doseRates: 
		oxygen = oxygenCurve(times, dDot, baseO2, REFoxDep, REFoxRec)
		oxCurves.append(oxygen)

	fig = plt.figure(figsize=(5,4))
	ax = plt.gca()
	for n,row in enumerate(oxCurves):
		roundedRate = round(doseRates[n],-(int(np.floor(np.log10(abs(doseRates[n])))-1)))
		if roundedRate%1<0.0001: roundedRate = int(roundedRate)
		lab = str(roundedRate)+" Gy/s"
		ax.plot(times, row, label=lab)
	ax.set_xscale('log')
	ax.set_xlim(0.001,1000)
	ax.set_ylim(0,0.05)
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Oxygen concentration (atmospheric fraction)")
	ax.legend()

	plt.tight_layout()
	plt.show()

# Plot of OER values as a function of dose-rate
def meanOERCurve(dose=20):
	doseRates = [pow(10,x) for x in np.arange(-1,2.1,1)]
	o2Levels = [0.25/pow(10,x) for x in np.arange(0,3.5,0.005)]
	print(o2Levels)
	print(doseRates)
	oerOut = []
	for dDot in doseRates:
		oers = [3*cumulativeOER(dose, dDot, baseO2, REFoxDep, REFoxRec, REFOERCenter) for baseO2 in o2Levels]
		oerOut.append(oers)

	fig = plt.figure(figsize=(5,4))
	ax = plt.gca()
	for n,row in enumerate(oerOut):
		lab = str(round(doseRates[n],2))+" Gy/s"
		ax.plot(o2Levels,row,label=lab)

	ax.set_xscale('log')
	ax.set_xlim(0.0002,0.2)
	ax.set_ylim(1,3)
	ax.set_xlabel(r'O$_2$ Concentration (atmospheric fraction)')
	ax.set_ylabel(r'$\overline{OER}$')
	ax.legend()

	plt.tight_layout()
	plt.show()

# Plot of survivals for oxic and hypoxic cells, at normal and FLASH dose rates
def survivalPlot():
	doses = [x for x in np.arange(0,25,0.1)]
	o2Levels = [0.2, 0.016]
	o2Labels = ["Normoxic","Hypoxic"]
	doseRates = [0.2333, 600]
	doseLabels = ["Conventional","FLASH"]

	alpha = 0.12
	beta = 0.027

	fig = plt.figure(figsize=(5,4))
	ax = plt.gca()
	for n,o2 in enumerate(o2Levels):
		for m,dr in enumerate(doseRates):
			exposures = [[d,dr,o2] for d in doses]
			survs = predictSurvival(exposures, alpha, beta)
			label = o2Labels[n]+", "+doseLabels[n]
			ax.plot(doses, survs, label=label)

	ax.set_yscale('log')
	ax.set_xlabel("Dose (Gy)")
	ax.set_ylabel("Surviving Fraction")
	ax.set_ylim([8E-5,1])
	ax.set_xlim([0,24])
	ax.legend()

	plt.tight_layout()
	plt.show()

# Plot of survivals as a function of dose rates
def doseRateSurvivalPlot():
	doses = [x for x in np.arange(0,25,0.1)]
	o2Level = 0.016
	o2Label = str(o2Level*100)+"% O2"
	doseRates = [0.1, 1, 10, 100, 1000, 1E10]
	doseLabels = [str(dr)+" Gy/s" for dr in doseRates]

	alpha = 0.12
	beta = 0.027

	fig = plt.figure(figsize=(5,4))
	ax = plt.gca()
	for n,dr in enumerate(doseRates):
		exposures = [[d,dr,o2Level] for d in doses]
		survs = predictSurvival(exposures, alpha, beta)
		label = o2Label+", "+doseLabels[n]
		ax.plot(doses, survs, label=label)

	ax.set_yscale('log')
	ax.set_xlabel("Dose (Gy)")
	ax.set_ylabel("Surviving Fraction")
	ax.set_ylim([8E-5,1])
	ax.set_xlim([0,24])
	ax.legend()

	plt.tight_layout()
	plt.show()

if __name__ == "__main__":
	oxygenDepletion()
	meanOERCurve()
	survivalPlot()
	doseRateSurvivalPlot()