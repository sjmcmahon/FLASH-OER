
# FLASH oxygen depletion model

This is a simple implementation of a model of oxygen depletion by FLASH radiotherapy. This model simulates O2 concentration-dependent or -independent depletion of oxygen by radiation exposure, and the impact this has on OER. It is assumed the OER for a given dose of radiation depends on the O2 level present when that radiation is delivered, and that these effects can be integrated across an entire radiation exposure.

A full description of this model can be found in its original publication [Petersson 2020]. 

The core model is implemented in `flashOERModel.py`, providing three main functions: 

- oxygenCurve, which provides the level of oxygen during radiation exposures;
- cumulativeOER, which calculates the average OER of an exposure incorporating oxygen depletion by radiation; and
- predictSurvival, which calculates the predicted survival level for a given set of LQ parameters and radiation exposure.

In addition to standard radiobiological parameters alpha, beta and the mid-point oxygen concentration of the OER curve (`REFOERCenter`) there are two other parameters, the oxygen depletion rate and the oxygen recovery rate. Reference parameters used in the model publication are stored in the main model file (`REFoxDep` and `REFoxRec` respectively), or alternative values can be passed into functions to explore their impact.

An illustrative file showing how these different models can be used is also provided, showing the use of each of these functions (`flashAnalysis.py`).

## Requirements

This code is written in python3, and requires the following libraries:

- numpy
- matplotlib (analysis file only)

## Contacts

For questions/comments/bug reports, please contact stephen.mcmahon (at) qub.ac.uk

## References

[Petersson2020] Petersson, K., Adrian, G., Butterworth, K., & McMahon, S. J. (2020). A quantitative analysis of the role of oxygen tension in FLASH radiotherapy. International Journal of Radiation Oncology* Biology* Physics, 107(3), 539-547. https://doi.org/10.1016/j.ijrobp.2020.02.634
