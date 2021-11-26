# OnlineCausal

# Title: Online Causal Inference with Application to NearReal-Time Post-Market Vaccine Safety Surveillance
# Version: 1.0
# Date: 2021-11-26
# Author: Lan Luo
# Maintainer: Lan Luo <luo-lan@uiowa.edu>

# Description: Online estimation and inference of average treatment effect. 
* [datagenerator.R] generating the simulated data

* [renew_ATE.R] core functions for evaluating bias and coverage probability)
* [renew_ATE_seq.R] core functions for evaluating sequential testing
* [offline_ATE.R] main function for offline estimation and inference
* [increATE.cpp] built-in functions for online estimation and inference (with an interaction term)
* [increATE_no_interaction.cpp] built-in functions for online estimation and inference (without interaction term)

* [run_causal.R] execute file for simulations (evaluating bias and coverage probability)
* [run_causal_seq.R] execute file for simulations (evaluating performance in sequential testing)
