# Functional-Connectivity-Analysis-Toolbox

@Authors: Ioannis Ziogas & Charalampos Lamprou Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh

The toolbox and code presented here is part of and was used for the implementation of the MSc thesis project that Charalampos Lamprou and Ioannis Ziogas conducted in order to receive their MSc in Electrical & Computer Engineering from the Faculty of Engineering, Aristotle University of Thessaloniki (ECE AUTh). This work was implemented with the support of the Signal Processing and Biomedical Technology Unit research laboratory of AUTh, under the supervision of Prof. Leontios Hadjileontiadis and Prof. Herbert Jelinek.

The manuscript corresponding to the work presented here can be found at https://ikee.lib.auth.gr/record/338217/?ln=en in Greek, whereas an abstract is also available in English.

The authors would like to acknowledge the use of the EEGLAB package for the implementation of this work, and would like to thank the authors of EEGLAB for their contribution to this work through EEGLAB functions and packages: Delorme A & Makeig S (2004) EEGLAB: an open-source toolbox for analysis of single-trial EEG dynamics, Journal of Neuroscience Methods 134:9-21.

# DESCRIPTION

The present toolbox constitutes a set of functions and processes designed to carry out a brain connectivity pipeline, from the very initial stage of the raw EEG data, to the final stage of a dataset with variables that represent the brain activity in terms of brain functional connectivity. 

Included are functions to implement: EEG preprocessing, Decomposition of brain activity through the Swarm Decomposition algorithm, Feature extraction.

# EEG Preprocessing

Several EEG preprocessing steps are available, including: Filtering, Line-noise removal, Re-referencing, Artifact and Channel Removal, Artifact removal through wavelet ICA algorithm, Time-Windowing of EEG. For the preprocessing step, EEGLAB functions are used. 

# Brain Functional Connectivity Feature Extraction

5 different functional connectivity measures are available: Spectral Coherence, Phase Locking Value, Higher-Order-Spectra (Antisymmetric Cross-Bispectrum & Cross-Bicoherence), Phase-Amplitude Coupling through Modulation Index & Phase-Amplitude Coupling through a Time-Frequency algorithm.

Moreover, 4 different brain regional scenaria are available, for analyzing brain activity at a regional level, rather than an electrode/sensor level. 

# Demonstration 

A demonstration of the operation of this toolbox can be found in the script: example_script.m
