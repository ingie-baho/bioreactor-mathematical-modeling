# Cyanobacteria Mathematical Growth Modeling

This repository contains all the mathematical modeling and simulation code used to model the growth of cyanobacteria under various environmental conditions. The models aim to predict biomass accumulation and carbon uptake based on nutrient availability and light conditions, supporting carbon fixation analysis in a bioreactor environment.

## Repository Structure

### `main_model.py`
Implements a modified Monod model to simulate biomass growth as a function of:
- Nitrogen concentration  
- Phosphorus concentration  
- Light intensity (with absorption modeled using empirically derived coefficients)  
- Available carbon content  

Optical density (OD) data from bioreactor experiments is also included to inform parameter estimation and analysis of the cyanobacteria absorption coefficient, which plays a key role in the light intensity equation.

### `sensitivity_analysis/`
Contains scripts to explore how changes in substrate concentrations affect growth rate. This includes one-variable-at-a-time sensitivity analyses to identify key drivers of the system.

### `carbon_capture/`
Models carbon consumption using a mass balance approach. Equations and simulations in this folder estimate how much carbon is fixed over time as the cyanobacteria grow.

## Purpose

These models were developed to better understand and predict cyanobacterial growth behavior in controlled environments, enabling more efficient design of photobioreactor systems for carbon sequestration.

## Contributions

This code was written and developed by **Ingie Baho**, with support from the MIT Bioinstrumentation Lab.  
The work was carried out under the supervision of **Professor Ian Hunter**, who provided guidance and lab resources throughout the project.

##
