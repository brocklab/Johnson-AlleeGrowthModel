% This script is a master script to run functions that give the necessary
% outputs from each model structure we are investigating. We want to both
% run simulations of the Gillespie algorithm for a given N and parameter
% set. We want to extract the summary statistics from each simulated
% data set, as well as the expected moments derived from the CME. Then, we
% can run all of the functions and check to see how they vary depending on
% both the parameters and the assumptions about the Allee effect on either
% the birth or death for the three deterministic model structures:
% 1. dN/dt = gN
% 2. dN/dt = g(N-A) * Allee on birth, death, or both
% 3. dN/dt = gN(1-(A+tau)/(N+tau)) * Allee on birth, death, or both

% Leads to 7 total models to look at.
% We want to probably make comprehensive sheet of mean and variances for
% each of these models, but in the fitting we need up to the 4th moment and
% 4th order variance.

% Function structure:
% [ Nsamp, tsamp,
