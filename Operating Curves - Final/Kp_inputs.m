%% Run this file to start the simulation
tic
clear
close all

%% Simulation Inputs
input.psi = 0.0;                                                                                   % Target liquid per well (uL)
input.min_Vr_range = 0.1;                                                                           % Minimum simulated resin volume (uL)
input.max_Vr_range = 50;                                                                            % Maximum simulated resin volume (uL)
input.intra_ep = 0.6;

%% Feasibility Parameter Inputs
input.sep_fact_thresh = 0.9;                                                                        % Refernce contour for separation factor threshold
input.Kp_thresh = 100;                                                                               % Reference contours for retention time
input.conc_thresh = 0.02;                                                                           % Reference contour for equilibrium concentration

%% Reference Isotherm Parameters
ref_isotherm = import_iso();                                                                        % Import reference isotherm data
ref_isotherm.start_salt = 1;                                                                        % Lowest salt concentration for isotherm fitting (mM)
ref_isotherm.end_salt = 1000;                                                                       % Highest salt concentration for isotherm fitting (mM)
ref_isotherm.n_salt = 200;                                                                           % Number of salt concentration grid points
ref_isotherm.start_c = 10^-10;                                                                      % Starting protein concentration for isotherm fitting (mg/mL)
ref_isotherm.end_c = 30;                                                                            % Starting protein concentration for isotherm fitting (mg/mL)
ref_isotherm.n_c = 100;                                                                             % Number of salt concentration grid points
ref_isotherm.Rsqlim = 0.9;                                                                          % Minimumn correlation coefficient for Langmuir fitting

h = waitbar(0, {'Initializing'; '0%'});
if isempty(dir('*mat')) == 1
%% Calculate Equilibrium Concentrations, Separation Factors, and Column Retention Times
[Ke, qm, Rsq, h] = iso(ref_isotherm, input);

%% Calculate Reference Isotherm Frontier Points
frnt_pts = lead_edge(qm, Ke, Rsq, ref_isotherm);
else 
     file = dir('*mat');
     load(file.name);
end

%% Calculate point of intersection between Kp_max and Reference Isotherm
[qm_star, Kl_star] = intersect_point(frnt_pts, Ke, qm, input);

%% Calculate Top Operating Curve
[b, Load_top] = top(input, qm_star, Ke, qm, frnt_pts, ref_isotherm, h);

%% Calculate Bottom Operating Curve
[Load_bottom] = bottom(input, b);

%% Plot Operating Curves
plotting(input, b, Load_top, Load_bottom);

toc

