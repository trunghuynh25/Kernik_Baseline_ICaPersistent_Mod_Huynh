% Main File to generate Baseline Model Figures 10-11
close all; clear; clc

load ICs_baseline
load baseline_parameter_inputs

% ========================================================================
%  START: PARAMETER AND TITLE CONFIGURATION BLOCK
% ========================================================================

% --- 1. Define which model to run ---
% Set ONLY ONE of these to `true` to select your experiment
run_baseline = false;
run_model_1 = false;  % <--- Run this model
run_model_2 = false; % <--- Ignore this model for now
run_model_4 = true; % <--- Ignore this model for now

% --- 2. Define Parameter Indices ---
g_Kr_index = 7;
g_Na_index = 40;
p_CaL_index = 30;

% --- 3. Apply Modifications Based on Selection ---
modified_params = baseline_parameter_inputs; % Start with the baseline

if run_model_1
    plot_title = 'Model 1: IKr Defect Only';
    modified_params(g_Kr_index) = modified_params(g_Kr_index) * (1 - 0.725);
    
elseif run_model_2
    plot_title = 'Model 2: IKr + INa,L';
    modified_params(g_Kr_index) = modified_params(g_Kr_index) * (1 - 0.725);
    modified_params(g_Na_index) = modified_params(g_Na_index) * 2.31;

elseif run_model_4
    plot_title = 'Model 4: All Three Defects';
    modified_params(g_Kr_index) = modified_params(g_Kr_index) * (1 - 0.725);
    modified_params(g_Na_index) = modified_params(g_Na_index) * 2.31;
    modified_params(p_CaL_index) = modified_params(p_CaL_index) * 2.29;
    
else % This runs if all flags above are 'false'
    plot_title = 'Baseline Model (Isogenic Control)';
end

% Overwrite the original variable with our modified one for the simulation
baseline_parameter_inputs = modified_params;
% ========================================================================
%  END: CONFIGURATION BLOCK
% ========================================================================


%% Run iPSC_function
options = odeset('MaxStep',1,'InitialStep',2e-2);
run_time=3e3;
[Time, values] = ode15s(@ipsc_function,[0, run_time],Y_init, options, baseline_parameter_inputs);
Cai=values(:,3);
Vm=values(:,1);

%% Calculate select current traces:
INaCa = zeros(size(Time));
IpCa = zeros(size(Time));
Iup = zeros(size(Time));
for i= 1:size(values,1)
    [~, update_step_i] =  ipsc_function(Time(i), values(i,:),  baseline_parameter_inputs);    
    INaCa(i) = update_step_i(8);
    IpCa(i) = update_step_i(9);
    Iup(i) = update_step_i(14);
end

% ========================================================================
%% Figure 10A & 10C: Calcium Flux analysis and Calcium Transient Trace
ca_analysis( Time, Iup, INaCa, IpCa, Cai, plot_title )
% ========================================================================


%% Figure 11A: action potential trace for baseline model 
figure,set(gcf,'color','w')
plot(Time, Vm,'Color', [.8 0 .18]);
set(gca,'box','off','tickdir','out')
ylabel('Voltage (mV)');
xlabel('Time (ms)')
title(plot_title) % Use the dynamic title variable here

