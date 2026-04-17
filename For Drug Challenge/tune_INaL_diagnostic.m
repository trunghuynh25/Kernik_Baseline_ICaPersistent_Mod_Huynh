% =========================================================================
% find_INaL_range.m
% =========================================================================
% This is a DIAGNOSTIC script.
% It loops through a set of guesses for C_persist_NaL and
% prints a table of the measured I_Na,L / I_Na,P ratio for each.
%
% GOAL: Use this table to find the two C_persist values that
% "bracket" your target ratio of 0.0060.
%
% =========================================================================

clear; clc; close all;

disp('Loading baseline files...');
load ICs_baseline.mat;
load baseline_parameter_inputs.mat;

% --- Define a wide range of guesses to test ---
% We'll test from 0.1% up to 20%
guesses_to_try = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1];

% --- Store results here ---
results_table = zeros(length(guesses_to_try), 2);

disp('--- Starting Manual Search for INa,L Range ---');
fprintf('Target Ratio is: 0.0026 (0.26%%)\n\n');
fprintf('%-20s | %-20s\n', 'C_persist_Guess', 'Measured_Ratio');
fprintf('------------------------------------------\n');

for k = 1:length(guesses_to_try)
    
    C_persist_guess = guesses_to_try(k);
    
    % Call the 'get_ratio' function (nested below)
    measured_ratio = get_measured_ratio(C_persist_guess, baseline_parameter_inputs, Y_init);
    
    % Store and print the result
    results_table(k, 1) = C_persist_guess;
    results_table(k, 2) = measured_ratio;
    
    fprintf('%-20.6f | %-20.6f (%.4f%%)\n', ...
            C_persist_guess, measured_ratio, measured_ratio*100);
end

disp('--- Manual Search Complete ---');
disp('Look at the table above.');
disp('Find the two C_persist_Guess values where the Measured_Ratio crosses 0.0026.');
disp('Example: If you see:');
disp('   0.020000 | 0.002100');
disp('   0.030000 | 0.002800');
disp('Then your new search range for tune_INaL.m is [0.02, 0.03]');


% =========================================================================
% NESTED FUNCTION TO RUN SIM
% =========================================================================
function measured_ratio = get_measured_ratio(C_persist_guess, baseline_parameter_inputs, Y_init)
    % This is the same logic as the 'calculate_error' function
    % It runs one VC sim and returns the measured ratio.

    sim_params = baseline_parameter_inputs;
    sim_params(87) = C_persist_guess; % Set INaL
    sim_params(86) = 0;               % Set ICaL to 0
    sim_params(83) = 0;               % stim_flag OFF
    sim_params(85) = 1;               % voltageclamp ON
    
    t_span = [0 1000];
    options = odeset('MaxStep', 0.2);
    
    [T, V] = ode15s(@ipsc_function, t_span, Y_init, options, sim_params);

    i_Na = zeros(length(T), 1);
    for i = 1:length(T)
        [~, all_currents] = ipsc_function(T(i), V(i,:), sim_params);
        i_Na(i) = all_currents(7); % i_Na is 7th
    end

    t_step_indices = find(T >= 500 & T <= 1000);
    if isempty(t_step_indices)
        measured_ratio = -1; % Error flag
        return;
    end
    
    i_Na_trace = i_Na(t_step_indices);
    I_Na_P = min(i_Na_trace);
    I_Na_L = mean(i_Na_trace(end-100:end));
    
    measured_ratio = abs(I_Na_L / I_Na_P);
end