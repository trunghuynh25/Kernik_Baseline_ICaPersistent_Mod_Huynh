% =========================================================================
% tune_INaL.m (v3 - Corrected for ALL scoping errors)
% =========================================================================
% This script tunes the 'C_persist_NaL' parameter (index 87) to find
% the value that accurately produces a target I_Na,L / I_Na,P ratio.
%
% =========================================================================

clear; clc; close all;

disp('Loading baseline files...');
load ICs_baseline.mat;
load baseline_parameter_inputs.mat; % <-- Variable is loaded here

% --- Define Your Target ---
target_ratio = 0.0026; % <-- Variable is defined here, aiming for 0.0026 or 0.0060

% --- Solver Settings ---
% This is the correct bracket you found in the last step
initial_search_range = [0.9, 1]; 

disp('--- Starting INa,L Tuning Process ---');
fprintf('Target Ratio: %.6f (%.4f%%)\n', target_ratio, target_ratio*100);
fprintf('Search Range: [%.2f, %.2f]\n', initial_search_range(1), initial_search_range(2));

% --- Call the Solver ---
options = optimset('TolX',1e-8); % Set a high precision

% ========================================================================
%  *** THIS IS THE FIX ***
%
% The anonymous function handle 'fHandle' now captures and passes
% ALL required variables:
% 1. baseline_parameter_inputs
% 2. Y_init
% 3. target_ratio
%
fHandle = @(C_persist) calculate_error(C_persist, baseline_parameter_inputs, Y_init, target_ratio);
% ========================================================================

[C_persist_final, final_error] = fzero(fHandle, initial_search_range, options);

disp('--- Tuning Complete ---');
fprintf('Target Ratio: %.6f (%.4f%%)\n', target_ratio, target_ratio*100);
fprintf('Optimized C_persist_NaL value: %.8f\n', C_persist_final);
fprintf('Error at final value (should be ~0): %e\n', final_error);

disp('*******************************************************************');
disp('ACTION REQUIRED:');
fprintf('Your new "baseline_C_persist_NaL" is %.8f\n', C_persist_final);
disp('Update your `baseline_parameter_inputs.mat` file to use this value.');
disp('*******************************************************************');


% =========================================================================
% NESTED ERROR FUNCTION
% =========================================================================
%
% *** THIS IS THE OTHER PART OF THE FIX ***
% The function must now be defined to *accept* target_ratio.
%
function error = calculate_error(C_persist_guess, baseline_parameter_inputs, Y_init, target_ratio)
    % This function is called by 'fHandle' and now receives all the
    % variables it was missing.

    % 1. Create a local copy of parameters for this simulation
    sim_params = baseline_parameter_inputs;

    % 2. Set the C_persist_NaL (param 87) for this guess
    sim_params(87) = C_persist_guess;
    sim_params(86) = 0; % Set persistent CaL (param 86) to 0

    % 3. Set up the VOLTAGE CLAMP protocol
    sim_params(83) = 0; % stim_flag OFF
    sim_params(85) = 1; % voltageclamp ON
    
    % 4. Run the simulation
    t_span = [0 1000];
    options = odeset('MaxStep', 0.2);
    
    [T, V] = ode15s(@ipsc_function, t_span, Y_init, options, sim_params);

    % 5. Get the sodium current trace
    i_Na = zeros(length(T), 1);
    for i = 1:length(T)
        [~, all_currents] = ipsc_function(T(i), V(i,:), sim_params);
        i_Na(i) = all_currents(7); % i_Na is the 7th current
    end

    % 6. Analyze the currents
    t_step_indices = find(T >= 500 & T <= 1000);
    
    if isempty(t_step_indices)
        error = 1e9; 
        return;
    end
    
    i_Na_trace = i_Na(t_step_indices);
    I_Na_P = min(i_Na_trace);
    I_Na_L = mean(i_Na_trace(end-100:end)); 
    
    % 7. Calculate the measured ratio and the error
    measured_ratio = abs(I_Na_L / I_Na_P);
    
    % This line will now work, as target_ratio is an argument
    error = measured_ratio - target_ratio; 

    fprintf('Guess: C_persist = %.8f | Measured Ratio = %.6f (%.4f%%) | Error = %e\n', ...
            C_persist_guess, measured_ratio, measured_ratio*100, error);
end