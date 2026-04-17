% Main File to generate Baseline Model Figures 10-11
close all; clear; clc

load ICs_baseline
% <<< MODIFICATION: Load parameters into a 'clean' variable
load baseline_parameter_inputs
original_baseline_params = baseline_parameter_inputs;
clear baseline_parameter_inputs; % Clear to avoid confusion

% ========================================================================
%  START: EXPERIMENTAL CONFIGURATION BLOCK
% ========================================================================

% <<< MODIFICATION: Define the list of models to run
models_to_run_list = 0:4; % Run all models from 0 (Baseline) to 4
model_names = {'Baseline', 'Model 1 (IKr)', 'Model 2 (IKr+INa,L)', 'Model 3 (IKr+ICa,P)', 'Model 4 (All)'};

% --- Parameter Indices (moved from config block) ---
g_Kr_index = 2;   % Conductance of IKr
g_Na_index = 7;  % Conductance of INa (proxy for late current)
p_CaL_index = 5; % Permeability of ICaL (proxy for persistent current)
g_K1_index = 1;   % Conductance of IK1 (Added for stabilization)
C_persist_CaL_index = 86; % Our new persistent ICaL parameter
C_persist_NaL_index = 87; % Persistent INaL parameter (baseline = 0.0026 or 0.26%)


% <<< MODIFICATION: Pre-allocate result structures
% Use cell arrays to accommodate potentially dissimilar structs
all_ca_results = cell(length(models_to_run_list), 1);
all_ap_results = cell(length(models_to_run_list), 1);

% ========================================================================
%  START: MAIN SIMULATION LOOP
% ========================================================================
fprintf('Starting simulations for all models...\n');

for model_idx = 1:length(models_to_run_list)
    
    model_to_run = models_to_run_list(model_idx);
    
    % --- Apply Modifications Based on Selection ---
    % <<< MODIFICATION: Reset to baseline parameters for *each* iteration
    modified_params = original_baseline_params; 
    modified_params(83) = 1; % Off = 0 | ON = 1 |% This is the master "on/off" switch for pacing.

    % This is what we will scale by 2.29-fold OR 2.31-fold
    baseline_C_persist_CaL = modified_params(C_persist_CaL_index);
    baseline_C_persist_NaL = modified_params(C_persist_NaL_index);

        modified_params(g_K1_index) = modified_params(g_K1_index) * 5; % Kernik-Clancy Mature iPSC-CM
        modified_params(g_Na_index) = modified_params(g_Na_index) * 1.45; % Kernik-Clancy Mature iPSC-CM

    if model_to_run == 1
        plot_title = 'Model 1: IKr Defect Only';
        % Apply 72.5% IKr reduction
        modified_params(g_Kr_index) = modified_params(g_Kr_index) * (1 - 0.725);
    elseif model_to_run == 2
        plot_title = 'Model 2: IKr + INa,L';
        % Apply IKr reduction
        modified_params(g_Kr_index) = modified_params(g_Kr_index) * (1 - 0.725);
        % Apply 2.31-fold INa,L increase
        modified_params(C_persist_NaL_index) = baseline_C_persist_NaL * 2.31;
        
    elseif model_to_run == 3
        plot_title = 'Model 3: IKr + ICa,P';
        % Apply IKr reduction
        modified_params(g_Kr_index) = modified_params(g_Kr_index) * (1 - 0.725);
        % Apply 2.29-fold ICa,P increase (from your previous query)
        modified_params(C_persist_CaL_index) = baseline_C_persist_CaL * 2.29;

    elseif model_to_run == 4
        plot_title = 'Model 4: All Three Defects';
        % Apply all three modifications   
        modified_params(g_Kr_index) = modified_params(g_Kr_index) * (1 - 0.725);
        modified_params(C_persist_NaL_index) = baseline_C_persist_NaL * 2.31;
        modified_params(C_persist_CaL_index) = baseline_C_persist_CaL * 2.29;

    else % Default to baseline if model_to_run is 0 or any other number
        plot_title = 'Baseline Model (Isogenic Control)';
    end

    % <<< MODIFICATION: Use a new variable for the simulation
    current_sim_params = modified_params;
    
    fprintf('--- Running: %s ---\n', plot_title);

    % ========================================================================
    %  RUN SIMULATION & ANALYSIS (INSIDE LOOP)
    % ========================================================================
    
    % <<< MODIFICATION: Removed the erroneous 'Run iPSC_function' line
    options = odeset('MaxStep',1,'InitialStep',2e-2);
    run_time=3e3;
    
    % <<< MODIFICATION: Pass the 'current_sim_params' to the solver
    [Time, values] = ode15s(@ipsc_function,[0, run_time],Y_init, options, current_sim_params);
    Cai=values(:,3);
    Vm=values(:,1);
    
    % --- Calculate select current traces ---
    INaCa = zeros(size(Time));
    IpCa = zeros(size(Time));
    Iup = zeros(size(Time));
    for i= 1:size(values,1)
        % <<< MODIFICATION: Pass 'current_sim_params' here too
        [~, update_step_i] =  ipsc_function(Time(i), values(i,:),  current_sim_params);
        INaCa(i) = update_step_i(8);
        IpCa(i) = update_step_i(9);
        Iup(i) = update_step_i(14);
    end
    
    % --- Run Calcium Analysis ---
    % <<< MODIFICATION: Pass plot_title, but we will suppress the plots
    [ca_results, validation_data] = ca_analysis_v6( Time, Iup, INaCa, IpCa, Cai, plot_title );
    
    % --- Run AP Morphology Analysis ---
    voltage = Vm;
    time = Time;
    ap_results = CalculateAPMorphology(time, voltage); % Call the local function
    
    % <<< MODIFICATION: STORE the results for this model
    all_ca_results{model_idx} = ca_results;
    all_ap_results{model_idx} = ap_results;

    % ========================================================================
    %  PLOTTING (INSIDE LOOP) - NOW COMMENTED OUT
    % ========================================================================
    % Running the loop will create 5x3=15 figures, which is messy.
    % You can re-enable these if you want all plots, or move this
    % logic outside the loop to create overlay plots.
    % ========================================================================
    
    % % Figure 11A: action potential trace for baseline model
    % figure,set(gcf,'color','w')
    % plot(Time, Vm,'Color', [.8 0 .18]);
    % set(gca,'box','off','tickdir','out')
    % ylabel('Voltage (mV)');
    % xlabel('Time (ms)')
    % title(plot_title)
    % 
    % % Figure 10A & 10C: TTP VALIDATION PLOT
    % figure, set(gcf,'color','w');
    % hold on;
    % % ... (rest of TTP validation plot code) ...
    % hold off;
    % title({plot_title; 'Time to Peak (TTP) Validation'});
    % 
    % % AP Morphology Plot
    % DisplayAPMorphology(time, voltage, ap_results);

end % <<< MODIFICATION: END OF MAIN SIMULATION LOOP

fprintf('... All simulations complete. Generating summary tables.\n');

% ========================================================================
%  START: CALCULATE AND DISPLAY DATA TABLES (OUTSIDE LOOP)
% ========================================================================
% <<< MODIFICATION: This entire block is new.
% It iterates through the *stored results* to build comparative tables.

% --- Helper: Generate the header row ---
fprintf('\n\n');
fprintf('%-30s |', 'Parameter');
for i = 1:length(model_names)
    fprintf(' %-18s |', model_names{i});
end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 32 + 21 * length(model_names)));

% --- Display Table 1: Calcium Transient Morphology ---
fprintf('\n## Table 1: Calcium Transient Morphology ##\n');
fprintf('%-30s |', 'Peak [Ca²⁺] (nM)');
for i = 1:length(all_ca_results), fprintf(' %-18.2f |', all_ca_results{i}.peak_ca_nM); end
fprintf('\n');

fprintf('%-30s |', 'Diastolic [Ca²⁺] (nM)');
for i = 1:length(all_ca_results), fprintf(' %-18.2f |', all_ca_results{i}.diastolic_ca_nM); end
fprintf('\n');

fprintf('%-30s |', 'Time to Peak (ms)');
for i = 1:length(all_ca_results), fprintf(' %-18.2f |', all_ca_results{i}.ttp_ms); end
fprintf('\n');

fprintf('%-30s |', 'Tau Decay Fast (ms)');
for i = 1:length(all_ca_results), fprintf(' %-18.2f |', all_ca_results{i}.tau_fast_ms); end
fprintf('\n');

fprintf('%-30s |', 'Tau Decay Slow (ms)');
for i = 1:length(all_ca_results), fprintf(' %-18.2f |', all_ca_results{i}.tau_slow_ms); end
fprintf('\n');

fprintf('%-30s |', 'R^2 Tau Decay');
for i = 1:length(all_ca_results), fprintf(' %-18.2f |', all_ca_results{i}.tau_decay_r_squared); end
fprintf('\n');


% --- Display Table 2: % Contribution of Calcium Flux ---
fprintf('\n## Table 2: %% Contribution of Calcium Flux ##\n');
fprintf('%-30s |', 'SERCA');
for i = 1:length(all_ca_results), fprintf(' %-18.1f%% |', all_ca_results{i}.pct_serca); end
fprintf('\n');

fprintf('%-30s |', 'NCX');
for i = 1:length(all_ca_results), fprintf(' %-18.1f%% |', all_ca_results{i}.pct_ncx); end
fprintf('\n');

fprintf('%-30s |', 'Non-NCX (I_pCa)');
for i = 1:length(all_ca_results), fprintf(' %-18.1f%% |', all_ca_results{i}.pct_ipca); end
fprintf('\n');


% --- Display Table 3: Action Potential Morphology ---
fprintf('\n## Table 3: Action Potential Morphology ##\n');
fprintf('%-30s |', 'Resting Membrane Potential (mV)');
for i = 1:length(all_ap_results), fprintf(' %-18.2f |', all_ap_results{i}.RMP); end
fprintf('\n');

fprintf('%-30s |', 'Amplitude (mV)');
for i = 1:length(all_ap_results), fprintf(' %-18.2f |', all_ap_results{i}.APA); end
fprintf('\n');

fprintf('%-30s |', 'Vmax (mV)');
for i = 1:length(all_ap_results), fprintf(' %-18.2f |', all_ap_results{i}.V_max); end
fprintf('\n');

fprintf('%-30s |', 'dV/dt max (V/s)');
for i = 1:length(all_ap_results), fprintf(' %-18.2f |', all_ap_results{i}.dVdt_max / 1000); end
fprintf('\n');

fprintf('%-30s |', 'APD30 (ms)');
for i = 1:length(all_ap_results), fprintf(' %-18.2f |', all_ap_results{i}.APD30); end
fprintf('\n');

fprintf('%-30s |', 'APD50 (ms)');
for i = 1:length(all_ap_results), fprintf(' %-18.2f |', all_ap_results{i}.APD50); end
fprintf('\n');

fprintf('%-30s |', 'APD90 (ms)');
for i = 1:length(all_ap_results), fprintf(' %-18.2f |', all_ap_results{i}.APD90); end
fprintf('\n');


% ========================================================================
%  LOCAL FUNCTIONS (Cleaned and properly indented)
% ========================================================================

%%
function results = CalculateAPMorphology(time, voltage)
    % Calculates action potential morphology parameters using linear interpolation.
    % CORRECTED VERSION: Assumes the 'time' vector is in microseconds (µs).

    % --- NEW: Define a search window for the AP peak (in microseconds) ---
    % 500 ms is a safe window, which equals 500,000 µs.
    peak_search_window_us = 500000; 

    % --- Find Max Upstroke Velocity (dV/dt) to anchor the analysis ---
    % dVdt will be in mV/µs
    dVdt = diff(voltage) ./ diff(time);
    [dVdt_max, idx_dVdt_max] = max(dVdt);
    t_depol = time(idx_dVdt_max); % This is time of dVdt_max (reference for APD)

    % --- MODIFIED: Search for the AP peak ONLY after the upstroke ---
    % 1. Determine the average time step (sampling interval) in µs
    sampling_interval_us = mean(diff(time));
    
    % 2. Calculate the window size in array indices
    window_indices = round(peak_search_window_us / sampling_interval_us);
    search_start_idx = idx_dVdt_max;
    search_end_idx = min(search_start_idx + window_indices, length(voltage));
    
    % 3. Find the peak voltage and its index within the defined window
    [V_max, local_idx_max] = max(voltage(search_start_idx:search_end_idx));
    
    % 4. Convert the local index back to the global index of the original trace
    idx_max = search_start_idx + local_idx_max - 1;
    t_peak = time(idx_max);

    % --- Standard Calculations ---
    RMP = mean(voltage(1:10)); 
    APA = V_max - RMP;

    % ***************************************************************
    % --- NEW: Find the true AP Start time (t_AP_start) using RMP + 1mV ---
    V_threshold = RMP + 6; % 6mV above RMP 

    % Search before the max upstroke point
    depol_phase_voltage = voltage(1:idx_dVdt_max);
    depol_phase_time = time(1:idx_dVdt_max);

    % Find the first index where voltage crosses the threshold
    idx1_start_cross = find(depol_phase_voltage > V_threshold, 1, 'first');

    if ~isempty(idx1_start_cross) && idx1_start_cross > 1
        % Interpolate the exact time of the threshold crossing
        t1 = depol_phase_time(idx1_start_cross - 1);
        v1 = depol_phase_voltage(idx1_start_cross - 1);
        t2 = depol_phase_time(idx1_start_cross);
        v2 = depol_phase_voltage(idx1_start_cross);
        t_AP_start = interpolate_time(t1, v1, t2, v2, V_threshold);
        V_AP_start = V_threshold;
    else
        % Fallback: Use the first data point if crossing is not found
        t_AP_start = time(1);
        V_AP_start = RMP;
    end
    % ***************************************************************

    % --- APD Calculations ---
    repol_phase_time = time(idx_max:end);
    repol_phase_voltage = voltage(idx_max:end);
    
    apd_levels = [30, 50, 90]; 
    
    results.t_rep = zeros(length(apd_levels), 1);
    results.v_rep = zeros(length(apd_levels), 1);
    
    for i = 1:length(apd_levels)
        level = apd_levels(i);
        V_rep_target = V_max - (level/100) * APA;
        idx2_repol = find(repol_phase_voltage < V_rep_target, 1, 'first');
        
        if ~isempty(idx2_repol) && idx2_repol > 1
            t1 = repol_phase_time(idx2_repol - 1);
            v1 = repol_phase_voltage(idx2_repol - 1);
            t2 = repol_phase_time(idx2_repol);
            v2 = repol_phase_voltage(idx2_repol);
            t_interpolated = interpolate_time(t1, v1, t2, v2, V_rep_target);
            
            % --- Store results with CORRECT unit conversion ---
            % Duration (t_interpolated - t_depol) 
            results.(['APD' num2str(level)]) = (t_interpolated - t_depol); % APD in ms
            results.t_rep(i) = t_interpolated;
            results.v_rep(i) = V_rep_target;
        else
            results.(['APD' num2str(level)]) = NaN;
            results.t_rep(i) = NaN;
            results.v_rep(i) = NaN;
        end
    end
    
    % --- Store other results with CORRECT unit conversion ---
    results.RMP = RMP;
    results.APA = APA;
    results.V_max = V_max;
    results.t_peak = t_peak;
    results.t_depol = t_depol; % Time of dVdt_max (used as APD reference)
    results.t_AP_start = t_AP_start; % Time of RMP + 1mV (used for plotting/onset)
    results.V_AP_start = V_AP_start; % The RMP + 1mV threshold voltage
    % dVdt_max is in mV/µs. (1 mV/µs = 1000 V/s).
    results.dVdt_max = dVdt_max * 1000; % dV/dt in V/s

end

% --- Helper function for linear interpolation (unchanged) ---
function t_interp = interpolate_time(t1, v1, t2, v2, v_target)
    t_interp = t1 + (v_target - v1) * (t2 - t1) / (v2 - v1);
end

function DisplayAPMorphology(time, voltage, results)
    % Displays the action potential and overlays the calculated morphology parameters.
    % NOTE: This will be called 5 times if left inside the loop
    
    figure;
    hold on;
    
    % Plot the main AP trace
    plot(time * 1000, voltage, 'k-', 'LineWidth', 1.5);
    title('Action Potential Morphology');
    xlabel('Time (s)');
    ylabel('Voltage (mV)');
    grid on;
    
    % --- Plot Interpolated APD Markers and Lines ---
    apd_levels = [30, 50, 90];
    colors = lines(length(apd_levels));
    
    for i = 1:length(apd_levels)
        t_rep_val = results.t_rep(i);
        v_rep_val = results.v_rep(i);
        
        if ~isnan(t_rep_val)
            % Plot a marker at the exact interpolated point
            plot(t_rep_val * 1000, v_rep_val, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
            
            % --- CHANGED: Plot horizontal line from t_depol to t_rep ---
            line([results.t_depol, t_rep_val] * 1000, [v_rep_val, v_rep_val], ...
                 'Color', colors(i,:), 'LineStyle', '--');
            
            % Add a text label
            text(t_rep_val * 1000 + 5, v_rep_val, ['APD' num2str(apd_levels(i))], 'Color', colors(i,:));
        end
    end
    
    % --- Plot other features ---

    % *** MODIFIED CODE for Depolarization Start Time (RMP + 1mV) ***
    t_plot = results.t_AP_start;
    V_plot = results.V_AP_start;
    
    % 2. Mark the AP start time (t_AP_start) with a red dot
    plot(t_plot * 1000, V_plot, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
    % **********************************************

    % Mark the peak voltage
    plot(results.t_peak * 1000, results.V_max, 'r^', 'MarkerFaceColor', 'r');
    
    % Mark the RMP
    line([time(1), time(end)]*1000, [results.RMP, results.RMP], 'Color', 'b', 'LineStyle', ':');
    text(time(5)*1000, results.RMP + 2, 'RMP');
    
    hold off;
    legend('AP Trace', 'APD30','', 'APD50','', 'APD90','', 'Depol. Start', 'V_{max}', 'Location', 'best');
end



function DisplayAPMorphologyTable(results)
    % Displays the calculated AP morphology parameters in a formatted table.
    % NOTE: This function is no longer called by the main script,
    % as we now use the new summary table generator.
    
    % Create a header for the table
    fprintf('\n--- Action Potential Morphology ---\n');
    fprintf('===================================\n');
    fprintf('%-12s %10s %-8s\n', 'Parameter', 'Value', 'Units');
    fprintf('-----------------------------------\n');
    
    % Populate the table with data from the results struct
    fprintf('%-12s %10.2f %-8s\n', 'Resting Membrane Potential', results.RMP, 'mV');
    fprintf('%-12s %10.2f %-8s\n', 'Amplitude', results.APA, 'mV');
    fprintf('%-12s %10.2f %-8s\n', 'Vmax', results.V_max, 'mV');
    fprintf('%-12s %10.2f %-8s\n', 'dV/dt max', results.dVdt_max / 1000, 'V/s');
    fprintf('%-12s %10.2f %-8s\n', 'APD30', results.APD30, 'ms');
    fprintf('%-12s %10.2f %-8s\n', 'APD50', results.APD50, 'ms');
    fprintf('%-12s %10.2f %-8s\n', 'APD90', results.APD90, 'ms');
    
    fprintf('===================================\n\n');
end