function [ ca_results, validation_data ]= ca_analysis_v6( Time, Iup, INaCa, IpCa, Ca, plot_title )
%calcium transient analysis
%output: plots for Figure 10A&C

% Constants (copied from ipsc_function)
V_tot=3960; %um^3 from hwang et al.
Vc_tenT=16404; VSR_tenT=1094; V_tot_tenT=Vc_tenT+VSR_tenT;
Vc=V_tot*(Vc_tenT/V_tot_tenT);
Cm = 60; %pF
F = 96.4853415;   % coulomb_per_mmole (in model_parameters)

%% Find first beat to analyze 1200 and 2400 originally | tried 1200 and 2400
inds_time_1200=find(Time>1200);inds_time_1200=inds_time_1200(1);
inds_time_2400=find(Time>2400); inds_time_2400=inds_time_2400(1);
[~,inds1]=min(Ca(1:inds_time_1200));
[~,inds2]=min(Ca(inds_time_1200:inds_time_2400)); inds2=inds_time_1200+inds2;



%% --- DEFINE BEAT VECTORS ---
% This is the critical step that defines the 'ca_beat' variable
ca_beat = Ca(inds1:inds2);
time_beat = Time(inds1:inds2);



%% Calculate Normalized Ca2+ flux 
%take integral
intJserca = cumtrapz(Time,Iup);
intIncx_ca = cumtrapz(Time,-INaCa*2*Cm/(2.0*Vc*F));
intIpca = cumtrapz(Time,IpCa*Cm/(2.0*Vc*F));

%integral for first beat
fluxJserca = intJserca(inds1:inds2)-intJserca(inds1);
fluxIncx_ca = intIncx_ca(inds1:inds2)-intIncx_ca(inds1);
fluxIpca = intIpca(inds1:inds2)-intIpca(inds1);

%Normalize flux
flux_total=fluxJserca+fluxIncx_ca+fluxIpca;
ref=max(flux_total);
fluxJserca_norm=fluxJserca./ref;
fluxIncx_ca_norm=fluxIncx_ca./ref;
fluxIpca_norm=fluxIpca./ref;
Time_flux=Time(inds1:inds2);

%% plot figure 10A 
figure, set(gcf,'color','w')
plot((Time(inds1:inds2)-Time(inds1))./1000, Ca(inds1:inds2).*1e6, 'Color', [.8 0 .18])
set(gca,'box','off','tickdir','out')
legend boxoff
ylabel('[Ca^{2+}] (nM)')
xlabel('Time (s)')
title(plot_title); % <--- ADD THIS LINE

%% plot figure 10C 
figure,set(gcf,'color','w')
plot((Time_flux-Time(inds1)),fluxJserca_norm,(Time_flux-Time(inds1)),fluxIncx_ca_norm,(Time_flux-Time(inds1)),fluxIpca_norm);
set(gca,'box','off','tickdir','out')
legend('SERCA', 'NCX', 'non-NCX')
legend boxoff
ylabel('Ca flux normalized')
xlabel('Time (ms)')
title(plot_title) % <--- ADD THIS LINE






%% =======================================================================
%  START: CALCULATE AND STORE RESULTS
%  =======================================================================

% --- Table 1: Calcium Transient Morphology ---
ca_results.diastolic_ca_nM = min(ca_beat) * 1e6;
[peak_val, peak_idx] = max(ca_beat);
ca_results.peak_ca_nM = peak_val * 1e6;

% --- CORRECTED Time to Peak (TTP) Calculation ---
% 1. Calculate the first derivative of the calcium transient (dCa/dt)
dCa_dt = diff(ca_beat) ./ diff(time_beat);
% 2. Find the index of the maximum derivative (point of fastest upstroke)
[~, upstroke_idx] = max(dCa_dt); 
% 3. Calculate TTP from the start of the upstroke to the peak
ca_results.ttp_ms = time_beat(peak_idx) - time_beat(upstroke_idx);

% --- Table 2: % Contribution of Calcium Flux ---
% This section reports the final contribution of each pathway to the
% total calcium removal over the entire analyzed beat. The values are
% taken from the last timepoint of the normalized cumulative flux vectors.

ca_results.pct_serca = fluxJserca_norm(end) * 100;
ca_results.pct_ncx   = fluxIncx_ca_norm(end) * 100;
ca_results.pct_ipca  = fluxIpca_norm(end) * 100;

% --- Tau Decay (Bi-exponential Fit) ---
% 1. Prepare full decay data
decay_time_full = time_beat(peak_idx:end) - time_beat(peak_idx);
decay_cai_full = ca_beat(peak_idx:end);
normalized_decay_full = (decay_cai_full - min(ca_beat)) / (peak_val - min(ca_beat));

% 2. Find the index where decay has dropped to 0% to bypass the "shoulder".
% If wanting to change to 90%, change "1, 1" to "0.9, 1"
fit_start_index = find(normalized_decay_full <= 1, 1, 'first');
if isempty(fit_start_index), fit_start_index = 1; end

% 3. Truncate the data to the true exponential portion
decay_time_trunc = decay_time_full(fit_start_index:end);
normalized_decay_trunc = normalized_decay_full(fit_start_index:end);

% --- ✅ FINAL CORRECTION: Re-zero and Re-normalize the truncated data ---
% This ensures the data passed to the fit function starts at (t=0, y=1)
decay_time_final = decay_time_trunc - decay_time_trunc(1); % Re-zero time
normalized_decay_final = normalized_decay_trunc / normalized_decay_trunc(1); % Re-normalize to 1.0

% 4. Define the BI-EXPONENTIAL model
bi_exp_model = fittype('(a1 * exp(-x/b1)) + ((1-a1) * exp(-x/b2))', 'independent', 'x', 'dependent', 'y');

% 5. Define fitting parameters
start_points = [0.8, 40, 300];     % [Amp_fast, Tau_fast, Tau_slow]
lower_bounds = [0, 5, 5];
upper_bounds = [1, 2000, 4000];

% 6. Perform the fit on the FINAL, correctly formatted data
[fit_object, gof] = fit(decay_time_final, normalized_decay_final, bi_exp_model, ...
    'StartPoint', start_points, ...
    'Lower', lower_bounds, ...
    'Upper', upper_bounds);

% 7. Extract, sort, and store results (code remains the same)
taus = sort([fit_object.b1, fit_object.b2]);
ca_results.tau_fast_ms = taus(1);
ca_results.tau_slow_ms = taus(2);

if fit_object.b1 == taus(1)
    ca_results.amp_fast = fit_object.a1;
    ca_results.amp_slow = 1 - fit_object.a1;
else
    ca_results.amp_fast = 1 - fit_object.a1;
    ca_results.amp_slow = fit_object.a1;
end
ca_results.tau_decay_r_squared = gof.rsquare;

%% =======================================================================
%  BUNDLE AND RETURN VALIDATION DATA
%  =======================================================================
validation_data.time_beat = time_beat;
validation_data.ca_beat = ca_beat;
validation_data.peak_idx = peak_idx;
validation_data.upstroke_idx = upstroke_idx;

%% =======================================================================
%  UPDATED: PLOT TO VISUALIZE THE BI-EXPONENTIAL FIT
%  =======================================================================
figure, set(gcf,'color','w');
hold on;

% --- Plot the data used for the fit
plot(decay_time_trunc, normalized_decay_trunc, 'bo', ...
    'MarkerFaceColor', 'b', ...
    'DisplayName', 'Data for Fit');

% --- Evaluate the fit object and scale it back to the original data's coordinates
% 1. Evaluate the fit using the re-zeroed time vector
y_fit_on_normed_scale = fit_object(decay_time_final);

% 2. "Un-normalize" the fitted y-values by multiplying by the original start value
y_fit_on_original_scale = y_fit_on_normed_scale * normalized_decay_trunc(1);

% 3. Create the legend text using sprintf for cleaner syntax
fit_legend_text = sprintf('Fit: \\tau_{fast} = %.1f ms (%.0f%%)', ...
    ca_results.tau_fast_ms, ca_results.amp_fast * 100);

% 4. Plot the correctly scaled fit line against the original truncated time vector
plot(decay_time_trunc, y_fit_on_original_scale, 'r-', ...
    'LineWidth', 2, ...
    'DisplayName', fit_legend_text);

hold off;

% --- Add titles, labels, and legend
title({plot_title; sprintf('Bi-Exponential Tau Decay Fit (R^2 = %.4f)', gof.rsquare)});
xlabel('Time from Peak (ms)');
ylabel('Normalized Amplitude');
legend('show', 'Location', 'northeast'); % Display the legend
legend boxoff;
set(gca, 'box', 'off', 'tickdir', 'out');


end