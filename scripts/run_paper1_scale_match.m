%% run_paper1_scale_match.m
%% Scale-matched run of AMD_Electrodialysis_MetalRecovery against Paper #1
%% in References/electrodialysis_metal_recovery_papers.xlsx:
%%   Isaksson et al., "Nickel and Cobalt Recovery from Spent Lithium-Ion
%%   Batteries via Electrodialysis Metathesis", Membranes, 2025.
%%   https://pmc.ncbi.nlm.nih.gov/articles/PMC12028818/
%%
%% Matches their STACK HARDWARE + FEED CONCENTRATION only (N_cp, membrane
%% area, current density, flow, feed molarity). Does NOT reproduce their
%% transport MECHANISM: their Ni/Co are EDTA-chelated and migrate as
%% anionic complexes through an AEM (electrodialysis metathesis), whereas
%% this model still treats Ni2+/Co2+ as free cations migrating through a
%% CEM. Treat results as an order-of-magnitude/scale sanity check, not a
%% mechanistic reproduction.
%%
%% Their run is also a BATCH recirculation of a fixed 250 cm^3 donor
%% volume to near-total depletion (conductivity 40 -> 0.3 mS/cm, ~99.25%
%% depleted), vs. this model's single-pass steady-state removal. See the
%% single-pass -> equivalent-pass-count note printed at the end.

clear; clc;
cd(fileparts(fileparts(mfilename('fullpath'))));   % repo root

feed_ov.Q_in    = 20;                  % L/h  (paper: 20 dm3/h per stream)
feed_ov.Ni_mgL  = 0.10*58.69*1000;     % mg/L, from 0.10 M NiSO4
feed_ov.Co_mgL  = 0.10*58.93*1000;     % mg/L, from 0.10 M CoSO4
feed_ov.SO4_mgL = 0.20*96.06*1000;     % mg/L, SO4 counter-ion from both sulfates

stack_ov.N_cp  = 5;         % paper: 5 repeating units
stack_ov.A_mem = 0.0064;    % paper: 64 cm^2 active membrane area per cell pair
stack_ov.i_app = 200;       % paper: initial current density, A/m^2 (they run
                             % constant-VOLTAGE at 10 V; this is only the
                             % initial current density, not a fixed setpoint)

R = AMD_Electrodialysis_MetalRecovery(feed_ov, stack_ov, true);

%% Single-pass -> equivalent recirculating-batch pass count
%% Approximates their fixed-volume recirculation as repeated single passes
%% through this stack, decaying geometrically by the per-pass removal
%% fraction, to reach their reported ~99.25% depletion endpoint.
target_depletion = 1 - 0.3/40;             % conductivity 40 -> 0.3 mS/cm
n_passes = log(1 - target_depletion) / log(1 - R.rem_frac.Ni);
donor_vol_L = 0.250;                       % paper: 250 cm^3 donor volume
t_equiv_h   = n_passes * donor_vol_L / feed_ov.Q_in;

fprintf('\n══════════════════════════════════════════════\n');
fprintf('  Single-pass -> batch-recirculation equivalence\n');
fprintf('══════════════════════════════════════════════\n');
fprintf('  Ni removed per pass  : %.2f%%\n', R.rem_frac.Ni*100);
fprintf('  Target depletion     : %.2f%% (paper conductivity endpoint)\n', target_depletion*100);
fprintf('  Equivalent passes    : %.0f\n', n_passes);
fprintf('  Equivalent run time  : %.0f min  (donor vol %.0f cm^3 @ %.0f L/h)\n', ...
    t_equiv_h*60, donor_vol_L*1000, feed_ov.Q_in);
fprintf('══════════════════════════════════════════════\n\n');
