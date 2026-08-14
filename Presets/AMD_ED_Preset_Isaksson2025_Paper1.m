function [feed_ov, stack_ov, meta] = AMD_ED_Preset_Isaksson2025_Paper1()
%% ╔══════════════════════════════════════════════════════════════════╗
%% ║   ED Preset — Isaksson et al. 2025 (Paper #1)                     ║
%% ║   Nickel and Cobalt Recovery from Spent Lithium-Ion Batteries via ║
%% ║   Electrodialysis Metathesis, Membranes, 2025.                    ║
%% ║   DOI: 10.3390/membranes15040097                                  ║
%% ║   https://pmc.ncbi.nlm.nih.gov/articles/PMC12028818/              ║
%% ╠══════════════════════════════════════════════════════════════════╣
%% ║  Returns feed/stack overrides for AMD_Electrodialysis_MetalRecovery║
%% ║  matching this paper's STACK HARDWARE + FEED CONCENTRATION only.  ║
%% ║  Does NOT reproduce their transport MECHANISM: their Ni/Co are    ║
%% ║  EDTA-chelated and migrate as anionic complexes through an AEM    ║
%% ║  (electrodialysis metathesis, swapped against SO4 from a          ║
%% ║  secondary Na2SO4 stream), whereas AMD_Electrodialysis_MetalRecovery║
%% ║  still treats Ni2+/Co2+ as free cations migrating through a CEM.  ║
%% ║  Use for order-of-magnitude/scale comparison, not a mechanistic   ║
%% ║  reproduction. meta carries their batch-run details for that      ║
%% ║  single-pass -> equivalent-batch-passes translation.              ║
%% ╚══════════════════════════════════════════════════════════════════╝

feed_ov.Q_in    = 20;                  % L/h  (paper: 20 dm3/h per stream)
feed_ov.Ni_mgL  = 0.10*58.69*1000;     % mg/L, from 0.10 M NiSO4
feed_ov.Co_mgL  = 0.10*58.93*1000;     % mg/L, from 0.10 M CoSO4
feed_ov.SO4_mgL = 0.20*96.06*1000;     % mg/L, SO4 counter-ion from both sulfates

stack_ov.N_cp  = 5;         % paper: 5 repeating units
stack_ov.A_mem = 0.0064;    % paper: 64 cm^2 active membrane area per cell pair
stack_ov.i_app = 200;       % paper: initial current density, A/m^2 (they run
                             % constant-VOLTAGE at 10 V; this is only the
                             % initial current density, not a fixed setpoint)

meta.citation         = 'Isaksson et al., "Nickel and Cobalt Recovery from Spent Lithium-Ion Batteries via Electrodialysis Metathesis", Membranes, 2025';
meta.doi               = '10.3390/membranes15040097';
meta.url                = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12028818/';
meta.stack_voltage_V    = 10;          % constant-voltage control mode (not modeled -- this preset is constant-i_app)
meta.donor_vol_L        = 0.250;       % their batch donor volume, L
meta.donor_conductivity_start_mScm = 40;
meta.donor_conductivity_end_mScm   = 0.3;
meta.target_depletion   = 1 - meta.donor_conductivity_end_mScm/meta.donor_conductivity_start_mScm;
meta.Ni_separated_pct   = 97.9;        % their reported batch-endpoint result
meta.Co_separated_pct   = 96.6;

end
