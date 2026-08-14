%% run_paper1_scale_match.m
%% Scale-matched run of AMD_Electrodialysis_MetalRecovery against Paper #1,
%% via the Presets/AMD_ED_Preset_Isaksson2025_Paper1 preset (single source
%% of truth for that paper's feed/stack numbers -- see that file's header
%% for the model-fidelity caveats, e.g. CEM/cation vs. their AEM/chelate
%% transport mechanism).

clear; clc;
repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
addpath(fullfile(repo_root, 'Presets'));

[feed_ov, stack_ov, meta] = AMD_ED_Preset_Isaksson2025_Paper1();

fprintf('\n══════════════════════════════════════════════\n');
fprintf('  Preset: %s\n', meta.citation);
fprintf('  %s\n', meta.url);
fprintf('══════════════════════════════════════════════\n\n');

R = AMD_Electrodialysis_MetalRecovery(feed_ov, stack_ov, true);

%% Single-pass -> equivalent recirculating-batch pass count
%% Approximates their fixed-volume recirculation as repeated single passes
%% through this stack, decaying geometrically by the per-pass removal
%% fraction, to reach their reported depletion endpoint.
n_passes  = log(1 - meta.target_depletion) / log(1 - R.rem_frac.Ni);
t_equiv_h = n_passes * meta.donor_vol_L / feed_ov.Q_in;

fprintf('\n══════════════════════════════════════════════\n');
fprintf('  Single-pass -> batch-recirculation equivalence\n');
fprintf('══════════════════════════════════════════════\n');
fprintf('  Ni removed per pass  : %.2f%%\n', R.rem_frac.Ni*100);
fprintf('  Target depletion     : %.2f%% (paper conductivity endpoint: %.1f -> %.1f mS/cm)\n', ...
    meta.target_depletion*100, meta.donor_conductivity_start_mScm, meta.donor_conductivity_end_mScm);
fprintf('  Equivalent passes    : %.0f\n', n_passes);
fprintf('  Equivalent run time  : %.0f min  (donor vol %.0f cm^3 @ %.0f L/h)\n', ...
    t_equiv_h*60, meta.donor_vol_L*1000, feed_ov.Q_in);
fprintf('  Paper reported result: Ni %.1f%% / Co %.1f%% separated at batch endpoint\n', ...
    meta.Ni_separated_pct, meta.Co_separated_pct);
fprintf('══════════════════════════════════════════════\n\n');
