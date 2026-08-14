%% run_chelator_pH_secondaryfeed_sensitivity.m
%% Sensitivity of Ni/Co recovery to the chemistry-proxy inputs added to
%% AMD_Electrodialysis_MetalRecovery: p.chelator, p.secondary_feed, p.pH.
%%
%% IMPORTANT: these are ASSUMPTION-FLAGGED proxy multipliers on eta_Ni/
%% eta_Co (see AMD_Electrodialysis_MetalRecovery.m header + chelatorFactor_ED/
%% secondaryFeedFactor_ED/pHFactor_ED), calibrated so EDTA + Na2SO4 + pH 7
%% (Isaksson et al. 2025's exact recipe) reduces to factor=1.0. Results here
%% are illustrative sensitivity, NOT a validated quantitative correlation --
%% no pH-sweep or multi-chelator Ni/Co dataset exists in the two papers
%% reviewed so far (References/electrodialysis_metal_recovery_papers.xlsx).
%%
%% Base case: Article-1 (Isaksson 2025) feed/stack preset, so results are
%% directly comparable to scripts/run_paper1_scale_match.m.

clear; clc;
repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);
addpath(fullfile(repo_root, 'Presets'));

[feed_ov, stack_ov_base, meta] = AMD_ED_Preset_Isaksson2025_Paper1();

chelators = {'EDTA','DTPA','HEDTA','GLDA','none'};
sec_feeds = {'Na2SO4','H2SO4'};
pH_sweep  = [3 5 7 9 11];

%% ── Sweep 1: chelator x secondary feed, at baseline pH 7 ────────────────
n1 = numel(chelators)*numel(sec_feeds);
T1_chelator = cell(n1,1); T1_secfeed = cell(n1,1);
T1_Ni_pct = zeros(n1,1); T1_Co_pct = zeros(n1,1);
row = 1;
for ic = 1:numel(chelators)
    for is = 1:numel(sec_feeds)
        stack_ov = stack_ov_base;
        stack_ov.chelator = chelators{ic};
        stack_ov.secondary_feed = sec_feeds{is};
        stack_ov.pH = 7;
        R = AMD_Electrodialysis_MetalRecovery(feed_ov, stack_ov, false);
        T1_chelator{row} = chelators{ic};
        T1_secfeed{row}  = sec_feeds{is};
        T1_Ni_pct(row) = R.rem_frac.Ni*100;
        T1_Co_pct(row) = R.rem_frac.Co*100;
        row = row + 1;
    end
end
T1 = table(T1_chelator, T1_secfeed, T1_Ni_pct, T1_Co_pct, ...
    'VariableNames', {'Chelator','SecondaryFeed','Ni_removed_pct','Co_removed_pct'});

%% ── Sweep 2: pH, for EDTA at each secondary feed ─────────────────────────
n2 = numel(pH_sweep)*numel(sec_feeds);
T2_pH = zeros(n2,1); T2_secfeed = cell(n2,1);
T2_Ni_pct = zeros(n2,1); T2_Co_pct = zeros(n2,1);
row = 1;
for ip = 1:numel(pH_sweep)
    for is = 1:numel(sec_feeds)
        stack_ov = stack_ov_base;
        stack_ov.chelator = 'EDTA';
        stack_ov.secondary_feed = sec_feeds{is};
        stack_ov.pH = pH_sweep(ip);
        R = AMD_Electrodialysis_MetalRecovery(feed_ov, stack_ov, false);
        T2_pH(row) = pH_sweep(ip);
        T2_secfeed{row} = sec_feeds{is};
        T2_Ni_pct(row) = R.rem_frac.Ni*100;
        T2_Co_pct(row) = R.rem_frac.Co*100;
        row = row + 1;
    end
end
T2 = table(T2_pH, T2_secfeed, T2_Ni_pct, T2_Co_pct, ...
    'VariableNames', {'pH','SecondaryFeed','Ni_removed_pct','Co_removed_pct'});

%% ── Save tables ───────────────────────────────────────────────────────
save_dir = fullfile(repo_root, 'ChemistrySensitivity_Plots');
if ~exist(save_dir,'dir'); mkdir(save_dir); end
writetable(T1, fullfile(save_dir,'chelator_secondaryfeed_sensitivity.csv'));
writetable(T2, fullfile(save_dir,'pH_secondaryfeed_sensitivity.csv'));

disp('=== Chelator x Secondary Feed (pH 7) ===');
disp(T1);
disp('=== pH x Secondary Feed (EDTA) ===');
disp(T2);

%% ── Dashboard: 2 panels ──────────────────────────────────────────────────
hDash = figure('Name','Chemistry Sensitivity — Chelator/pH/Secondary Feed', 'Position',[120 80 1400 560]);
tl = tiledlayout(hDash, 1, 2, 'TileSpacing','compact', 'Padding','compact');

% Panel 1: grouped bar, chelator x secondary feed, Ni removal
ax1 = nexttile(tl);
Ni_matrix = reshape(T1_Ni_pct, numel(sec_feeds), numel(chelators))';  % chelator rows, secfeed cols
b = bar(ax1, Ni_matrix, 'grouped');
b(1).FaceColor = [0.10 0.30 0.80]; b(2).FaceColor = [0.85 0.55 0.10];
set(ax1, 'XTickLabel', chelators);
ylabel(ax1, 'Ni removed per pass (%)');
legend(ax1, sec_feeds, 'Location','best');
title(ax1, {'Chelator x Secondary Feed', 'Ni removal, pH 7 (illustrative proxy)'}, 'FontSize',10);
grid(ax1,'on');

% Panel 2: line plot, pH sweep, Ni removal, one line per secondary feed
ax2 = nexttile(tl);
hold(ax2,'on');
colors = [0.10 0.30 0.80; 0.85 0.55 0.10];
for is = 1:numel(sec_feeds)
    mask = strcmp(T2_secfeed, sec_feeds{is});
    plot(ax2, T2_pH(mask), T2_Ni_pct(mask), '-o', 'LineWidth',2, 'Color',colors(is,:), ...
        'DisplayName', sec_feeds{is});
end
xline(ax2, 7, '--k', 'pH 7 (paper baseline)', 'FontSize',7, 'HandleVisibility','off');
xlabel(ax2, 'pH'); ylabel(ax2, 'Ni removed per pass (%)');
legend(ax2, 'Location','best');
title(ax2, {'pH Sweep (EDTA)', 'Ni removal vs. pH (illustrative proxy)'}, 'FontSize',10);
grid(ax2,'on');

title(tl, 'Chemistry Sensitivity — ASSUMPTION-FLAGGED proxy, not a validated correlation', ...
    'FontSize',11, 'FontWeight','bold');

exportgraphics(hDash, fullfile(save_dir,'ChemistrySensitivity_Dashboard.png'), 'Resolution',150);
fprintf('\nSaved tables + dashboard to: %s\n', save_dir);
