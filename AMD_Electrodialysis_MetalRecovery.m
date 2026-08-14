function R = AMD_Electrodialysis_MetalRecovery(feed_ov, stack_ov, make_plots)
%% ╔══════════════════════════════════════════════════════════════════╗
%% ║   Electrodialysis (ED) — METAL RECOVERY STACK  [standalone]     ║
%% ║   Steady-state ED stack: concentrates Ni/Co from a dilute feed  ║
%% ║   into a recirculating concentrate loop via an applied current, ║
%% ║   with SO4 tracked as the balancing counter-ion (AEM side).     ║
%% ╠══════════════════════════════════════════════════════════════════╣
%% ║  INDEPENDENT DESIGN — not wired into AMD_TreatmentTrain_Seq4.    ║
%% ║  Feed is synthetic/user-specified, not read from another module. ║
%% ╠══════════════════════════════════════════════════════════════════╣
%% ║  MODEL (steady-state, algebraic — same style as Lime/Thickener): ║
%% ║   I        = i_app * A_mem                    [A, per cell pair] ║
%% ║   eta_eff_i= eta_i * chelatorFactor*secFeedFactor*pHFactor(Ni/Co only) [chemistry proxy, see below] ║
%% ║   Ndot_i   = N_cp*eta_eff_i*I/(z_i*F) * 3600   [mol/h transferred]║
%% ║   C_i,dil_out = C_i,in - Ndot_i/Q_d            [diluate outlet]   ║
%% ║   Q_c      sized so Ni hits CF_target          [concentrate flow] ║
%% ║   C_i,conc_out = Ndot_i/Q_c                    [concentrate, fresh-water makeup] ║
%% ║   i_lim_i  = z_i*F*D_i*C_i,dil_out/delta       [limiting current, worst-case outlet] ║
%% ║   V_stack  = N_cp*(V_junction + i_app*R_area)  [stack voltage]    ║
%% ╠══════════════════════════════════════════════════════════════════╣
%% ║  CHEMISTRY PROXY (chelator / secondary feed / pH):                ║
%% ║  This model still treats Ni2+/Co2+ as free cations through a CEM  ║
%% ║  (see Presets/AMD_ED_Preset_Isaksson2025_Paper1.m header) — it does║
%% ║  NOT simulate chelation equilibria or AEM anionic-complex transport║
%% ║  mechanistically. p.chelator/p.secondary_feed/p.pH instead apply   ║
%% ║  literature-informed, ASSUMPTION-FLAGGED multiplier factors onto   ║
%% ║  eta_Ni/eta_Co, calibrated to factor=1.0 at Isaksson et al.'s exact║
%% ║  recipe (EDTA, Na2SO4, pH~7). Treat swept results as illustrative  ║
%% ║  sensitivity, not a validated quantitative correlation — see       ║
%% ║  chelatorFactor_ED/secondaryFeedFactor_ED/pHFactor_ED for sources. ║
%% ╚══════════════════════════════════════════════════════════════════╝

fprintf('\n╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║   Electrodialysis — METAL RECOVERY STACK  [standalone design]     ║\n');
fprintf('╚══════════════════════════════════════════════════════════════════╝\n\n');

if nargin < 1; feed_ov  = []; end
if nargin < 2; stack_ov = []; end
if nargin < 3 || isempty(make_plots); make_plots = true; end

feed = defineFeed_ED(feed_ov);
p    = defineParameters_ED(feed, stack_ov);
R    = runED(p);
if make_plots; plotDashboard_ED(R, p); end
R.p = p;

end


%% ═══════════════════════════════════════════════════════════════════════
%%  FEED DEFINITION — synthetic/user-specified (no upstream module dependency)
%% ═══════════════════════════════════════════════════════════════════════
function feed = defineFeed_ED(feed_ov)

feed.Q_in    = 80;      % L/h   diluate feed flow
feed.Ni_mgL  = 250;     % mg/L  Ni2+  (concentrated stream, e.g. IX eluate/bleed -- ED needs
feed.Co_mgL  = 60;      % mg/L  Co2+   this range to keep i_app under the limiting current;
feed.SO4_mgL = 5000;    % mg/L  SO4^2- (dominant counter-anion)  dilute AMD-effluent-level
                         % feeds (tens of mg/L) push i_lim below any practical i_app -- see EP3.

if ~isempty(feed_ov)
    f = fieldnames(feed_ov);
    for k = 1:numel(f); feed.(f{k}) = feed_ov.(f{k}); end
end

fprintf('══════════════════════════════════════════════\n');
fprintf('  Electrodialysis — synthetic diluate FEED\n');
fprintf('══════════════════════════════════════════════\n');
fprintf('  Q_in : %.1f L/h\n', feed.Q_in);
fprintf('  Ni   : %.2f mg/L\n', feed.Ni_mgL);
fprintf('  Co   : %.2f mg/L\n', feed.Co_mgL);
fprintf('  SO4  : %.0f mg/L\n', feed.SO4_mgL);
fprintf('══════════════════════════════════════════════\n\n');

end


%% ═══════════════════════════════════════════════════════════════════════
%%  PARAMETER DEFINITION — stack design + derived molar feed
%% ═══════════════════════════════════════════════════════════════════════
function p = defineParameters_ED(feed, stack_ov)

p.MW.Ni = 58.69; p.MW.Co = 58.93; p.MW.SO4 = 96.06;   % g/mol
p.z.Ni  = 2;      p.z.Co  = 2;     p.z.SO4  = 2;        % charge
p.F     = 96485;                                        % C/mol

p.in.Q   = feed.Q_in;
p.in.Ni  = feed.Ni_mgL  / 1000 / p.MW.Ni;   % mol/L
p.in.Co  = feed.Co_mgL  / 1000 / p.MW.Co;   % mol/L
p.in.SO4 = feed.SO4_mgL / 1000 / p.MW.SO4;  % mol/L

%── Stack design (defaults — typical industrial ED ranges) ───────────────
p.N_cp     = 13;       % number of cell pairs
p.A_mem    = 1.0;      % m^2  effective membrane area per cell pair -- sized for total
                        % current, NOT i_app: i_app is capped by the limiting current
                        % density (concentration-polarization physics, see EP3), which
                        % falls as the diluate depletes toward the outlet -- so more
                        % recovery must come from more membrane area, not higher i_app.
p.i_app    = 1.5;       % A/m^2  applied current density (dilute feed -> low i_lim, see below)
% eta_Ni + eta_Co MUST stay <= 1: both compete for the SAME cation-side current
% budget through the CEM (charge conservation/electroneutrality). eta_SO4 is an
% independent anion-side budget through the AEM. Any remainder of either budget
% (here 0.30 cation-side, 0.40 anion-side) is carried by untracked ions (Ca,
% Mg, Na, H+ / Cl, HCO3) — do not size eta_Ni+eta_Co above 1, it silently
% overstates transferable Ni/Co beyond what the current can physically carry.
p.eta.Ni   = 0.55;     % current efficiency, Ni2+ transport (CEM side)
p.eta.Co   = 0.15;     % current efficiency, Co2+ transport (CEM side)
p.eta.SO4  = 0.60;     % current efficiency, SO4^2- transport (AEM side, independent budget)
p.D.Ni     = 6.6e-10;  % m^2/s  diffusivity of Ni2+ in water
p.D.Co     = 7.3e-10;  % m^2/s  diffusivity of Co2+ in water
p.delta    = 20e-6;    % m      diffusion boundary layer thickness (turbulence-promoting spacer)
p.R_area   = 8e-4;     % ohm*m^2  areal resistance (membranes + solution), i.e. 8 ohm*cm^2
p.V_junct  = 0.4;      % V  membrane potential + electrode overpotential, per cell pair
p.CF_target = 8;       % -  target concentration factor for Ni in the concentrate loop
p.i_margin  = 0.8;     % -  design margin vs limiting current density (avoid water splitting)

%── Chemistry proxy inputs — see file-header CHEMISTRY PROXY note ────────
p.chelator       = 'EDTA';    % 'EDTA' | 'HEDTA' | 'DTPA' | 'GLDA' | 'none'
p.secondary_feed = 'Na2SO4';  % 'Na2SO4' | 'H2SO4'
p.pH             = 7;         % feed pH, dimensionless

if ~isempty(stack_ov)
    f = fieldnames(stack_ov);
    for k = 1:numel(f); p.(f{k}) = stack_ov.(f{k}); end
end

%── Derived: total stack current (series circuit — same I through every cell pair) ─
p.I = p.i_app * p.A_mem;   % A

%── Chemistry proxy factors -> effective (chemistry-adjusted) Ni/Co efficiency ──
p.chem.chelatorFactor  = chelatorFactor_ED(p.chelator);
p.chem.secFeedFactor   = secondaryFeedFactor_ED(p.secondary_feed);
p.chem.pHFactor        = pHFactor_ED(p.pH);
p.chem.combinedFactor  = p.chem.chelatorFactor * p.chem.secFeedFactor * p.chem.pHFactor;
p.eta_eff.Ni  = p.eta.Ni  * p.chem.combinedFactor;
p.eta_eff.Co  = p.eta.Co  * p.chem.combinedFactor;
p.eta_eff.SO4 = p.eta.SO4;   % SO4 (AEM/counter-ion side) not adjusted -- see header

%── Design summary ────────────────────────────────────────────────────────
fprintf('══════════════════════════════════════════════\n');
fprintf('  Electrodialysis Stack — Design\n');
fprintf('══════════════════════════════════════════════\n');
fprintf('  Cell pairs N_cp    : %d\n', p.N_cp);
fprintf('  Membrane area/pair : %.3f m^2  (total %.1f m^2)\n', p.A_mem, p.N_cp*p.A_mem*2);
fprintf('  Applied i_app      : %.0f A/m^2  ->  I = %.2f A\n', p.i_app, p.I);
fprintf('  Chemistry (proxy)  : chelator=%s (x%.2f)  secondary feed=%s (x%.2f)  pH=%.1f (x%.2f)\n', ...
    p.chelator, p.chem.chelatorFactor, p.secondary_feed, p.chem.secFeedFactor, p.pH, p.chem.pHFactor);
fprintf('  Current efficiency : Ni %.0f%% -> eff %.0f%%   Co %.0f%% -> eff %.0f%%   SO4 %.0f%%\n', ...
    p.eta.Ni*100, p.eta_eff.Ni*100, p.eta.Co*100, p.eta_eff.Co*100, p.eta_eff.SO4*100);
fprintf('  Target CF (Ni)     : %.1fx\n', p.CF_target);
fprintf('══════════════════════════════════════════════\n\n');

end


%% ═══════════════════════════════════════════════════════════════════════
%%  CHEMISTRY PROXY FACTORS — chelator / secondary feed / pH
%%  ASSUMPTION-FLAGGED multipliers on eta_Ni/eta_Co (see file header). All
%%  factors = 1.0 at Isaksson et al. 2025's exact recipe (EDTA, Na2SO4,
%%  pH~7) so the default run is unchanged from before this feature existed.
%% ═══════════════════════════════════════════════════════════════════════
function factor = chelatorFactor_ED(chelator)
% Isaksson et al. 2025 used EDTA (+ Na2SO4) to reach 97.9%/96.6% Ni/Co
% separation -- taken as baseline (1.00). Xing & Srinivasan 2023 (Table 1,
% LCO leachate) found DTPA > HEDTA > GLDA > EDTA for THEIR target ion (Li)
% under bipolar-membrane ED, a different mechanism/target metal; the
% relative ordering is carried over here only as a qualitative proxy for
% "how readily this chelator supports metal-complex transport", NOT a
% validated Ni/Co number for any agent besides EDTA. 'none' (no chelator)
% has no reported EDM/BMED Ni-Co result in either paper we reviewed --
% treated as a substantial illustrative penalty.
switch chelator
    case 'EDTA';  factor = 1.00;   % literature-anchored (Isaksson et al. 2025)
    case 'DTPA';  factor = 0.90;   % ASSUMPTION, proxy from Xing & Srinivasan Table 1 ordering
    case 'HEDTA'; factor = 0.80;   % ASSUMPTION, proxy
    case 'GLDA';  factor = 0.75;   % ASSUMPTION, proxy
    case 'none';  factor = 0.20;   % ASSUMPTION, illustrative penalty -- no reported EDM data
    otherwise
        error('AMD_ED:chelator', 'Unknown p.chelator "%s" (expected EDTA/DTPA/HEDTA/GLDA/none)', chelator);
end
end

function factor = secondaryFeedFactor_ED(secondary_feed)
% Isaksson et al. 2025 report that Na2SO4 (vs. H2SO4) mitigates membrane
% fouling from EDTA precipitation in the acidic H2SO4 secondary feed --
% baseline Na2SO4 = 1.00; H2SO4 penalty is an ILLUSTRATIVE proxy for that
% qualitative fouling effect, not a measured efficiency loss.
switch secondary_feed
    case 'Na2SO4'; factor = 1.00;  % literature-anchored (fouling-mitigating choice)
    case 'H2SO4';  factor = 0.75;  % ASSUMPTION, illustrative fouling penalty
    otherwise
        error('AMD_ED:secondary_feed', 'Unknown p.secondary_feed "%s" (expected Na2SO4/H2SO4)', secondary_feed);
end
end

function factor = pHFactor_ED(pH)
% ASSUMPTION: EDTA-metal chelate stability (and thus effective transport)
% is strongest near neutral pH -- Xing & Srinivasan's protocol adjusts
% leachate to pH 7 before electrodialysis. Modeled as a Gaussian centered
% at pH0=7 or width sigma=3 (illustrative shape, not fit to a reported
% pH-sweep dataset for Ni/Co specifically -- neither paper we reviewed ran
% one). factor=1.0 exactly at pH 7.
pH0 = 7; sigma = 3;
factor = exp(-0.5*((pH-pH0)/sigma)^2);
end


%% ═══════════════════════════════════════════════════════════════════════
%%  STEADY-STATE STACK SOLUTION
%% ═══════════════════════════════════════════════════════════════════════
function R = runED(p)

%% Faraday's-law molar transfer rates (mol/h), summed over N_cp cell pairs.
%% Ni/Co use the chemistry-adjusted eta_eff (chelator/secondary-feed/pH
%% proxy factors, see file header + chelatorFactor_ED/secondaryFeedFactor_ED/
%% pHFactor_ED); SO4 uses the raw AEM-side eta (unadjusted).
R.Ndot.Ni  = p.N_cp * p.eta_eff.Ni  * p.I / (p.z.Ni  * p.F) * 3600;
R.Ndot.Co  = p.N_cp * p.eta_eff.Co  * p.I / (p.z.Co  * p.F) * 3600;
R.Ndot.SO4 = p.N_cp * p.eta_eff.SO4 * p.I / (p.z.SO4 * p.F) * 3600;

%% Diluate outlet (steady-state mass balance, single pass)
R.dil.Ni  = max(p.in.Ni  - R.Ndot.Ni /p.in.Q, 0);
R.dil.Co  = max(p.in.Co  - R.Ndot.Co /p.in.Q, 0);
R.dil.SO4 = max(p.in.SO4 - R.Ndot.SO4/p.in.Q, 0);

%% Removal fractions (feasibility check — Ndot should not exceed feed load)
R.rem_frac.Ni  = min(R.Ndot.Ni /(p.in.Q*p.in.Ni),  1);
R.rem_frac.Co  = min(R.Ndot.Co /(p.in.Q*p.in.Co),  1);
R.rem_frac.SO4 = min(R.Ndot.SO4/(p.in.Q*p.in.SO4), 1);
if R.Ndot.Ni > p.in.Q*p.in.Ni
    warning('ED:OverCurrent', ['Applied current transfers more Ni than the feed ', ...
        'supplies (Ndot_Ni=%.4f mol/h > feed load=%.4f mol/h). i_app is oversized ', ...
        'for this feed/flow — diluate Ni floors at 0, real removal is feed-limited ', ...
        'not current-limited.'], R.Ndot.Ni, p.in.Q*p.in.Ni);
end

%% Concentrate loop — flow sized so Ni hits the target concentration factor,
%% fresh-water makeup assumed (concentrate inlet salt ~0 at steady state)
p.Q_c = R.Ndot.Ni / (p.CF_target * p.in.Ni);
R.Q_c = p.Q_c;
R.conc.Ni  = R.Ndot.Ni  / p.Q_c;
R.conc.Co  = R.Ndot.Co  / p.Q_c;
R.conc.SO4 = R.Ndot.SO4 / p.Q_c;
R.CF.Ni  = R.conc.Ni  / p.in.Ni;
R.CF.Co  = R.conc.Co  / p.in.Co;
R.CF.SO4 = R.conc.SO4 / p.in.SO4;

%% Limiting current density — evaluated at the diluate OUTLET (worst-case,
%% most depleted point in the stack) for the two depleting cations
R.i_lim.Ni = p.z.Ni * p.F * p.D.Ni * R.dil.Ni*1000 / p.delta;   % C_i in mol/m^3 -> *1000
R.i_lim.Co = p.z.Co * p.F * p.D.Co * R.dil.Co*1000 / p.delta;
R.i_lim_binding = min(R.i_lim.Ni, R.i_lim.Co);
R.i_app_over_lim = p.i_app / R.i_lim_binding;
R.ilim_ok = p.i_app <= p.i_margin * R.i_lim_binding;

%% Electroneutrality diagnostic — cation-side vs anion-side charge accounted
%% for by the TRACKED ions only (Ca/Mg/Na/H+ etc. are untracked and carry
%% the remainder of the current; a gap here is expected, not a leak).
R.eq_cation_h = p.z.Ni*R.Ndot.Ni + p.z.Co*R.Ndot.Co;   % eq/h
R.eq_anion_h  = p.z.SO4*R.Ndot.SO4;                    % eq/h
% Each of the N_cp cell pairs carries the same current I and has its own
% charge budget (I*3600/F eq/h) -- total available scales with N_cp exactly
% like Ndot does, since Ndot was built as N_cp * eta*I/(z*F)*3600.
R.eq_total_available_h = p.N_cp * p.I*3600/p.F;        % eq/h total charge passed, all cell pairs
R.frac_current_tracked_cation = R.eq_cation_h / R.eq_total_available_h;
R.frac_current_tracked_anion  = R.eq_anion_h  / R.eq_total_available_h;

%% Stack voltage / power / specific energy consumption
R.V_stack = p.N_cp * (p.V_junct + p.i_app*p.R_area);   % V
R.P_stack = R.V_stack * p.I;                            % W
R.metal_kg_h = (R.Ndot.Ni*p.MW.Ni + R.Ndot.Co*p.MW.Co) / 1000;  % kg/h
R.SEC_kWh_per_kgMetal = (R.P_stack/1000) / max(R.metal_kg_h, 1e-9);
R.SEC_kWh_per_m3      = (R.P_stack/1000) / (p.in.Q/1000);

%% Print
fprintf('══════════════════════════════════════════════\n');
fprintf('  Electrodialysis Stack — Results\n');
fprintf('══════════════════════════════════════════════\n');
fprintf('  Diluate  Ni : %.2f -> %.2f mg/L  (%.1f%% removed)\n', ...
    p.in.Ni*p.MW.Ni*1000, R.dil.Ni*p.MW.Ni*1000, R.rem_frac.Ni*100);
fprintf('  Diluate  Co : %.2f -> %.2f mg/L  (%.1f%% removed)\n', ...
    p.in.Co*p.MW.Co*1000, R.dil.Co*p.MW.Co*1000, R.rem_frac.Co*100);
fprintf('  Diluate SO4 : %.0f -> %.0f mg/L  (%.1f%% removed)\n', ...
    p.in.SO4*p.MW.SO4*1000, R.dil.SO4*p.MW.SO4*1000, R.rem_frac.SO4*100);
fprintf('  Concentrate loop Q_c : %.2f L/h\n', p.Q_c);
fprintf('  Concentrate Ni/Co/SO4: %.0f / %.0f / %.0f mg/L  (CF %.1fx / %.1fx / %.1fx)\n', ...
    R.conc.Ni*p.MW.Ni*1000, R.conc.Co*p.MW.Co*1000, R.conc.SO4*p.MW.SO4*1000, ...
    R.CF.Ni, R.CF.Co, R.CF.SO4);
fprintf('  Limiting current    : i_app %.0f A/m^2 vs i_lim(binding) %.0f A/m^2  (%.0f%% of limit)  %s\n', ...
    p.i_app, R.i_lim_binding, R.i_app_over_lim*100, ternary_str(R.ilim_ok,'OK','** EXCEEDS MARGIN **'));
fprintf('  Stack voltage/power : %.1f V   %.2f kW   SEC %.2f kWh/kg-metal (%.2f kWh/m^3)\n', ...
    R.V_stack, R.P_stack/1000, R.SEC_kWh_per_kgMetal, R.SEC_kWh_per_m3);
fprintf('  Current accounted   : cations(Ni+Co) %.1f%%, anion(SO4) %.1f%% of total charge\n', ...
    R.frac_current_tracked_cation*100, R.frac_current_tracked_anion*100);
fprintf('  (remainder carried by untracked ions — Ca/Mg/Na/H+ — not a mass-balance error)\n');
fprintf('══════════════════════════════════════════════\n\n');

end

function s = ternary_str(cond, a, b)
if cond; s = a; else; s = b; end
end


%% ═══════════════════════════════════════════════════════════════════════
%%  DASHBOARD  (2×2)
%% ═══════════════════════════════════════════════════════════════════════
function plotDashboard_ED(R, p)

save_dir = fullfile(pwd, 'ElectrodialysisMetalRecovery_Plots');
if ~exist(save_dir,'dir'); mkdir(save_dir); end

hDash = figure('Name','Electrodialysis — Metal Recovery Stack Dashboard','Position',[120 80 1500 940]);
tl = tiledlayout(hDash, 2, 2, 'TileSpacing','compact', 'Padding','compact');

draw_EP1_diluate(nexttile(tl), R, p);
draw_EP2_concentrate(nexttile(tl), R, p);
draw_EP3_ilim(nexttile(tl), R, p);
draw_EP4_energy(nexttile(tl), R, p);

title(tl, { ...
    sprintf('Electrodialysis — Metal Recovery Stack   |   N_{cp}=%d   |   i_{app}=%.0f A/m^2   |   I=%.1f A', ...
        p.N_cp, p.i_app, p.I), ...
    sprintf('Ni %.0f%% recovered, Co %.0f%% recovered   |   SEC %.2f kWh/kg-metal', ...
        R.rem_frac.Ni*100, R.rem_frac.Co*100, R.SEC_kWh_per_kgMetal) }, ...
    'FontSize',12, 'FontWeight','bold');

exportgraphics(hDash, fullfile(save_dir,'ED_Dashboard_Full.png'),'Resolution',150);

panels = {
    'EP1_Diluate',    @(ax) draw_EP1_diluate(ax, R, p);
    'EP2_Concentrate',@(ax) draw_EP2_concentrate(ax, R, p);
    'EP3_ILim',       @(ax) draw_EP3_ilim(ax, R, p);
    'EP4_Energy',     @(ax) draw_EP4_energy(ax, R, p);
};
for ii = 1:4
    hf = figure('Visible','off','Position',[0 0 720 540]);
    panels{ii,2}(axes());
    exportgraphics(hf, fullfile(save_dir, sprintf('ED_%s.png', panels{ii,1})),'Resolution',300);
    close(hf);
end
fprintf('  All Electrodialysis panels saved to: %s\n', save_dir);

end


%% ─── Panel functions ─────────────────────────────────────────────────────
function draw_EP1_diluate(ax, R, p)
axes(ax);
vals_in  = [p.in.Ni*p.MW.Ni, p.in.Co*p.MW.Co, p.in.SO4*p.MW.SO4/10] * 1000;  % SO4/10 for scale
vals_out = [R.dil.Ni*p.MW.Ni, R.dil.Co*p.MW.Co, R.dil.SO4*p.MW.SO4/10] * 1000;
b = bar([vals_in; vals_out]','grouped');
b(1).FaceColor = [0.55 0.65 0.85]; b(2).FaceColor = [0.10 0.30 0.80];
set(ax,'XTickLabel',{'Ni','Co','SO_4 (/10)'})
ylabel('mg/L'); legend({'in','diluate out'},'Location','best','FontSize',7)
title({'EP1: Diluate In vs Out', ...
    sprintf('Ni %.1f%%, Co %.1f%%, SO_4 %.1f%% removed', R.rem_frac.Ni*100, R.rem_frac.Co*100, R.rem_frac.SO4*100)}, ...
    'FontSize',9)
grid on
end

function draw_EP2_concentrate(ax, R, p)
axes(ax);
CFs = [R.CF.Ni, R.CF.Co, R.CF.SO4];
b = bar(CFs,'FaceColor','flat');
b.CData = [0.10 0.30 0.80; 0.85 0.10 0.10; 0.20 0.60 0.20];
set(ax,'XTickLabel',{'Ni','Co','SO_4'})
yline(p.CF_target, '--k', sprintf('target %.1fx',p.CF_target), 'FontSize',7);
ylabel('Concentration factor (x feed)')
title({'EP2: Concentrate Loop Buildup', sprintf('Q_c = %.1f L/h', R.Q_c)}, 'FontSize',9)
for i=1:3
    text(i, CFs(i)+max(CFs)*0.03, sprintf('%.1fx',CFs(i)), ...
        'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Parent',ax)
end
grid on
end

function draw_EP3_ilim(ax, R, p)
axes(ax);
vals = [p.i_app, R.i_lim.Ni, R.i_lim.Co];
b = bar(vals,'FaceColor','flat');
b.CData = [0.85 0.55 0.10; 0.10 0.30 0.80; 0.85 0.10 0.10];
set(ax,'XTickLabel',{'i_{app}','i_{lim}(Ni)','i_{lim}(Co)'})
ylabel('A/m^2')
statusTxt = ternary_str(R.ilim_ok, 'within margin', '** EXCEEDS MARGIN **');
title({'EP3: Limiting Current Check', ...
    sprintf('i_{app} = %.0f%% of binding limit — %s', R.i_app_over_lim*100, statusTxt)}, 'FontSize',9)
grid on
end

function draw_EP4_energy(ax, R, p)
axes(ax);
V_ohmic  = p.N_cp * p.i_app * p.R_area;
V_junct  = p.N_cp * p.V_junct;
b = bar([V_junct, V_ohmic],'FaceColor','flat');
b.CData = [0.40 0.20 0.60; 0.85 0.55 0.10];
set(ax,'XTickLabel',{'Membrane/junction','Ohmic (i \times R)'})
ylabel('V (of stack total)')
title({'EP4: Stack Voltage Breakdown', ...
    sprintf('V=%.1fV  P=%.2fkW  SEC=%.2f kWh/kg-metal', R.V_stack, R.P_stack/1000, R.SEC_kWh_per_kgMetal)}, ...
    'FontSize',9)
grid on
end
