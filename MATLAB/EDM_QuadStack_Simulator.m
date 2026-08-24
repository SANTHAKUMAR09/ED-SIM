function R = EDM_QuadStack_Simulator(cfg_ov, make_plots)
%% ╔══════════════════════════════════════════════════════════════════╗
%% ║   EDM Quad Stack — MATLAB EQUIVALENT of WebUI/edm_stack_simulator.html ║
%% ║   Single-unit (N=1), corrected termination: D2/CEM at the cathode  ║
%% ║   end, C1/CEM at the anode end -- keeps metal cations out of the   ║
%% ║   cathode rinse and Cl- out of the anode rinse.                    ║
%% ╠══════════════════════════════════════════════════════════════════╣
%% ║  This is a PORT, not a re-derivation: same explicit-Euler stepping ║
%% ║  algorithm, same 8 equations, same default numbers as the web tool ║
%% ║  -- given identical inputs, results match the web tool's history   ║
%% ║  to within floating-point precision. If you change the physics in  ║
%% ║  one, change it in the other; they are meant to stay in lockstep.  ║
%% ╠══════════════════════════════════════════════════════════════════╣
%% ║  COMPARTMENTS (order): cathode, D2, C2, D1(Feed), C1, anode         ║
%% ║  MEMBRANES (order, between consecutive compartments): CEM,AEM,CEM,AEM,CEM ║
%% ╠══════════════════════════════════════════════════════════════════╣
%% ║  EQUATIONS (see WebUI page's own "Model equations" section for the ║
%% ║  full narrative explanation of each):                              ║
%% ║   1  eta_i    = |z_i|*lambda_i*C_i / sum_j(|z_j|*lambda_j*C_j)      ║
%% ║   2  J_i      = eta_i*I / (|z_i|*F)                                 ║
%% ║   3  dC_i     = (J_in - J_out)*dt / V_compartment                   ║
%% ║   4  d[OH-]/dt|cathode = I/F ,  d[H+]/dt|anode = I/F                ║
%% ║   5  kappa = 0.001*sum(lambda_i*|z_i|*C_i) ; V = I*R_total + Vover + Vwatersplit ║
%% ║   6  pH = -log10(C_H,eff)                                           ║
%% ║   7  eta_Co(t) = min(n_Co,C2(t)/(I*t/(z_Co*F)), 1)*100%             ║
%% ║   8  I_lim = A*sum(|z_i|*F*D_i*C_i/(1000*delta)) ; D_i = R*T*lambda_i/(|z_i|*F^2) ║
%% ╚══════════════════════════════════════════════════════════════════╝

if nargin < 1; cfg_ov = []; end
if nargin < 2 || isempty(make_plots); make_plots = true; end

p = defineParameters_EDM(cfg_ov);
R = runEDM(p);
if make_plots; plotDashboard_EDM(R, p); end
R.p = p;

end


%% ═══════════════════════════════════════════════════════════════════════
%%  PARAMETERS — species table, compartments, stack design, all defaults
%%  matching WebUI/edm_stack_simulator.html exactly
%% ═══════════════════════════════════════════════════════════════════════
function p = defineParameters_EDM(cfg_ov)

p.F    = 96485;    % C/mol
p.Kw   = 1e-14;     % water ion product, 25 degC
p.Rgas = 8.314;     % J/mol/K
p.T    = 298.15;    % K (25 degC -- matches the 25 degC lambda values below)

%── Species: id, charge z, limiting equivalent ionic conductivity lambda ──
%   (S*cm^2/eq, 25degC literature values). D derived via Nernst-Einstein,
%   D_i = R*T*lambda_i / (|z_i|*F^2) -- verified against literature: gives
%   D_Na+ = 1.33e-5 and D_Ca2+ = 7.92e-6 cm^2/s, matching tabulated values
%   to 3 sig figs. Order: Co Ca Na SO4 Cl H OH.
p.sp_id     = {'Co','Ca','Na','SO4','Cl','H','OH'};
p.sp_label  = {'Co2+','Ca2+','Na+','SO4^2-','Cl-','H+','OH-'};
p.sp_z      = [ 2,    2,    1,   -2,    -1,   1,    -1  ];
p.sp_lambda = [53.0, 59.5, 50.1, 80.0, 76.3, 349.8, 198.0];
p.sp_D      = (p.Rgas*p.T .* p.sp_lambda) ./ (abs(p.sp_z) * p.F^2);  % cm^2/s
p.nSp       = numel(p.sp_id);
p.iCo=1; p.iCa=2; p.iNa=3; p.iSO4=4; p.iCl=5; p.iH=6; p.iOH=7;
p.conservedIdx = [p.iCo p.iCa p.iNa p.iSO4 p.iCl];   % H/OH are source/sink species, excluded

%── Compartments: physical left-to-right order (cathode -> anode) ────────
p.comp_id    = {'cathode','d2','c2','d1','c1','anode'};
p.comp_label = {'Cathode rinse','Sol.1 -> D2','Sol.2 -> C2','Feed -> D1','Sol.3 -> C1','Anode rinse'};
p.nComp = 6;
% Default initial concentrations (mol/L), rows=compartments (order above),
% cols=species (Co Ca Na SO4 Cl H OH) -- matches WebUI defaults exactly.
p.conc0 = [ ...
    0,      0,      0.500,  0.250,  0,      1e-7, 1e-7;   % cathode rinse
    0,      0,      0.500,  0,      0.500,  1e-7, 1e-7;   % D2 (0.5 M NaCl donor)
    0,      0,      0.050,  0,      0.050,  1e-7, 1e-7;   % C2 (0.05 M NaCl seed)
    0.0500, 0.0050, 0.1082, 0.0541, 0.1100, 1e-7, 1e-7;   % D1 Feed (CoCl2/Na2SO4/CaCl2)
    0,      0,      0.100,  0.050,  0,      1e-7, 1e-7;   % C1 (0.05 M Na2SO4 seed)
    0,      0,      0.500,  0.250,  0,      1e-7, 1e-7 ]; % anode rinse
p.vol_L = 0.250 * ones(p.nComp,1);   % L, 250 mL each

%── Membranes: fixed by the corrected topology, between consecutive comps ─
p.membrane_types = {'CEM','AEM','CEM','AEM','CEM'};
p.nMem = 5;

%── Stack design / operating point ────────────────────────────────────────
p.area_cm2   = 10;     % cm^2, membrane area
p.I_A        = 0.5;    % A, applied constant current
p.duration_h = 4;      % h, run duration
p.dt_min     = 1;      % min, timestep

%── Resistance / overpotential placeholders (flagged, not measured) ──────
p.AREAL_R_MEMBRANE = 5;      % ohm*cm^2, per membrane
p.GAP_CM           = 0.5;    % cm, compartment spacer thickness
p.V_OVERPOTENTIAL  = 1.5;    % V, fixed electrode overpotential
p.KAPPA_WATER_FLOOR = 5.5e-8; % S/cm, pure water's own conductivity at 25degC --
                               % physical floor so a depleted compartment gets a
                               % very HIGH but finite resistance, not zero.
p.V_STACK_CEILING  = 200;    % V, safety-net ceiling (not the primary voltage-cap
                               % mechanism -- eq. 8's water-splitting overpotential
                               % is what should realistically bound voltage).

%── Over-limiting current / water-splitting at interior membranes (eq. 8) ─
p.DELTA_CM = 50e-4;   % cm (50 micron), Nernst diffusion boundary layer thickness --
                        % typical for a turbulence-promoting spacer, placeholder.
p.V_WATERSPLIT_OVERPOTENTIAL = 0.8;   % V, per over-limiting membrane

%── Apply user overrides (any field above may be overridden) ─────────────
if ~isempty(cfg_ov)
    fns = fieldnames(cfg_ov);
    for k = 1:numel(fns); p.(fns{k}) = cfg_ov.(fns{k}); end
end

p.dt_sec   = p.dt_min * 60;
p.tEnd_sec = p.duration_h * 3600;

%── Design summary ────────────────────────────────────────────────────────
fprintf('══════════════════════════════════════════════\n');
fprintf('  EDM Quad Stack — Design\n');
fprintf('══════════════════════════════════════════════\n');
fprintf('  Membrane area      : %.1f cm^2\n', p.area_cm2);
fprintf('  Current (constant) : %.3f A\n', p.I_A);
fprintf('  Duration / dt      : %.1f h / %.1f min\n', p.duration_h, p.dt_min);
fprintf('  Compartments       : %s\n', strjoin(p.comp_label, ' | '));
fprintf('  Membranes          : %s\n', strjoin(p.membrane_types, ' - '));
fprintf('══════════════════════════════════════════════\n\n');

end


%% ═══════════════════════════════════════════════════════════════════════
%%  ONE EXPLICIT-EULER STEP (mirrors WebUI stepOnce() exactly)
%% ═══════════════════════════════════════════════════════════════════════
function [conc, overLimitFlags] = stepOnce_EDM(conc, p, I, dt)
% conc: nComp x nSp matrix (mol/L), mutated and returned
molBefore = conc .* p.vol_L;   % nComp x nSp, mol
overLimitFlags = false(1, p.nMem);

for k = 1:p.nMem
    type = p.membrane_types{k};
    left = k; right = k+1;
    if strcmp(type, 'CEM')
        source = right; destIdx = left;    % cations: right -> left (toward cathode)
        wantSign = 1;
    else
        source = left; destIdx = right;     % AEM, anions: left -> right (toward anode)
        wantSign = -1;
    end

    carriers = find(sign(p.sp_z) == wantSign);
    srcConc  = max(conc(source, carriers), 0);
    denom    = sum(abs(p.sp_z(carriers)) .* p.sp_lambda(carriers) .* srcConc);

    % Limiting current (eq. 8): max current this membrane's own ion
    % inventory can sustain via diffusion to the surface. denom and I_lim
    % are zero under exactly the same condition (all carrier concentrations
    % zero), so I_ionic naturally goes to 0 and the full current routes
    % through the water-splitting branch below -- no separate guard needed.
    I_lim = p.area_cm2 * sum(abs(p.sp_z(carriers)) .* p.F .* p.sp_D(carriers) .* srcConc / 1000 / p.DELTA_CM);

    overLimit = I > I_lim;
    overLimitFlags(k) = overLimit;
    I_ionic = min(I, I_lim);

    if denom > 1e-30 && I_ionic > 1e-30
        share       = (abs(p.sp_z(carriers)) .* p.sp_lambda(carriers) .* srcConc) / denom;
        fluxMolPerS = share * I_ionic ./ (abs(p.sp_z(carriers)) * p.F);
        dMol        = fluxMolPerS * dt;
        % Cap to what's actually available in the source THIS step (avoids
        % an Euler-overshoot mass-conservation leak when a species fully
        % depletes mid-step).
        dMol = min(dMol, max(molBefore(source, carriers), 0));
        molBefore(source, carriers) = molBefore(source, carriers) - dMol;
        molBefore(destIdx, carriers) = molBefore(destIdx, carriers) + dMol;
    end

    if overLimit
        % Excess current beyond I_lim is carried by water splitting right at
        % this membrane. CEM: H+ (a cation) crosses through to dest, OH-
        % stays in source. AEM: OH- crosses to dest, H+ stays in source.
        dMolWS = (I - I_lim) / p.F * dt;
        if strcmp(type, 'CEM')
            molBefore(destIdx, p.iH)  = molBefore(destIdx, p.iH)  + dMolWS;
            molBefore(source,  p.iOH) = molBefore(source,  p.iOH) + dMolWS;
        else
            molBefore(destIdx, p.iOH) = molBefore(destIdx, p.iOH) + dMolWS;
            molBefore(source,  p.iH)  = molBefore(source,  p.iH)  + dMolWS;
        end
    end
end

% Electrode reactions (Faradaic source terms):
% cathode (compartment 1): 2H2O + 2e- -> H2 + 2OH-  =>  d(OH-)/dt = I/F
molBefore(1, p.iOH) = molBefore(1, p.iOH) + (I/p.F)*dt;
% anode (compartment nComp): 2H2O -> O2 + 4H+ + 4e-  =>  d(H+)/dt = I/F
molBefore(p.nComp, p.iH) = molBefore(p.nComp, p.iH) + (I/p.F)*dt;

% write back concentrations, clamp at 0
conc = max(molBefore ./ p.vol_L, 0);

end


%% ═══════════════════════════════════════════════════════════════════════
%%  DERIVED QUANTITIES
%% ═══════════════════════════════════════════════════════════════════════
function kappa = conductivity_EDM(concRow, p)
% kappa [S/cm] = 0.001 * sum( lambda_i * |z_i| * C_i[mol/L] )
kappa = 0.001 * sum(p.sp_lambda .* abs(p.sp_z) .* max(concRow,0));
end

function V = stackVoltage_EDM(conc, overLimitCount, I, p)
Rtot = 0;
for i = 1:p.nComp
    kappa = max(conductivity_EDM(conc(i,:), p), p.KAPPA_WATER_FLOOR);
    Rtot = Rtot + p.GAP_CM/(kappa*p.area_cm2);
end
Rtot = Rtot + p.nMem * (p.AREAL_R_MEMBRANE/p.area_cm2);
V_watersplit = overLimitCount * p.V_WATERSPLIT_OVERPOTENTIAL;
V = min(I*Rtot + p.V_OVERPOTENTIAL + V_watersplit, p.V_STACK_CEILING);
end

function pH = phOf_EDM(cH, cOH, Kw)
net = cH - cOH;
if net > 1e-14
    cHeff = net;
elseif net < -1e-14
    cHeff = Kw / (-net);
else
    cHeff = 1e-7;   % no net H+/OH- signal -- default to neutral, not the
                     % misleading pH 0 a naive Kw/1e-14 fallback would give
end
pH = -log10(max(cHeff, 1e-14));
end


%% ═══════════════════════════════════════════════════════════════════════
%%  RUN — time-steps from t=0 to duration, records full history
%% ═══════════════════════════════════════════════════════════════════════
function R = runEDM(p)

nSteps = round(p.tEnd_sec / p.dt_sec) + 1;   % +1 for the t=0 snapshot
conc = p.conc0;

R.t_sec = zeros(nSteps,1);
R.conc  = zeros(nSteps, p.nComp, p.nSp);
R.pH    = zeros(nSteps, p.nComp);
R.V     = zeros(nSteps,1);
R.overLimitFlags = false(nSteps, p.nMem);
R.overLimitCount = zeros(nSteps,1);

initialTotals = squeeze(sum(p.conc0 .* p.vol_L, 1));   % 1 x nSp, mol

% t=0 snapshot (before any stepping -- matches WebUI's initial history push)
R.t_sec(1) = 0;
R.conc(1,:,:) = conc;
for c = 1:p.nComp; R.pH(1,c) = phOf_EDM(conc(c,p.iH), conc(c,p.iOH), p.Kw); end
R.V(1) = stackVoltage_EDM(conc, 0, p.I_A, p);

t = 0;
for step = 2:nSteps
    [conc, overLimitFlags] = stepOnce_EDM(conc, p, p.I_A, p.dt_sec);
    t = t + p.dt_sec;

    R.t_sec(step) = t;
    R.conc(step,:,:) = conc;
    for c = 1:p.nComp; R.pH(step,c) = phOf_EDM(conc(c,p.iH), conc(c,p.iOH), p.Kw); end
    R.overLimitFlags(step,:) = overLimitFlags;
    R.overLimitCount(step) = sum(overLimitFlags);
    R.V(step) = stackVoltage_EDM(conc, R.overLimitCount(step), p.I_A, p);
end

R.t_min = R.t_sec / 60;

%── Cumulative Co2+ current efficiency vs time (eq. 7) ────────────────────
c2Idx = find(strcmp(p.comp_id,'c2'));
R.coEfficiency_pct = zeros(nSteps,1);
for step = 1:nSteps
    coInC2 = R.conc(step, c2Idx, p.iCo) * p.vol_L(c2Idx);
    theoreticalMax = (p.I_A * R.t_sec(step)) / (p.sp_z(p.iCo) * p.F);
    if theoreticalMax > 1e-12
        R.coEfficiency_pct(step) = min(coInC2/theoreticalMax*100, 100);
    end
end

%── Mass balance closure vs time (conserved species only) ────────────────
R.massBalance_pct = zeros(nSteps, numel(p.conservedIdx));
for step = 1:nSteps
    for j = 1:numel(p.conservedIdx)
        spIdx = p.conservedIdx(j);
        tot = sum(squeeze(R.conc(step,:,spIdx))' .* p.vol_L);
        if initialTotals(spIdx) > 1e-12
            R.massBalance_pct(step,j) = tot/initialTotals(spIdx)*100;
        else
            R.massBalance_pct(step,j) = 100;
        end
    end
end
R.massBalanceClosure_worst_pct = R.massBalance_pct(end, ...
    find(max(abs(R.massBalance_pct(end,:)-100)) == abs(R.massBalance_pct(end,:)-100), 1));

%── Final KPIs ─────────────────────────────────────────────────────────
d1Idx = find(strcmp(p.comp_id,'d1'));
c1Idx = find(strcmp(p.comp_id,'c1'));
cathodeIdx = find(strcmp(p.comp_id,'cathode'));

R.KPI.massBalanceClosure_pct = R.massBalanceClosure_worst_pct;
R.KPI.coRecovered_pct = min(R.conc(end,c2Idx,p.iCo)*p.vol_L(c2Idx) / max(p.conc0(d1Idx,p.iCo)*p.vol_L(d1Idx),1e-12) * 100, 100);
c1 = squeeze(R.conc(end,c1Idx,:))';
wantedEq = c1(p.iNa)*1 + c1(p.iSO4)*2;
contamEq = c1(p.iH)*1 + c1(p.iCl)*1;
R.KPI.c1Purity_pct = 100 * wantedEq / max(wantedEq+contamEq, 1e-12);
R.KPI.cathodeNaGained_mmol = (R.conc(end,cathodeIdx,p.iNa) - p.conc0(cathodeIdx,p.iNa)) * p.vol_L(cathodeIdx) * 1000;
R.KPI.currentEfficiencyCo_pct = R.coEfficiency_pct(end);
R.KPI.membranesOverLimiting = R.overLimitCount(end);
R.KPI.membranesTotal = p.nMem;

fprintf('══════════════════════════════════════════════\n');
fprintf('  EDM Quad Stack — Results (t = %.0f h)\n', p.duration_h);
fprintf('══════════════════════════════════════════════\n');
fprintf('  Mass balance closure   : %.2f%%\n', R.KPI.massBalanceClosure_pct);
fprintf('  Co2+ recovered -> C2   : %.1f%%\n', R.KPI.coRecovered_pct);
fprintf('  C1 product purity      : %.1f%%  (Na2SO4 vs H+/Cl- leak)\n', R.KPI.c1Purity_pct);
fprintf('  Cathode Na+ gained     : %.2f mmol (sacrificial, by design)\n', R.KPI.cathodeNaGained_mmol);
fprintf('  Current efficiency, Co : %.1f%%\n', R.KPI.currentEfficiencyCo_pct);
fprintf('  Membranes over-limiting: %d / %d  (water-splitting active, eq. 8)\n', R.KPI.membranesOverLimiting, R.KPI.membranesTotal);
fprintf('  Final stack voltage    : %.2f V\n', R.V(end));
fprintf('══════════════════════════════════════════════\n\n');

end


%% ═══════════════════════════════════════════════════════════════════════
%%  DASHBOARD (2×3)
%% ═══════════════════════════════════════════════════════════════════════
function plotDashboard_EDM(R, p)

save_dir = fullfile(pwd, 'EDM_QuadStack_Plots');
if ~exist(save_dir,'dir'); mkdir(save_dir); end

hDash = figure('Name','EDM Quad Stack Dashboard','Position',[80 60 1600 900]);
tl = tiledlayout(hDash, 2, 3, 'TileSpacing','compact', 'Padding','compact');

draw_EQ1_conc(nexttile(tl), R, p);
draw_EQ2_balance(nexttile(tl), R, p);
draw_EQ3_voltage(nexttile(tl), R, p);
draw_EQ4_ph(nexttile(tl), R, p);
draw_EQ5_efficiency(nexttile(tl), R, p);
draw_EQ6_overlimit(nexttile(tl), R, p);

title(tl, { ...
    sprintf('EDM Quad Stack — I=%.3f A, A=%.0f cm^2, %.0f h run', p.I_A, p.area_cm2, p.duration_h), ...
    sprintf('Mass balance %.2f%%  |  Co recovered %.1f%%  |  C1 purity %.1f%%  |  %d/%d membranes over-limiting', ...
        R.KPI.massBalanceClosure_pct, R.KPI.coRecovered_pct, R.KPI.c1Purity_pct, R.KPI.membranesOverLimiting, R.KPI.membranesTotal) }, ...
    'FontSize',12, 'FontWeight','bold');

exportgraphics(hDash, fullfile(save_dir,'EDM_Dashboard_Full.png'),'Resolution',150);
fprintf('  Dashboard saved to: %s\n', save_dir);

end

function draw_EQ1_conc(ax, R, p)
axes(ax); hold(ax,'on');
d1Idx = find(strcmp(p.comp_id,'d1')); c2Idx = find(strcmp(p.comp_id,'c2')); c1Idx = find(strcmp(p.comp_id,'c1'));
plot(ax, R.t_min, R.conc(:,d1Idx,p.iCo), '-', 'Color',[0.02 0.42 0.63],'LineWidth',2, 'DisplayName','D1 (Feed) Co^{2+}');
plot(ax, R.t_min, R.conc(:,c2Idx,p.iCo), '-', 'Color',[0.12 0.53 0.90],'LineWidth',2, 'DisplayName','C2 Co^{2+}');
plot(ax, R.t_min, R.conc(:,d1Idx,p.iSO4),'-', 'Color',[0.49 0.23 0.93],'LineWidth',2, 'DisplayName','D1 (Feed) SO_4^{2-}');
plot(ax, R.t_min, R.conc(:,c1Idx,p.iSO4),'-', 'Color',[0.31 0.13 0.60],'LineWidth',2, 'DisplayName','C1 SO_4^{2-}');
xlabel(ax,'Time (min)'); ylabel(ax,'Concentration (mol/L)');
title(ax,'EQ1: Concentrations (key species)','FontSize',10);
legend(ax,'Location','best','FontSize',7); grid(ax,'on');
end

function draw_EQ2_balance(ax, R, p)
axes(ax); hold(ax,'on');
cols = lines(numel(p.conservedIdx));
for j = 1:numel(p.conservedIdx)
    plot(ax, R.t_min, R.massBalance_pct(:,j), '-', 'Color', cols(j,:), 'LineWidth',1.8, ...
        'DisplayName', p.sp_label{p.conservedIdx(j)});
end
yline(ax,100,'--k','HandleVisibility','off');
xlabel(ax,'Time (min)'); ylabel(ax,'Mass balance closure (%)');
ylim(ax,[95 105]);
title(ax,'EQ2: Mass Balance Closure','FontSize',10);
legend(ax,'Location','best','FontSize',7); grid(ax,'on');
end

function draw_EQ3_voltage(ax, R, p)
axes(ax);
plot(ax, R.t_min, R.V, '-', 'Color',[0.27 0.22 0.79],'LineWidth',2.2);
xlabel(ax,'Time (min)'); ylabel(ax,'Stack voltage (V)');
title(ax, {'EQ3: Voltage & Current (eq. 5, 8)', sprintf('final %.2f V @ I=%.3f A const.', R.V(end), p.I_A)}, 'FontSize',10);
grid(ax,'on');
end

function draw_EQ4_ph(ax, R, p)
axes(ax); hold(ax,'on');
cols = lines(p.nComp);
for c = 1:p.nComp
    plot(ax, R.t_min, R.pH(:,c), '-', 'Color', cols(c,:), 'LineWidth',1.6, 'DisplayName', p.comp_label{c});
end
ylim(ax,[0 14]);
xlabel(ax,'Time (min)'); ylabel(ax,'pH');
title(ax,'EQ4: pH (eq. 6)','FontSize',10);
legend(ax,'Location','best','FontSize',6); grid(ax,'on');
end

function draw_EQ5_efficiency(ax, R, p)
axes(ax);
plot(ax, R.t_min, R.coEfficiency_pct, '-', 'Color',[0.02 0.42 0.63],'LineWidth',2.2);
ylim(ax,[0 100]);
xlabel(ax,'Time (min)'); ylabel(ax,'Co^{2+} current efficiency, cumulative (%)');
title(ax, {'EQ5: Current Efficiency (eq. 7)', sprintf('final %.1f%%', R.coEfficiency_pct(end))}, 'FontSize',10);
grid(ax,'on');
end

function draw_EQ6_overlimit(ax, R, p)
axes(ax);
plot(ax, R.t_min, R.overLimitCount, '-o', 'Color',[0.71 0.46 0.04],'LineWidth',1.8,'MarkerSize',3,'MarkerFaceColor',[0.71 0.46 0.04]);
ylim(ax,[-0.3 p.nMem+0.3]); yticks(ax, 0:p.nMem);
xlabel(ax,'Time (min)'); ylabel(ax,'Membranes over-limiting (of 5)');
title(ax,'EQ6: Over-Limiting Membranes (eq. 8)','FontSize',10);
grid(ax,'on');
end
