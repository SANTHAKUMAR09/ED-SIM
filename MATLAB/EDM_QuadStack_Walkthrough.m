%% ═══════════════════════════════════════════════════════════════════════
%%  EDM QUAD STACK — BEGINNER WALKTHROUGH SCRIPT (full 27-equation model)
%% ═══════════════════════════════════════════════════════════════════════
%
%  This is the SAME simulation as EDM_QuadStack_Simulator.m (same physics,
%  same numbers, same results for the default single-unit run) but written
%  as ONE FLAT SCRIPT instead of a file full of functions. Nothing here is
%  hidden inside a function call -- every line runs top to bottom, in
%  order, and every variable it creates stays visible in your Workspace
%  panel afterwards so you can click on it and look at the actual numbers.
%
%  WHAT CHANGED SINCE THE FIRST VERSION OF THIS FILE: the reference model
%  grew from 8 equations / 2 metals to 27 equations / 15 species (10
%  metals). This version has grown to match it -- same 15 species, same
%  precipitation/scaling chemistry, same water transport, same membrane
%  fouling, same configurable initial pH. The two things this file does
%  NOT reproduce from the reference model are (a) N > 1 repeating units
%  (this walkthrough always runs N = 1, i.e. a single quad-stack cell) and
%  (b) the generic "pass in any config struct" override mechanism -- this
%  file has ONE fixed configuration, the same default the reference model
%  itself ships with, because the point here is to see every number, not
%  to build a reusable tool. For deeper rationale behind any constant
%  below, search for it by name in EDM_QuadStack_Simulator.m -- its own
%  comments go into much more depth than this file repeats.
%
%  HOW TO USE THIS FILE AS A LEARNING TOOL:
%   1. Open it in the MATLAB Editor.
%   2. Each "%%" line below starts a new SECTION (MATLAB highlights the
%      section you're in with a pale yellow background).
%   3. Click anywhere inside a section, then press Ctrl+Enter (Cmd+Enter on
%      Mac). MATLAB runs just THAT section and stops -- it does not run the
%      whole file. Look at the Workspace panel to see what got created.
%   4. Work your way down, section by section. By the time you reach the
%      big "MAIN LOOP" section you'll already recognize every variable it
%      uses, because you built them yourself in the sections above.
%   5. To run the ENTIRE script at once (like the function version), press
%      the green "Run" button, or type its filename in the Command Window.
%
%  THE PHYSICAL PICTURE (read this before touching any code):
%   Picture 6 small tanks of salty water in a row, with a membrane (a thin
%   sheet that only lets certain ions pass through) between each pair of
%   neighboring tanks. A DC electrical current is pushed through the whole
%   row from one end to the other. Positively-charged ions (metals, Na+)
%   drift toward one end (the cathode); negatively-charged ions (Cl-,
%   SO4^2-) drift toward the other end (the anode). Each membrane only lets
%   ONE type of charge through, so ions pile up in some tanks and get
%   stripped out of others -- that's how this stack separates and
%   concentrates the metals.
%
%   Tank layout, left (cathode, negative electrode) to right (anode, positive):
%     Tank 1: Cathode rinse   (just a rinse solution, sits at the electrode)
%     Tank 2: D2               (0.5 M NaCl "donor" solution)
%     Tank 3: C2                (starts dilute -- this is where metals end up)
%     Tank 4: D1 = Feed          (the actual metal-bearing solution)
%     Tank 5: C1                 (starts dilute -- this is where SO4^2- ends up)
%     Tank 6: Anode rinse        (just a rinse solution, sits at the electrode)
%
%   Membrane layout (5 membranes, one between each pair of neighboring tanks):
%     Membrane 1 (Tank 1 | Tank 2): CEM  (passes + ions only)
%     Membrane 2 (Tank 2 | Tank 3): AEM  (passes - ions only)
%     Membrane 3 (Tank 3 | Tank 4): CEM
%     Membrane 4 (Tank 4 | Tank 5): AEM
%     Membrane 5 (Tank 5 | Tank 6): CEM
%
%   This specific left-right arrangement keeps metals away from Tank 1 (so
%   they never plate onto the cathode) and keeps chloride away from Tank 6
%   (so no chlorine gas forms at the anode).
%
%   ON TOP OF THAT SIMPLE PICTURE, this version adds: real membranes are not
%   perfectly selective (some current leaks the "wrong" way); ions have to
%   physically diffuse to a membrane's surface, so there's a speed limit on
%   how much current any membrane can carry before it starts splitting
%   water instead; concentrated solutions conduct worse than dilute theory
%   predicts; water itself moves between tanks (dragged along with ions,
%   and drawn osmotically toward the concentrated side); suspended solids
%   deposit on membranes and foul them; and if a tank gets supersaturated
%   in some metal hydroxide (which happens here -- this stack drives some
%   tanks to pH 1 and others to pH 13+), that metal can precipitate out as
%   solid scale, taking it out of solution and off the current's books.

clear;      % wipe the Workspace so you start from a clean slate
close all;  % close any old figure windows
clc;        % clear the Command Window text


%% SECTION 1 — Physical constants
% These four numbers never change; they come from physics/chemistry, not
% from anything about OUR stack.

F    = 96485;     % Faraday's constant, Coulombs per mole of charge --
                   % "how much electric charge one mole of singly-charged
                   % ions carries." Converts between "current" and "moles
                   % of ions moved" everywhere in this script.
Kw   = 1e-14;      % water self-ionization constant, C_H * C_OH at equilibrium
Rgas = 8.314;      % gas constant, J/mol/K
T    = 298.15;     % temperature, Kelvin (25 degC -- what every lambda/Ksp
                   % value below was measured at)


%% SECTION 2 — The 15 ions we're tracking (10 metals + Na, SO4, Cl, H, OH)
% Every ion in this simulation is one of these 15. We give each one:
%   - a short id and a display label
%   - z = its electric charge
%   - lambda = its "limiting equivalent ionic conductivity" (S*cm^2/eq,
%     25 degC literature value) -- bigger means it moves faster under the
%     same electric field. This is the EQUIVALENT form (per unit charge),
%     which is why the di/trivalent metals sit in the same 50-70 range as
%     Na+ -- do not substitute molar conductivities here.
%
% Oxidation states assumed (the dominant forms in an acidic leachate): Fe
% as Fe2+, Al and Cr as the trivalent hydrated cations. In reality, at the
% pH extremes this stack produces, Fe3+/Al3+/Cr3+ would hydrolyse and
% precipitate -- that specific chemistry isn't modelled, so treat trivalent
% transport results as an upper bound.
sp_id     = {'Co','Ni','Fe','Al','Cr','Mn','Mg','Cu','Zn','Ca','Na','SO4','Cl','H','OH'};
sp_label  = {'Co2+','Ni2+','Fe2+','Al3+','Cr3+','Mn2+','Mg2+','Cu2+','Zn2+','Ca2+','Na+','SO4^2-','Cl-','H+','OH-'};
sp_z      = [ 2,     2,     2,     3,     3,     2,     2,     2,     2,     2,     1,   -2,      -1,   1,    -1  ];
sp_lambda = [53.0,  49.6,  54.0,  61.0,  67.0,  53.5,  53.1,  53.6,  52.8,  59.5,  50.1, 80.0,   76.3, 349.8, 198.0];
n_species = numel(sp_id);

% Named indices -- position of each ion in the arrays above. We'll use
% these everywhere instead of remembering "species number 3 is iron".
iCo=1; iNi=2; iFe=3; iAl=4; iCr=5; iMn=6; iMg=7; iCu=8; iZn=9; iCa=10;
iNa=11; iSO4=12; iCl=13; iH=14; iOH=15;

metal_idx    = [iCo iNi iFe iAl iCr iMn iMg iCu iZn iCa];   % the 10 metals
% H+ and OH- are deliberately created/destroyed by the electrode and
% water-splitting reactions, so they're never expected to "balance" --
% everything else should. Both the running history and the final check use
% this same definition of "conserved."
conserved_idx = setdiff(1:n_species, [iH iOH]);

% Diffusivity D -- how fast each ion spreads out on its own, no electric
% field needed. Related to lambda by the Nernst-Einstein relation, computed
% once for all 15 ions in one line (units cm^2/s). This reproduces the
% textbook diffusivity of Na+ to 3 significant figures: 1.33e-5 cm^2/s.
sp_D = (Rgas .* T .* sp_lambda) ./ (abs(sp_z) .* F.^2);


%% SECTION 3 — The 6 tanks and what's in them at the start
tank_id    = {'cathode','d2','c2','d1','c1','anode'};
tank_name  = {'Cathode rinse','Sol.1 -> D2','Sol.2 -> C2','Feed -> D1','Sol.3 -> C1','Anode rinse'};
n_tanks    = 6;
i_cathode=1; i_d2=2; i_c2=3; i_d1=4; i_c1=5; i_anode=6;

% "conc" is a table: 6 tanks (rows) x 15 ions (columns), mol/L. This is the
% state that changes every timestep. H+/OH- are NOT set here -- Section 3b
% below sets them from an "Initial pH" per tank, same as the reference
% model's cfg.pH0.
conc = zeros(n_tanks, n_species);
conc(i_cathode, [iNa iSO4])            = [0.500 0.250];   % 0.25 M Na2SO4 rinse
conc(i_d2,      [iNa iCl])             = [0.500 0.500];   % 0.5 M NaCl donor
conc(i_c2,      [iNa iCl])             = [0.050 0.050];   % 0.05 M NaCl seed
conc(i_d1,      [iCa iNa iSO4 iCl])    = [0.0050 0.1082 0.0541 0.1100];  % background salts
conc(i_c1,      [iNa iSO4])            = [0.100 0.050];   % 0.05 M Na2SO4 seed
conc(i_anode,   [iNa iSO4])            = [0.500 0.250];   % 0.25 M Na2SO4 rinse

% Feed metals: dropped into D1 (Tank 4) on top of the background salts
% above. This placeholder ten-metal mixed-leachate composition is the same
% one the reference model and the web app both default to -- replace it
% with your own assay before drawing conclusions.
conc(i_d1, iCo) = 0.0500;
conc(i_d1, iNi) = 0.0200;
conc(i_d1, iMn) = 0.0100;
conc(i_d1, iCu) = 0.0080;
conc(i_d1, iZn) = 0.0060;
conc(i_d1, iMg) = 0.0050;
conc(i_d1, iCa) = 0.0050;
conc(i_d1, iFe) = 0.0040;
conc(i_d1, iAl) = 0.0020;
conc(i_d1, iCr) = 0.0010;

% Charge-balance the feed: ten metal cations at various charges were just
% added, so D1 now has a cation excess unless something balances it. Real
% leachates come with an anion attached to whatever acid dissolved the ore
% in the first place -- here that's chloride. Top up Cl- by however much
% is needed to bring D1's net charge back to zero.
netEq_d1 = sum(sp_z .* conc(i_d1,:));       % >0 = excess cation, eq/L
if netEq_d1 > 0
    conc(i_d1, iCl) = conc(i_d1, iCl) + netEq_d1 / abs(sp_z(iCl));
end

tank_volume_L = 0.250 * ones(n_tanks,1);   % L, all six start at 250 mL


%% SECTION 3b — Initial pH per tank (and the acid/base it takes to get there)
% Every tank previously started at a hard-wired pH 7. Now each tank can
% start at its own acidity -- set pH0 below, one value per tank, in the
% same cathode..anode order as everything else. [OH-] then follows from
% water equilibrium, [OH-] = Kw/[H+].
pH0 = 7 * ones(1, n_tanks);
% Try, e.g.:  pH0(i_d1) = 2;   (an acidic feed -- see what it does to Co/Ni
%                                recovery and to which metal hydroxides
%                                start precipitating in Section 13)

H0  = 10 .^ (-pH0);      % mol/L
OH0 = Kw ./ H0;          % mol/L, exactly consistent with Kw by construction
conc(:, iH)  = H0(:);
conc(:, iOH) = OH0(:);

% Acid/base is a dosed REAGENT, so it brings a counter-ion along: H2SO4
% adds SO4^2- (2 mol H+ per mol SO4^2-), NaOH adds Na+. Without this, a
% tank set to low pH would have a pile of "bare" H+ with nothing to
% balance its charge -- chemically impossible. This is not tied to one
% specific acid: change acid_anion below to 'Cl' to dose HCl instead.
acid_anion  = iSO4;   % which anion is added with excess H+ (H2SO4 here)
base_cation = iNa;    % which cation is added with excess OH- (NaOH)
pH_dosed_eqL = zeros(1, n_tanks);
for i = 1:n_tanks
    net = conc(i,iH) - conc(i,iOH);      % eq/L of excess acid(+) or base(-)
    if net > 1e-9
        conc(i,acid_anion) = conc(i,acid_anion) + net/abs(sp_z(acid_anion));
        pH_dosed_eqL(i) = net;
    elseif net < -1e-9
        conc(i,base_cation) = conc(i,base_cation) + (-net)/abs(sp_z(base_cation));
        pH_dosed_eqL(i) = net;
    end
end

% Remember the STARTING total of each CONSERVED ion (everything but H+/OH-),
% added up across all 6 tanks. We compare against this throughout the run
% as a mass-balance sanity check.
starting_total_moles = sum(conc .* tank_volume_L, 1);   % 1 x 15


%% SECTION 4 — The 5 membranes between the tanks
% membrane_type{m} sits between Tank m and Tank (m+1).
%   'CEM'      standard cation-exchange: passes all cations, no selectivity.
%   'AEM'      standard anion-exchange, the mirror image.
%   'Mono-CEM' monovalent-SELECTIVE cation-exchange: passes Na+, holds back
%              Co2+/Ca2+/etc (a thin surface layer electrostatically repels
%              multivalent cations). Use it to reject metals and bleed Na+.
%   'Mono-AEM' monovalent-selective anion-exchange: Cl- passes, SO4^2- held
%              back.
%   'BPM'      bipolar: carries NO salt. Under reverse bias it splits water
%              at its internal junction -- H+ toward the cathode side, OH-
%              toward the anode side -- so it makes acid/base in situ
%              instead of you dosing them. Costs an extra ~0.9 V per copy.
% This script's default is the plain CEM/AEM chain below -- try changing
% one entry to 'Mono-CEM' or 'BPM' and re-running to see what changes.
membrane_type = {'CEM','AEM','CEM','AEM','CEM'};
n_membranes   = 5;

% membranes(k,:) = [leftTank, rightTank, code, multiplicity]. code: 1 =
% cation-conducting (CEM/Mono-CEM), 0 = anion-conducting (AEM/Mono-AEM),
% 2 = bipolar. multiplicity is how many physical copies this row
% represents -- always 1 here (N = 1, a single quad-stack cell; the
% reference model's N_units > 1 tiles this same table with higher
% multiplicities for a multi-cell stack, not reproduced in this file).
membranes = zeros(n_membranes, 4);
mem_is_mono = false(1, n_membranes);
for k = 1:n_membranes
    switch membrane_type{k}
        case 'CEM';      code = 1;
        case 'AEM';      code = 0;
        case 'Mono-CEM'; code = 1; mem_is_mono(k) = true;
        case 'Mono-AEM'; code = 0; mem_is_mono(k) = true;
        case 'BPM';      code = 2;
    end
    membranes(k,:) = [k, k+1, code, 1];
end
mem_mult = membranes(:,4)';
n_mem_copies = sum(mem_mult);
n_bpm_copies = sum(mem_mult(membranes(:,3)==2));

% Which tank feeds each membrane's boundary layer -- the SOURCE side, since
% that's the side being depleted and where polarisation happens.
mem_source = zeros(1, n_membranes);
for k = 1:n_membranes
    if membranes(k,3) == 1;      mem_source(k) = membranes(k,2);   % cations: from the right
    elseif membranes(k,3) == 0;  mem_source(k) = membranes(k,1);   % anions: from the left
    else;                        mem_source(k) = membranes(k,1);   % BPM: no salt source
    end
end

% Monovalent selectivity (eq. 22): a mono-selective membrane's rejection is
% electrostatic, so it strengthens with counter-ion charge. Modelled as a
% permeability factor S_i = MONO_SELECTIVITY^-(|z_i|-1) dividing each ion's
% transport-number weight: monovalent unaffected (S=1), divalent suppressed
% by MONO_SELECTIVITY, trivalent by its square. 20 is mid-range for
% commercial mono-selective grades (published ratios span ~5-50) --
% PLACEHOLDER, measure it for the membrane you actually buy.
MONO_SELECTIVITY = 20;

% Permselectivity / co-ion leakage (eq. 23): no real membrane is perfectly
% selective. A fraction (1-alpha) of the current leaks the WRONG way as
% co-ions -- anions through a CEM, cations through an AEM -- doing no
% separation work. This is the main reason measured current efficiency
% falls short of the ideal-membrane prediction. Commercial grades run
% about 0.90-0.99; AEMs typically a little below CEMs. PLACEHOLDERS --
% datasheet figures are usually measured in dilute KCl and overstate what
% you get in a real leachate.
PERMSEL_CEM = 0.95; PERMSEL_AEM = 0.92; PERMSEL_MONO_CEM = 0.95; PERMSEL_MONO_AEM = 0.92;
permselectivity = zeros(1, n_membranes);
for k = 1:n_membranes
    switch membrane_type{k}
        case 'CEM';      permselectivity(k) = PERMSEL_CEM;
        case 'AEM';      permselectivity(k) = PERMSEL_AEM;
        case 'Mono-CEM'; permselectivity(k) = PERMSEL_MONO_CEM;
        case 'Mono-AEM'; permselectivity(k) = PERMSEL_MONO_AEM;
        case 'BPM';      permselectivity(k) = 1;    % carries no salt either way
    end
end

% Bipolar membrane junction potential and drag (only relevant if you change
% a membrane_type above to 'BPM' -- with none in the default chain these
% simply never get used).
V_BPM      = 0.9;   % V, extra junction potential per BPM copy (thermodynamic
                     % minimum to dissociate water is 0.83 V; real ones run
                     % a bit above it -- PLACEHOLDER)
N_DRAG_BPM = 0;      % mol H2O/mol charge -- H+ and OH- leave a BPM junction
                      % in opposite directions so their drag largely cancels

% Per-membrane clean areal resistance (ohm*cm^2) -- same value on every
% membrane here; give each row its own if you have per-membrane data.
AREAL_R_MEMBRANE = 5;
areal_R = AREAL_R_MEMBRANE * ones(1, n_membranes);


%% SECTION 5 — How we're running the stack (the "knobs" you can turn)
area_cm2          = 10;    % membrane area, cm^2
applied_current_A = 0.30;  % constant current, Amps (30 mA/cm^2). Kept below
                            % this stack's limiting current -- see Section 8
                            % and the main loop for what happens if you push
                            % past it.
duration_hours    = 4;
timestep_minutes  = 1;

timestep_seconds = timestep_minutes * 60;
total_seconds    = duration_hours * 3600;
n_steps          = round(total_seconds / timestep_seconds);   % ticks after t=0


%% SECTION 6 — Voltage/resistance placeholders, over-limiting, polarisation
% Hardware numbers not yet measured for your specific stack -- reasonable
% engineering estimates, flagged so you know exactly what to replace.
GAP_CM               = 0.1;    % cm, compartment spacer thickness (1 mm)
V_OVERPOTENTIAL      = 1.5;    % V, fixed electrode overpotential
KAPPA_WATER_FLOOR    = 5.5e-8; % S/cm, pure water's own conductivity -- a
                                % tank can never have less conductivity
                                % than plain water
V_STACK_CEILING      = 200;    % V, safety cap (should basically never be
                                % hit if the water-splitting term below is
                                % working correctly)
DELTA_CM             = 50e-4;  % cm (50 micron), Nernst diffusion boundary
                                % layer -- the fallback used when a tank has
                                % no forced flow (Section 8 below)
V_WATERSPLIT_OVERPOTENTIAL = 0.8;   % V, added per over-limiting membrane

% Ionic-strength correction to conductivity (eq. 5b). Used raw, the lambda
% values above overestimate a concentrated solution's real conductivity --
% 0.5 M NaCl comes out ~40% too high, which understates stack resistance
% and so understates voltage. This Kohlrausch-shaped correction (fitted to
% NaCl at 25 degC) brings it back in line. It scales every lambda alike, so
% it cancels out of the eq. 1 transport-number ratio -- it only moves
% voltage, never which ion wins a membrane's current.
KAPPA_IONIC_CORRECTION = true;
KAPPA_B = 0.67;

% Concentration polarisation (eq. 26): as current nears a membrane's
% limiting current, the counter-ion concentration AT its surface collapses
% toward zero, and that depleted film is highly resistive -- its resistance
% diverges smoothly rather than jumping in a fixed step. This is what makes
% the model's voltage plateau look like a real polarisation curve.
ENABLE_CONC_POLARISATION = true;
CS_FLOOR = 1e-3;   % floor on c_surface/c_bulk, bounding the divergence


%% SECTION 7 — Water transport through the membranes (eq. 19-20)
% Tank volumes are NOT constant in a real stack. Water crosses membranes by
% two mechanisms -- the common lab observation of one tank steadily gaining
% while its neighbor drains (looking like a leak) is usually these two, not
% a failed seal:
%   (19) ELECTRO-OSMOTIC DRAG -- every ion crossing drags its hydration
%        shell with it, so water moves with the counter-ion: toward the
%        cathode through a CEM, toward the anode through an AEM.
%   (20) OSMOSIS -- water moves from the dilute side toward the
%        concentrated side, i.e. it works to UNDO the separation.
% These move water BETWEEN tanks; total stack volume is conserved (we check
% this in Section 14). Evaporation and real leaks are not modelled.
ENABLE_WATER_TRANSPORT = true;
N_DRAG_CEM   = 6.0;   % mol H2O per mol charge through a CEM (lit. ~4-8 for
                       % a hydrated Na+ form) -- MEASURE IT, nothing here
                       % calibrates this term
N_DRAG_AEM   = 4.0;   % mol H2O per mol charge through an AEM
LP_WATER     = 2e-6;  % cm/s/bar, membrane osmotic permeability (spans
                       % 5e-7 to 1e-5 across commercial membranes)
V_WATER_MOLAR_CM3 = 18.0;   % cm^3/mol
MIN_VOL_FRACTION  = 0.05;   % floor on tank volume (fraction of initial),
                             % keeps things finite if a tank is pumped dry


%% SECTION 8 — Per-tank flow and the resulting diffusion boundary layer (eq. 18)
% Flow is set PER TANK (real stacks run separate hydraulic loops at
% different rates). What flow actually changes is the Nernst boundary
% layer at the membrane surface: crossflow thins it, raising the limiting
% current and delaying water-splitting. Default is ZERO flow everywhere,
% so every tank falls back to the fixed DELTA_CM above -- set an entry
% below to activate the correlation for that tank.
flow_Lpm = zeros(n_tanks, 1);     % L/min through each tank
% e.g. flow_Lpm = [0.2; 1.0; 1.0; 1.5; 1.0; 0.2];   % cathode .. anode

CHANNEL_LENGTH_CM  = sqrt(area_cm2);   % cm, assume a square cell
SPACER_POROSITY    = 0.75;             % open fraction of the channel
KIN_VISC_CM2_S     = 0.00893;          % cm^2/s, water at 25 degC
D_REF_CM2_S        = sp_D(iNa);        % reference diffusivity for Schmidt number

chan_width_cm = area_cm2 / CHANNEL_LENGTH_CM;
Axs_cm2       = chan_width_cm * GAP_CM * SPACER_POROSITY;   % open cross-section
dh_cm         = 2 * GAP_CM * SPACER_POROSITY;               % hydraulic diameter
Sc            = KIN_VISC_CM2_S / D_REF_CM2_S;                % Schmidt number

velocity_cm_s = (flow_Lpm * 1000/60) / Axs_cm2;
delta_cm      = zeros(n_tanks,1);
for i = 1:n_tanks
    u = velocity_cm_s(i);
    if u <= 0
        delta_cm(i) = DELTA_CM;      % no forced flow: fixed placeholder
        continue
    end
    % Leveque correlation, laminar developing flow in a thin channel:
    % Sh = 1.85*(Re*Sc*d_h/L)^(1/3),  delta = d_h/Sh.
    Re = u * dh_cm / KIN_VISC_CM2_S;
    Sh = 1.85 * (Re * Sc * dh_cm / CHANNEL_LENGTH_CM)^(1/3);
    delta_cm(i) = min(dh_cm/max(Sh, eps), 0.5*GAP_CM);   % can't exceed half the gap
end

% Crossflow shears suspended particles off a membrane, so their deposition
% velocity falls with tank velocity (Section 9).
DEPOSITION_VEL_CM_S  = 1e-5;   % cm/s at zero flow -- PLACEHOLDER, lumps
                                % electrophoresis + settling + polarisation
SHEAR_DEP_REF_CM_S   = 5.0;    % cm/s at which crossflow halves deposition
dep_vel_cm_s = DEPOSITION_VEL_CM_S ./ (1 + velocity_cm_s/SHEAR_DEP_REF_CM_S);


%% SECTION 9 — Suspended solids and membrane fouling (eq. 9-11, 27)
% TSS is NOT an ionic species: particles carry no Faradaic current and
% never cross a membrane. What they DO is deposit on the membrane they're
% driven against, forming a cake that (a) adds resistance and (b) thickens
% the diffusion boundary layer, lowering the limiting current -- fouling
% shows up as rising voltage AND earlier water-splitting together.
tss0_mgL = zeros(n_tanks, 1);   % mg/L suspended solids at t=0. Default zero
                                 % (clean-water run); try tss0_mgL(i_d1)=200
                                 % to see fouling develop.
PARTICLE_DENSITY_G_CM3 = 2.5;   % g/cm^3, typical mineral floc/silt
TSS_CHARGE_SIGN        = -1;    % sign of the particle zeta potential; -1
                                 % (typical for mineral colloids) drives
                                 % particles toward the anode, fouling the
                                 % membrane on a tank's ANODE side

% Deposit resistance/thickness derived from loading, density and porosity
% (Bruggeman relation) rather than fitted constants -- a deposit of loading
% w (mg/cm^2) and porosity eps is a layer of thickness w/(rho*(1-eps)) whose
% pores conduct at kappa*eps^1.5.
CAKE_POROSITY  = 0.50;   % loose particulate cake; published range 0.3-0.7
SCALE_POROSITY = 0.10;   % dense crystalline scale; published range 0.05-0.2


%% SECTION 10 — Mineral scaling table (eq. 12-16)
% This stack drives some tanks to pH 1 and others to pH 13+ while
% concentrating metals -- exactly the conditions that precipitate gypsum
% and metal hydroxides. Ksp values are 25 degC thermodynamic solubility
% products; treat them as order-of-magnitude here since real scaling also
% depends on nucleation, seed availability and crystal habit (none of
% which are captured beyond the simple induction-time barrier below).
ENABLE_SCALING     = true;
ENABLE_DISSOLUTION = true;
USE_DAVIES         = true;   % Davies activity coefficients (valid to ~0.5 M;
                              % this stack's ionic strength can exceed that,
                              % so treat results near the extremes as an
                              % extrapolation)

%           name          cation  anion  nu_cat nu_an   Ksp        MW (g/mol)
scale_name  = {'CaSO4.2H2O','Ca(OH)2','Co(OH)2','Ni(OH)2','Fe(OH)2', ...
               'Al(OH)3','Cr(OH)3','Mn(OH)2','Mg(OH)2','Cu(OH)2','Zn(OH)2'};
scale_iCat  = [iCa, iCa, iCo, iNi, iFe, iAl, iCr, iMn, iMg, iCu, iZn];
scale_iAn   = [iSO4,iOH, iOH, iOH, iOH, iOH, iOH, iOH, iOH, iOH, iOH];
scale_nuCat = [ 1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1  ];
scale_nuAn  = [ 1,    2,    2,    2,    2,    3,    3,    2,    2,    2,    2  ];
scale_Ksp   = [3.14e-5, 5.5e-6, 5.9e-15, 5.5e-16, 4.9e-17, 3.0e-34, ...
               6.3e-31, 1.9e-13, 5.6e-12, 2.2e-20, 3.0e-17];
scale_MW    = [172.17, 74.09, 92.95, 92.71, 89.86, 78.00, 103.02, 88.95, 58.32, 97.56, 99.42];
n_scale     = numel(scale_name);

K_PRECIP_PER_S        = 1e-3;   % 1/s, relaxation rate toward equilibrium --
                                 % ~17 min time constant, fast vs. a 4 h run
                                 % but not instantaneous
K_DISSOL_PER_S        = 5e-4;   % 1/s, dissolution relaxation (slower --
                                 % dissolving a formed crystal is surface-
                                 % limited, precipitation can proceed in the bulk)
SCALE_SURFACE_FRACTION = 0.30;  % fraction of new solid that nucleates
                                 % directly on the membrane rather than the bulk

% Nucleation barrier (eq. 15): crystals don't appear the instant SI crosses
% zero. No NEW solid forms until SI exceeds SI_NUCLEATION AND the tank has
% spent t_ind = TAU_NUC_REF_S/SI^2 above that threshold. Once ANY seed of
% that solid exists (suspended or on a membrane), the barrier is gone.
SI_NUCLEATION = 0.5;     % log10 units of supersaturation needed to nucleate
TAU_NUC_REF_S = 1800;    % s, induction time at SI = 1 (falls as 1/SI^2)


%% SECTION 11 — Set up empty storage for the whole run's history
n_pts = n_steps + 1;   % +1 for the t=0 snapshot

time_min                = zeros(n_pts, 1);
conc_history            = zeros(n_pts, n_tanks, n_species);
pH_history              = zeros(n_pts, n_tanks);
voltage_history         = zeros(n_pts, 1);
Vfoul_history           = zeros(n_pts, 1);
Vscale_history          = zeros(n_pts, 1);
overlimit_count_history = zeros(n_pts, 1);
co_efficiency_history   = zeros(n_pts, 1);
mass_balance_history    = zeros(n_pts, numel(conserved_idx));
vol_history             = zeros(n_pts, n_tanks);
tss_history             = zeros(n_pts, n_tanks);
cake_history            = zeros(n_pts, n_membranes);    % mg/cm^2, total per membrane
scale_history           = zeros(n_pts, n_membranes);    % mg/cm^2, total per membrane
scale_mass_history      = zeros(n_pts, n_scale);        % mg, summed over tanks
SI_history              = zeros(n_pts, n_tanks, n_scale);
Ilim_history            = zeros(n_pts, n_membranes);
Iionic_history          = zeros(n_pts, n_membranes);

% Solids/scaling STATE, carried from one step to the next (this is the
% "st" struct's fields in the reference model, kept here as plain arrays):
tss_inert_mgL = tss0_mgL(:);           % fed solids, never dissolve
solid_bulk_mg = zeros(n_tanks, n_scale);   % precipitate still suspended
cake_inert    = zeros(1, n_membranes);     % mg/cm^2, deposited inert particulate
cake_solid    = zeros(n_membranes, n_scale);   % mg/cm^2, deposited precipitate
scale_mem     = zeros(n_membranes, n_scale);   % mg/cm^2, scale grown in place
ind_clock_s   = zeros(n_tanks, n_scale);       % s, time above the nucleation SI
precip_mol    = zeros(1, n_species);           % mol, ions currently locked in solid (stack-wide)
formed_mg     = zeros(n_tanks, n_scale);       % mg, gross ever precipitated
dissolved_mg  = zeros(n_tanks, n_scale);       % mg, gross ever redissolved
R_bl          = zeros(1, n_membranes);         % ohm, concentration-polarisation resistance


%% SECTION 12 — Record the very first point (t = 0), before we run anything
time_min(1) = 0;

% Every tank starts at its own dosed pH from Section 3b, so re-speciate
% once onto the water equilibrium for a clean, fully consistent snapshot
% (closed-form solve, eq. 17 -- see the main loop for the full explanation).
for i = 1:n_tanks
    net  = conc(i,iH) - conc(i,iOH);
    disc = sqrt(net*net + 4*Kw);
    if net >= 0
        cH = 0.5*(net+disc); cOH = Kw/max(cH,realmin);
    else
        cOH = 0.5*(-net+disc); cH = Kw/max(cOH,realmin);
    end
    conc(i,iH) = cH; conc(i,iOH) = cOH;
end
conc_history(1,:,:) = conc;
vol_history(1,:)    = tank_volume_L;
tss_history(1,:)    = tss_inert_mgL' + (sum(solid_bulk_mg,2) ./ tank_volume_L)';

for i = 1:n_tanks
    pH_history(1,i) = -log10(max(conc(i,iH), 1e-30));
end

% Stack voltage at t=0 (no membrane is over-limiting yet, no cake/scale
% has formed yet, so this reduces to the plain solution+membrane resistance).
total_R = 0;
for i = 1:n_tanks
    kap = 0.001 * sum(sp_lambda .* abs(sp_z) .* max(conc(i,:),0));
    kap = max(kap, KAPPA_WATER_FLOOR);
    total_R = total_R + GAP_CM/(kap*area_cm2);
end
total_R = total_R + sum(mem_mult .* areal_R / area_cm2);
voltage_history(1) = min(applied_current_A*total_R + V_OVERPOTENTIAL, V_STACK_CEILING);
Vfoul_history(1)  = 0;
Vscale_history(1) = 0;

overlimit_count_history(1) = 0;
co_efficiency_history(1)   = 0;
current_totals = sum(conc .* tank_volume_L, 1);
mass_balance_history(1,:) = current_totals(conserved_idx) ./ starting_total_moles(conserved_idx) * 100;
SI_history(1,:,:) = -Inf(n_tanks, n_scale);   % nothing dissolved yet in most tanks; recomputed properly inside the loop from step 2 on
for i = 1:n_tanks
    for k = 1:n_scale
        cC = max(conc(i,scale_iCat(k)),0); cA = max(conc(i,scale_iAn(k)),0);
        if cC>0 && cA>0
            SI_history(1,i,k) = log10((cC^scale_nuCat(k))*(cA^scale_nuAn(k)) / scale_Ksp(k));
        end
    end
end


%% SECTION 13 — THE MAIN LOOP: step the simulation forward in time
% One pass through this loop = one timestep. Each pass:
%   (a) moles = conc x volume, for every tank/ion right now
%   (b) go through the 5 membranes: transport-number share (eq.1), monovalent
%       selectivity (eq.22), limiting current with the fouling-thickened
%       boundary layer (eq.8, 11, 18, 27), permselectivity/co-ion leakage
%       (eq.23), concentration-polarisation resistance (eq.26), and
%       water-splitting for whatever current the ions can't carry (eq.8 cont.)
%   (c) the two electrode reactions (eq.4)
%   (d) water transport between tanks (eq.19-20)
%   (e) moles -> concentrations, using the (possibly changed) volumes
%   (f) water self-ionization re-speciation (eq.6, 17)
%   (g) precipitation / nucleation / dissolution (eq.12-16)
%   (h) suspended-solids deposition onto membranes (eq.9-10)
%   (i) stack voltage this step (eq.5, 5b, 11, 26, 27, water-split, BPM)
%   (j) save everything into the history arrays from Section 11

for step = 1:n_steps

    %% --- (a) moles right now ---
    moles = conc .* tank_volume_L;
    overlimit_flag_this_step = false(1, n_membranes);
    Ilim_this_step   = zeros(1, n_membranes);
    Iionic_this_step = zeros(1, n_membranes);

    % Effective boundary layer per membrane: this tank's own flow-derived
    % delta (eq.18), thickened by whatever cake and scale (eq.9-16) have
    % deposited on it, via their derived thickness (eq.27).
    cake_tot  = cake_inert + sum(cake_solid, 2)';     % mg/cm^2, per membrane
    scale_tot = sum(scale_mem, 2)';                    % mg/cm^2, per membrane
    t_cake_cm  = cake_tot  / (PARTICLE_DENSITY_G_CM3 * (1 - CAKE_POROSITY )) / 1000;
    t_scale_cm = scale_tot / (PARTICLE_DENSITY_G_CM3 * (1 - SCALE_POROSITY)) / 1000;
    delta_eff = delta_cm(mem_source)' + t_cake_cm + t_scale_cm;

    %% --- (b) go through the 5 membranes ---
    for k = 1:n_membranes
        left = membranes(k,1); right = membranes(k,2); code = membranes(k,3); mult = membranes(k,4);

        if code == 2
            %── BIPOLAR: carries no salt, splits water at the junction ──
            dMolBPM = mult * applied_current_A/F * timestep_seconds;
            moles(left,  iH)  = moles(left,  iH)  + dMolBPM;
            moles(right, iOH) = moles(right, iOH) + dMolBPM;
            R_bl(k) = 0;
            continue
        end

        if code == 1
            source = right; destination = left;   wantSign = 1;   % cations: right -> left
        else
            source = left;  destination = right;  wantSign = -1;  % anions: left -> right
        end

        carriers = find(sign(sp_z) == wantSign);
        source_conc = max(conc(source, carriers), 0);

        % --- EQUATION 1 (+ 22): transport-number share, with monovalent
        % selectivity dividing a multivalent ion's weight if this membrane
        % is Mono-CEM/Mono-AEM (sel=1 for every ion on a standard membrane).
        if mem_is_mono(k)
            sel = MONO_SELECTIVITY .^ -(abs(sp_z(carriers)) - 1);
        else
            sel = ones(1, numel(carriers));
        end
        weight = sel .* abs(sp_z(carriers)) .* sp_lambda(carriers) .* source_conc;
        denom  = sum(weight);

        % --- EQUATION 8 (+ 18, 27): limiting current from diffusion to the
        % membrane surface through its (fouling-thickened) boundary layer.
        I_lim = area_cm2 * sum(sel .* abs(sp_z(carriers)) .* F .* sp_D(carriers) ...
                    .* source_conc / 1000 / delta_eff(k));

        % --- EQUATION 23: only the COUNTER-ion share of the current (alpha
        % fraction) has to cross as counter-ions, so that's what I_lim caps.
        alpha = permselectivity(k);
        I_counter = alpha * applied_current_A;
        is_over_limit = I_counter > I_lim;
        overlimit_flag_this_step(k) = is_over_limit;
        I_ionic = min(I_counter, I_lim);
        Ilim_this_step(k)   = I_lim;
        Iionic_this_step(k) = I_ionic;

        % --- EQUATION 26: concentration-polarisation resistance from the
        % depleted film's log-mean conductivity between bulk and surface.
        if ENABLE_CONC_POLARISATION && I_lim > 1e-30
            ratio = min(I_counter/I_lim, 1 - CS_FLOOR);
            x = max(1 - ratio, CS_FLOOR);             % c_surface / c_bulk
            if ratio < 1e-9
                fac = 1;
            else
                fac = (1-x) / log(1/x);
            end
            kap_src = 0.001 * sum(sp_lambda .* abs(sp_z) .* max(conc(source,:),0));
            kap_src = max(kap_src, KAPPA_WATER_FLOOR);
            R_bl(k) = delta_eff(k)/(area_cm2*kap_src) * (1/max(fac,CS_FLOOR) - 1);
        else
            R_bl(k) = 0;
        end

        % --- EQUATION 2: move each carrier ion in proportion to its share.
        if denom > 1e-30 && I_ionic > 1e-30
            share = weight / denom;
            flux_mol_per_s = mult * share * I_ionic ./ (abs(sp_z(carriers)) * F);
            dMol = flux_mol_per_s * timestep_seconds;
            dMol = min(dMol, max(moles(source, carriers), 0));   % availability cap
            moles(source, carriers)      = moles(source, carriers)      - dMol;
            moles(destination, carriers) = moles(destination, carriers) + dMol;
        end

        % --- EQUATION 23 (continued): co-ion leakage. The remaining
        % (1-alpha) fraction of the current crosses as CO-ions going the
        % OTHER way, drawn from the opposite tank by their own transport
        % numbers. This moves salt backwards, undoing separation.
        if alpha < 1
            I_co = (1-alpha) * applied_current_A;
            src_co = destination; dst_co = source;
            carriers_co = find(sign(sp_z) == -wantSign);
            src_conc_co = max(conc(src_co, carriers_co), 0);
            weight_co = abs(sp_z(carriers_co)) .* sp_lambda(carriers_co) .* src_conc_co;
            denom_co  = sum(weight_co);
            if denom_co > 1e-30 && I_co > 1e-30
                share_co = weight_co / denom_co;
                flux_co  = mult * share_co * I_co ./ (abs(sp_z(carriers_co)) * F);
                dMol_co  = min(flux_co * timestep_seconds, max(moles(src_co, carriers_co), 0));
                moles(src_co, carriers_co) = moles(src_co, carriers_co) - dMol_co;
                moles(dst_co, carriers_co) = moles(dst_co, carriers_co) + dMol_co;
            end
        end

        % --- EQUATION 8 (continued): water-splitting for the shortfall.
        % CEM: H+ crosses to the destination, OH- stays in the source.
        % AEM: the mirror image. Branch on CODE (cation/anion), not the
        % type string, so Mono-CEM splits water the same way a plain CEM
        % does.
        if is_over_limit
            dMolWS = mult * (I_counter - I_lim) / F * timestep_seconds;
            if code == 1
                moles(destination, iH)  = moles(destination, iH)  + dMolWS;
                moles(source,      iOH) = moles(source,      iOH) + dMolWS;
            else
                moles(destination, iOH) = moles(destination, iOH) + dMolWS;
                moles(source,      iH)  = moles(source,      iH)  + dMolWS;
            end
        end
    end   % end of the membrane loop


    %% --- (c) the two electrode reactions (eq. 4) ---
    % Cathode (Tank 1): 2H2O + 2e- -> H2 + 2OH-   =>  d(OH-)/dt = I/F
    moles(i_cathode, iOH) = moles(i_cathode, iOH) + (applied_current_A/F)*timestep_seconds;
    % Anode (Tank 6):   2H2O -> O2 + 4H+ + 4e-    =>  d(H+)/dt  = I/F
    moles(i_anode,   iH)  = moles(i_anode,   iH)  + (applied_current_A/F)*timestep_seconds;


    %% --- (d) water transport between tanks (eq. 19-20) ---
    % Applied to volumes BEFORE converting moles back to concentrations, so
    % the concentrating/diluting effect shows up this same step.
    if ENABLE_WATER_TRANSPORT
        dVol = zeros(n_tanks,1);
        for k = 1:n_membranes
            left = membranes(k,1); right = membranes(k,2); mult = membranes(k,4);
            switch membranes(k,3)
                case 1;   src = right; dst = left;  nDrag = N_DRAG_CEM;   % cations drag water right->left
                case 0;   src = left;  dst = right; nDrag = N_DRAG_AEM;   % anions drag water left->right
                otherwise; src = left; dst = right; nDrag = N_DRAG_BPM;
            end
            % (19) electro-osmotic drag: proportional to charge passed.
            molCharge = mult * applied_current_A/F * timestep_seconds;
            dV_drag = nDrag * molCharge * V_WATER_MOLAR_CM3 / 1000;   % L
            % (20) osmosis, van 't Hoff: water toward the concentrated side.
            osmL = sum(max(conc(left,:),0));  osmR = sum(max(conc(right,:),0));
            dPi = Rgas/100 * T * (osmR - osmL);   % bar
            dV_osm = LP_WATER * mult * area_cm2 * dPi * timestep_seconds / 1000;   % L
            dVol(src) = dVol(src) - dV_drag;  dVol(dst) = dVol(dst) + dV_drag;
            dVol(left) = dVol(left) - dV_osm; dVol(right) = dVol(right) + dV_osm;
        end
        % Floor: scale back the WHOLE step's transfer by one factor so no
        % tank is driven below its minimum volume, while keeping
        % sum(dVol)=0 exactly (so total volume stays conserved).
        volMin = MIN_VOL_FRACTION * tank_volume_L;
        f = 1;
        for i = 1:n_tanks
            if dVol(i) < 0 && (tank_volume_L(i)+dVol(i)) < volMin(i)
                f = min(f, max(volMin(i)-tank_volume_L(i), 0) / dVol(i));
            end
        end
        volBefore = tank_volume_L;
        tank_volume_L = max(tank_volume_L + f*dVol, volMin);
        % Water moves; suspended particles it leaves behind do not, so their
        % CONCENTRATION has to be rescaled or mass is silently created/lost.
        tss_inert_mgL = tss_inert_mgL .* volBefore ./ tank_volume_L;
    end


    %% --- (e) moles -> concentrations ---
    conc = max(moles ./ tank_volume_L, 0);


    %% --- (f) water self-ionization re-speciation (eq. 6, 17) ---
    % H+ and OH- were just transported/generated independently, so nothing
    % above stops a tank ending up with high H+ AND high OH- at once --
    % chemically impossible (they neutralize on contact). Re-speciate onto
    % C_H*C_OH = Kw, holding the NET (C_H - C_OH) fixed, since net is what
    % transport/electrodes actually determined and holding it fixed
    % conserves charge exactly. Solved in closed form:
    %   cH = (net + sqrt(net^2 + 4 Kw)) / 2      (evaluated to avoid
    % cancellation for whichever sign net has).
    for i = 1:n_tanks
        net  = conc(i,iH) - conc(i,iOH);
        disc = sqrt(net*net + 4*Kw);
        if net >= 0
            cH = 0.5*(net+disc); cOH = Kw/max(cH,realmin);
        else
            cOH = 0.5*(-net+disc); cH = Kw/max(cOH,realmin);
        end
        conc(i,iH) = cH; conc(i,iOH) = cOH;
    end


    %% --- (g) precipitation / nucleation / dissolution (eq. 12-16) ---
    if ENABLE_SCALING
        relaxP = 1 - exp(-K_PRECIP_PER_S * timestep_seconds);
        relaxD = 1 - exp(-K_DISSOL_PER_S * timestep_seconds);
        for i = 1:n_tanks
            mems_here = find(membranes(:,1)==i | membranes(:,2)==i)';
            for k = 1:n_scale
                ic = scale_iCat(k); ia = scale_iAn(k);
                a = scale_nuCat(k); b = scale_nuAn(k);
                cC = max(conc(i,ic),0); cA = max(conc(i,ia),0);

                % Davies activity coefficients (this tank's ionic strength).
                if USE_DAVIES
                    Ion = 0.5*sum(sp_z.^2 .* max(conc(i,:),0));
                    sI = sqrt(Ion); Ad = 0.509;
                    gC = 10^(-Ad*sp_z(ic)^2*(sI/(1+sI)-0.3*Ion));
                    gA = 10^(-Ad*sp_z(ia)^2*(sI/(1+sI)-0.3*Ion));
                else
                    gC = 1; gA = 1;
                end
                Ksp = scale_Ksp(k);

                rowArea = mem_mult(mems_here) * area_cm2;   % cm^2 per bounding-membrane row
                seed_mg = solid_bulk_mg(i,k) + ...
                    sum(cake_solid(mems_here,k)'.*rowArea) + sum(scale_mem(mems_here,k)'.*rowArea);

                IAP = (gC*cC)^a * (gA*cA)^b;
                if cC<=0 || cA<=0; IAP = 0; end

                if IAP > Ksp
                    %── SUPERSATURATED: grow, subject to the nucleation barrier (eq.15) ──
                    SI = log10(IAP/Ksp);
                    if seed_mg <= 1e-12
                        if SI < SI_NUCLEATION
                            ind_clock_s(i,k) = max(ind_clock_s(i,k) - timestep_seconds, 0);
                            continue
                        end
                        ind_clock_s(i,k) = ind_clock_s(i,k) + timestep_seconds;
                        t_ind = TAU_NUC_REF_S / SI^2;
                        if ind_clock_s(i,k) < t_ind; continue; end
                    end
                    % Extent x that restores equilibrium, by bisection: IAP
                    % decreases monotonically as x grows (ions are being
                    % removed), so no derivative is needed.
                    hi = min(cC/a, cA/b);
                    f_precip = @(xx) (gC*max(cC-a*xx,0))^a * (gA*max(cA-b*xx,0))^b - Ksp;
                    lo_b = 0; hi_b = hi;
                    if hi_b <= lo_b
                        x = 0;
                    elseif f_precip(hi_b) > 0
                        x = hi_b;
                    else
                        for bi = 1:60
                            mid_b = 0.5*(lo_b+hi_b);
                            if f_precip(mid_b) > 0; lo_b = mid_b; else; hi_b = mid_b; end
                        end
                        x = 0.5*(lo_b+hi_b);
                    end
                    x = x * relaxP;
                    if x <= 0; continue; end
                    conc(i,ic) = max(conc(i,ic)-a*x,0);
                    conc(i,ia) = max(conc(i,ia)-b*x,0);
                    precip_mol(ic) = precip_mol(ic) + a*x*tank_volume_L(i);
                    precip_mol(ia) = precip_mol(ia) + b*x*tank_volume_L(i);
                    newMass = x * tank_volume_L(i) * scale_MW(k) * 1000;   % mol -> mg
                    formed_mg(i,k) = formed_mg(i,k) + newMass;
                    if isempty(mems_here)
                        solid_bulk_mg(i,k) = solid_bulk_mg(i,k) + newMass;
                    else
                        surfMass = newMass * SCALE_SURFACE_FRACTION;
                        shareA = rowArea / sum(rowArea);
                        scale_mem(mems_here,k) = scale_mem(mems_here,k) + (surfMass*shareA ./ rowArea)';
                        solid_bulk_mg(i,k) = solid_bulk_mg(i,k) + (newMass - surfMass);
                    end
                elseif ENABLE_DISSOLUTION && seed_mg > 1e-12
                    %── UNDERSATURATED with solid present: redissolve (eq.16) ──
                    ind_clock_s(i,k) = max(ind_clock_s(i,k) - timestep_seconds, 0);
                    maxX = seed_mg / (scale_MW(k)*1000) / tank_volume_L(i);
                    % Same equilibrium condition, opposite direction: IAP
                    % increases monotonically as x grows (ions are being
                    % added back into solution).
                    g_dissol = @(xx) Ksp - (gC*(cC+a*xx))^a * (gA*(cA+b*xx))^b;
                    lo_b = 0; hi_b = maxX;
                    if hi_b <= lo_b
                        x = 0;
                    elseif g_dissol(hi_b) > 0
                        x = hi_b;
                    else
                        for bi = 1:60
                            mid_b = 0.5*(lo_b+hi_b);
                            if g_dissol(mid_b) > 0; lo_b = mid_b; else; hi_b = mid_b; end
                        end
                        x = 0.5*(lo_b+hi_b);
                    end
                    x = x * relaxD;
                    if x <= 0; continue; end
                    conc(i,ic) = conc(i,ic) + a*x;
                    conc(i,ia) = conc(i,ia) + b*x;
                    precip_mol(ic) = precip_mol(ic) - a*x*tank_volume_L(i);
                    precip_mol(ia) = precip_mol(ia) - b*x*tank_volume_L(i);
                    need = x * tank_volume_L(i) * scale_MW(k) * 1000;   % mg
                    dissolved_mg(i,k) = dissolved_mg(i,k) + need;
                    take = min(need, solid_bulk_mg(i,k));
                    solid_bulk_mg(i,k) = solid_bulk_mg(i,k) - take; need = need - take;
                    % Consume loosest-first: suspended crystal, then cake, then in-situ scale.
                    for mi = 1:numel(mems_here)
                        if need<=0; break; end
                        m = mems_here(mi); aM = rowArea(mi);
                        take = min(need, cake_solid(m,k)*aM);
                        cake_solid(m,k) = cake_solid(m,k) - take/aM; need = need - take;
                    end
                    for mi = 1:numel(mems_here)
                        if need<=0; break; end
                        m = mems_here(mi); aM = rowArea(mi);
                        take = min(need, scale_mem(m,k)*aM);
                        scale_mem(m,k) = scale_mem(m,k) - take/aM; need = need - take;
                    end
                else
                    ind_clock_s(i,k) = max(ind_clock_s(i,k) - timestep_seconds, 0);
                end
            end
        end
        % Re-impose water equilibrium: hydroxide precipitation/dissolution
        % just consumed/released OH-, so re-speciate again (operator
        % splitting -- one pass per step, not iterated to joint convergence).
        for i = 1:n_tanks
            net  = conc(i,iH) - conc(i,iOH);
            disc = sqrt(net*net + 4*Kw);
            if net >= 0
                cH = 0.5*(net+disc); cOH = Kw/max(cH,realmin);
            else
                cOH = 0.5*(-net+disc); cH = Kw/max(cOH,realmin);
            end
            conc(i,iH) = cH; conc(i,iOH) = cOH;
        end
    end


    %% --- (h) suspended-solids deposition onto membranes (eq. 9-11) ---
    tss_tot = tss_inert_mgL + (sum(solid_bulk_mg,2) ./ tank_volume_L);
    for i = 1:n_tanks
        if tss_tot(i) <= 0; continue; end
        if TSS_CHARGE_SIGN < 0
            memList = find(membranes(:,1)==i)';   % this tank's anode-side membrane(s)
        else
            memList = find(membranes(:,2)==i)';   % cathode-side
        end
        if isempty(memList); continue; end        % electrode-facing tank: no membrane
        depArea = sum(membranes(memList,4)) * area_cm2;
        dMass = dep_vel_cm_s(i) * (tss_tot(i)/1000) * timestep_seconds * depArea;   % mg
        dMass = min(dMass, max(tss_tot(i)*tank_volume_L(i), 0));
        if dMass <= 0; continue; end
        frac = dMass / (tss_tot(i)*tank_volume_L(i));
        wgt = membranes(memList,4)' * area_cm2 / depArea;

        dInert = tss_inert_mgL(i)*tank_volume_L(i) * frac;
        tss_inert_mgL(i) = max(tss_inert_mgL(i) - dInert/tank_volume_L(i), 0);
        cake_inert(memList) = cake_inert(memList) + (dInert*wgt) ./ (membranes(memList,4)'*area_cm2);

        dEach = solid_bulk_mg(i,:) * frac;
        solid_bulk_mg(i,:) = max(solid_bulk_mg(i,:) - dEach, 0);
        for mm = 1:numel(memList)
            cake_solid(memList(mm),:) = cake_solid(memList(mm),:) + ...
                dEach*wgt(mm) / (membranes(memList(mm),4)*area_cm2);
        end
    end


    %% --- (i) stack voltage this step (eq. 5, 5b, 11, 26, 27, water-split, BPM) ---
    tss_tot = tss_inert_mgL + (sum(solid_bulk_mg,2) ./ tank_volume_L);
    Rcomp = zeros(1, n_tanks); kappa_now = zeros(1, n_tanks);
    for i = 1:n_tanks
        lam = sp_lambda;
        if KAPPA_IONIC_CORRECTION
            Ion = 0.5*sum(sp_z.^2 .* max(conc(i,:),0));
            sI = sqrt(Ion);
            lam = lam * max(1 - KAPPA_B*sI/(1+sI), 0.05);
        end
        kap = 0.001 * sum(lam .* abs(sp_z) .* max(conc(i,:),0));
        if tss_tot(i) > 0
            phi = min(tss_tot(i)/(PARTICLE_DENSITY_G_CM3*1e6), 0.6);
            kap = kap*(1-phi)/(1+phi/2);
        end
        kap = max(kap, KAPPA_WATER_FLOOR);
        kappa_now(i) = kap;
        Rcomp(i) = GAP_CM/(kap*area_cm2);
    end
    [cake_tot, scale_tot] = deal(cake_inert + sum(cake_solid,2)', sum(scale_mem,2)');
    Rmem_clean = mem_mult .* areal_R / area_cm2;
    kapAt = max(kappa_now(mem_source), KAPPA_WATER_FLOOR);
    % Areal resistance (ohm*cm^2) of a deposit of loading w (mg/cm^2),
    % porosity eps and particle density rho: thickness t = w/(rho*(1-eps)),
    % conducting only through its pores at kappa*eps^1.5 (Bruggeman).
    t_cake_layer  = cake_tot  / (PARTICLE_DENSITY_G_CM3 * (1-CAKE_POROSITY )) / 1000;
    t_scale_layer = scale_tot / (PARTICLE_DENSITY_G_CM3 * (1-SCALE_POROSITY)) / 1000;
    rCake  = t_cake_layer  ./ max(kapAt * CAKE_POROSITY^1.5,  1e-12);
    rScale = t_scale_layer ./ max(kapAt * SCALE_POROSITY^1.5, 1e-12);
    Rmem_cake  = mem_mult .* rCake  / area_cm2;
    Rmem_scale = mem_mult .* rScale / area_cm2;
    Rmem_polar = mem_mult .* R_bl;
    Rmem_total = Rmem_clean + Rmem_cake + Rmem_scale + Rmem_polar;
    Rtotal = sum(Rcomp) + sum(Rmem_total);
    overlimit_count_history(step+1) = sum(overlimit_flag_this_step);
    Vfoul_this  = applied_current_A * sum(Rmem_cake);
    Vscale_this = applied_current_A * sum(Rmem_scale);
    V_watersplit = overlimit_count_history(step+1)*V_WATERSPLIT_OVERPOTENTIAL + n_bpm_copies*V_BPM;
    V_this = min(applied_current_A*Rtotal + V_OVERPOTENTIAL + V_watersplit, V_STACK_CEILING);


    %% --- (j) save this step into the history arrays ---
    time_min(step+1) = time_min(step) + timestep_minutes;
    conc_history(step+1,:,:) = conc;
    vol_history(step+1,:)    = tank_volume_L;
    tss_history(step+1,:)    = tss_tot';
    cake_history(step+1,:)   = cake_tot;
    scale_history(step+1,:)  = scale_tot;
    scale_mass_history(step+1,:) = sum(solid_bulk_mg,1) + ...
        (mem_mult*cake_solid + mem_mult*scale_mem) * area_cm2;
    voltage_history(step+1) = V_this;
    Vfoul_history(step+1)   = Vfoul_this;
    Vscale_history(step+1)  = Vscale_this;
    Ilim_history(step+1,:)   = Ilim_this_step;
    Iionic_history(step+1,:) = Iionic_this_step;

    for i = 1:n_tanks
        pH_history(step+1,i) = -log10(max(conc(i,iH), 1e-30));
        for k = 1:n_scale
            cC = max(conc(i,scale_iCat(k)),0); cA = max(conc(i,scale_iAn(k)),0);
            if cC>0 && cA>0
                SI_history(step+1,i,k) = log10((cC^scale_nuCat(k))*(cA^scale_nuAn(k)) / scale_Ksp(k));
            else
                SI_history(step+1,i,k) = -Inf;
            end
        end
    end

    % Cumulative Co2+ current efficiency vs time (eq. 7): of all the charge
    % pushed through so far, what fraction actually moved Co2+ into Tank 3?
    co_gain = conc(i_c2,iCo)*tank_volume_L(i_c2) - conc_history(1,i_c2,iCo)*vol_history(1,i_c2);
    theoretical_max = (applied_current_A * time_min(step+1)*60) / (sp_z(iCo)*F);
    if theoretical_max > 1e-12
        co_efficiency_history(step+1) = co_gain/theoretical_max*100;
    end

    % Mass balance closure: dissolved inventory PLUS whatever's locked in
    % solids (precipitation removes ions from solution but not from the
    % system, so both have to be counted or scaling would look like a leak).
    current_totals = sum(conc .* tank_volume_L, 1) + precip_mol;
    mass_balance_history(step+1,:) = current_totals(conserved_idx) ./ starting_total_moles(conserved_idx) * 100;

end   % end of the main time-stepping loop


%% SECTION 14 — Check the results: did we conserve mass? What did we get?
mass_balance_pct = mass_balance_history(end,:);
worst_mass_balance_pct = mass_balance_pct(find(abs(mass_balance_pct-100) == max(abs(mass_balance_pct-100)), 1));

co_recovered_pct = (conc(i_c2,iCo)*tank_volume_L(i_c2) - conc_history(1,i_c2,iCo)*vol_history(1,i_c2)) ...
    / max(conc_history(1,i_d1,iCo)*vol_history(1,i_d1), 1e-12) * 100;

wanted_eq = conc(i_c1,iNa)*1 + conc(i_c1,iSO4)*2;
contam_eq = conc(i_c1,iH)*1 + conc(i_c1,iCl)*1 + sum(abs(sp_z(metal_idx)) .* conc(i_c1,metal_idx));
c1_purity_pct = 100 * wanted_eq / max(wanted_eq+contam_eq, 1e-12);

cathode_Na_gained_mmol = (conc(i_cathode,iNa) - conc_history(1,i_cathode,iNa)) * tank_volume_L(i_cathode) * 1000;

initial_solids_mg = sum(tss0_mgL .* (0.250*ones(n_tanks,1)));
precip_total_mg   = sum(scale_mass_history(end,:));
expected_solids_mg = initial_solids_mg + precip_total_mg;
if expected_solids_mg > 1e-12
    final_solids_mg = sum(tss_history(end,:)' .* vol_history(end,:)') ...
        + sum(cake_history(end,:)  .* mem_mult) * area_cm2 ...
        + sum(scale_history(end,:) .* mem_mult) * area_cm2;
    solids_balance_pct = final_solids_mg / expected_solids_mg * 100;
else
    solids_balance_pct = 100;
end

fprintf('\n══════════════════════════════════════════════\n');
fprintf('  EDM Quad Stack (walkthrough script) — Results\n');
fprintf('══════════════════════════════════════════════\n');
fprintf('  Mass balance closure    : %.2f%%\n', worst_mass_balance_pct);
fprintf('  Solids balance          : %.2f%%\n', solids_balance_pct);
fprintf('  Co2+ recovered -> C2    : %.1f%%\n', co_recovered_pct);
fprintf('  C1 product purity       : %.1f%%\n', c1_purity_pct);
fprintf('  Cathode Na+ gained      : %.2f mmol\n', cathode_Na_gained_mmol);
fprintf('  Current efficiency, Co  : %.1f%%\n', co_efficiency_history(end));
fprintf('  Membranes over-limiting : %d / %d\n', overlimit_count_history(end), n_membranes);
fprintf('  Final stack voltage     : %.2f V  (%.3f V fouling, %.3f V scaling)\n', ...
    voltage_history(end), Vfoul_history(end), Vscale_history(end));
fprintf('  Total scale/precipitate : %.2f mg\n', sum(scale_mass_history(end,:)));
fprintf('──────────────────────────────────────────────\n');
fprintf('  Per-metal recovery (feed D1 -> product C2):\n');
for j = 1:numel(metal_idx)
    s = metal_idx(j);
    m0  = conc_history(1,i_d1,s) * vol_history(1,i_d1);
    if m0 <= 1e-12; continue; end
    mC2 = conc(i_c2,s)*tank_volume_L(i_c2) - conc_history(1,i_c2,s)*vol_history(1,i_c2);
    fprintf('    %-6s recovered %6.1f%%  (%.3f mmol of %.3f mmol fed)\n', ...
        sp_label{s}, mC2/m0*100, mC2*1000, m0*1000);
end
fprintf('══════════════════════════════════════════════\n\n');
if ~isempty(find(scale_mass_history(end,:) > 1e-9, 1))
    fprintf('  Scale actually formed: %s\n\n', strjoin(scale_name(scale_mass_history(end,:) > 1e-9), ', '));
end


%% SECTION 15 — Plot the results
figure('Name','EDM Quad Stack — Walkthrough Results','Position',[60 40 1700 950]);

subplot(3,3,1);
plot(time_min, conc_history(:,i_d1,iCo), 'LineWidth',2); hold on;
plot(time_min, conc_history(:,i_c2,iCo), 'LineWidth',2);
plot(time_min, conc_history(:,i_d1,iNi), '--', 'LineWidth',1.5);
plot(time_min, conc_history(:,i_c2,iNi), '--', 'LineWidth',1.5);
xlabel('Time (min)'); ylabel('Concentration (mol/L)');
legend('D1 Co2+','C2 Co2+','D1 Ni2+','C2 Ni2+','Location','best');
title('Key metal concentrations'); grid on;

subplot(3,3,2);
plot(time_min, mass_balance_history, 'LineWidth',1.5);
legend(sp_label(conserved_idx), 'Location','eastoutside','FontSize',6);
yline(100,'--k','HandleVisibility','off');
xlabel('Time (min)'); ylabel('Mass balance closure (%)'); ylim([95 105]);
title('Mass balance closure'); grid on;

subplot(3,3,3);
plot(time_min, voltage_history, 'LineWidth',2); hold on;
plot(time_min, Vfoul_history, 'LineWidth',1); plot(time_min, Vscale_history, 'LineWidth',1);
legend('Total','Fouling','Scaling','Location','best');
xlabel('Time (min)'); ylabel('Voltage (V)');
title(sprintf('Voltage (final %.2f V)', voltage_history(end))); grid on;

subplot(3,3,4);
plot(time_min, pH_history, 'LineWidth',1.5);
xlabel('Time (min)'); ylabel('pH'); ylim([0 14]);
legend(tank_name, 'Location','eastoutside','FontSize',6);
title('pH of every tank'); grid on;

subplot(3,3,5);
plot(time_min, co_efficiency_history, 'LineWidth',2);
xlabel('Time (min)'); ylabel('Co^{2+} current efficiency (%)'); ylim([0 100]);
title(sprintf('Current efficiency (final %.1f%%)', co_efficiency_history(end))); grid on;

subplot(3,3,6);
plot(time_min, overlimit_count_history, '-o','LineWidth',1.5,'MarkerSize',3);
xlabel('Time (min)'); ylabel('Membranes over-limiting (of 5)'); ylim([-0.3 5.3]);
title('Over-limiting membranes'); grid on;

subplot(3,3,7);
plot(time_min, vol_history, 'LineWidth',1.5);
legend(tank_name, 'Location','eastoutside','FontSize',6);
xlabel('Time (min)'); ylabel('Volume (L)');
title('Tank volumes (water transport)'); grid on;

subplot(3,3,8);
plot(time_min, cake_history, 'LineWidth',1.5); hold on;
plot(time_min, scale_history, '--', 'LineWidth',1.5);
xlabel('Time (min)'); ylabel('Deposit (mg/cm^2)');
title('Membrane fouling (solid) / scaling (dashed)'); grid on;

subplot(3,3,9);
plot(time_min, sum(scale_mass_history,2), 'LineWidth',2);
xlabel('Time (min)'); ylabel('Total scale (mg)');
title('Precipitate/scale formed, stack-wide'); grid on;


%% SECTION 16 — Things to try next
% Now that you've seen every line, here are some easy experiments. Change
% ONE number above, then run the WHOLE script again (green Run button) and
% compare the printed results / charts to before.
%
%   - applied_current_A = 0.5   (higher current -> pushes membrane 3 over
%                                 its limiting current partway through the
%                                 run -- watch the "over-limiting" chart
%                                 and the voltage plateau that follows)
%   - pH0(i_d1) = 2               (an acidic feed -- see Section 3b; changes
%                                   which metal hydroxides can precipitate,
%                                   and how much SO4^2- gets dosed in)
%   - tss0_mgL(i_d1) = 200         (turbid feed -- watch cake build on
%                                   membrane 3 and the fouling voltage rise)
%   - membrane_type{3} = 'Mono-CEM'  (reject metals instead of passing them
%                                      -- see Section 4)
%   - duration_hours = 8            (run twice as long -- does D1 ever
%                                     fully drain? does scale keep growing?)
