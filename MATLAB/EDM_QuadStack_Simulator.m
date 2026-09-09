function R = EDM_QuadStack_Simulator(cfg_ov, make_plots)
%% ╔══════════════════════════════════════════════════════════════════╗
%% ║   EDM Quad Stack — EXTENDED MULTI-METAL / MULTI-PHYSICS MODEL       ║
%% ║   Single-unit by default (N tileable), corrected termination:      ║
%% ║   D2/CEM at the cathode end, C1/CEM at the anode end -- keeps      ║
%% ║   metal cations out of the cathode rinse and Cl- out of the anode. ║
%% ╠══════════════════════════════════════════════════════════════════╣
%% ║  HISTORY: originally ported 1:1 from WebUI/edm_stack_simulator.html║
%% ║  (8 equations, 2 metals, N=1, fixed volumes). This version has     ║
%% ║  SUBSTANTIALLY DIVERGED from that web tool -- 15 species (10       ║
%% ║  metals), N repeating units, monovalent-selective and bipolar      ║
%% ║  membranes, real precipitation/scaling chemistry, water transport, ║
%% ║  membrane fouling, and 19 additional equations (9-27 below). It no ║
%% ║  longer reproduces the web tool's numbers, even at N=1 with the    ║
%% ║  Co/Ca-only feed -- eq. 5 (conductivity), eq. 6 (pH), and the      ║
%% ║  default I_A/GAP_CM values have all changed. If the web tool is    ║
%% ║  updated to match, note that fact there; until then treat the two  ║
%% ║  as independent models built on a shared lineage, not equivalent.  ║
%% ╠══════════════════════════════════════════════════════════════════╣
%% ║  COMPARTMENTS (order): cathode, D2, C2, D1(Feed), C1, anode         ║
%% ║  MEMBRANES (order, between consecutive compartments): CEM,AEM,CEM,AEM,CEM ║
%% ╠══════════════════════════════════════════════════════════════════╣
%% ║  EQUATIONS 1-8 (original core, unchanged in form -- see each       ║
%% ║  parameter/function's own comments for 9-27, listed here by name): ║
%% ║   1  eta_i    = |z_i|*lambda_i*C_i / sum_j(|z_j|*lambda_j*C_j)      ║
%% ║   2  J_i      = eta_i*I / (|z_i|*F)                                 ║
%% ║   3  dC_i     = (J_in - J_out)*dt / V_compartment                   ║
%% ║   4  d[OH-]/dt|cathode = I/F ,  d[H+]/dt|anode = I/F                ║
%% ║   5  kappa = 0.001*sum(lambda_i*|z_i|*C_i) ; V = I*R_total + Vover + Vwatersplit ║
%% ║      5b  ionic-strength (Davies-type) correction to lambda_i        ║
%% ║   6  pH = -log10(activity(C_H)), via full water-equilibrium solve   ║
%% ║   7  eta_Co(t) = min(n_Co,C2(t)/(I*t/(z_Co*F)), 1)*100%             ║
%% ║   8  I_lim = A*sum(|z_i|*F*D_i*C_i/(1000*delta)) ; D_i = R*T*lambda_i/(|z_i|*F^2) ║
%% ║   9-11   suspended solids, deposition, fouling resistance/delta     ║
%% ║   12-14  precipitation: saturation index, relaxation, bulk/surface split ║
%% ║   15     nucleation induction-time barrier                          ║
%% ║   16     dissolution of previously-formed solid                     ║
%% ║   17     water self-ionization re-speciation (closed form)          ║
%% ║   18     per-compartment flow -> boundary layer (Leveque correlation)║
%% ║   19-20  electro-osmotic drag and osmotic water transport           ║
%% ║   21     N repeating units, tiled/lumped by compartment type        ║
%% ║   22     monovalent-selective membrane permeability weighting       ║
%% ║   23     permselectivity / co-ion leakage                           ║
%% ║   24-25  flux decline and conductivity decline reporting            ║
%% ║   26     concentration-polarisation boundary-layer resistance       ║
%% ║   27     fouling-layer resistance derived from porosity             ║
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
%   (S*cm^2/eq, 25degC literature values). These are EQUIVALENT
%   conductivities, lambda(1/|z| * M^z+), which is why the di/trivalent
%   metals all sit in the same 50-70 range as Na+ -- do NOT substitute molar
%   conductivities here, eqs. 1/5/8 all assume the per-equivalent form.
%   D derived via Nernst-Einstein, D_i = R*T*lambda_i / (|z_i|*F^2) --
%   verified against literature: gives D_Na+ = 1.33e-5 and
%   D_Ca2+ = 7.92e-6 cm^2/s, matching tabulated values to 3 sig figs.
%
%   Oxidation states assumed (the dominant aqueous forms in an acidic
%   leachate): Fe as Fe2+, Al and Cr as the trivalent hydrated cations. For
%   a ferric stream override sp_z/sp_lambda for Fe (z=3, lambda=68.0); note
%   that at the pH excursions this model produces, Fe3+/Al3+/Cr3+ would in
%   reality hydrolyse and precipitate as hydroxides -- that chemistry is NOT
%   modelled here, so treat trivalent results as an upper bound on transport.
%
%   Order: metals (Co Ni Fe Al Cr Mn Mg Cu Zn Ca), Na, anions, H/OH.
p.sp_id     = {'Co','Ni','Fe','Al','Cr','Mn','Mg','Cu','Zn','Ca','Na','SO4','Cl','H','OH'};
p.sp_label  = {'Co2+','Ni2+','Fe2+','Al3+','Cr3+','Mn2+','Mg2+','Cu2+','Zn2+','Ca2+','Na+','SO4^2-','Cl-','H+','OH-'};
p.sp_z      = [ 2,     2,     2,     3,     3,     2,     2,     2,     2,     2,     1,   -2,      -1,   1,    -1  ];
p.sp_lambda = [53.0,  49.6,  54.0,  61.0,  67.0,  53.5,  53.1,  53.6,  52.8,  59.5,  50.1, 80.0,   76.3, 349.8, 198.0];
p.sp_D      = (p.Rgas*p.T .* p.sp_lambda) ./ (abs(p.sp_z) * p.F^2);  % cm^2/s
p.nSp       = numel(p.sp_id);

% Named indices, resolved from sp_id so the table above can be reordered or
% extended without hand-editing hardcoded positions.
for k = 1:p.nSp; p.(['i' p.sp_id{k}]) = k; end
p.metalIdx     = find(strcmp(p.sp_id,'Co') | strcmp(p.sp_id,'Ni') | strcmp(p.sp_id,'Fe') | ...
                      strcmp(p.sp_id,'Al') | strcmp(p.sp_id,'Cr') | strcmp(p.sp_id,'Mn') | ...
                      strcmp(p.sp_id,'Mg') | strcmp(p.sp_id,'Cu') | strcmp(p.sp_id,'Zn') | ...
                      strcmp(p.sp_id,'Ca'));
p.conservedIdx = setdiff(1:p.nSp, [p.iH p.iOH]);   % H/OH are source/sink species, excluded

%── Compartments: physical left-to-right order (cathode -> anode) ────────
p.comp_id    = {'cathode','d2','c2','d1','c1','anode'};
p.comp_label = {'Cathode rinse','Sol.1 -> D2','Sol.2 -> C2','Feed -> D1','Sol.3 -> C1','Anode rinse'};
p.nComp = 6;
% Default initial concentrations (mol/L). Assigned by NAMED species index
% rather than as a positional matrix, so adding a species cannot silently
% shift a column. Anything not named below starts at 0.
% Values are unchanged from the WebUI defaults.
p.conc0 = zeros(p.nComp, p.nSp);
p.conc0(:, p.iH)  = 1e-7;                                    % neutral water
p.conc0(:, p.iOH) = 1e-7;
p.conc0(1, [p.iNa p.iSO4])                   = [0.500 0.250];  % cathode rinse, 0.25 M Na2SO4
p.conc0(2, [p.iNa p.iCl])                    = [0.500 0.500];  % D2, 0.5 M NaCl donor
p.conc0(3, [p.iNa p.iCl])                    = [0.050 0.050];  % C2, 0.05 M NaCl seed
p.conc0(4, [p.iCo p.iCa p.iNa p.iSO4 p.iCl]) = ...             % D1 Feed (CoCl2/Na2SO4/CaCl2)
                                               [0.0500 0.0050 0.1082 0.0541 0.1100];
p.conc0(5, [p.iNa p.iSO4])                   = [0.100 0.050];  % C1, 0.05 M Na2SO4 seed
p.conc0(6, [p.iNa p.iSO4])                   = [0.500 0.250];  % anode rinse, 0.25 M Na2SO4

% DEFAULT FEED METALS, written into D1 after the overrides below. The
% balancing anion is then topped up automatically so the feed starts
% electroneutral. The default run therefore exercises all ten metals and
% they all appear on the dashboard.
%   - Original Co/Ca-only WebUI baseline:
%       R = EDM_QuadStack_Simulator(struct('feed_metals', struct()));
%   - Your own assay (replaces this table entirely):
%       R = EDM_QuadStack_Simulator(struct('feed_metals', ...
%               struct('Co',0.05,'Ni',0.02,'Mn',0.01)));
% A metal named here overrides whatever conc0 holds for it; metals not named
% keep their conc0 value. These numbers are a PLACEHOLDER mixed-leachate
% composition -- replace them with your own before drawing conclusions.
p.feed_metals = struct('Co',0.0500, 'Ni',0.0200, 'Mn',0.0100, 'Cu',0.0080, ...
                       'Zn',0.0060, 'Mg',0.0050, 'Ca',0.0050, 'Fe',0.0040, ...
                       'Al',0.0020, 'Cr',0.0010);
p.FEED_BALANCE_ANION = 'Cl';   % anion topped up to charge-balance D1. Set to
                                % 'SO4' for a sulfate-matrix leachate, or ''
                                % to skip balancing (the electroneutrality
                                % check below then warns instead).
p.vol_L = 0.250 * ones(p.nComp,1);   % L, 250 mL each

%── Membranes: fixed by the corrected topology, between consecutive comps ─
%  Five selectable membrane types. Pick them to suit what the stack is FOR:
%    'CEM'       standard cation-exchange. Passes all cations, sharing the
%                current by |z|*lambda*C. No chemical selectivity.
%    'AEM'       standard anion-exchange, the mirror image.
%    'Mono-CEM'  monovalent-SELECTIVE cation-exchange. A thin, like-charged
%                surface layer electrostatically repels multivalent cations,
%                so Na+ passes and Co2+/Ca2+ are largely held back. This is
%                the commercial route to separating monovalent from divalent
%                cations -- use it when the objective is to REJECT the metals
%                at a membrane and take Na+ out instead.
%    'Mono-AEM'  monovalent-selective anion-exchange: Cl- passes, SO4^2- is
%                held back. Use it to split chloride from sulfate.
%    'BPM'       bipolar membrane. Does NOT transport salt at all: under
%                reverse bias it dissociates water at its internal junction,
%                sending H+ toward the cathode side and OH- toward the anode
%                side. This is how you make acid and base in situ instead of
%                dosing them (BMED). It carries the full current as water
%                splitting and costs an extra ~0.9 V junction potential.
%
%  Selecting by objective:
%    maximise metal recovery into C2  -> CEM at the C2|D1 boundary
%    reject metals, bleed off Na+     -> Mono-CEM there instead
%    split Cl- from SO4^2-            -> Mono-AEM at the D1|C1 boundary
%    regenerate acid/base on site     -> BPM against the relevant compartment
p.membrane_types = {'CEM','AEM','CEM','AEM','CEM'};
p.nMem = 5;

%── Monovalent selectivity and bipolar membranes ─────────────────────────
%  A mono-selective membrane's rejection is electrostatic, so it strengthens
%  with counter-ion charge. Modelled as a permeability factor dividing the
%  transport-number weight of eq. 1:
%       S_i = MONO_SELECTIVITY^-(|z_i| - 1)
%  so monovalent ions are unaffected (S = 1), divalent are suppressed by
%  MONO_SELECTIVITY, trivalent by its square. A value of 20 is mid-range for
%  commercial mono-selective grades (published permselectivity ratios span
%  roughly 5-50) -- PLACEHOLDER, measure it for the membrane you actually buy.
%
%  Note this affects only the SHARE of current each ion takes, not the
%  limiting current: selectivity is a property of the membrane, whereas
%  I_lim is set by diffusion through the boundary layer up to its surface,
%  which the membrane cannot influence.
p.MONO_SELECTIVITY = 20;

%── Permselectivity / co-ion leakage (eq. 23) ────────────────────────────
%  No real membrane is perfectly permselective. A fraction (1 - alpha) of
%  the current is carried by CO-IONS travelling the other way: anions
%  leaking through a CEM, cations through an AEM. That leak does no useful
%  separation work -- it is the main reason measured current efficiency
%  falls short of the ideal-membrane prediction, and it worsens as the
%  concentrate builds up (Donnan exclusion weakens with concentration).
%
%  alpha is defined here as the fraction of current carried by counter-ions,
%  so counter-ions carry alpha*I and co-ions carry (1-alpha)*I in the
%  opposite direction. Total charge transfer is still exactly I: a co-ion
%  moving one way is electrically the same as a counter-ion moving the other.
%
%  Set per MEMBRANE TYPE here, or override p.permselectivity with one value
%  per membrane row to give individual membranes different values.
%  Commercial grades run about 0.90-0.99; AEMs typically sit a little below
%  CEMs. These are PLACEHOLDERS -- a membrane's datasheet figure is usually
%  measured in dilute KCl and overstates what you get in a real leachate.
p.PERMSEL_CEM      = 0.95;
p.PERMSEL_AEM      = 0.92;
p.PERMSEL_MONO_CEM = 0.95;
p.PERMSEL_MONO_AEM = 0.92;
p.permselectivity  = [];   % [] = build from the type values above;
                            % or supply a vector, one entry per membrane row.
                            % Set to 1 everywhere for ideal membranes (this
                            % is what the model assumed before eq. 23).

%  Bipolar membrane: the thermodynamic minimum to dissociate water is 0.83 V;
%  real BPMs run somewhat above it. PLACEHOLDER.
p.V_BPM       = 0.9;    % V, junction potential per bipolar membrane
p.N_DRAG_BPM  = 0;      % mol H2O per mol charge. H+ and OH- leave the
                         % junction in OPPOSITE directions, so their drag
                         % largely cancels; 0 is the sensible default.

%── Repeating units (eq. 21) ──────────────────────────────────────────────
%  N_units tiles the quad-stack cell N times between the two electrode
%  rinses. Compartments are LUMPED BY TYPE, exactly as a real stack is
%  manifolded: all N of the D1 chambers share one feed tank, all N of the C2
%  chambers share one product tank, and so on. So the compartment list stays
%  six long however large N is; what grows is the number of membrane COPIES
%  between those tanks.
%
%  All unit cells are electrically in series, so every membrane copy passes
%  the same current I. A tank bounded by N copies therefore sees N times the
%  ionic flux of a single cell -- which is the whole point of a stack, and
%  why voltage rises with N while treatment rate per tank rises too.
%
%  N = 1 reproduces the single-cell model exactly.
%
%  Tiling (from the corrected D2-start topology):
%     cathode |CEM| D2 |AEM| C2 |CEM| D1 |AEM| C1 |CEM| D2(next unit) ...
%     ... and the last C1 closes onto the anode rinse through a CEM.
%  The C1 -> D2 handoff is why the topology cannot be a simple linear chain:
%  it links the last compartment of one unit back to the second of the next.
p.N_units = 1;

%  Optional explicit override: an M x 4 table [leftComp, rightComp,
%  isCEM(1)/isAEM(0), multiplicity]. Leave empty to have it built from
%  membrane_types (N = 1) or from the tiling above (N > 1). Supply your own
%  for a topology neither of those covers.
p.membranes = [];

%── Stack design / operating point ────────────────────────────────────────
p.area_cm2   = 10;     % cm^2, membrane area
p.I_A        = 0.30;   % A, applied constant current (30 mA/cm^2).
                        % Was 0.5 A (50 mA/cm^2), which sits close enough to
                        % the stack's limiting current that modest depletion
                        % of the feed tipped membrane 3 over it at t = 73 min.
                        % Past that point the current can no longer be carried
                        % by metal ions, so it splits water instead, the feed
                        % turns alkaline, and the metals drop out as hydroxide
                        % sludge -- 19% of the feed metal, lost as scale.
                        % At 0.30 A no membrane exceeds its limit and that
                        % whole failure mode disappears.
p.duration_h = 4;      % h, run duration
p.dt_min     = 1;      % min, timestep

%── Resistance / overpotential placeholders (flagged, not measured) ──────
p.AREAL_R_MEMBRANE = 5;      % ohm*cm^2. Scalar applies to every membrane;
                              % supply a VECTOR (one entry per membrane row)
                              % to give individual membranes their own value,
                              % e.g. a thin low-resistance CEM beside a
                              % thicker AEM. Resolved into p.areal_R below.
p.GAP_CM           = 0.1;    % cm, compartment spacer thickness (1 mm).
                              % Was 0.5 cm (5 mm), which sits well outside the
                              % 0.5-2 mm used in commercial ED equipment and
                              % inflated the solution resistance roughly
                              % fivefold. Measure your own cell and set this:
                              % voltage scales almost linearly with it.
p.V_OVERPOTENTIAL  = 1.5;    % V, fixed electrode overpotential
%── Ionic-strength correction for conductivity (eq. 5b) ──────────────────
%  sp_lambda holds LIMITING (infinite-dilution) conductivities. Used raw,
%  eq. 5 overestimates conductivity badly once the solution is concentrated:
%  0.5 M NaCl comes out at 0.0632 S/cm against a measured 0.0452, a 40%
%  error, which understates compartment resistance and so understates stack
%  voltage. Ion-ion interaction slows every ion down as concentration rises.
%
%  Corrected with a Kohlrausch-type empirical factor of the same shape as
%  the Davies expression:
%       lambda_i(I) = lambda_i(0) * (1 - b*sqrt(I)/(1+sqrt(I)))
%  Fitted to NaCl at 25 degC, this holds to ~1% from 1 mM to 1 M:
%       I (M)          0.001   0.01    0.1     0.5     1.0
%       measured       123.7   118.5   106.7    90.4    85.0   (S cm^2/eq)
%       corrected      123.8   118.7   106.0    91.4    84.1
%
%  Fitted on a 1:1 electrolyte; multivalent solutions deviate more, so treat
%  it as an engineering correction, not a first-principles result. It cancels
%  out of the transport numbers (eq. 1 is a ratio and every lambda is scaled
%  alike), so it moves voltage only. The diffusivities behind I_lim still
%  come from the limiting values.
p.KAPPA_IONIC_CORRECTION = true;
p.KAPPA_B = 0.67;      % fitted coefficient of the correction above

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

%── Concentration polarisation resistance (eq. 26) ───────────────────────
%  As the current approaches I_lim, the counter-ion concentration AT the
%  membrane surface falls toward zero, and that depleted film is highly
%  resistive. Film theory gives the surface concentration as
%       c_s / c_bulk = 1 - I / I_lim
%  so the boundary layer's conductivity collapses as the limit is neared and
%  its resistance diverges. That divergence is what produces the plateau in
%  a measured polarisation curve, and hence the minimum in a Cowan plot.
%
%  Without this term the model treats over-limiting as a fixed overpotential
%  per membrane, which makes V rise in steps rather than through a plateau --
%  it cannot reproduce the shape of the standard ED characterisation
%  measurement, and a Cowan construction on it finds no minimum.
%
%  Only the EXCESS over the undepleted film is added, so the term vanishes
%  at low current and does not double-count the compartment resistance the
%  gap already carries.
p.ENABLE_CONC_POLARISATION = true;
p.CS_FLOOR = 1e-3;    % floor on c_surface/c_bulk, bounding the divergence

%── Water transport through the membranes (eq. 19-20) ────────────────────
%  Compartment volumes are NOT constant in a real stack. Water crosses the
%  membranes by two mechanisms, and the usual lab observation -- one
%  compartment steadily gaining while its neighbour loses, looking like a
%  leak -- is these two, not a failed seal:
%
%   (19) ELECTRO-OSMOTIC DRAG. Every ion crossing carries its hydration
%        shell with it, so water is pulled in the direction the current
%        carries the counter-ion: toward the cathode through a CEM, toward
%        the anode through an AEM. Scales with charge passed, so it is
%        unavoidable and grows with current and run time.
%
%   (20) OSMOSIS. Water moves from the dilute side toward the concentrated
%        side, i.e. it works to UNDO the separation you are performing.
%        Driven by the osmotic pressure difference (van 't Hoff).
%
%  These move water between compartments; total stack volume is conserved
%  (checked in runEDM). Evaporation and actual leaks are NOT modelled -- if
%  your measured total volume falls, that part is a genuine leak.
p.ENABLE_WATER_TRANSPORT = true;

p.N_DRAG_CEM = 6.0;    % mol H2O per mol charge through a CEM. Literature is
                        % ~4-8 for a hydrated Na+ form; higher for divalent
                        % and for looser membranes. MEASURE IT: this sets the
                        % whole drag term and nothing here calibrates it.
p.N_DRAG_AEM = 4.0;    % mol H2O per mol charge through an AEM (Cl- carries
                        % a smaller shell than Na+, hence the lower value)
p.LP_WATER   = 2e-6;   % cm/s/bar, membrane hydraulic/osmotic permeability.
                        % Spans 5e-7 to 1e-5 across commercial membranes --
                        % a decade of spread that maps to 1.6 vs 32 mL over a
                        % 4 h run, so it is worth measuring rather than guessing.
p.V_WATER_MOLAR_CM3 = 18.0;   % cm^3/mol, molar volume of liquid water
p.MIN_VOL_FRACTION  = 0.05;   % floor on compartment volume, as a fraction of
                               % its initial value, to keep the ODEs finite if
                               % a compartment is being pumped dry

%── Per-compartment flow / hydrodynamics (eq. 18) ────────────────────────
%  Flow rate is set PER COMPARTMENT, because in a real stack the chambers
%  are on separate hydraulic loops and rarely run at the same rate: the
%  dilute loops are usually pushed hard to hold up mass transfer, while the
%  rinses are often much slower.
%
%  What the flow actually changes is the Nernst boundary layer. Crossflow
%  thins it, which raises the limiting current, delays water splitting, and
%  shears deposited particles off the membrane. A single global DELTA_CM
%  forces every chamber to share one boundary layer regardless of how hard
%  it is being pumped; this replaces that with a per-chamber delta computed
%  from the chamber's own velocity.
%
%  delta comes from the Leveque correlation for developing laminar flow in a
%  thin channel,  Sh = 1.85 (Re Sc d_h / L)^(1/3),  delta = d_h / Sh, so
%  delta ~ velocity^(-1/3). Valid for laminar, hydrodynamically developing
%  flow -- the normal regime in a spacer-filled ED cell, though a real
%  turbulence-promoting spacer does better than this (it is a conservative
%  estimate of mass transfer, i.e. a pessimistic I_lim).
%
%  DEFAULT IS ZERO = no forced flow: every chamber falls back to the fixed
%  DELTA_CM, reproducing the previous behaviour exactly. Set a flow to
%  activate the correlation for that chamber, e.g.
%     cfg.flow_Lpm = [0.2; 1.0; 1.0; 1.5; 1.0; 0.2];   % cathode..anode
p.flow_Lpm = zeros(p.nComp,1);    % L/min through each compartment

p.CHANNEL_LENGTH_CM   = [];       % cm, flow-path length. [] = sqrt(area),
                                   % i.e. assume a square cell.
p.SPACER_POROSITY     = 0.75;     % open fraction of the channel cross-section
p.KIN_VISC_CM2_S      = 0.00893;  % cm^2/s, water at 25 degC
p.D_REF_CM2_S         = [];       % cm^2/s reference diffusivity for the Schmidt
                                   % number. [] = use Na+. delta depends on it
                                   % only as Sc^(-1/3), so a single reference
                                   % across species is a mild approximation.
p.SHEAR_DEP_REF_CM_S  = 5.0;      % cm/s at which crossflow halves the particle
                                   % deposition velocity: v_dep is scaled by
                                   % 1/(1+u/u_ref). Placeholder -- the real
                                   % dependence is on wall shear and particle
                                   % size, and nothing here calibrates it.

%── Suspended solids & turbidity (eq. 9-11) ───────────────────────────────
%  TSS is NOT an ionic species: particles carry no Faradaic current and do
%  not cross a membrane, so they are held as their own per-compartment state
%  rather than as a column of conc. What they DO is deposit on the membrane
%  they are driven against, forming a cake that (a) adds areal resistance
%  and (b) thickens the diffusion boundary layer, lowering I_lim. Both
%  effects push the stack toward water splitting -- fouling shows up as a
%  rising voltage and earlier over-limiting, which is the whole reason to
%  track feed turbidity in the first place.
%
%  Turbidity is treated as an OBSERVABLE of TSS, not an independent state:
%  it is converted to TSS on input and back out for reporting. The
%  correlation is site-specific (1-3 mg/L per NTU is the usual range for
%  natural waters and is strongly dependent on particle size and colour) --
%  MEASURE IT for your own stream rather than trusting this default.
p.MGL_PER_NTU = 2.0;    % mg/L of TSS per NTU -- placeholder correlation

%── Water self-ionization and pH (eq. 6, 17) ─────────────────────────────
%  H+ and OH- are transported independently by the membranes and generated
%  independently at the electrodes, so nothing in the transport step stops a
%  compartment ending up with high H+ AND high OH- at the same time. That is
%  chemically impossible: they neutralize on contact. Without this step the
%  C1 compartment finishes at [H+] = 0.16 M and [OH-] = 0.11 M, an ion
%  product 1e12 times Kw.
%
%  So after every transport step each compartment is re-speciated onto the
%  water equilibrium C_H * C_OH = Kw, holding the NET (C_H - C_OH) fixed --
%  net is what transport and the electrodes actually determine, and holding
%  it fixed conserves charge exactly. pH is then a real quantity read
%  straight off C_H, not inferred from a net-charge workaround.
p.ENFORCE_WATER_EQUILIBRIUM = true;

%  pH is defined on activity, not concentration. With Davies coefficients
%  and ionic strengths around 0.5 M, gamma_H is ~0.7, i.e. ~0.15 pH units --
%  small but not nothing. Set false to report -log10(concentration).
p.PH_USE_ACTIVITY = true;

% Initial TSS per compartment (mg/L), same row order as conc0. Defaults to
% zero everywhere, so a run with no solids reproduces the clean-water model
% exactly. Set p.turbidity0_NTU instead if you measured NTU (see below).
p.tss0_mgL       = zeros(p.nComp,1);
p.turbidity0_NTU = [];   % optional alternative input; converted via MGL_PER_NTU

p.PARTICLE_DENSITY_G_CM3 = 2.5;   % g/cm^3, typical mineral floc/silt
p.DEPOSITION_VEL_CM_S    = 1e-5;  % cm/s, particle drift velocity onto the
                                    % membrane it is driven against (eq. 10).
                                    % Lumps electrophoresis + settling +
                                    % concentration polarisation -- a single
                                    % fitted placeholder, NOT a measurement.
p.TSS_CHARGE_SIGN        = -1;     % sign of the particle zeta potential.
                                    % -1 (typical for mineral colloids and
                                    % NOM) drives particles toward the anode,
                                    % so they foul the membrane on a
                                    % compartment's ANODE side. +1 flips it.
%── Deposit layer properties, DERIVED rather than assumed (eq. 27) ───────
%  R_CAKE_AREAL and R_SCALE_AREAL were lumped constants with no published
%  counterpart, so they looked unfittable. They are not: a deposit of
%  loading w (mg/cm^2) and porosity eps is simply a layer
%       thickness  t   = w / (rho_particle * (1 - eps))
%  whose pores hold electrolyte, so by the Bruggeman relation it conducts at
%       kappa_eff      = kappa_local * eps^1.5
%  and its areal resistance is t/kappa_eff. The same thickness also adds
%  directly to the diffusion boundary layer. So all three fouling-resistance
%  constants follow from ONE physically meaningful quantity per deposit --
%  its porosity -- which does have published ranges: roughly 0.3-0.7 for a
%  loose particle cake and 0.05-0.2 for dense crystalline scale.
%
%  This matters because the assumed constants implied porosities of 0.6% and
%  0.3%, i.e. essentially impermeable solid, and they dominated the predicted
%  stack voltage. Deriving them removes two arbitrary numbers and replaces
%  them with a quantity that can be measured or looked up.
%
%  Set DERIVE_FOULING_RESISTANCE = false to fall back on the lumped
%  constants below, e.g. when fitting them directly to a measured voltage rise.
p.DERIVE_FOULING_RESISTANCE = true;
p.CAKE_POROSITY  = 0.50;   % loose particulate cake; published range 0.3-0.7
p.SCALE_POROSITY = 0.10;   % dense crystalline scale; published range 0.05-0.2

p.R_CAKE_AREAL           = 20;     % ohm*cm^2 per (mg/cm^2) of cake -- the
                                    % specific cake resistance (eq. 11).
p.CAKE_DELTA_FACTOR      = 0.05;   % cm^2/mg, boundary-layer thickening per
                                    % unit cake load: delta_eff = delta*(1+f*w)

%── Mineral scaling (eq. 12-14) ───────────────────────────────────────────
%  This stack drives compartments to pH 1 and pH 13.5 while concentrating
%  metals -- exactly the conditions that precipitate gypsum and metal
%  hydroxides. Scaling is modelled as: (12) an ion-activity product against
%  Ksp, (13) first-order relaxation of the supersaturation toward
%  equilibrium, and (14) a split of the new solid between the bulk (where it
%  becomes TSS and then feeds the existing deposition/fouling model) and
%  direct heterogeneous growth on the adjacent membrane surfaces.
%
%  Ksp values are 25 degC thermodynamic solubility products from standard
%  tables. Treat them as order-of-magnitude for this purpose: real scaling
%  is governed by nucleation induction time, seed availability and
%  crystal habit, none of which are modelled. This tells you WHERE and WHEN
%  a stream goes supersaturated -- not how fast a real crystal grows.
p.ENABLE_SCALING = true;

%  Solid | cation | nu_cat | anion | nu_an | Ksp | molar mass (g/mol)
p.scale_name   = {'CaSO4.2H2O','Ca(OH)2','Co(OH)2','Ni(OH)2','Fe(OH)2', ...
                  'Al(OH)3','Cr(OH)3','Mn(OH)2','Mg(OH)2','Cu(OH)2','Zn(OH)2'};
p.scale_cat    = {'Ca','Ca','Co','Ni','Fe','Al','Cr','Mn','Mg','Cu','Zn'};
p.scale_anion  = {'SO4','OH','OH','OH','OH','OH','OH','OH','OH','OH','OH'};
p.scale_nuCat  = [ 1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1  ];
p.scale_nuAn   = [ 1,    2,    2,    2,    2,    3,    3,    2,    2,    2,    2  ];
p.scale_Ksp    = [3.14e-5, 5.5e-6, 5.9e-15, 5.5e-16, 4.9e-17, 3.0e-34, ...
                  6.3e-31, 1.9e-13, 5.6e-12, 2.2e-20, 3.0e-17];
p.scale_MW     = [172.17, 74.09, 92.95, 92.71, 89.86, 78.00, 103.02, 88.95, 58.32, 97.56, 99.42];

p.K_PRECIP_PER_S        = 1e-3;  % 1/s, relaxation rate toward equilibrium.
                                  % 1e-3 gives a ~17 min time constant --
                                  % fast vs. a 4 h run but not instantaneous.
                                  % Raise toward 1 for equilibrium-limited
                                  % behaviour, lower for kinetically hindered.
p.SCALE_SURFACE_FRACTION = 0.30;  % fraction of new solid that nucleates
                                   % directly on the membrane rather than in
                                   % the bulk. Membrane surfaces are the
                                   % preferred nucleation site (highest local
                                   % supersaturation), hence >0.
p.R_SCALE_AREAL          = 50;    % ohm*cm^2 per (mg/cm^2) -- dense crystalline
                                   % scale is more resistive than a loose
                                   % particle cake (R_CAKE_AREAL).
%── Nucleation barrier (eq. 15) ───────────────────────────────────────────
%  Crystals do not appear the instant SI crosses zero. Primary nucleation
%  needs a finite driving force AND a finite induction time; classical
%  nucleation theory gives t_ind ~ exp(B/(ln S)^2), i.e. it falls steeply as
%  supersaturation rises. Modelled here as: no NEW solid forms in a
%  compartment until SI exceeds SI_NUCLEATION and the compartment has spent
%  t_ind = TAU_NUC_REF_S / SI^2 above that threshold. Once any seed of that
%  solid is present (suspended or on an adjacent membrane), the barrier is
%  gone -- secondary growth on existing crystal needs no nucleation -- and
%  growth resumes whenever SI > 0.
%
%  This is what stops the model precipitating everything the moment the pH
%  moves. Gypsum in particular tolerates substantial supersaturation before
%  it nucleates, which is why scale-control chemistry works at all.
p.SI_NUCLEATION  = 0.5;      % log10 units of supersaturation needed to nucleate
p.TAU_NUC_REF_S  = 1800;     % s, induction time at SI = 1 (scales as 1/SI^2)

%── Dissolution (eq. 16) ──────────────────────────────────────────────────
%  Solid formed earlier redissolves when its compartment later becomes
%  undersaturated -- which happens constantly here, because the membranes
%  keep moving ions and the electrodes keep moving pH. Without this, scale
%  only ever accumulates and the model badly over-predicts fouling.
%  Dissolution draws on solid in contact with that compartment's solution:
%  suspended crystals first (highest specific area), then the particulate
%  cake, then dense in-situ scale last.
p.ENABLE_DISSOLUTION = true;
p.K_DISSOL_PER_S     = 5e-4;  % 1/s, dissolution relaxation rate. Slower than
                               % K_PRECIP_PER_S: dissolving a formed crystal
                               % is transport-limited at its surface, whereas
                               % precipitation can proceed throughout the bulk.

p.USE_DAVIES             = true;  % Davies activity coefficients in eq. 12.
                                   % Ionic strength here reaches ~0.5 M, the
                                   % upper edge of Davies validity -- beyond
                                   % that it is an extrapolation. Set false
                                   % to fall back on concentrations.

%── Apply user overrides (any field above may be overridden) ─────────────
if ~isempty(cfg_ov)
    fns = fieldnames(cfg_ov);
    for k = 1:numel(fns); p.(fns{k}) = cfg_ov.(fns{k}); end
end

%── Re-derive anything the overrides may have invalidated ────────────────
if numel(p.sp_z)~=p.nSp || numel(p.sp_lambda)~=p.nSp || size(p.conc0,2)~=p.nSp
    error('EDM:speciesSizeMismatch', ...
        ['Species arrays are inconsistent (nSp=%d, z=%d, lambda=%d, conc0 cols=%d). ' ...
         'If you override sp_id/sp_z/sp_lambda you must override conc0 to match.'], ...
        p.nSp, numel(p.sp_z), numel(p.sp_lambda), size(p.conc0,2));
end
p.sp_D = (p.Rgas*p.T .* p.sp_lambda) ./ (abs(p.sp_z) * p.F^2);   % overrides may change lambda/z/T

%── Resolve the TSS / turbidity input pair ───────────────────────────────
% Accept either; turbidity is converted to TSS, which is the state variable.
p.tss0_mgL = p.tss0_mgL(:);
if ~isempty(p.turbidity0_NTU)
    ntu = p.turbidity0_NTU(:);
    if numel(ntu) ~= p.nComp
        error('EDM:turbiditySize','turbidity0_NTU must have %d entries (one per compartment).', p.nComp);
    end
    if any(p.tss0_mgL > 0)
        warning('EDM:tssAndTurbidity', ...
            ['Both tss0_mgL and turbidity0_NTU were supplied. Using turbidity0_NTU ' ...
             '(x %.2f mg/L per NTU) and discarding tss0_mgL.'], p.MGL_PER_NTU);
    end
    p.tss0_mgL = ntu * p.MGL_PER_NTU;
end
if numel(p.tss0_mgL) ~= p.nComp
    error('EDM:tssSize','tss0_mgL must have %d entries (one per compartment).', p.nComp);
end
if any(p.tss0_mgL < 0)
    error('EDM:tssNegative','tss0_mgL cannot be negative.');
end
p.hasSolids = any(p.tss0_mgL > 0);

%── Structural consistency of the stack configuration ────────────────────
% The transport engine is generic in the number of compartments and the
% membrane sequence, so a different stack can be simulated by overriding
% comp_id/nComp/nMem/membrane_types/conc0/vol_L together. Catch a mismatched
% set here with a clear message rather than deep inside the reporting layer.
if isempty(p.membranes) && p.N_units <= 1 && p.nMem ~= p.nComp - 1
    error('EDM:topology', ...
        ['For the default linear chain, nMem must equal nComp-1 (%d compartments ' ...
         'in a line have %d membranes between them); got nComp=%d, nMem=%d. ' ...
         'Supply p.membranes explicitly for a non-chain topology.'], ...
        p.nComp, p.nComp-1, p.nComp, p.nMem);
end
chk = {'comp_id',numel(p.comp_id); 'comp_label',numel(p.comp_label); ...
       'vol_L',numel(p.vol_L); 'conc0 rows',size(p.conc0,1)};
for k = 1:size(chk,1)
    if chk{k,2} ~= p.nComp
        error('EDM:topology','%s has %d entries but nComp = %d.', chk{k,1}, chk{k,2}, p.nComp);
    end
end
if isempty(p.membranes)
    if p.N_units <= 1
        if numel(p.membrane_types) ~= p.nMem
            error('EDM:topology','membrane_types has %d entries but nMem = %d.', ...
                numel(p.membrane_types), p.nMem);
        end
    elseif ~ismember(numel(p.membrane_types), [5 6])
        % The N-unit tiling has six rows: cathode|D2, D2|C2, C2|D1, D1|C1,
        % C1|D2(next), C1|anode. Give six to choose a type at each position,
        % or leave the default five and get the standard arrangement.
        error('EDM:topology', ...
            ['With N_units > 1, membrane_types must have 6 entries (one per ' ...
             'tiling row: cathode|D2, D2|C2, C2|D1, D1|C1, C1|D2next, C1|anode) ' ...
             'or be left at the 5-entry default; got %d.'], numel(p.membrane_types));
    end
end
validTypes = {'CEM','AEM','Mono-CEM','Mono-AEM','BPM'};
bad = ~ismember(p.membrane_types, validTypes);
if any(bad)
    error('EDM:membraneType','membrane_types must be one of %s; found "%s".', ...
        strjoin(validTypes,'/'), p.membrane_types{find(bad,1)});
end

%── Resolve per-compartment hydrodynamics (eq. 18) ───────────────────────
p.flow_Lpm = p.flow_Lpm(:);
if numel(p.flow_Lpm) ~= p.nComp
    error('EDM:flowSize','flow_Lpm must have %d entries (one per compartment).', p.nComp);
end
if any(p.flow_Lpm < 0)
    error('EDM:flowNegative','flow_Lpm cannot be negative.');
end
if isempty(p.CHANNEL_LENGTH_CM); p.CHANNEL_LENGTH_CM = sqrt(p.area_cm2); end
if isempty(p.D_REF_CM2_S);       p.D_REF_CM2_S = p.sp_D(p.iNa); end

chanWidth = p.area_cm2 / p.CHANNEL_LENGTH_CM;                 % cm
Axs       = chanWidth * p.GAP_CM * p.SPACER_POROSITY;          % cm^2 open cross-section
dh        = 2 * p.GAP_CM * p.SPACER_POROSITY;                  % cm, thin-slit hydraulic diameter
Sc        = p.KIN_VISC_CM2_S / p.D_REF_CM2_S;

p.velocity_cm_s = (p.flow_Lpm * 1000/60) / Axs;                % cm/s in each chamber
p.delta_cm      = zeros(p.nComp,1);
p.Re            = zeros(p.nComp,1);
for i = 1:p.nComp
    u = p.velocity_cm_s(i);
    if u <= 0
        p.delta_cm(i) = p.DELTA_CM;         % no forced flow: fixed placeholder
        continue
    end
    p.Re(i) = u * dh / p.KIN_VISC_CM2_S;
    Sh = 1.85 * (p.Re(i) * Sc * dh / p.CHANNEL_LENGTH_CM)^(1/3);
    % delta cannot exceed half the channel -- past that the "boundary layer"
    % has consumed the bulk and the film model no longer applies.
    p.delta_cm(i) = min(dh/max(Sh,eps), 0.5*p.GAP_CM);
end
% Crossflow shears particles off, so the deposition velocity falls with u.
p.dep_vel_cm_s = p.DEPOSITION_VEL_CM_S ./ (1 + p.velocity_cm_s/p.SHEAR_DEP_REF_CM_S);

% Which compartment feeds each membrane's boundary layer: the SOURCE side,
% since that is the side being depleted and where polarisation happens.
% Membrane topology table: [left, right, isCEM, multiplicity].
% Type code in column 3: 1 = cation-conducting (CEM or Mono-CEM),
% 0 = anion-conducting (AEM or Mono-AEM), 2 = bipolar (no salt transport).
typeCode = @(t) 2*strcmp(t,'BPM') + 1*(strcmp(t,'CEM') || strcmp(t,'Mono-CEM'));
if isempty(p.membranes)
    if p.N_units <= 1
        % Linear chain: membrane k joins compartments k and k+1, one copy each.
        p.membranes = zeros(p.nMem,4);
        for k = 1:p.nMem
            p.membranes(k,:) = [k, k+1, typeCode(p.membrane_types{k}), 1];
        end
    else
        % N-unit tiling. Requires the standard six lumped compartments.
        need = {'cathode','d2','c2','d1','c1','anode'};
        if ~isequal(p.comp_id(:)', need)
            error('EDM:tiling', ...
                ['N_units > 1 uses the standard lumped compartments %s. For a ' ...
                 'different arrangement, supply p.membranes explicitly as an ' ...
                 'M x 4 table [left right isCEM multiplicity].'], strjoin(need,'/'));
        end
        N = p.N_units;
        iCat=1; iD2=2; iC2=3; iD1=4; iC1=5; iAn=6;
        % Row order: cathode|D2, D2|C2, C2|D1, D1|C1, C1|D2(next), C1|anode.
        % Supply membrane_types with one entry per row to choose the type at
        % each position; otherwise the standard all-CEM/AEM set is used.
        tileTypes = {'CEM','AEM','CEM','AEM','CEM','CEM'};
        if numel(p.membrane_types) == 6
            tileTypes = p.membrane_types;
        end
        rows = [ ...
            iCat, iD2, typeCode(tileTypes{1}), 1;      % terminal: cathode rinse | D2
            iD2,  iC2, typeCode(tileTypes{2}), N;      % D2 | C2
            iC2,  iD1, typeCode(tileTypes{3}), N;      % C2 | D1
            iD1,  iC1, typeCode(tileTypes{4}), N;      % D1 | C1
            iC1,  iD2, typeCode(tileTypes{5}), N-1;    % C1 | D2 of the NEXT unit
            iC1,  iAn, typeCode(tileTypes{6}), 1];     % terminal: C1 | anode rinse
        keep = rows(:,4) > 0;                           % drop N=1's empty handoff
        p.membranes = rows(keep,:);
        p.membrane_types = tileTypes(keep);
    end
end
p.nMem     = size(p.membranes,1);      % number of membrane ROWS (not copies)
p.mem_mult = p.membranes(:,4)';        % copies of each row
p.nMemCopies = sum(p.mem_mult);        % total membranes in the stack

% Which compartment feeds each membrane's boundary layer: the SOURCE side,
% since that is the side being depleted and where polarisation happens.
p.mem_source = zeros(1, p.nMem);
p.mem_isMono = false(1, p.nMem);
for k = 1:p.nMem
    p.mem_isMono(k) = startsWith(p.membrane_types{k}, 'Mono-');
    switch p.membranes(k,3)
        case 1;  p.mem_source(k) = p.membranes(k,2);   % cations come from the right
        case 0;  p.mem_source(k) = p.membranes(k,1);   % anions come from the left
        otherwise, p.mem_source(k) = p.membranes(k,1); % BPM: no salt source; index kept valid
    end
end
p.nBPMCopies = sum(p.membranes(p.membranes(:,3)==2, 4));

% Per-membrane clean areal resistance (ohm*cm^2).
if isscalar(p.AREAL_R_MEMBRANE)
    p.areal_R = p.AREAL_R_MEMBRANE * ones(1, p.nMem);
else
    p.areal_R = p.AREAL_R_MEMBRANE(:)';
    if numel(p.areal_R) ~= p.nMem
        error('EDM:arealR', ...
            'AREAL_R_MEMBRANE must be scalar or have one entry per membrane row (%d); got %d.', ...
            p.nMem, numel(p.areal_R));
    end
end

% Per-membrane permselectivity (eq. 23), from type unless overridden.
if isempty(p.permselectivity)
    p.permselectivity = zeros(1, p.nMem);
    for k = 1:p.nMem
        switch p.membrane_types{k}
            case 'CEM';      p.permselectivity(k) = p.PERMSEL_CEM;
            case 'AEM';      p.permselectivity(k) = p.PERMSEL_AEM;
            case 'Mono-CEM'; p.permselectivity(k) = p.PERMSEL_MONO_CEM;
            case 'Mono-AEM'; p.permselectivity(k) = p.PERMSEL_MONO_AEM;
            case 'BPM';      p.permselectivity(k) = 1;   % carries no salt either way
        end
    end
else
    p.permselectivity = p.permselectivity(:)';
    if numel(p.permselectivity) == 1
        p.permselectivity = p.permselectivity * ones(1, p.nMem);
    end
    if numel(p.permselectivity) ~= p.nMem
        error('EDM:permsel', ...
            'permselectivity must be scalar or have one entry per membrane row (%d); got %d.', ...
            p.nMem, numel(p.permselectivity));
    end
end
if any(p.permselectivity < 0 | p.permselectivity > 1)
    error('EDM:permsel','permselectivity must lie in [0,1].');
end
p.hasFlow = any(p.flow_Lpm > 0);

%── Single-pass conversion: does the well-mixed-tank assumption hold? ─────
%  The model treats each compartment as ONE well-mixed tank. Physically the
%  solution is pumped tank -> cell -> tank, so that is only valid if a
%  parcel of liquid changes little during one pass through the cell. If the
%  pump is too slow the cell outlet differs sharply from the tank, a
%  concentration profile develops along the channel, and a single lumped
%  tank no longer represents it.
%
%  Fraction of the ionic equivalents in a parcel converted in one pass:
%      X = I / (F * Q * C_eq)
%  (the cell hold-up volume cancels: longer residence in a bigger cell is
%  exactly offset by there being more solution to convert). So it depends
%  only on current, pump rate and concentration -- nothing about geometry.
p.cellVolume_L    = p.area_cm2 * p.GAP_CM * p.SPACER_POROSITY / 1000;
p.residenceTime_s = zeros(p.nComp,1);
p.turnoverTime_s  = zeros(p.nComp,1);
p.singlePassConv  = zeros(p.nComp,1);
for i = 1:p.nComp
    Q_Ls = p.flow_Lpm(i)/60;
    if Q_Ls <= 0
        p.residenceTime_s(i) = Inf; p.turnoverTime_s(i) = Inf;
        p.singlePassConv(i)  = NaN;      % undefined without a pump rate
        continue
    end
    Ceq = sum(abs(p.sp_z) .* max(p.conc0(i,:),0));    % equivalents/L
    p.residenceTime_s(i) = p.cellVolume_L / Q_Ls;
    p.turnoverTime_s(i)  = p.vol_L(i)     / Q_Ls;
    if Ceq > 0
        p.singlePassConv(i) = p.I_A / (p.F * Q_Ls * Ceq);
    else
        p.singlePassConv(i) = Inf;        % no carriers: any current is total conversion
    end
end

%── Resolve the scaling table to species indices ─────────────────────────
p.nScale = numel(p.scale_name);
p.scale_iCat = zeros(1,p.nScale); p.scale_iAn = zeros(1,p.nScale);
for k = 1:p.nScale
    ic = find(strcmp(p.sp_id, p.scale_cat{k}),   1);
    ia = find(strcmp(p.sp_id, p.scale_anion{k}), 1);
    if isempty(ic) || isempty(ia)
        error('EDM:scaleSpecies','Scale "%s" references a species not in sp_id.', p.scale_name{k});
    end
    p.scale_iCat(k) = ic; p.scale_iAn(k) = ia;
end
if ~(numel(p.scale_Ksp)==p.nScale && numel(p.scale_MW)==p.nScale && ...
     numel(p.scale_nuCat)==p.nScale && numel(p.scale_nuAn)==p.nScale)
    error('EDM:scaleTableSize','Scaling table columns have inconsistent lengths.');
end

%── feed_metals shortcut: drop metals straight into D1 (compartment 4) ────
d1 = find(strcmp(p.comp_id,'d1'));
mfns = fieldnames(p.feed_metals);
for k = 1:numel(mfns)
    idx = find(strcmp(p.sp_id, mfns{k}), 1);
    if isempty(idx)
        error('EDM:unknownSpecies','feed_metals: "%s" is not in sp_id.', mfns{k});
    end
    p.conc0(d1, idx) = p.feed_metals.(mfns{k});
end

% Charge-balance the feed. Changing the metal loading changes the cation
% equivalents, so the balancing anion has to follow or D1 starts with a net
% charge -- which the transport model would happily propagate through the
% whole run.
if ~isempty(mfns) && ~isempty(p.FEED_BALANCE_ANION)
    jb = find(strcmp(p.sp_id, p.FEED_BALANCE_ANION), 1);
    if isempty(jb)
        error('EDM:balanceAnion','FEED_BALANCE_ANION "%s" is not in sp_id.', p.FEED_BALANCE_ANION);
    end
    netEq = sum(p.sp_z .* p.conc0(d1,:));      % >0 means excess cation
    if netEq > 0
        p.conc0(d1,jb) = p.conc0(d1,jb) + netEq/abs(p.sp_z(jb));
        p.feedBalanceAdded_M = netEq/abs(p.sp_z(jb));
    elseif netEq < 0
        warning('EDM:anionExcess', ...
            ['Feed has %.4f eq/L of ANION excess; %s cannot be reduced below zero ' ...
             'to balance it. Lower the background anions in conc0.'], -netEq, p.sp_id{jb});
        p.feedBalanceAdded_M = 0;
    else
        p.feedBalanceAdded_M = 0;
    end
else
    p.feedBalanceAdded_M = 0;
end

p.dt_sec   = p.dt_min * 60;
p.tEnd_sec = p.duration_h * 3600;

%── Electroneutrality check on the initial condition (warn, don't block) ──
for i = 1:p.nComp
    netEq = sum(p.sp_z .* p.conc0(i,:));
    totEq = sum(abs(p.sp_z) .* p.conc0(i,:));
    if totEq > 1e-9 && abs(netEq)/totEq > 0.01
        warning('EDM:chargeImbalance', ...
            '%s initial charge imbalance %+.4f eq/L (%.1f%% of total) -- add counter-ions.', ...
            p.comp_label{i}, netEq, 100*abs(netEq)/totEq);
    end
end

%── Design summary ────────────────────────────────────────────────────────
fprintf('══════════════════════════════════════════════\n');
fprintf('  EDM Quad Stack — Design\n');
fprintf('══════════════════════════════════════════════\n');
fprintf('  Membrane area      : %.1f cm^2\n', p.area_cm2);
fprintf('  Current (constant) : %.3f A\n', p.I_A);
fprintf('  Duration / dt      : %.1f h / %.1f min\n', p.duration_h, p.dt_min);
fprintf('  Compartments       : %s\n', strjoin(p.comp_label, ' | '));
fprintf('  Membranes          : %s\n', strjoin(p.membrane_types, ' - '));
if p.N_units > 1
    fprintf('  Repeating units    : N = %d  (%d membrane copies, compartments lumped by type)\n', ...
        p.N_units, p.nMemCopies);
else
    fprintf('  Repeating units    : N = 1  (%d membranes)\n', p.nMemCopies);
end
fprintf('  Species tracked    : %d (%s)\n', p.nSp, strjoin(p.sp_label, ', '));
inFeed = p.metalIdx(p.conc0(d1, p.metalIdx) > 0);
if isempty(inFeed)
    fprintf('  Metals in feed     : (none)\n');
else
    fprintf('  Metals in feed     : %s\n', strjoin(arrayfun(@(j) ...
        sprintf('%s %.4g M', p.sp_label{j}, p.conc0(d1,j)), inFeed, 'UniformOutput', false), ', '));
end
if p.feedBalanceAdded_M > 0
    fprintf('  Feed balanced with : %.4f M %s (auto)\n', p.feedBalanceAdded_M, ...
        p.sp_label{strcmp(p.sp_id, p.FEED_BALANCE_ANION)});
end
if p.hasFlow
    fprintf('  Flow per chamber   : %s L/min\n', strjoin(arrayfun(@(v) sprintf('%.3g',v), ...
        p.flow_Lpm', 'UniformOutput', false), ' | '));
    fprintf('  Velocity           : %s cm/s\n', strjoin(arrayfun(@(v) sprintf('%.3g',v), ...
        p.velocity_cm_s', 'UniformOutput', false), ' | '));
    fprintf('  Boundary layer     : %s um  (Leveque; %s = fixed DELTA_CM)\n', ...
        strjoin(arrayfun(@(v) sprintf('%.1f',v*1e4), p.delta_cm', 'UniformOutput', false), ' | '), ...
        'stagnant chambers');
    fprintf('  Cell hold-up       : %.1f mL; tank turnover %s s\n', p.cellVolume_L*1000, ...
        strjoin(arrayfun(@(v) sprintf('%.0f',v), p.turnoverTime_s', 'UniformOutput', false), ' | '));
    fprintf('  Single-pass conv.  : %s %%  (well-mixed tank needs to stay small)\n', ...
        strjoin(arrayfun(@(v) sprintf('%.2f',v*100), p.singlePassConv', 'UniformOutput', false), ' | '));
    bad = find(p.singlePassConv > 0.05 & p.flow_Lpm > 0);
    for bb = bad'
        warning('EDM:singlePass', ...
            ['%s converts %.1f%% of its ions in ONE pass at %.3g L/min. The ', ...
             'model lumps it into a single well-mixed tank, which assumes this ', ...
             'is small -- a real channel would develop an axial concentration ', ...
             'profile and deplete at the outlet. Raise the pump rate or treat ', ...
             'this compartment result as indicative only.'], ...
            p.comp_label{bb}, p.singlePassConv(bb)*100, p.flow_Lpm(bb));
    end
    if any(p.Re > 2000)
        warning('EDM:turbulent', ...
            ['Re = %.0f exceeds ~2000 in at least one chamber. The Leveque ', ...
             'correlation used for delta assumes LAMINAR developing flow and ', ...
             'will under-predict mass transfer here.'], max(p.Re));
    end
else
    fprintf('  Flow per chamber   : none set -- fixed DELTA_CM = %.0f um everywhere\n', p.DELTA_CM*1e4);
end
if p.hasSolids
    fprintf('  Suspended solids   : %s mg/L\n', strjoin(arrayfun(@(v) sprintf('%.4g',v), ...
        p.tss0_mgL', 'UniformOutput', false), ' | '));
    fprintf('  Turbidity (derived): %s NTU  (@ %.2f mg/L per NTU)\n', strjoin(arrayfun(@(v) ...
        sprintf('%.4g', v/p.MGL_PER_NTU), p.tss0_mgL', 'UniformOutput', false), ' | '), p.MGL_PER_NTU);
else
    fprintf('  Suspended solids   : none (clean-water run)\n');
end
fprintf('══════════════════════════════════════════════\n\n');

end


%% ═══════════════════════════════════════════════════════════════════════
%%  ONE EXPLICIT-EULER STEP (mirrors WebUI stepOnce() exactly)
%% ═══════════════════════════════════════════════════════════════════════
function [st, overLimitFlags, diagn] = stepOnce_EDM(st, p, I, dt)
% st holds the full stack state. Solids are tracked BY COMPOSITION, not as
% lumped mass: dissolution needs to know which mineral is sitting where.
%   .conc       nComp x nSp     (mol/L)
%   .tssInert   nComp x 1       (mg/L)    fed solids; never dissolve
%   .solidBulk  nComp x nScale  (mg)      precipitate still in suspension
%   .cakeInert  1 x nMem        (mg/cm^2) deposited inert particulate
%   .cakeSolid  nMem x nScale   (mg/cm^2) deposited precipitate particles
%   .scaleMem   nMem x nScale   (mg/cm^2) scale nucleated in place on the film
%   .indClock   nComp x nScale  (s)       time spent above the nucleation SI
%   .precipMol  1 x nSp         (mol)     NET ions locked in solids, STACK-WIDE
%                                          (global, not per-compartment: solid
%                                          migrates between compartments via the
%                                          membranes, so only the total is
%                                          physically meaningful)
%   .formed_mg  nComp x nScale  (mg)      gross mass ever precipitated
%   .dissolved_mg nComp x nScale (mg)     gross mass ever redissolved
conc = st.conc;
molBefore = conc .* st.vol_L;   % nComp x nSp, mol
overLimitFlags = false(1, p.nMem);
% Per-membrane diagnostics, so the run can report what each membrane
% actually passed rather than only what the stack did in aggregate.
diagn.I_lim    = zeros(1, p.nMem);   % A, diffusion-limited counter-ion current
diagn.I_ionic  = zeros(1, p.nMem);   % A, counter-ion current actually carried
diagn.I_split  = zeros(1, p.nMem);   % A, current diverted to water splitting
% Charge accounting, coulombs carried across each membrane this step, split
% by mechanism. Every membrane must total I*dt*mult -- that is Faraday's law
% and it is the strongest single check on the transport stage.
diagn.q_counter = zeros(1, p.nMem);   % C, carried by counter-ions
diagn.q_coion   = zeros(1, p.nMem);   % C, carried by leaking co-ions
diagn.q_split   = zeros(1, p.nMem);   % C, carried by split water
diagn.molMem    = zeros(p.nMem, p.nSp);  % mol of each species moved across

% Effective boundary layer per membrane, thickened by cake AND scale
% (eq. 11). Applied to I_lim below: a fouled membrane goes over-limiting
% sooner, which is the mechanism by which fouling and scaling feed back
% into the electrochemistry rather than just adding a resistor.
[cakeTot, scaleTot] = surfaceLoads_EDM(st);
% Each membrane inherits the boundary layer of the compartment feeding it
% (eq. 18), then that layer is thickened by its own cake and scale (eq. 11).
if p.DERIVE_FOULING_RESISTANCE
    % The deposit's own thickness adds directly to the diffusion path, so
    % this too follows from loading, density and porosity (eq. 27) rather
    % than from a fitted multiplier.
    tCake  = cakeTot  / (p.PARTICLE_DENSITY_G_CM3 * (1 - p.CAKE_POROSITY )) / 1000;
    tScale = scaleTot / (p.PARTICLE_DENSITY_G_CM3 * (1 - p.SCALE_POROSITY)) / 1000;
    deltaEff = p.delta_cm(p.mem_source)' + tCake + tScale;
else
    deltaEff = p.delta_cm(p.mem_source)' .* (1 + p.CAKE_DELTA_FACTOR * (cakeTot + scaleTot));
end

for k = 1:p.nMem
    type  = p.membrane_types{k};
    left  = p.membranes(k,1); right = p.membranes(k,2);
    code  = p.membranes(k,3);
    mult  = p.membranes(k,4);          % copies of this membrane in the stack

    if code == 2
        %── BIPOLAR MEMBRANE ────────────────────────────────────────────
        % Carries NO salt. Water dissociates at the internal junction and
        % the products migrate apart: H+ toward the cathode (left), OH-
        % toward the anode (right). The whole current goes this way, so
        % there is no limiting current to exceed and no over-limiting flag.
        dMolBPM = mult * I/p.F * dt;
        molBefore(left,  p.iH)  = molBefore(left,  p.iH)  + dMolBPM;
        molBefore(right, p.iOH) = molBefore(right, p.iOH) + dMolBPM;
        overLimitFlags(k) = false;
        diagn.I_lim(k)   = Inf;      % not salt-limited: it splits water by design
        diagn.I_ionic(k) = 0;        % carries no salt at all
        diagn.I_split(k) = I;
        diagn.q_split(k) = mult * I * dt;
        diagn.molMem(k,p.iH)  = dMolBPM;
        diagn.molMem(k,p.iOH) = dMolBPM;
        st.R_bl(k) = 0;          % a BPM has no salt film to deplete
        continue
    end

    if code == 1
        source = right; destIdx = left;    % cations: right -> left (toward cathode)
        wantSign = 1;
    else
        source = left; destIdx = right;     % anions: left -> right (toward anode)
        wantSign = -1;
    end

    carriers = find(sign(p.sp_z) == wantSign);
    srcConc  = max(conc(source, carriers), 0);

    % Monovalent selectivity (eq. 22): a mono-selective membrane's surface
    % layer repels multivalent counter-ions electrostatically, so their
    % share of the current is divided by MONO_SELECTIVITY^(|z|-1). Standard
    % grades have sel = 1 for every ion, leaving eq. 1 untouched.
    if p.mem_isMono(k)
        sel = p.MONO_SELECTIVITY .^ -(abs(p.sp_z(carriers)) - 1);
    else
        sel = ones(1, numel(carriers));
    end
    wgtSrc = sel .* abs(p.sp_z(carriers)) .* p.sp_lambda(carriers) .* srcConc;
    denom  = sum(wgtSrc);

    % Limiting current (eq. 8): max current this membrane's own ion
    % inventory can sustain via diffusion to the surface. denom and I_lim
    % are zero under exactly the same condition (all carrier concentrations
    % zero), so I_ionic naturally goes to 0 and the full current routes
    % through the water-splitting branch below -- no separate guard needed.
    % Limiting current (eq. 8). The selectivity weights belong here too: an
    % ion the membrane rejects cannot carry current no matter how fast it
    % diffuses to the surface, so counting it would overstate what this
    % membrane can pass. Omitting it lets a mono-selective membrane appear
    % able to sustain a current its permitted ions cannot actually supply --
    % the deficit then gets silently absorbed by the availability cap below,
    % which breaks charge conservation in the source compartment.
    I_lim = p.area_cm2 * sum(sel .* abs(p.sp_z(carriers)) .* p.F .* p.sp_D(carriers) .* srcConc / 1000 / deltaEff(k));

    % Only the counter-ion share of the current has to cross as counter-ions,
    % so it is that share which the limiting current constrains (eq. 23).
    alpha     = p.permselectivity(k);
    I_counter = alpha * I;
    overLimit = I_counter > I_lim;
    overLimitFlags(k) = overLimit;
    I_ionic = min(I_counter, I_lim);
    diagn.I_lim(k)   = I_lim;
    diagn.I_ionic(k) = I_ionic;
    diagn.I_split(k) = max(I_counter - I_lim, 0);

    % Concentration polarisation (eq. 26): the depleted film's excess
    % resistance, from its log-mean conductivity between bulk and surface.
    if p.ENABLE_CONC_POLARISATION && I_lim > 1e-30
        ratio = min(I_counter/I_lim, 1 - p.CS_FLOOR);
        x     = max(1 - ratio, p.CS_FLOOR);        % c_surface / c_bulk
        if ratio < 1e-9
            fac = 1;                                % undepleted; no excess
        else
            fac = (1 - x) / log(1/x);               % log-mean of the film
        end
        kap = max(conductivity_EDM(conc(source,:), 0, p), p.KAPPA_WATER_FLOOR);
        st.R_bl(k) = deltaEff(k)/(p.area_cm2*kap) * (1/max(fac,p.CS_FLOOR) - 1);
    else
        st.R_bl(k) = 0;
    end

    if denom > 1e-30 && I_ionic > 1e-30
        share       = wgtSrc / denom;
        % Each of the `mult` copies passes the same current I, so a lumped
        % tank bounded by mult copies sees mult times the single-cell flux.
        fluxMolPerS = mult * share * I_ionic ./ (abs(p.sp_z(carriers)) * p.F);
        dMol        = fluxMolPerS * dt;
        % Cap to what's actually available in the source THIS step (avoids
        % an Euler-overshoot mass-conservation leak when a species fully
        % depletes mid-step).
        dMol = min(dMol, max(molBefore(source, carriers), 0));
        molBefore(source, carriers) = molBefore(source, carriers) - dMol;
        molBefore(destIdx, carriers) = molBefore(destIdx, carriers) + dMol;
        diagn.q_counter(k) = sum(abs(p.sp_z(carriers)) .* dMol) * p.F;
        diagn.molMem(k,carriers) = diagn.molMem(k,carriers) + dMol;
    end

    %── Co-ion leakage (eq. 23) ─────────────────────────────────────────
    % The remaining (1-alpha)*I crosses as co-ions going the OTHER way:
    % anions through a CEM, cations through an AEM. They are drawn from the
    % opposite compartment and share that current by their own transport
    % numbers. This moves salt backwards -- undoing separation -- which is
    % exactly what imperfect permselectivity costs you.
    if alpha < 1
        I_co  = (1 - alpha) * I;
        srcCo = destIdx; dstCo = source;          % opposite direction
        carriersCo = find(sign(p.sp_z) == -wantSign);
        srcConcCo  = max(conc(srcCo, carriersCo), 0);
        wgtCo      = abs(p.sp_z(carriersCo)) .* p.sp_lambda(carriersCo) .* srcConcCo;
        denomCo    = sum(wgtCo);
        if denomCo > 1e-30 && I_co > 1e-30
            shareCo = wgtCo / denomCo;
            fluxCo  = mult * shareCo * I_co ./ (abs(p.sp_z(carriersCo)) * p.F);
            dMolCo  = min(fluxCo * dt, max(molBefore(srcCo, carriersCo), 0));
            molBefore(srcCo, carriersCo) = molBefore(srcCo, carriersCo) - dMolCo;
            molBefore(dstCo, carriersCo) = molBefore(dstCo, carriersCo) + dMolCo;
            diagn.q_coion(k) = sum(abs(p.sp_z(carriersCo)) .* dMolCo) * p.F;
            diagn.molMem(k,carriersCo) = diagn.molMem(k,carriersCo) + dMolCo;
        end
    end

    if overLimit
        % Excess current beyond I_lim is carried by water splitting right at
        % this membrane. CEM: H+ (a cation) crosses through to dest, OH-
        % stays in source. AEM: OH- crosses to dest, H+ stays in source.
        dMolWS = mult * (I_counter - I_lim) / p.F * dt;
        diagn.q_split(k) = dMolWS * p.F;
        % Branch on the TYPE CODE, not the type string: 'Mono-CEM' is a
        % cation-exchange membrane and must split water the same way a plain
        % CEM does. Testing strcmp(type,'CEM') silently sent Mono-CEM down
        % the anion branch, putting H+ and OH- on the wrong sides and
        % breaking charge conservation in the source compartment.
        if code == 1
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

% Transport stage done -- write back before the chemistry stage, which works
% in concentrations. Order matters: precipitation must see the post-electrode
% pH, since that is what drives hydroxide scaling.
%── Water transport across the membranes (eq. 19-20) ─────────────────────
% Applied to volumes BEFORE converting moles back to concentrations, so the
% concentrating/diluting effect of the water movement is captured this step
% rather than lagging one behind.
if p.ENABLE_WATER_TRANSPORT
    volBefore = st.vol_L;
    st.vol_L  = waterTransport_EDM(st.vol_L, conc, p, I, dt);
    % Water moves; the particles it leaves behind do not. tssInert is a
    % CONCENTRATION, so it has to be rescaled or suspended-solid mass is
    % silently created/destroyed whenever a compartment gains or loses
    % volume. (solidBulk is stored as mass, so it needs no rescaling.)
    st.tssInert = st.tssInert .* volBefore ./ st.vol_L;
end

conc = max(molBefore ./ st.vol_L, 0);

conc = max(conc, 0);

%── Water self-ionization (eq. 17) ───────────────────────────────────────
% Before the solids chemistry, because it is the equilibrated OH- that sets
% hydroxide saturation -- running precipitation against a non-equilibrium
% OH- inventory would scale on hydroxide that cannot physically be there.
if p.ENFORCE_WATER_EQUILIBRIUM
    conc = waterEquilibrium_EDM(conc, p);
end
st.conc = conc;

%── Solids chemistry: nucleation, growth, dissolution (eq. 12-16) ────────
if p.ENABLE_SCALING
    st = solidsChemistry_EDM(st, p, dt);
    % Hydroxide precipitation consumes OH- and dissolution releases it, so
    % the compartment leaves the chemistry stage off the water equilibrium
    % again. Re-impose it: water re-dissociates to partly buffer what the
    % precipitation just took, which is real chemistry, not bookkeeping.
    % (Operator splitting -- one pass per timestep, not iterated to joint
    % convergence. With dt = 1 min against a ~17 min precipitation time
    % constant the splitting error is small, but shorten dt if you raise
    % K_PRECIP_PER_S toward equilibrium-limited behaviour.)
    if p.ENFORCE_WATER_EQUILIBRIUM
        st.conc = waterEquilibrium_EDM(st.conc, p);
    end
end

%── Suspended-solids deposition onto membranes (eq. 10) ──────────────────
% Charged particles drift toward one electrode and pile up on the first
% membrane in their path. With the default negative zeta potential they
% move toward the anode (rightward), so compartment i fouls membrane i --
% the membrane on its anode side. A positive zeta flips this to membrane
% i-1. Particles never CROSS a membrane; they only leave the bulk by
% depositing, so bulk solids decay as cake grows.
%
% Deposition is composition-blind: whatever is suspended deposits in
% proportion to what it is made of. Tracking WHICH solid ends up on the
% membrane is what lets it dissolve again later.
tssTot = totalTSS_EDM(st, p);
for i = 1:p.nComp
    if tssTot(i) <= 0; continue; end
    if p.TSS_CHARGE_SIGN < 0
        memList = find(p.membranes(:,1) == i)';   % membranes on this tank's anode side
    else
        memList = find(p.membranes(:,2) == i)';   % on its cathode side
    end
    if isempty(memList); continue; end            % electrode-facing: no membrane
    mem = memList(1);
    % Area available to deposit on: every copy of every qualifying row.
    depArea = sum(p.membranes(memList,4)) * p.area_cm2;
    % flux [mg/cm^2/s] = v[cm/s] * C[mg/cm^3];  C = tss/1000
    dMass = p.dep_vel_cm_s(i) * (tssTot(i)/1000) * dt * depArea;      % mg
    % Cap at the solids actually present this step (same Euler-overshoot
    % guard as the ionic transport above).
    dMass = min(dMass, max(tssTot(i)*st.vol_L(i), 0));
    if dMass <= 0; continue; end
    frac = dMass / (tssTot(i)*st.vol_L(i));      % same fraction of every phase

    % Spread the deposited mass over the qualifying rows in proportion to
    % their area, then express it per unit area on each row.
    wgt = p.membranes(memList,4)' * p.area_cm2 / depArea;

    dInert = st.tssInert(i)*st.vol_L(i) * frac;
    st.tssInert(i) = max(st.tssInert(i) - dInert/st.vol_L(i), 0);
    st.cakeInert(memList) = st.cakeInert(memList) + ...
        (dInert*wgt) ./ (p.membranes(memList,4)'*p.area_cm2);

    dEach = st.solidBulk(i,:) * frac;             % mg of each solid
    st.solidBulk(i,:) = max(st.solidBulk(i,:) - dEach, 0);
    for mm = 1:numel(memList)
        st.cakeSolid(memList(mm),:) = st.cakeSolid(memList(mm),:) + ...
            dEach*wgt(mm) / (p.membranes(memList(mm),4)*p.area_cm2);
    end
end

end


function tssTot = totalTSS_EDM(st, p)
% Bulk suspended solids = inert (fed) + precipitate still in suspension.
tssTot = st.tssInert + (sum(st.solidBulk, 2) ./ st.vol_L);
end


%% ═══════════════════════════════════════════════════════════════════════
%%  DERIVED QUANTITIES
%% ═══════════════════════════════════════════════════════════════════════
%% ═══════════════════════════════════════════════════════════════════════
%%  PRECIPITATION / SCALING (eq. 12-13)
%% ═══════════════════════════════════════════════════════════════════════
function gamma = activityCoeffs_EDM(concRow, p)
% Davies equation: log10(gamma_i) = -A z_i^2 (sqrt(I)/(1+sqrt(I)) - 0.3 I)
% Valid to about I = 0.5 M; beyond that it is an extrapolation (flagged).
if ~p.USE_DAVIES
    gamma = ones(1, p.nSp);
    return
end
Ion = 0.5 * sum(p.sp_z.^2 .* max(concRow,0));       % ionic strength, mol/L
A   = 0.509;                                          % water, 25 degC
sI  = sqrt(Ion);
gamma = 10.^(-A * p.sp_z.^2 * (sI/(1+sI) - 0.3*Ion));
end

function SI = saturationIndex_EDM(conc, p)
% SI = log10(IAP/Ksp): 0 = saturated, >0 supersaturated (scaling risk),
% <0 undersaturated. Reported for every compartment/solid pair whether or
% not precipitation is enabled, so you can see the risk even in a run with
% ENABLE_SCALING = false.
SI = zeros(p.nComp, p.nScale);
for i = 1:p.nComp
    gamma = activityCoeffs_EDM(conc(i,:), p);
    for k = 1:p.nScale
        ic = p.scale_iCat(k); ia = p.scale_iAn(k);
        cC = max(conc(i,ic),0); cA = max(conc(i,ia),0);
        if cC <= 0 || cA <= 0
            SI(i,k) = -Inf;
            continue
        end
        IAP = (gamma(ic)*cC)^p.scale_nuCat(k) * (gamma(ia)*cA)^p.scale_nuAn(k);
        SI(i,k) = log10(IAP / p.scale_Ksp(k));
    end
end
end

function vol = waterTransport_EDM(vol, conc, p, I, dt)
% Electro-osmotic drag (eq. 19) and osmosis (eq. 20) across every membrane.
% Water is MOVED between compartments, never created or destroyed, so the
% stack total is conserved exactly -- runEDM checks this. A real rig also
% loses water to evaporation and genuine leaks; neither is modelled, so if
% your measured TOTAL falls, that difference is the real leak.
dVol = zeros(p.nComp,1);   % L, net change per compartment this step

for k = 1:p.nMem
    left = p.membranes(k,1); right = p.membranes(k,2);
    mult = p.membranes(k,4);
    switch p.membranes(k,3)
        case 1   % cation-conducting: drag carries water right -> left
            src = right; dst = left; nDrag = p.N_DRAG_CEM;
        case 0   % anion-conducting: drag carries water left -> right
            src = left; dst = right; nDrag = p.N_DRAG_AEM;
        otherwise % bipolar: H+ and OH- leave in opposite directions, so the
                  % two drags oppose and largely cancel
            src = left; dst = right; nDrag = p.N_DRAG_BPM;
    end

    % (19) Electro-osmotic drag: proportional to charge passed.
    molCharge = mult * I/p.F * dt;                           % mol of charge (all copies)
    dV_drag   = nDrag * molCharge * p.V_WATER_MOLAR_CM3 / 1000;   % L

    % (20) Osmosis: van 't Hoff, water toward the more concentrated side.
    % Osmolarity counted over all ions -- a dissociated salt contributes one
    % osmotic unit per ion, which is why total ion molarity is the right
    % measure here rather than salt molarity.
    osmL = sum(max(conc(left,:), 0));
    osmR = sum(max(conc(right,:),0));
    dPi  = p.Rgas/100 * p.T * (osmR - osmL);   % bar (R/100 = L*bar/mol/K)
    % Positive dPi means the right side is more concentrated, so water is
    % drawn left -> right.
    dV_osm = p.LP_WATER * mult * p.area_cm2 * dPi * dt / 1000;   % L

    dVol(src) = dVol(src) - dV_drag;
    dVol(dst) = dVol(dst) + dV_drag;
    dVol(left)  = dVol(left)  - dV_osm;
    dVol(right) = dVol(right) + dV_osm;
end

% Floor: keep a compartment from being driven to zero (or negative) volume,
% which would make every concentration infinite.
%
% Clamping each compartment independently would CREATE water: dVol sums to
% zero by construction, but max(vol,volMin) does not preserve that sum, so a
% clamped run silently gains volume (135% over a 96 h run before this was
% fixed). Instead the whole step's water transfer is scaled back by the
% single factor that just keeps the tightest compartment on its floor. That
% preserves sum(dVol) = 0 exactly, so volume is still conserved, and it is
% also the more physical statement: a compartment that has run dry cannot
% keep donating water.
volMin = p.MIN_VOL_FRACTION * p.vol_L;
f = 1;
for i = 1:p.nComp
    if dVol(i) < 0 && (vol(i) + dVol(i)) < volMin(i)
        f = min(f, max(volMin(i) - vol(i), 0) / dVol(i));
    end
end
vol = vol + f*dVol;
vol = max(vol, volMin);          % belt and braces against round-off
end

function R = recordResistance_(R, br, step)
R.Rmem_ohm(step,:)       = br.Rmem;
R.Rmem_clean_ohm(step,:) = br.Rmem_clean;
R.Rmem_cake_ohm(step,:)  = br.Rmem_cake;
R.Rmem_scale_ohm(step,:) = br.Rmem_scale;
R.Rcomp_ohm(step,:)      = br.Rcomp;
R.Rtotal_ohm(step)       = br.Rtotal;
end

function r = layerAreal_(w, eps, rho, kappa)
% Areal resistance (ohm cm^2) of a deposit of loading w (mg/cm^2), porosity
% eps and particle density rho, sitting in solution of conductivity kappa.
% Bruggeman for the porous layer; the deposit conducts only through the
% electrolyte held in its pores.
t = w / (rho * (1 - eps)) / 1000;          % cm of layer
r = t ./ max(kappa * eps^1.5, 1e-12);
end

function [cakeTot, scaleTot] = surfaceLoads_EDM(st)
% Total surface loading per membrane (mg/cm^2), summed over composition.
cakeTot  = st.cakeInert + sum(st.cakeSolid, 2)';
scaleTot = sum(st.scaleMem, 2)';
end

function st = solidsChemistry_EDM(st, p, dt)
% Precipitation (with a nucleation barrier) and dissolution, per compartment
% per solid. Both directions solve the SAME equilibrium condition -- the
% extent x that brings IAP to Ksp -- and differ only in sign and in what
% limits them: precipitation by the dissolved ions available, dissolution by
% the solid actually in contact with that compartment.
relaxP = 1 - exp(-p.K_PRECIP_PER_S * dt);
relaxD = 1 - exp(-p.K_DISSOL_PER_S * dt);

for i = 1:p.nComp
    % Membranes bounding this compartment: their deposits are in contact
    % with this solution and can dissolve back into it.
    mems = find(p.membranes(:,1) == i | p.membranes(:,2) == i)';

    for k = 1:p.nScale
        ic = p.scale_iCat(k); ia = p.scale_iAn(k);
        a  = p.scale_nuCat(k); b = p.scale_nuAn(k);
        cC = max(st.conc(i,ic), 0); cA = max(st.conc(i,ia), 0);

        gamma = activityCoeffs_EDM(st.conc(i,:), p);
        gC = gamma(ic); gA = gamma(ia);
        Ksp = p.scale_Ksp(k);

        % Solid of this type in contact with this compartment's solution
        rowArea = p.membranes(mems,4)' * p.area_cm2;      % cm^2 per row (all copies)
        seed_mg = st.solidBulk(i,k) + ...
                  sum(st.cakeSolid(mems,k)'.*rowArea) + sum(st.scaleMem(mems,k)'.*rowArea);

        IAP = (gC*cC)^a * (gA*cA)^b;
        if cC <= 0 || cA <= 0; IAP = 0; end

        if IAP > Ksp
            %── SUPERSATURATED: grow, subject to the nucleation barrier ──
            SI = log10(IAP/Ksp);
            if seed_mg <= 1e-12
                % No seed crystal: primary nucleation must be cleared first
                % (eq. 15). Below the threshold driving force nothing forms
                % and the induction clock unwinds.
                if SI < p.SI_NUCLEATION
                    st.indClock(i,k) = max(st.indClock(i,k) - dt, 0);
                    continue
                end
                st.indClock(i,k) = st.indClock(i,k) + dt;
                t_ind = p.TAU_NUC_REF_S / SI^2;      % steep fall with driving force
                if st.indClock(i,k) < t_ind
                    continue                          % still in the induction period
                end
            end

            % Extent that restores equilibrium, by bisection: IAP is
            % monotonically decreasing in x, so no derivative is needed and
            % it is stable at any timestep.
            hi = min(cC/a, cA/b);
            f  = @(x) (gC*max(cC-a*x,0))^a * (gA*max(cA-b*x,0))^b - Ksp;
            x  = bisectExtent_EDM(f, 0, hi);
            x  = x * relaxP;                          % kinetic damping (eq. 13)
            if x <= 0; continue; end

            st.conc(i,ic) = max(st.conc(i,ic) - a*x, 0);
            st.conc(i,ia) = max(st.conc(i,ia) - b*x, 0);
            st.precipMol(ic) = st.precipMol(ic) + a*x*st.vol_L(i);
            st.precipMol(ia) = st.precipMol(ia) + b*x*st.vol_L(i);
            % Per-compartment copy, for auditing only. It can go negative in
            % a compartment when solid formed elsewhere dissolves into it --
            % solid migrates via the membranes it deposits on -- which is
            % exactly what makes a per-compartment balance close.
            st.precipMolComp(i,ic) = st.precipMolComp(i,ic) + a*x*st.vol_L(i);
            st.precipMolComp(i,ia) = st.precipMolComp(i,ia) + b*x*st.vol_L(i);

            newMass = x * st.vol_L(i) * p.scale_MW(k) * 1000;   % mol -> mg
            st.formed_mg(i,k) = st.formed_mg(i,k) + newMass;
            % Split between in-situ growth on the bounding membranes and the
            % bulk, where it becomes suspended solid (eq. 14).
            if isempty(mems)
                st.solidBulk(i,k) = st.solidBulk(i,k) + newMass;
            else
                surfMass = newMass * p.SCALE_SURFACE_FRACTION;
                % Share by area across the bounding rows, then convert to a
                % per-unit-area loading on each.
                shareA = rowArea / sum(rowArea);
                st.scaleMem(mems,k) = st.scaleMem(mems,k) + ...
                                      (surfMass*shareA ./ rowArea)';
                st.solidBulk(i,k)   = st.solidBulk(i,k) + (newMass - surfMass);
            end

        elseif p.ENABLE_DISSOLUTION && seed_mg > 1e-12
            %── UNDERSATURATED with solid present: redissolve (eq. 16) ──
            % Same equilibrium condition, solved in the opposite direction:
            % x is now mol/L of solid returning to solution.
            st.indClock(i,k) = max(st.indClock(i,k) - dt, 0);
            maxX = seed_mg / (p.scale_MW(k)*1000) / st.vol_L(i);  % mg -> mol/L
            g  = @(x) Ksp - (gC*(cC+a*x))^a * (gA*(cA+b*x))^b;    % increasing in x
            x  = bisectExtent_EDM(g, 0, maxX);
            x  = x * relaxD;
            if x <= 0; continue; end

            st.conc(i,ic) = st.conc(i,ic) + a*x;
            st.conc(i,ia) = st.conc(i,ia) + b*x;
            st.precipMol(ic) = st.precipMol(ic) - a*x*st.vol_L(i);
            st.precipMol(ia) = st.precipMol(ia) - b*x*st.vol_L(i);
            st.precipMolComp(i,ic) = st.precipMolComp(i,ic) - a*x*st.vol_L(i);
            st.precipMolComp(i,ia) = st.precipMolComp(i,ia) - b*x*st.vol_L(i);

            % Consume solid loosest-first: suspended crystals have the most
            % exposed surface, dense in-situ scale the least.
            need = x * st.vol_L(i) * p.scale_MW(k) * 1000;        % mg
            st.dissolved_mg(i,k) = st.dissolved_mg(i,k) + need;
            take = min(need, st.solidBulk(i,k));
            st.solidBulk(i,k) = st.solidBulk(i,k) - take;  need = need - take;
            for src = 1:2
                if need <= 0; break; end
                for mi = 1:numel(mems)
                    if need <= 0; break; end
                    m = mems(mi); aM = rowArea(mi);
                    if src == 1
                        take = min(need, st.cakeSolid(m,k)*aM);
                        st.cakeSolid(m,k) = st.cakeSolid(m,k) - take/aM;
                    else
                        take = min(need, st.scaleMem(m,k)*aM);
                        st.scaleMem(m,k)  = st.scaleMem(m,k)  - take/aM;
                    end
                    need = need - take;
                end
            end
        else
            st.indClock(i,k) = max(st.indClock(i,k) - dt, 0);
        end
    end
end
end

function x = bisectExtent_EDM(f, lo, hi)
% Bisection for a monotonic f on [lo,hi] with f(lo) > 0 >= f(hi) expected.
% If the whole interval stays positive the extent is capped at hi (the
% reaction runs out of reagent before reaching equilibrium).
if hi <= lo; x = 0; return; end
if f(hi) > 0; x = hi; return; end
for it = 1:60                     % ~1e-18 relative precision
    mid = 0.5*(lo+hi);
    if f(mid) > 0; lo = mid; else; hi = mid; end
end
x = 0.5*(lo+hi);
end

function kappa = conductivity_EDM(concRow, tssRow, p)
% kappa [S/cm] = 0.001 * sum( lambda_i * |z_i| * C_i[mol/L] )   (eq. 5)
lam = p.sp_lambda;
if p.KAPPA_IONIC_CORRECTION
    % Ion-ion interaction slows every ion as the solution concentrates
    % (eq. 5b). Applied to all ions alike, so transport numbers are
    % unaffected -- only conductivity, and hence voltage, changes.
    Ion = 0.5 * sum(p.sp_z.^2 .* max(concRow,0));
    sI  = sqrt(Ion);
    lam = lam * max(1 - p.KAPPA_B * sI/(1+sI), 0.05);
end
kappa = 0.001 * sum(lam .* abs(p.sp_z) .* max(concRow,0));
% Suspended solids are non-conducting inclusions: Maxwell's relation for a
% dilute suspension of insulating spheres, kappa_eff = kappa*(1-phi)/(1+phi/2).
% phi is tiny at realistic TSS (100 mg/L of rho=2.5 solids is phi=4e-5), so
% this is a completeness term, not a driver -- fouling dominates by orders
% of magnitude.
if tssRow > 0
    phi   = min(tssRow / (p.PARTICLE_DENSITY_G_CM3 * 1e6), 0.6);   % mg/L -> volume fraction
    kappa = kappa * (1-phi) / (1 + phi/2);
end
end

function [V, Rtot, V_foul, V_scale, br] = stackVoltage_EDM(st, overLimitCount, I, p)
% br breaks the resistance out element by element, so the stack can be
% audited: which membrane is costing what, and how much of that is the
% membrane itself versus what the reactions have deposited on it.
tssTot = totalTSS_EDM(st, p);
br.Rcomp = zeros(1, p.nComp);        % solution resistance of each compartment
br.kappa = zeros(1, p.nComp);        % S/cm, reported so its decline can be tracked
for i = 1:p.nComp
    kappa = max(conductivity_EDM(st.conc(i,:), tssTot(i), p), p.KAPPA_WATER_FLOOR);
    br.kappa(i) = kappa;
    br.Rcomp(i) = p.GAP_CM/(kappa*p.area_cm2);
end
% Membrane term, per row: clean areal resistance + particulate cake +
% crystalline scale (eq. 11, 14). Scale is the more resistive of the two per
% unit mass. Every copy of a row is in series, hence mem_mult.
[cakeTot, scaleTot] = surfaceLoads_EDM(st);
br.Rmem_clean = p.mem_mult .* p.areal_R / p.area_cm2;
if p.DERIVE_FOULING_RESISTANCE
    % Areal resistance of each deposit from its own thickness and porosity
    % (eq. 27), using the conductivity of the solution the layer sits in.
    kapAt = zeros(1,p.nMem);
    for k = 1:p.nMem
        kapAt(k) = max(br.kappa(p.mem_source(k)), p.KAPPA_WATER_FLOOR);
    end
    rCake  = layerAreal_(cakeTot,  p.CAKE_POROSITY,  p.PARTICLE_DENSITY_G_CM3, kapAt);
    rScale = layerAreal_(scaleTot, p.SCALE_POROSITY, p.PARTICLE_DENSITY_G_CM3, kapAt);
else
    rCake  = p.R_CAKE_AREAL  * cakeTot;
    rScale = p.R_SCALE_AREAL * scaleTot;
end
br.Rmem_cake  = p.mem_mult .* rCake  / p.area_cm2;
br.Rmem_scale = p.mem_mult .* rScale / p.area_cm2;
% Concentration polarisation sits electrically in series with the membrane
% it belongs to, so it is reported as part of that membrane's resistance.
if isfield(st,'R_bl'); br.Rmem_polar = p.mem_mult .* st.R_bl; else; br.Rmem_polar = zeros(1,p.nMem); end
br.Rmem       = br.Rmem_clean + br.Rmem_cake + br.Rmem_scale + br.Rmem_polar;
R_clean = sum(br.Rmem_clean);
R_cake  = sum(br.Rmem_cake);
R_scale = sum(br.Rmem_scale);
R_polar = sum(br.Rmem_polar);
Rtot    = sum(br.Rcomp) + R_clean + R_cake + R_scale + R_polar;
br.Rpolarisation = R_polar;
br.Rsolution = sum(br.Rcomp);
br.Rtotal    = Rtot;
V_foul  = I * R_cake;                        % volts attributable to fouling
V_scale = I * R_scale;                       % volts attributable to scaling
% Bipolar membranes each add their water-dissociation junction potential,
% which is a standing cost of the stack, not a consequence of polarisation.
V_watersplit = overLimitCount * p.V_WATERSPLIT_OVERPOTENTIAL + p.nBPMCopies * p.V_BPM;
V = min(I*Rtot + p.V_OVERPOTENTIAL + V_watersplit, p.V_STACK_CEILING);
end

function conc = waterEquilibrium_EDM(conc, p)
% Re-speciate H+/OH- onto C_H*C_OH = Kw in every compartment, holding the net
% (C_H - C_OH) fixed (eq. 17).
%
% Net is the physically meaningful quantity: membrane transport and the
% electrode reactions determine how much EXCESS acid or base a compartment
% holds, not the individual H+ and OH- inventories. Holding net fixed also
% conserves charge exactly, since neutralization removes one positive and
% one negative charge together.
%
% Solved in closed form from  cH - cOH = net,  cH*cOH = Kw:
%   cH = (net + sqrt(net^2 + 4Kw))/2
% evaluated in whichever algebraic form avoids cancellation for the sign of
% net -- the naive form loses all precision for net << -sqrt(Kw).
for i = 1:p.nComp
    net  = conc(i,p.iH) - conc(i,p.iOH);
    disc = sqrt(net*net + 4*p.Kw);
    if net >= 0
        cH  = 0.5*(net + disc);
        cOH = p.Kw / max(cH, realmin);
    else
        cOH = 0.5*(-net + disc);
        cH  = p.Kw / max(cOH, realmin);
    end
    conc(i,p.iH)  = cH;
    conc(i,p.iOH) = cOH;
end
end

function [pH, pOH] = phOf_EDM(cH, cOH, p, concRow)
% pH from the equilibrated H+ concentration (eq. 6). By definition pH is
% -log10 of ACTIVITY; with p.PH_USE_ACTIVITY the Davies coefficient already
% used by the scaling block is applied here too, for consistency between
% what sets pH and what sets hydroxide saturation.
gH = 1; gOH = 1;
if p.PH_USE_ACTIVITY && nargin >= 4
    gamma = activityCoeffs_EDM(concRow, p);
    gH = gamma(p.iH); gOH = gamma(p.iOH);
end
pH  = -log10(max(gH*cH,   1e-30));
pOH = -log10(max(gOH*cOH, 1e-30));
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
R.pOH   = zeros(nSteps, p.nComp);   % reported alongside pH: in a stack that
R.cH_M  = zeros(nSteps, p.nComp);   % splits water at the membranes, the base
R.cOH_M = zeros(nSteps, p.nComp);   % side matters as much as the acid side
R.V     = zeros(nSteps,1);
R.overLimitFlags = false(nSteps, p.nMem);
R.overLimitCount = zeros(nSteps,1);
R.tss_mgL        = zeros(nSteps, p.nComp);
R.turbidity_NTU  = zeros(nSteps, p.nComp);
R.cake_mgcm2     = zeros(nSteps, p.nMem);
R.V_fouling      = zeros(nSteps,1);
R.scale_mgcm2    = zeros(nSteps, p.nMem);
R.V_scaling      = zeros(nSteps,1);
R.scaleMass_mg   = zeros(nSteps, p.nScale);   % summed over compartments
R.SI             = zeros(nSteps, p.nComp, p.nScale);
R.precipMol_hist = zeros(nSteps, p.nSp);   % stack-wide ions held as solid

st.conc         = conc;
st.vol_L        = p.vol_L;      % volume is now STATE, not a constant
st.tssInert     = p.tss0_mgL;
st.solidBulk    = zeros(p.nComp, p.nScale);
st.cakeInert    = zeros(1, p.nMem);
st.cakeSolid    = zeros(p.nMem, p.nScale);
st.scaleMem     = zeros(p.nMem, p.nScale);
st.indClock     = zeros(p.nComp, p.nScale);
st.precipMol     = zeros(1, p.nSp);
st.precipMolComp = zeros(p.nComp, p.nSp);   % audit diagnostic, see above
st.R_bl          = zeros(1, p.nMem);        % concentration-polarisation resistance
st.formed_mg    = zeros(p.nComp, p.nScale);
st.dissolved_mg = zeros(p.nComp, p.nScale);
initialSolids_mg = sum(p.tss0_mgL .* p.vol_L);
R.dissolved_mg   = zeros(nSteps, p.nScale);
R.vol_L          = zeros(nSteps, p.nComp);   % L, per compartment over time
% Resistance broken out per element (ohm), recorded every step.
R.Rmem_ohm       = zeros(nSteps, p.nMem);   % total, per membrane row
R.Rmem_clean_ohm = zeros(nSteps, p.nMem);   % the membrane itself
R.Rmem_cake_ohm  = zeros(nSteps, p.nMem);   % particulate deposited on it
R.Rmem_scale_ohm = zeros(nSteps, p.nMem);   % crystalline scale grown on it
R.Rcomp_ohm      = zeros(nSteps, p.nComp);  % solution resistance per compartment
R.Rtotal_ohm     = zeros(nSteps, 1);
% Flux decline (eq. 24): the counter-ion flux each membrane actually carries.
% At constant current the TOTAL charge flux is fixed by the power supply, so
% what declines is the share carried by salt: as a membrane depletes or fouls,
% I_lim falls, the counter-ion current is capped, and the remainder is diverted
% into water splitting. That diverted fraction is wasted energy, and tracking
% it is how you see a membrane failing before the voltage says so.
R.flux_eq_m2h    = zeros(nSteps, p.nMem);   % eq/(m^2 h) of counter-ions
R.fluxNorm       = zeros(nSteps, p.nMem);   % flux / initial flux
R.Ilim_A         = zeros(nSteps, p.nMem);   % limiting current per membrane
R.saltFraction   = zeros(nSteps, p.nMem);   % fraction of I carried by salt
% Conductivity decline (eq. 25): per compartment, and normalised to t = 0.
R.kappa_S_cm     = zeros(nSteps, p.nComp);
R.kappaNorm      = zeros(nSteps, p.nComp);
% Cumulative charge accounting per membrane (C), by mechanism, and the moles
% of each species each membrane has moved. Used by Audit_EDM to check
% Faraday's law membrane by membrane.
R.q_counter_C    = zeros(1, p.nMem);
R.q_coion_C      = zeros(1, p.nMem);
R.q_split_C      = zeros(1, p.nMem);
R.molMem_mol     = zeros(p.nMem, p.nSp);
R.vol_L(1,:)     = p.vol_L;

initialTotals = squeeze(sum(p.conc0 .* p.vol_L, 1));   % 1 x nSp, mol

% t=0 snapshot (before any stepping -- matches WebUI's initial history push)
R.t_sec(1) = 0;
R.conc(1,:,:) = conc;
if p.ENFORCE_WATER_EQUILIBRIUM
    conc = waterEquilibrium_EDM(conc, p);   % t=0 snapshot on the equilibrium too
    st.conc = conc; R.conc(1,:,:) = conc;
end
for c = 1:p.nComp
    [R.pH(1,c), R.pOH(1,c)] = phOf_EDM(conc(c,p.iH), conc(c,p.iOH), p, conc(c,:));
    R.cH_M(1,c) = conc(c,p.iH); R.cOH_M(1,c) = conc(c,p.iOH);
end
[R.V(1), ~, R.V_fouling(1), R.V_scaling(1), br] = stackVoltage_EDM(st, 0, p.I_A, p);
R = recordResistance_(R, br, 1);
R.kappa_S_cm(1,:) = br.kappa;
% t = 0 flux: evaluated on the initial state, before anything has depleted,
% so it is the reference every later flux is measured against.
[~, ~, d0] = stepOnce_EDM(st, p, p.I_A, 0);
R.Ilim_A(1,:)       = d0.I_lim;
R.flux_eq_m2h(1,:)  = d0.I_ionic / p.F * 3600 / (p.area_cm2/1e4);
R.saltFraction(1,:) = d0.I_ionic / p.I_A;
flux0 = R.flux_eq_m2h(1,:);
R.tss_mgL(1,:)       = totalTSS_EDM(st, p);
R.turbidity_NTU(1,:) = R.tss_mgL(1,:) / p.MGL_PER_NTU;
R.SI(1,:,:)          = saturationIndex_EDM(st.conc, p);

t = 0;
for step = 2:nSteps
    [st, overLimitFlags, diagn] = stepOnce_EDM(st, p, p.I_A, p.dt_sec);
    conc = st.conc;
    t = t + p.dt_sec;

    R.t_sec(step) = t;
    R.conc(step,:,:) = conc;
    for c = 1:p.nComp
        [R.pH(step,c), R.pOH(step,c)] = phOf_EDM(conc(c,p.iH), conc(c,p.iOH), p, conc(c,:));
        R.cH_M(step,c) = conc(c,p.iH); R.cOH_M(step,c) = conc(c,p.iOH);
    end
    R.overLimitFlags(step,:) = overLimitFlags;
    R.overLimitCount(step) = sum(overLimitFlags);
    [cakeTot, scaleTot]     = surfaceLoads_EDM(st);
    R.tss_mgL(step,:)       = totalTSS_EDM(st, p);
    R.turbidity_NTU(step,:) = R.tss_mgL(step,:) / p.MGL_PER_NTU;
    R.cake_mgcm2(step,:)    = cakeTot;
    R.scale_mgcm2(step,:)   = scaleTot;
    % Solid CURRENTLY present (formed minus redissolved), not gross formed.
    R.scaleMass_mg(step,:)  = sum(st.solidBulk,1) + ...
        (p.mem_mult*st.cakeSolid + p.mem_mult*st.scaleMem) * p.area_cm2;
    R.dissolved_mg(step,:)  = sum(st.dissolved_mg, 1);
    R.vol_L(step,:)         = st.vol_L;
    R.SI(step,:,:)          = saturationIndex_EDM(st.conc, p);
    R.precipMol_hist(step,:) = st.precipMol;
    [R.V(step), ~, R.V_fouling(step), R.V_scaling(step), br] = ...
        stackVoltage_EDM(st, R.overLimitCount(step), p.I_A, p);
    R = recordResistance_(R, br, step);
    R.kappa_S_cm(step,:)   = br.kappa;
    R.Ilim_A(step,:)       = diagn.I_lim;
    R.flux_eq_m2h(step,:)  = diagn.I_ionic / p.F * 3600 / (p.area_cm2/1e4);
    R.saltFraction(step,:) = diagn.I_ionic / p.I_A;
    R.q_counter_C = R.q_counter_C + diagn.q_counter;
    R.q_coion_C   = R.q_coion_C   + diagn.q_coion;
    R.q_split_C   = R.q_split_C   + diagn.q_split;
    R.molMem_mol  = R.molMem_mol  + diagn.molMem;
end
R.precipMol          = st.precipMol;
R.precipMolComp      = st.precipMolComp;
R.formed_mg          = st.formed_mg;        % gross precipitated, per comp/solid
R.dissolvedByComp_mg = st.dissolved_mg;     % gross redissolved
R.solidBulk_mg       = st.solidBulk;
R.scaleMem_mgcm2     = st.scaleMem;

%── Solids balance: bulk + cake must equal what we started with ──────────
% Solids in = fed TSS + everything precipitated. Solids out = bulk TSS +
% particulate cake + crystalline scale.
precipTotal_mg = sum(R.scaleMass_mg(end,:));
expectedSolids_mg = initialSolids_mg + precipTotal_mg;
if expectedSolids_mg > 1e-12
    finalSolids_mg = sum(R.tss_mgL(end,:)' .* R.vol_L(end,:)') ...
                   + sum(R.cake_mgcm2(end,:)  .* p.mem_mult) * p.area_cm2 ...
                   + sum(R.scale_mgcm2(end,:) .* p.mem_mult) * p.area_cm2;
    R.solidsBalance_pct = finalSolids_mg / expectedSolids_mg * 100;
else
    R.solidsBalance_pct = 100;
end

R.t_min = R.t_sec / 60;

%── Cumulative Co2+ current efficiency vs time (eq. 7) ────────────────────
% Named-compartment KPIs. A re-configured stack may not have a 'c2' at all,
% so every lookup below degrades to NaN rather than erroring -- the physics
% does not depend on these names, only the reporting does.
c2Idx = find(strcmp(p.comp_id,'c2'), 1);
R.coEfficiency_pct = nan(nSteps,1);
if isempty(c2Idx)
    R.coEfficiency_pct(:) = NaN;
else
% Current efficiency counts what the CURRENT delivered, so it must use the
% GAIN in the product compartment, not the inventory sitting there. With a
% seeded product the two differ badly: 0.02 M of Co pre-loaded into C2 makes
% the total-based figure read 35% against a true 22%.
coInC2_0 = R.conc(1, c2Idx, p.iCo) * R.vol_L(1, c2Idx);
for step = 1:nSteps
    coGain = R.conc(step, c2Idx, p.iCo) * R.vol_L(step, c2Idx) - coInC2_0;
    theoreticalMax = (p.I_A * R.t_sec(step)) / (p.sp_z(p.iCo) * p.F);
    if theoreticalMax > 1e-12
        R.coEfficiency_pct(step) = coGain/theoreticalMax*100;
    else
        R.coEfficiency_pct(step) = 0;
    end
end
% Deliberately NOT clamped to 100%. A value above 100% is physically
% impossible and means something is wrong -- clamping would hide exactly the
% error worth seeing. runEDM warns instead.
if any(R.coEfficiency_pct > 100 + 1e-6)
    warning('EDM:currentEfficiency', ...
        ['Co current efficiency exceeded 100%% (peak %.1f%%), which is not ' ...
         'physically possible. Check the initial condition and the charge ' ...
         'accounting.'], max(R.coEfficiency_pct));
end
end

%── Mass balance closure vs time (conserved species only) ────────────────
R.massBalance_pct = zeros(nSteps, numel(p.conservedIdx));
for step = 1:nSteps
    for j = 1:numel(p.conservedIdx)
        spIdx = p.conservedIdx(j);
        % Dissolved inventory PLUS whatever has been locked into solids --
        % precipitation removes ions from solution but not from the system,
        % so the balance must count them or scaling would look like a leak.
        tot = sum(squeeze(R.conc(step,:,spIdx))' .* R.vol_L(step,:)') + R.precipMol_hist(step,spIdx);
        if initialTotals(spIdx) > 1e-12
            R.massBalance_pct(step,j) = tot/initialTotals(spIdx)*100;
        else
            R.massBalance_pct(step,j) = 100;
        end
    end
end
R.massBalanceClosure_worst_pct = R.massBalance_pct(end, ...
    find(max(abs(R.massBalance_pct(end,:)-100)) == abs(R.massBalance_pct(end,:)-100), 1));

%── Volume balance and water transport report (eq. 19-20) ───────────────
R.volumeChange_mL   = (R.vol_L(end,:) - p.vol_L') * 1000;
R.volumeBalance_pct = sum(R.vol_L(end,:)) / sum(p.vol_L) * 100;
R.KPI.worstVolumeChange_mL = max(abs(R.volumeChange_mL));
R.KPI.volumeBalance_pct    = R.volumeBalance_pct;
R.KPI.hitVolumeFloor = any(R.vol_L(end,:)' <= p.MIN_VOL_FRACTION*p.vol_L*1.001);

%── Flux and conductivity decline (eq. 24, 25) ──────────────────────────
% Normalised against t = 0, the convention for reporting fouling: a value of
% 1 means unchanged, 0.6 means 40% of the initial flux has been lost.
R.fluxNorm  = R.flux_eq_m2h ./ max(flux0, eps);
R.fluxNorm(:, flux0 <= 0) = NaN;                 % undefined where nothing flowed
kappa0      = R.kappa_S_cm(1,:);
R.kappaNorm = R.kappa_S_cm ./ max(kappa0, eps);

R.KPI.fluxDecline_pct  = 100 * (1 - min(R.fluxNorm(end,:)));
[~, R.KPI.worstFluxMembrane] = min(R.fluxNorm(end,:));
R.KPI.saltFractionEnd_pct = 100 * min(R.saltFraction(end,:));
% Conductivity moves BOTH ways here: dilute compartments lose it, concentrate
% compartments gain it, so report the worst decline and the largest rise.
R.KPI.kappaDecline_pct = 100 * (1 - min(R.kappaNorm(end,:)));
[~, R.KPI.worstKappaComp] = min(R.kappaNorm(end,:));
R.KPI.kappaRise_pct    = 100 * (max(R.kappaNorm(end,:)) - 1);
[~, R.KPI.bestKappaComp]  = max(R.kappaNorm(end,:));

%── Membrane resistance KPIs ────────────────────────────────────────────
R.KPI.Rtotal_ohm     = R.Rtotal_ohm(end);
R.KPI.Rmembranes_ohm = sum(R.Rmem_ohm(end,:));
R.KPI.Rsolution_ohm  = sum(R.Rcomp_ohm(end,:));
[R.KPI.worstMembraneR_ohm, R.KPI.worstMembraneR_idx] = max(R.Rmem_ohm(end,:));
R.KPI.membraneShare_pct = 100 * R.KPI.Rmembranes_ohm / max(R.KPI.Rtotal_ohm, eps);
% How much of each membrane's resistance the run itself created, as opposed
% to the membrane's own clean value.
R.KPI.foulingGrowth_pct = 100 * (sum(R.Rmem_ohm(end,:)) - sum(R.Rmem_clean_ohm(end,:))) ...
                              / max(sum(R.Rmem_clean_ohm(end,:)), eps);

%── Final KPIs ─────────────────────────────────────────────────────────
d1Idx      = find(strcmp(p.comp_id,'d1'), 1);
c1Idx      = find(strcmp(p.comp_id,'c1'), 1);
cathodeIdx = find(strcmp(p.comp_id,'cathode'), 1);
% A re-configured stack may lack any of these. Report NaN for the KPIs that
% depend on a missing compartment instead of failing the whole run.
haveC2 = ~isempty(c2Idx); haveD1 = ~isempty(d1Idx);
haveC1 = ~isempty(c1Idx); haveCat = ~isempty(cathodeIdx);

R.KPI.massBalanceClosure_pct = R.massBalanceClosure_worst_pct;
if haveC2 && haveD1
    % Recovery is what MOVED from the feed into the product, so the product's
    % starting inventory has to be subtracted. Not clamped, for the same
    % reason as above: a figure above 100% is a signal, not a nuisance.
    coGained = R.conc(end,c2Idx,p.iCo)*R.vol_L(end,c2Idx) - ...
               p.conc0(c2Idx,p.iCo)*p.vol_L(c2Idx);
    R.KPI.coRecovered_pct = coGained / max(p.conc0(d1Idx,p.iCo)*p.vol_L(d1Idx),1e-12) * 100;
else
    R.KPI.coRecovered_pct = NaN;
end
if haveC1
    c1 = squeeze(R.conc(end,c1Idx,:))';
else
    c1 = nan(1,p.nSp);
end
wantedEq = c1(p.iNa)*1 + c1(p.iSO4)*2;
% Any metal that reaches C1 is contamination too -- zero under the intended
% topology (C1's only cation-facing membrane draws from the anode rinse),
% but counted so a re-arranged stack reports it honestly.
contamEq = c1(p.iH)*1 + c1(p.iCl)*1 + sum(abs(p.sp_z(p.metalIdx)) .* c1(p.metalIdx));
R.KPI.c1Purity_pct = 100 * wantedEq / max(wantedEq+contamEq, 1e-12);
if haveCat
    R.KPI.cathodeNaGained_mmol = (R.conc(end,cathodeIdx,p.iNa)*R.vol_L(end,cathodeIdx) ...
        - p.conc0(cathodeIdx,p.iNa)*p.vol_L(cathodeIdx)) * 1000;
else
    R.KPI.cathodeNaGained_mmol = NaN;
end
R.KPI.currentEfficiencyCo_pct = R.coEfficiency_pct(end);
R.KPI.membranesOverLimiting = R.overLimitCount(end);
R.KPI.membranesTotal = p.nMem;

%── Solids / fouling KPIs ────────────────────────────────────────────────
R.KPI.solidsBalance_pct     = R.solidsBalance_pct;
R.KPI.feedTurbidity0_NTU    = R.turbidity_NTU(1,   d1Idx);
R.KPI.feedTurbidityEnd_NTU  = R.turbidity_NTU(end, d1Idx);
R.KPI.feedTSS0_mgL          = R.tss_mgL(1,   d1Idx);
R.KPI.feedTSSEnd_mgL        = R.tss_mgL(end, d1Idx);
if R.KPI.feedTSS0_mgL > 1e-12
    R.KPI.feedTSSRemoved_pct = (1 - R.KPI.feedTSSEnd_mgL/R.KPI.feedTSS0_mgL) * 100;
else
    R.KPI.feedTSSRemoved_pct = 0;
end
[R.KPI.worstCake_mgcm2, R.KPI.worstCakeMembrane] = max(R.cake_mgcm2(end,:));
R.KPI.foulingVoltage_V      = R.V_fouling(end);
R.KPI.totalCake_mg          = sum(R.cake_mgcm2(end,:) .* p.mem_mult) * p.area_cm2;

%── Scaling KPIs ─────────────────────────────────────────────────────────
[R.KPI.worstScale_mgcm2, R.KPI.worstScaleMembrane] = max(R.scale_mgcm2(end,:));
R.KPI.scalingVoltage_V = R.V_scaling(end);
R.KPI.totalScale_mg    = sum(R.scaleMass_mg(end,:));
% Which solids actually formed, and where the risk sits even if they didn't.
R.KPI.scaleFormed = p.scale_name(R.scaleMass_mg(end,:) > 1e-9);
SIend = squeeze(R.SI(end,:,:));
R.KPI.maxSI = max(SIend(:));
[iw, kw] = find(SIend == R.KPI.maxSI, 1);
if ~isempty(iw)
    R.KPI.maxSI_where = sprintf('%s in %s', p.scale_name{kw}, p.comp_label{iw});
else
    R.KPI.maxSI_where = '(none)';
end
R.KPI.totalDissolved_mg = sum(R.dissolved_mg(end,:));
R.scaleReport = struct('name',{},'mass_mg',{},'formed_mg',{},'dissolved_mg',{}, ...
                       'maxSI',{},'where',{});
for k = 1:p.nScale
    if R.scaleMass_mg(end,k) <= 1e-9 && max(SIend(:,k)) < 0; continue; end
    [si, ib] = max(SIend(:,k));
    n = numel(R.scaleReport) + 1;
    R.scaleReport(n).name         = p.scale_name{k};
    R.scaleReport(n).mass_mg      = R.scaleMass_mg(end,k);        % still present
    R.scaleReport(n).formed_mg    = sum(R.formed_mg(:,k));         % gross formed
    R.scaleReport(n).dissolved_mg = sum(R.dissolvedByComp_mg(:,k));
    R.scaleReport(n).maxSI        = si;
    R.scaleReport(n).where        = p.comp_label{ib};
end

% Supersaturated but held back by the nucleation barrier -- these are the
% pairs that would have precipitated under the old instantaneous model.
R.KPI.heldByNucleation = 0;
for i = 1:p.nComp
    for k = 1:p.nScale
        present = R.solidBulk_mg(i,k) > 1e-12 || any(R.scaleMem_mgcm2(:,k) > 1e-12);
        if isfinite(SIend(i,k)) && SIend(i,k) > 0 && ~present
            R.KPI.heldByNucleation = R.KPI.heldByNucleation + 1;
        end
    end
end

%── Per-metal transport summary (feed -> C2), for every metal actually fed ─
if haveD1
    R.metalsInFeed = p.metalIdx(p.conc0(d1Idx, p.metalIdx) > 0);
else
    R.metalsInFeed = [];
end
% The per-metal table tracks feed -> product, so it needs both a 'd1' and a
% 'c2'. Without them, skip it rather than fail.
if ~(haveD1 && haveC2); R.metalsInFeed = []; end
nM = numel(R.metalsInFeed);
R.metalReport = struct('label',{},'feed_mmol',{},'d1_final_mmol',{}, ...
                       'removed_pct',{},'c2_final_mmol',{},'recovered_pct',{}, ...
                       'currentEff_pct',{});
for j = 1:nM
    s   = R.metalsInFeed(j);
    m0  = p.conc0(d1Idx,s)      * p.vol_L(d1Idx);
    mD1 = R.conc(end,d1Idx,s)   * R.vol_L(end,d1Idx);
    mC2 = R.conc(end,c2Idx,s)*R.vol_L(end,c2Idx) - p.conc0(c2Idx,s)*p.vol_L(c2Idx);
    thMax = (p.I_A * R.t_sec(end)) / (abs(p.sp_z(s)) * p.F);   % mol if this ion carried ALL charge
    R.metalReport(j).label          = p.sp_label{s};
    R.metalReport(j).feed_mmol      = m0  * 1000;
    R.metalReport(j).d1_final_mmol  = mD1 * 1000;
    R.metalReport(j).removed_pct    = (1 - mD1/max(m0,1e-12)) * 100;
    R.metalReport(j).c2_final_mmol  = mC2 * 1000;
    R.metalReport(j).recovered_pct  = mC2/max(m0,1e-12)*100;
    R.metalReport(j).currentEff_pct = mC2/max(thMax,1e-12)*100;
end

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
fprintf('----------------------------------------------------------------\n');
fprintf('  Membrane  Type       Clean    Cake   Scale   TOTAL    dV      %% of\n');
fprintf('                       (ohm)   (ohm)   (ohm)   (ohm)    (V)   stack\n');
for k = 1:p.nMem
    fprintf('  M%-8d %-8s %7.3f %7.3f %7.3f %7.3f %6.2f %6.1f\n', k, p.membrane_types{k}, ...
        R.Rmem_clean_ohm(end,k), R.Rmem_cake_ohm(end,k), R.Rmem_scale_ohm(end,k), ...
        R.Rmem_ohm(end,k), p.I_A*R.Rmem_ohm(end,k), ...
        100*R.Rmem_ohm(end,k)/max(R.Rtotal_ohm(end),eps));
end
fprintf('  %-19s %23.3f %6.2f %6.1f   (all membranes)\n', '', ...
    R.KPI.Rmembranes_ohm, p.I_A*R.KPI.Rmembranes_ohm, R.KPI.membraneShare_pct);
fprintf('  %-19s %23.3f %6.2f %6.1f   (solutions)\n', '', ...
    R.KPI.Rsolution_ohm, p.I_A*R.KPI.Rsolution_ohm, 100-R.KPI.membraneShare_pct);
if R.KPI.foulingGrowth_pct > 0.05
    fprintf(['  Membrane resistance is %.0f%% above clean -- the excess is cake and\n' ...
             '  scale laid down by the run itself.\n'], R.KPI.foulingGrowth_pct);
end
fprintf('----------------------------------------------------------------\n');
fprintf('  Membrane  Flux eq/(m2 h)   of initial   salt-carried   I_lim\n');
fprintf('              start    end        (%%)        current(%%)     (A)\n');
for k = 1:p.nMem
    fprintf('  M%-8d %6.1f %6.1f %10.1f %14.1f %8.3f\n', k, ...
        R.flux_eq_m2h(1,k), R.flux_eq_m2h(end,k), 100*R.fluxNorm(end,k), ...
        100*R.saltFraction(end,k), R.Ilim_A(end,k));
end
fprintf('  (flux decline = counter-ion flux lost to depletion and fouling;\n');
fprintf('   the shortfall is diverted into water splitting, which does no\n');
fprintf('   separation work)\n');
fprintf('----------------------------------------------------------------\n');
fprintf('  Compartment        kappa start   kappa end   change     R (ohm)\n');
for c = 1:p.nComp
    fprintf('  %-18s %10.5f %11.5f %8.1f%% %9.3f\n', p.comp_label{c}, ...
        R.kappa_S_cm(1,c), R.kappa_S_cm(end,c), 100*(R.kappaNorm(end,c)-1), ...
        R.Rcomp_ohm(end,c));
end
fprintf('  (conductivity falls where a compartment is depleted and rises\n');
fprintf('   where it concentrates -- the dilute side sets the resistance)\n');
if p.ENABLE_WATER_TRANSPORT
    fprintf('----------------------------------------------------------------\n');
    fprintf('  Compartment          V0 (mL)   V_end (mL)    change      %%\n');
    for c = 1:p.nComp
        fprintf('  %-18s %9.1f %11.1f %+10.1f %+7.1f\n', p.comp_label{c}, ...
            p.vol_L(c)*1000, R.vol_L(end,c)*1000, R.volumeChange_mL(c), ...
            100*R.volumeChange_mL(c)/(p.vol_L(c)*1000));
    end
    fprintf('  Stack volume balance : %.4f%% (water moves between chambers,\n', R.volumeBalance_pct);
    fprintf('   it is not created or lost -- evaporation and real leaks are\n');
    fprintf('   NOT modelled, so a measured total loss is a genuine leak)\n');
    if R.KPI.hitVolumeFloor
        warning('EDM:volumeFloor', ...
            ['A compartment hit the minimum-volume floor (%.0f%% of initial). ', ...
             'The run past that point is not physical -- shorten it, lower the ', ...
             'current, or check N_DRAG/LP_WATER.'], p.MIN_VOL_FRACTION*100);
    end
end
fprintf('----------------------------------------------------------------\n');
fprintf('  Compartment          pH     pOH      [H+] M     [OH-] M   pH+pOH\n');
for c = 1:p.nComp
    fprintf('  %-18s %5.2f  %6.2f  %10.3e  %10.3e   %6.2f\n', p.comp_label{c}, ...
        R.pH(end,c), R.pOH(end,c), R.cH_M(end,c), R.cOH_M(end,c), R.pH(end,c)+R.pOH(end,c));
end
if p.ENFORCE_WATER_EQUILIBRIUM
    if p.PH_USE_ACTIVITY
        fprintf('  (pH+pOH sits above pKw = %.2f by -log10(gamma_H*gamma_OH); both are\n', -log10(p.Kw));
        fprintf('   activity-based, so the offset is the ionic-strength correction)\n');
    else
        fprintf('  (pH+pOH = pKw = %.2f confirms water equilibrium is enforced)\n', -log10(p.Kw));
    end
else
    fprintf('  (WARNING: water equilibrium NOT enforced -- H+ and OH- coexist\n');
    fprintf('   out of equilibrium, so these pH values are not meaningful)\n');
end
R.anySolids = p.hasSolids || any(R.tss_mgL(:) > 0) || any(R.scale_mgcm2(end,:) > 0);
if R.anySolids
    fprintf('  Solids balance         : %.2f%% (bulk + cake + scale)\n', R.KPI.solidsBalance_pct);
    fprintf('  Feed turbidity         : %.1f -> %.1f NTU  (%.1f%% of TSS deposited)\n', ...
        R.KPI.feedTurbidity0_NTU, R.KPI.feedTurbidityEnd_NTU, R.KPI.feedTSSRemoved_pct);
    fprintf('  Worst membrane cake    : %.4f mg/cm^2 on membrane %d (%s)\n', ...
        R.KPI.worstCake_mgcm2, R.KPI.worstCakeMembrane, p.membrane_types{R.KPI.worstCakeMembrane});
    fprintf('  Voltage from fouling   : %.3f V of the %.2f V total\n', ...
        R.KPI.foulingVoltage_V, R.V(end));
end
if ~isempty(R.scaleReport)
    fprintf('----------------------------------------------------------------\n');
    fprintf('  Scale            Present  Formed  Dissolved   max SI  Worst compartment\n');
    fprintf('                      (mg)    (mg)       (mg)\n');
    for k = 1:numel(R.scaleReport)
        s = R.scaleReport(k);
        fprintf('  %-14s %8.3f %7.3f %10.3f %8.2f  %s\n', ...
            s.name, s.mass_mg, s.formed_mg, s.dissolved_mg, s.maxSI, s.where);
    end
    fprintf('  Scale on membranes     : %.4f mg/cm^2 worst (membrane %d), +%.3f V\n', ...
        R.KPI.worstScale_mgcm2, R.KPI.worstScaleMembrane, R.KPI.scalingVoltage_V);
    if R.KPI.heldByNucleation > 0
        fprintf('  Held by nucleation     : %d supersaturated pair(s) with no solid yet\n', ...
            R.KPI.heldByNucleation);
    end
    fprintf('  (SI = log10(IAP/Ksp); >0 supersaturated. Nucleation needs SI > %.2f\n', p.SI_NUCLEATION);
    fprintf('   sustained for tau/SI^2; seeded growth and dissolution need no barrier.)\n');
end
if ~isempty(R.metalReport)
    fprintf('----------------------------------------------------------------\n');
    fprintf('  Metal            Feed    D1 left   Removed    C2 gain   Recov.\n');
    fprintf('                  (mmol)   (mmol)      (%%)      (mmol)     (%%)\n');
    for j = 1:numel(R.metalReport)
        m = R.metalReport(j);
        fprintf('  %-12s %8.3f %8.3f %9.1f %10.3f %8.1f\n', ...
            m.label, m.feed_mmol, m.d1_final_mmol, m.removed_pct, m.c2_final_mmol, m.recovered_pct);
    end
    fprintf('  (all metals share the CEM current by eq. 1 -- no chemical\n');
    fprintf('   selectivity is modelled, only |z|*lambda*C weighting)\n');
end
fprintf('══════════════════════════════════════════════\n\n');

end


%% ═══════════════════════════════════════════════════════════════════════
%%  DASHBOARD (2×3)
%% ═══════════════════════════════════════════════════════════════════════
function plotDashboard_EDM(R, p)

save_dir = fullfile(pwd, 'EDM_QuadStack_Plots');
if ~exist(save_dir,'dir'); mkdir(save_dir); end

hDash = figure('Name','EDM Quad Stack Dashboard','Position',[80 60 1600 900]);
% 2x4 when there are solids to show, otherwise the original 2x3.
showVol = p.ENABLE_WATER_TRANSPORT && R.KPI.worstVolumeChange_mL > 1e-6;
showRes = true;                      % resistance breakdown is always useful
nTiles  = 6 + 2*R.anySolids + showVol + showRes + 2;   % +2: flux, conductivity
nCols   = 4;
nRows   = ceil(nTiles/nCols);
tl = tiledlayout(hDash, nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

draw_EQ1_conc(nexttile(tl), R, p);
draw_EQ2_balance(nexttile(tl), R, p);
draw_EQ3_voltage(nexttile(tl), R, p);
if R.anySolids; draw_EQ7_turbidity(nexttile(tl), R, p); end
draw_EQ4_ph(nexttile(tl), R, p);
draw_EQ5_efficiency(nexttile(tl), R, p);
draw_EQ6_overlimit(nexttile(tl), R, p);
if R.anySolids; draw_EQ8_fouling(nexttile(tl), R, p); end
if showVol; draw_EQ9_volume(nexttile(tl), R, p); end
if showRes; draw_EQ10_resistance(nexttile(tl), R, p); end
draw_EQ11_flux(nexttile(tl), R, p);
draw_EQ12_kappa(nexttile(tl), R, p);

title(tl, { ...
    sprintf('EDM Quad Stack — I=%.3f A, A=%.0f cm^2, %.0f h run', p.I_A, p.area_cm2, p.duration_h), ...
    sprintf('Mass balance %.2f%%  |  Co recovered %.1f%%  |  C1 purity %.1f%%  |  %d/%d membranes over-limiting', ...
        R.KPI.massBalanceClosure_pct, R.KPI.coRecovered_pct, R.KPI.c1Purity_pct, R.KPI.membranesOverLimiting, R.KPI.membranesTotal) }, ...
    'FontSize',12, 'FontWeight','bold');

exportgraphics(hDash, fullfile(save_dir,'EDM_Dashboard_Full.png'),'Resolution',150);

% Scaling risk map gets its own figure -- an nComp x nScale grid does not
% compress into a dashboard tile legibly.
SIend = squeeze(R.SI(end,:,:));
if any(isfinite(SIend(:)) & SIend(:) > -6)
    hSI = figure('Name','EDM Scaling Risk','Position',[120 100 1000 520]);
    plotSaturationMap_EDM(hSI, R, p);
    exportgraphics(hSI, fullfile(save_dir,'EDM_Scaling_SI.png'),'Resolution',150);
end
fprintf('  Dashboard saved to: %s\n', save_dir);

end

function draw_EQ1_conc(ax, R, p)
axes(ax); hold(ax,'on');
d1Idx = find(strcmp(p.comp_id,'d1')); c2Idx = find(strcmp(p.comp_id,'c2')); c1Idx = find(strcmp(p.comp_id,'c1'));
% One colour per fed metal: solid = depletion in D1, dashed = build-up in C2.
mIdx = R.metalsInFeed;
cols = lines(max(numel(mIdx),1));
for j = 1:numel(mIdx)
    s = mIdx(j);
    plot(ax, R.t_min, R.conc(:,d1Idx,s), '-',  'Color',cols(j,:),'LineWidth',1.8, ...
        'DisplayName',sprintf('D1 %s', p.sp_label{s}));
    plot(ax, R.t_min, R.conc(:,c2Idx,s), '--', 'Color',cols(j,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('C2 %s', p.sp_label{s}));
end
plot(ax, R.t_min, R.conc(:,d1Idx,p.iSO4),'-', 'Color',[0.31 0.13 0.60],'LineWidth',2, 'DisplayName','D1 SO_4^{2-}');
plot(ax, R.t_min, R.conc(:,c1Idx,p.iSO4),'--','Color',[0.31 0.13 0.60],'LineWidth',1.6,'DisplayName','C1 SO_4^{2-}');
xlabel(ax,'Time (min)'); ylabel(ax,'Concentration (mol/L)');
title(ax,'EQ1: Concentrations (key species)','FontSize',10);
legend(ax,'Location','best','FontSize',7); grid(ax,'on');
end

function draw_EQ2_balance(ax, R, p)
axes(ax); hold(ax,'on');
% Only plot species actually present -- absent ones are flat 100% by
% construction and would just crowd the legend.
present = find(arrayfun(@(s) sum(p.conc0(:,s).*p.vol_L) > 1e-12, p.conservedIdx));
cols = lines(max(numel(present),1));
for k = 1:numel(present)
    j = present(k);
    plot(ax, R.t_min, R.massBalance_pct(:,j), '-', 'Color', cols(k,:), 'LineWidth',1.8, ...
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

function plotSaturationMap_EDM(hFig, R, p)
% Final-state saturation index for every compartment/solid pair.
% Clipped to [-6, +6]: beyond that the number stops being informative, and
% -Inf (an ion completely absent) would otherwise swamp the colour scale.
SIend = squeeze(R.SI(end,:,:));
SIplot = max(min(SIend, 6), -6);
SIplot(~isfinite(SIend)) = NaN;

ax = axes(hFig);
imagesc(ax, SIplot, 'AlphaData', ~isnan(SIplot));
set(ax, 'XTick', 1:p.nScale, 'XTickLabel', p.scale_name, 'XTickLabelRotation', 40, ...
        'YTick', 1:p.nComp,  'YTickLabel', p.comp_label, 'TickLabelInterpreter','none');
% Diverging map centred on SI = 0, the saturation boundary.
colormap(ax, flipud(brewerLikeRdBu_EDM(256)));
caxis(ax, [-6 6]);
cb = colorbar(ax); cb.Label.String = 'SI = log_{10}(IAP/K_{sp})   ( >0 = scaling risk )';
% Annotate the supersaturated cells -- those are the ones that matter.
for i = 1:p.nComp
    for k = 1:p.nScale
        if isfinite(SIend(i,k)) && SIend(i,k) > 0
            text(ax, k, i, sprintf('%+.1f', SIend(i,k)), 'HorizontalAlignment','center', ...
                'FontSize',7,'FontWeight','bold','Color',[0 0 0]);
        end
    end
end
title(ax, {sprintf('Saturation index at t = %.0f h  (blank = ion absent)', p.duration_h), ...
    sprintf('worst: %s, SI = %+.2f', R.KPI.maxSI_where, R.KPI.maxSI)}, 'FontSize',11);
end

function cmap = brewerLikeRdBu_EDM(n)
% Small blue-white-red diverging map, so no toolbox dependency is added.
half = floor(n/2);
b = [linspace(0.13,1,half)', linspace(0.40,1,half)', linspace(0.67,1,half)'];
r = [linspace(1,0.70,n-half)', linspace(1,0.09,n-half)', linspace(1,0.17,n-half)'];
cmap = [b; r];
end

function draw_EQ11_flux(ax, R, p)
% Flux decline, normalised to t = 0. At constant current the total charge
% flux is fixed, so a falling curve means salt transport is being replaced by
% water splitting -- current still flows, but it stops doing separation work.
axes(ax); hold(ax,'on');
cols = lines(p.nMem);
for k = 1:p.nMem
    if all(isnan(R.fluxNorm(:,k))); continue; end
    plot(ax, R.t_min, 100*R.fluxNorm(:,k), '-', 'Color', cols(k,:), 'LineWidth',1.7, ...
        'DisplayName', sprintf('M%d %s', k, p.membrane_types{k}));
end
yline(ax, 100, '--k', 'HandleVisibility','off');
ylim(ax,[0 110]);
xlabel(ax,'Time (min)'); ylabel(ax,'Counter-ion flux (% of initial)');
title(ax, {'EQ11: Flux Decline (eq. 24)', ...
    sprintf('worst %.0f%% lost on M%d', R.KPI.fluxDecline_pct, R.KPI.worstFluxMembrane)}, ...
    'FontSize',10);
legend(ax,'Location','best','FontSize',6); grid(ax,'on');
end

function draw_EQ12_kappa(ax, R, p)
% Conductivity per compartment. It falls where the stack is depleting a
% stream and rises where it is concentrating one; the dilute side is what
% sets the stack resistance.
axes(ax); hold(ax,'on');
cols = lines(p.nComp);
for c = 1:p.nComp
    plot(ax, R.t_min, R.kappa_S_cm(:,c)*1000, '-', 'Color', cols(c,:), 'LineWidth',1.6, ...
        'DisplayName', p.comp_label{c});
end
set(ax,'YScale','log');
xlabel(ax,'Time (min)'); ylabel(ax,'Conductivity (mS/cm)');
title(ax, {'EQ12: Conductivity Decline (eq. 25)', ...
    sprintf('worst -%.0f%% (%s)', R.KPI.kappaDecline_pct, p.comp_id{R.KPI.worstKappaComp})}, ...
    'FontSize',10);
legend(ax,'Location','best','FontSize',6); grid(ax,'on');
end

function draw_EQ10_resistance(ax, R, p)
% Per-membrane resistance over time. The clean value is flat by definition,
% so any rise is entirely cake and scale -- this is the plot that says which
% membrane the run is actually fouling.
axes(ax); hold(ax,'on');
cols = lines(p.nMem);
for k = 1:p.nMem
    plot(ax, R.t_min, R.Rmem_ohm(:,k), '-', 'Color', cols(k,:), 'LineWidth',1.7, ...
        'DisplayName', sprintf('M%d %s', k, p.membrane_types{k}));
end
plot(ax, R.t_min, sum(R.Rcomp_ohm,2), '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.4, ...
    'DisplayName','solutions (all)');
xlabel(ax,'Time (min)'); ylabel(ax,'Resistance (\Omega)');
title(ax, {'EQ10: Resistance per Membrane', ...
    sprintf('total %.1f ohm; membranes %.0f%%%%', R.KPI.Rtotal_ohm, R.KPI.membraneShare_pct)}, ...
    'FontSize',10);
legend(ax,'Location','best','FontSize',6); grid(ax,'on');
end

function draw_EQ9_volume(ax, R, p)
% Compartment volumes over time -- the electro-osmosis/osmosis result that
% shows up in the lab as one vessel filling while its neighbour drains.
axes(ax); hold(ax,'on');
cols = lines(p.nComp);
for c = 1:p.nComp
    plot(ax, R.t_min, R.vol_L(:,c)*1000, '-', 'Color', cols(c,:), 'LineWidth',1.7, ...
        'DisplayName', p.comp_label{c});
end
yline(ax, p.vol_L(1)*1000, '--k', 'HandleVisibility','off');
xlabel(ax,'Time (min)'); ylabel(ax,'Compartment volume (mL)');
title(ax, {'EQ9: Water Transport (eq. 19,20)', ...
    sprintf('worst %+.1f mL; stack total %.3f%%', ...
        max(R.volumeChange_mL(abs(R.volumeChange_mL)==max(abs(R.volumeChange_mL)))), ...
        R.volumeBalance_pct)}, 'FontSize',10);
legend(ax,'Location','best','FontSize',6); grid(ax,'on');
end

function draw_EQ7_turbidity(ax, R, p)
axes(ax); hold(ax,'on');
cols = lines(p.nComp);
shown = false;
for c = 1:p.nComp
    if max(R.turbidity_NTU(:,c)) <= 0; continue; end   % skip clean compartments
    plot(ax, R.t_min, R.turbidity_NTU(:,c), '-', 'Color', cols(c,:), 'LineWidth',1.8, ...
        'DisplayName', p.comp_label{c});
    shown = true;
end
xlabel(ax,'Time (min)'); ylabel(ax,'Turbidity (NTU)');
title(ax, {'EQ7: Turbidity / TSS (eq. 9,10)', ...
    sprintf('feed %.1f -> %.1f NTU', R.KPI.feedTurbidity0_NTU, R.KPI.feedTurbidityEnd_NTU)}, 'FontSize',10);
if shown; legend(ax,'Location','best','FontSize',6); end
grid(ax,'on');
% Right axis in mg/L -- same curve, the units people actually spec against.
yyaxis(ax,'right');
ylabel(ax,'TSS (mg/L)');
ylim(ax, ax.YAxis(1).Limits * p.MGL_PER_NTU);
ax.YAxis(2).Color = [0.35 0.35 0.35];
yyaxis(ax,'left');
end

function draw_EQ8_fouling(ax, R, p)
axes(ax); hold(ax,'on');
cols = lines(p.nMem);
for k = 1:p.nMem
    if max(R.cake_mgcm2(:,k)) > 0
        plot(ax, R.t_min, R.cake_mgcm2(:,k), '-', 'Color', cols(k,:), 'LineWidth',1.8, ...
            'DisplayName', sprintf('M%d %s cake', k, p.membrane_types{k}));
    end
    if max(R.scale_mgcm2(:,k)) > 0
        plot(ax, R.t_min, R.scale_mgcm2(:,k), '--', 'Color', cols(k,:), 'LineWidth',1.6, ...
            'DisplayName', sprintf('M%d %s scale', k, p.membrane_types{k}));
    end
end
xlabel(ax,'Time (min)'); ylabel(ax,'Surface loading (mg/cm^2)');
title(ax, {'EQ8: Fouling & Scaling (eq. 10,11,14)', ...
    sprintf('+%.3f V cake, +%.3f V scale', R.KPI.foulingVoltage_V, R.KPI.scalingVoltage_V)}, 'FontSize',10);
legend(ax,'Location','best','FontSize',6); grid(ax,'on');
end

function draw_EQ6_overlimit(ax, R, p)
axes(ax);
plot(ax, R.t_min, R.overLimitCount, '-o', 'Color',[0.71 0.46 0.04],'LineWidth',1.8,'MarkerSize',3,'MarkerFaceColor',[0.71 0.46 0.04]);
ylim(ax,[-0.3 p.nMem+0.3]); yticks(ax, 0:p.nMem);
xlabel(ax,'Time (min)'); ylabel(ax,'Membranes over-limiting (of 5)');
title(ax,'EQ6: Over-Limiting Membranes (eq. 8)','FontSize',10);
grid(ax,'on');
end
