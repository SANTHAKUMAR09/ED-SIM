%% ═══════════════════════════════════════════════════════════════════════
%%  EDM QUAD STACK — BEGINNER WALKTHROUGH SCRIPT
%% ═══════════════════════════════════════════════════════════════════════
%
%  This is the SAME simulation as EDM_QuadStack_Simulator.m (same physics,
%  same numbers, same results) but written as ONE FLAT SCRIPT instead of a
%  file full of functions. Nothing here is hidden inside a function call --
%  every line runs top to bottom, in order, and every variable it creates
%  stays visible in your Workspace panel afterwards so you can click on it
%  and look at the actual numbers.
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
%   row from one end to the other. Positively-charged ions (like Na+, Co2+)
%   drift toward one end (the cathode); negatively-charged ions (like Cl-,
%   SO4^2-) drift toward the other end (the anode). Each membrane only lets
%   ONE type of charge through, so ions pile up in some tanks and get
%   stripped out of others -- that's how this stack separates and
%   concentrates the cobalt.
%
%   Tank layout, left (cathode, negative electrode) to right (anode, positive):
%     Tank 1: Cathode rinse   (just a rinse solution, sits at the electrode)
%     Tank 2: D2              (0.5 M NaCl "donor" solution)
%     Tank 3: C2               (starts dilute -- this is where Co2+ ends up)
%     Tank 4: D1 = Feed         (the actual cobalt-bearing solution)
%     Tank 5: C1                (starts dilute -- this is where SO4^2- ends up)
%     Tank 6: Anode rinse       (just a rinse solution, sits at the electrode)
%
%   Membrane layout (5 membranes, one between each pair of neighboring tanks):
%     Membrane 1 (between Tank 1 and Tank 2): type CEM (passes + ions only)
%     Membrane 2 (between Tank 2 and Tank 3): type AEM (passes - ions only)
%     Membrane 3 (between Tank 3 and Tank 4): type CEM (passes + ions only)
%     Membrane 4 (between Tank 4 and Tank 5): type AEM (passes - ions only)
%     Membrane 5 (between Tank 5 and Tank 6): type CEM (passes + ions only)
%
%   This specific left-right arrangement was chosen deliberately: it keeps
%   cobalt away from Tank 1 (so it never plates onto the cathode) and keeps
%   chloride away from Tank 6 (so no chlorine gas forms at the anode).

clear;      % wipe the Workspace so you start from a clean slate
close all;  % close any old figure windows
clc;        % clear the Command Window text


%% SECTION 1 — Physical constants
% These three numbers never change; they come from physics/chemistry, not
% from anything about OUR stack.

F    = 96485;     % Faraday's constant, in Coulombs per mole of charge.
                   % This is "how much electric charge is carried by one
                   % mole of singly-charged ions." Every calculation that
                   % converts between "electric current" and "how many
                   % ions moved" goes through this number.
Kw   = 1e-14;      % the water equilibrium constant (used later for pH)
Rgas = 8.314;      % the gas constant, Joules per mole per Kelvin
T    = 298.15;     % temperature, in Kelvin (this is 25 degrees Celsius --
                   % room temperature, and the temperature the other
                   % constants below were measured at)


%% SECTION 2 — The 7 ions we're tracking
% Every ion in this simulation is one of these 7. We give each one:
%   - a short name
%   - z = its electric charge (how many + or - charges it carries)
%   - lambda = a measure of how "mobile" it is in water (bigger number =
%     moves faster under the same electric field). These are standard,
%     published values for each ion at 25 degrees C.
%
% We store all 7 as entries in matching lists (this is just 7 separate
% pieces of information about each of the 7 ions, kept in parallel arrays).

species_name   = {'Co2+', 'Ca2+', 'Na+', 'SO4^2-', 'Cl-',  'H+',   'OH-'};
species_charge = [ 2,      2,      1,     -2,       -1,     1,      -1  ];
species_mobility = [53.0,  59.5,   50.1,   80.0,     76.3,   349.8,  198.0];  % lambda, S*cm^2/eq

% We'll refer to each ion by its position in these lists all through the
% script. Rather than remembering "species number 1 means cobalt", we give
% each position a readable name right now:
iCo  = 1;   % Co2+  is species #1
iCa  = 2;   % Ca2+  is species #2
iNa  = 3;   % Na+   is species #3
iSO4 = 4;   % SO4^2- is species #4
iCl  = 5;   % Cl-   is species #5
iH   = 6;   % H+    is species #6
iOH  = 7;   % OH-   is species #7
n_species = 7;

% H+ and OH- are deliberately created/destroyed by the electrode and
% water-splitting reactions later in this script, so they're never
% expected to "balance" -- everything ELSE should. We list those 5 here so
% both the running history (Section 9) and the final check (Section 10)
% use the exact same definition of "conserved."
conserved_ions = [iCo, iCa, iNa, iSO4, iCl];

% Every ion also has a "diffusivity" D -- roughly, how fast it spreads out
% on its own (no electric field needed). We don't look this up separately;
% it's mathematically related to the mobility number above by a formula
% called the Nernst-Einstein relation. We compute it once, for all 7 ions
% in one line:
species_diffusivity = (Rgas .* T .* species_mobility) ./ (abs(species_charge) .* F.^2);
% (units: cm^2/s. If you're curious, this line alone reproduces the
%  textbook diffusivity of Na+ to 3 significant figures: 1.33e-5 cm^2/s.)


%% SECTION 3 — The 6 tanks (compartments) and what's in them at the start
% "conc" is a table (a matrix): 6 rows (one per tank), 7 columns (one per
% ion), holding the STARTING concentration of every ion in every tank, in
% mol/L (molarity). This table gets updated every timestep as the
% simulation runs -- think of it as "the current state of all 6 tanks."

tank_name = {'Cathode rinse', 'D2 (donor)', 'C2', 'D1 (Feed)', 'C1', 'Anode rinse'};
n_tanks = 6;

%              Co2+     Ca2+     Na+      SO4^2-   Cl-      H+     OH-
conc = [ ...
          0,       0,       0.500,   0.250,   0,       1e-7,  1e-7;   % Tank 1: Cathode rinse
          0,       0,       0.500,   0,       0.500,   1e-7,  1e-7;   % Tank 2: D2 (0.5 M NaCl)
          0,       0,       0.050,   0,       0.050,   1e-7,  1e-7;   % Tank 3: C2 (dilute seed)
          0.0500,  0.0050,  0.1082,  0.0541,  0.1100,  1e-7,  1e-7;   % Tank 4: D1 = Feed
          0,       0,       0.100,   0.050,   0,       1e-7,  1e-7;   % Tank 5: C1 (dilute seed)
          0,       0,       0.500,   0.250,   0,       1e-7,  1e-7];  % Tank 6: Anode rinse

% Every tank's volume, in liters. All six start at 250 mL = 0.25 L.
tank_volume_L = [0.25; 0.25; 0.25; 0.25; 0.25; 0.25];

% Remember the STARTING total of each ion, added up across all 6 tanks.
% We'll compare against this at the very end to check nothing was created
% or destroyed by mistake (a "mass balance" sanity check).
starting_total_moles = sum(conc .* tank_volume_L, 1);   % 1x7, one total per ion


%% SECTION 4 — The 5 membranes between the tanks
% membrane_type{m} tells you what kind of membrane sits between Tank m and
% Tank (m+1). 'CEM' = cation-exchange membrane (only + ions cross it).
% 'AEM' = anion-exchange membrane (only - ions cross it).

membrane_type = {'CEM', 'AEM', 'CEM', 'AEM', 'CEM'};
n_membranes = 5;


%% SECTION 5 — How we're running the stack (the "knobs" you can turn)
% These are the settings you'd type into a real power supply and a timer.
% Try changing these later and re-running the whole script to see what
% happens -- that's the easiest way to build intuition.

area_cm2          = 10;    % membrane area, cm^2
applied_current_A = 0.5;   % constant current, Amps
duration_hours    = 4;     % how long to run the simulation for
timestep_minutes  = 1;     % how big a "tick" of simulated time we take

timestep_seconds = timestep_minutes * 60;
total_seconds    = duration_hours * 3600;
n_steps          = round(total_seconds / timestep_seconds);   % number of ticks


%% SECTION 6 — Placeholder numbers for the voltage/resistance model
% These describe things about the PHYSICAL hardware (spacer thickness,
% membrane resistance, etc.) that we haven't measured for your specific
% stack -- they're reasonable engineering estimates, flagged here so you
% know exactly which numbers to replace once you have real measurements.

areal_resistance_per_membrane = 5;      % ohm*cm^2
compartment_gap_cm            = 0.5;    % cm (spacer thickness)
electrode_overpotential_V     = 1.5;    % V (fixed)
water_conductivity_floor      = 5.5e-8; % S/cm (pure water's own conductivity --
                                          % a tank can never have LESS
                                          % conductivity than plain water)
voltage_ceiling_V             = 200;    % V (a safety cap, should basically
                                          % never be hit if eq. 8 below is
                                          % working correctly)

% "Over-limiting current" settings (this is the part of the model that
% explains why real stack voltage plateaus instead of climbing forever --
% see the long comment inside the main loop below for the full story).
boundary_layer_cm       = 50e-4;  % cm (50 micron)
water_split_overpotential_V = 0.8; % V, added per membrane that's "maxed out"


%% SECTION 7 — Set up empty storage for the whole run's history
% We're about to simulate n_steps+1 points in time (the "+1" is for the
% very start, before anything has happened). For each of those points we
% want to remember: every tank's concentrations, every tank's pH, the
% stack voltage, and how many membranes were "maxed out." So we make empty
% containers now and fill them in as we go.

time_minutes           = zeros(n_steps+1, 1);
conc_history           = zeros(n_steps+1, n_tanks, n_species);
pH_history              = zeros(n_steps+1, n_tanks);
voltage_history         = zeros(n_steps+1, 1);
overlimit_count_history = zeros(n_steps+1, 1);
co_efficiency_history   = zeros(n_steps+1, 1);
mass_balance_history    = zeros(n_steps+1, 5);   % one column per conserved ion (Co,Ca,Na,SO4,Cl)


%% SECTION 8 — Record the very first point (t = 0), before we run anything
% This block computes "what is the pH and voltage of the stack RIGHT NOW,
% before any current has flowed" and saves it as row 1 of our history. The
% pH and voltage calculations here are the exact same calculations we'll
% repeat inside the main loop below -- seeing them once, standalone, here,
% should make them easier to recognize later.

time_minutes(1) = 0;
conc_history(1,:,:) = conc;

for tank = 1:n_tanks
    % --- pH of this tank, from its H+ and OH- concentrations ---
    c_H  = conc(tank, iH);
    c_OH = conc(tank, iOH);
    net = c_H - c_OH;
    if net > 1e-14
        c_H_effective = net;
    elseif net < -1e-14
        c_H_effective = Kw / (-net);
    else
        c_H_effective = 1e-7;   % no real signal either way -> call it neutral
    end
    pH_history(1, tank) = -log10(max(c_H_effective, 1e-14));
end

% --- stack voltage at t=0 (no membranes are over-limiting yet, so that
%     part of the formula contributes zero) ---
total_resistance = 0;
for tank = 1:n_tanks
    conductivity = 0.001 * sum(species_mobility .* abs(species_charge) .* max(conc(tank,:),0));
    conductivity = max(conductivity, water_conductivity_floor);
    total_resistance = total_resistance + compartment_gap_cm / (conductivity * area_cm2);
end
total_resistance = total_resistance + n_membranes * (areal_resistance_per_membrane / area_cm2);
voltage_history(1) = min(applied_current_A*total_resistance + electrode_overpotential_V, voltage_ceiling_V);

overlimit_count_history(1) = 0;
co_efficiency_history(1) = 0;

% Mass balance at t=0 is trivially 100% (nothing has moved yet), but we
% compute it with the real formula anyway so Section 9 can just repeat the
% exact same line every step without a special case for step 1.
current_totals = sum(conc .* tank_volume_L, 1);
mass_balance_history(1,:) = current_totals(conserved_ions) ./ starting_total_moles(conserved_ions) * 100;


%% SECTION 9 — THE MAIN LOOP: step the simulation forward in time
% This is the heart of the whole script. We repeat this block once per
% timestep (n_steps times total). Each pass through the loop:
%   (a) works out how many MOLES of each ion are in each tank right now
%   (b) for each of the 5 membranes, works out how many ions cross it
%       during this one timestep, and moves them
%   (c) applies the two electrode reactions (these always happen,
%       regardless of the membranes)
%   (d) converts moles back into concentrations, and saves this timestep
%       into the history arrays from Section 7

for step = 1:n_steps

    %% --- (a) how many moles of each ion are in each tank, right now ---
    % moles = concentration x volume. This is the same table shape as
    % "conc" (6 tanks x 7 ions), just in moles instead of mol/L.
    moles = conc .* tank_volume_L;

    overlimit_flag_this_step = false(1, n_membranes);

    %% --- (b) go through the 5 membranes one at a time ---
    for m = 1:n_membranes

        left_tank  = m;       % the tank on the cathode side of this membrane
        right_tank = m + 1;   % the tank on the anode side of this membrane

        if strcmp(membrane_type{m}, 'CEM')
            % CEM passes CATIONS (+ ions). Cations always drift TOWARD the
            % cathode, i.e. from the right-hand tank into the left-hand tank.
            source_tank      = right_tank;
            destination_tank = left_tank;
            carrier_sign     = 1;    % we only care about POSITIVE ions here
        else
            % AEM passes ANIONS (- ions). Anions always drift TOWARD the
            % anode, i.e. from the left-hand tank into the right-hand tank.
            source_tank      = left_tank;
            destination_tank = right_tank;
            carrier_sign     = -1;   % we only care about NEGATIVE ions here
        end

        % Which of the 7 ions are even ALLOWED through this membrane?
        % (the ones whose charge sign matches carrier_sign)
        carrier_ions = find(sign(species_charge) == carrier_sign);

        % How much of the SOURCE tank's current content is made up of each
        % of those carrier ions? (can't be negative, even from floating-
        % point rounding, hence the max(...,0))
        source_conc_of_carriers = max(conc(source_tank, carrier_ions), 0);

        % --- EQUATION 1: transport-number "share" ---
        % If several ion types are all allowed through this membrane, they
        % don't split the current 50/50 -- each one gets a SHARE of the
        % current proportional to (its charge) x (its mobility) x (its
        % concentration). "denom" is the sum of that quantity across every
        % carrier ion -- we divide by it in a moment to turn each ion's raw
        % number into a fraction that adds up to exactly 1 (100%).
        denom = sum(abs(species_charge(carrier_ions)) .* species_mobility(carrier_ions) .* source_conc_of_carriers);

        % --- EQUATION 8: how much current can this membrane ACTUALLY carry
        %     just by moving the ions that are already there? ---
        % This is called the "limiting current." Physically: ions have to
        % physically diffuse to the membrane surface to cross it, and they
        % can only do that so fast. boundary_layer_cm stands in for "how
        % thick the sluggish layer of water right next to the membrane is."
        limiting_current_A = 0;
        for k = 1:numel(carrier_ions)
            ion = carrier_ions(k);
            limiting_current_A = limiting_current_A + ...
                abs(species_charge(ion)) * F * species_diffusivity(ion) * source_conc_of_carriers(k) / 1000 / boundary_layer_cm;
        end
        limiting_current_A = limiting_current_A * area_cm2;

        % If the current we're TRYING to push through is bigger than what
        % the ions can carry, this membrane is "over-limiting" -- the ions
        % alone can't keep up, and (see part below) water-splitting has to
        % pick up the slack.
        is_over_limit = applied_current_A > limiting_current_A;
        overlimit_flag_this_step(m) = is_over_limit;

        % The current actually available to move real ions is capped at
        % whatever the membrane can sustain:
        ionic_current_A = min(applied_current_A, limiting_current_A);

        % --- EQUATION 2: move each carrier ion, in proportion to its share ---
        if denom > 1e-30 && ionic_current_A > 1e-30
            for k = 1:numel(carrier_ions)
                ion = carrier_ions(k);
                share = (abs(species_charge(ion)) * species_mobility(ion) * source_conc_of_carriers(k)) / denom;
                flux_mol_per_second = share * ionic_current_A / (abs(species_charge(ion)) * F);
                moles_to_move = flux_mol_per_second * timestep_seconds;

                % Safety cap: never try to move more of an ion than the
                % source tank actually has left (this prevents a rare
                % rounding glitch where a nearly-empty tank would otherwise
                % go slightly negative).
                moles_to_move = min(moles_to_move, max(moles(source_tank, ion), 0));

                moles(source_tank, ion)      = moles(source_tank, ion)      - moles_to_move;
                moles(destination_tank, ion) = moles(destination_tank, ion) + moles_to_move;
            end
        end

        % --- water-splitting for the "leftover" current (eq. 8, continued) ---
        % If this membrane IS over-limiting, the shortfall (applied_current
        % minus limiting_current) doesn't just vanish -- physically, water
        % (H2O) splits apart right at the membrane surface into H+ and OH-,
        % and those brand-new ions carry the rest of the current. A CEM
        % lets the new H+ through (so it shows up in the destination tank);
        % the new OH- is left behind in the source tank (OH- is a NEGATIVE
        % ion, so a CEM blocks it). An AEM does the mirror image.
        if is_over_limit
            leftover_current_A = applied_current_A - limiting_current_A;
            moles_from_water_split = (leftover_current_A / F) * timestep_seconds;
            if strcmp(membrane_type{m}, 'CEM')
                moles(destination_tank, iH)  = moles(destination_tank, iH)  + moles_from_water_split;
                moles(source_tank,      iOH) = moles(source_tank,      iOH) + moles_from_water_split;
            else
                moles(destination_tank, iOH) = moles(destination_tank, iOH) + moles_from_water_split;
                moles(source_tank,      iH)  = moles(source_tank,      iH)  + moles_from_water_split;
            end
        end

    end   % end of the "go through the 5 membranes" loop


    %% --- (c) the two electrode reactions ---
    % These happen at the very ends of the stack (Tank 1 and Tank 6),
    % regardless of what any membrane is doing. At the cathode (Tank 1),
    % water is reduced to hydrogen gas and hydroxide: this always adds
    % OH- to Tank 1, at a rate set directly by the applied current. At the
    % anode (Tank 6, the last one), water is oxidized to oxygen gas and
    % acid: this always adds H+ to Tank 6.
    moles(1,        iOH) = moles(1,        iOH) + (applied_current_A/F) * timestep_seconds;
    moles(n_tanks,   iH) = moles(n_tanks,   iH)  + (applied_current_A/F) * timestep_seconds;


    %% --- (d) turn moles back into concentrations, and save this step ---
    conc = max(moles ./ tank_volume_L, 0);

    time_minutes(step+1) = time_minutes(step) + timestep_minutes;
    conc_history(step+1,:,:) = conc;
    overlimit_count_history(step+1) = sum(overlimit_flag_this_step);

    % pH of every tank, same formula as Section 8
    for tank = 1:n_tanks
        c_H  = conc(tank, iH);
        c_OH = conc(tank, iOH);
        net = c_H - c_OH;
        if net > 1e-14
            c_H_effective = net;
        elseif net < -1e-14
            c_H_effective = Kw / (-net);
        else
            c_H_effective = 1e-7;
        end
        pH_history(step+1, tank) = -log10(max(c_H_effective, 1e-14));
    end

    % Stack voltage, same formula as Section 8, but now including the
    % water-splitting overpotential for however many membranes are
    % over-limiting THIS step.
    total_resistance = 0;
    for tank = 1:n_tanks
        conductivity = 0.001 * sum(species_mobility .* abs(species_charge) .* max(conc(tank,:),0));
        conductivity = max(conductivity, water_conductivity_floor);
        total_resistance = total_resistance + compartment_gap_cm / (conductivity * area_cm2);
    end
    total_resistance = total_resistance + n_membranes * (areal_resistance_per_membrane / area_cm2);
    water_split_voltage = overlimit_count_history(step+1) * water_split_overpotential_V;
    voltage_history(step+1) = min(applied_current_A*total_resistance + electrode_overpotential_V + water_split_voltage, voltage_ceiling_V);

    % Cumulative Co2+ current efficiency: of all the charge we've pushed
    % through so far, what fraction actually ended up moving Co2+ into
    % Tank 3 (C2)? "theoretical_max" is what we'd have gotten if every
    % single electron so far had gone to Co2+ and nothing else.
    co_in_C2 = conc(3, iCo) * tank_volume_L(3);
    theoretical_max_mol = (applied_current_A * time_minutes(step+1)*60) / (species_charge(iCo) * F);
    if theoretical_max_mol > 1e-12
        co_efficiency_history(step+1) = min(co_in_C2/theoretical_max_mol*100, 100);
    end

    % Mass balance closure: total moles of each conserved ion, summed
    % across all 6 tanks, as a percentage of what we started with. Should
    % stay at (very close to) 100% every single step -- if it drifts, that
    % means ions are being created or destroyed somewhere by mistake.
    current_totals = sum(conc .* tank_volume_L, 1);
    mass_balance_history(step+1,:) = current_totals(conserved_ions) ./ starting_total_moles(conserved_ions) * 100;

end   % end of the main time-stepping loop


%% SECTION 10 — Check the results: did we conserve mass? What did we get?
% "Conserved" ions (everything except H+ and OH-, which are deliberately
% created/destroyed by the electrode and water-splitting reactions) should
% sum to EXACTLY the same total at the end as they started with. This is
% the single best sanity check that the code above has no bugs in it.

ending_total_moles = sum(conc .* tank_volume_L, 1);
mass_balance_pct = ending_total_moles(conserved_ions) ./ starting_total_moles(conserved_ions) * 100;
worst_mass_balance_pct = mass_balance_pct(find(abs(mass_balance_pct-100) == max(abs(mass_balance_pct-100)), 1));

co_recovered_pct = min( conc(3,iCo)*tank_volume_L(3) / max(conc_history(1,4,iCo)*tank_volume_L(4), 1e-12) * 100, 100);

wanted_equivalents      = conc(5,iNa)*1 + conc(5,iSO4)*2;
contamination_equivalents = conc(5,iH)*1 + conc(5,iCl)*1;
c1_purity_pct = 100 * wanted_equivalents / max(wanted_equivalents+contamination_equivalents, 1e-12);

cathode_Na_gained_mmol = (conc(1,iNa) - conc_history(1,1,iNa)) * tank_volume_L(1) * 1000;

fprintf('\n══════════════════════════════════════════════\n');
fprintf('  EDM Quad Stack (walkthrough script) — Results\n');
fprintf('══════════════════════════════════════════════\n');
fprintf('  Mass balance closure    : %.2f%%\n', worst_mass_balance_pct);
fprintf('  Co2+ recovered -> C2    : %.1f%%\n', co_recovered_pct);
fprintf('  C1 product purity       : %.1f%%\n', c1_purity_pct);
fprintf('  Cathode Na+ gained      : %.2f mmol\n', cathode_Na_gained_mmol);
fprintf('  Current efficiency, Co  : %.1f%%\n', co_efficiency_history(end));
fprintf('  Membranes over-limiting : %d / %d\n', overlimit_count_history(end), n_membranes);
fprintf('  Final stack voltage     : %.2f V\n', voltage_history(end));
fprintf('══════════════════════════════════════════════\n\n');


%% SECTION 11 — Plot the results
% Six small charts, one per thing we tracked. Nothing fancy -- this is
% deliberately close to the simplest MATLAB plotting code there is, so you
% can see exactly how "numbers in an array" become "a line on a chart."

figure('Name','EDM Quad Stack — Walkthrough Results','Position',[80 60 1600 900]);

subplot(2,3,1);
plot(time_minutes, conc_history(:,4,iCo), 'LineWidth',2); hold on;
plot(time_minutes, conc_history(:,3,iCo), 'LineWidth',2);
plot(time_minutes, conc_history(:,4,iSO4), 'LineWidth',2);
plot(time_minutes, conc_history(:,5,iSO4), 'LineWidth',2);
xlabel('Time (min)'); ylabel('Concentration (mol/L)');
legend('D1 (Feed) Co2+','C2 Co2+','D1 (Feed) SO4^{2-}','C1 SO4^{2-}','Location','best');
title('Key concentrations'); grid on;

subplot(2,3,2);
plot(time_minutes, mass_balance_history, 'LineWidth',1.5);
legend(species_name(conserved_ions), 'Location','best','FontSize',7);
yline(100,'--k','HandleVisibility','off');
xlabel('Time (min)'); ylabel('Mass balance closure (%)');
ylim([95 105]);
title('Mass balance closure'); grid on;

subplot(2,3,3);
plot(time_minutes, voltage_history, 'LineWidth',2);
xlabel('Time (min)'); ylabel('Stack voltage (V)');
title(sprintf('Voltage (final %.2f V)', voltage_history(end))); grid on;

subplot(2,3,4);
plot(time_minutes, pH_history, 'LineWidth',1.5);
xlabel('Time (min)'); ylabel('pH'); ylim([0 14]);
legend(tank_name, 'Location','best','FontSize',7);
title('pH of every tank'); grid on;

subplot(2,3,5);
plot(time_minutes, co_efficiency_history, 'LineWidth',2);
xlabel('Time (min)'); ylabel('Co^{2+} current efficiency (%)'); ylim([0 100]);
title(sprintf('Current efficiency (final %.1f%%)', co_efficiency_history(end))); grid on;

subplot(2,3,6);
plot(time_minutes, overlimit_count_history, '-o','LineWidth',1.5,'MarkerSize',3);
xlabel('Time (min)'); ylabel('Membranes over-limiting (of 5)'); ylim([-0.3 5.3]);
title('Over-limiting membranes'); grid on;


%% SECTION 12 — Things to try next
% Now that you've seen every line, here are some easy experiments. Change
% ONE number in Section 5 or Section 6, then run the WHOLE script again
% (green Run button) and compare the printed results / charts to before.
%
%   - applied_current_A = 0.3   (lower current -> slower, gentler run)
%   - duration_hours = 8         (run twice as long -- does D1 ever fully
%                                  drain? does the voltage keep climbing?)
%   - boundary_layer_cm = 20e-4  (thinner boundary layer -> higher limiting
%                                  current -> fewer membranes go over-limit)
%   - Try changing the STARTING concentrations in Section 3, e.g. make the
%     Feed (Tank 4) more concentrated in Co2+ and see how that changes
%     Co2+ recovery.
