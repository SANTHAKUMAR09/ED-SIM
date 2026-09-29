# ED-SIM — Project Status

Portable session handoff for this repo. Update and commit this file at every
session checkpoint ("save this session") so any machine can `git pull` and
resume with full context — this repo has no dependency on ARWA or its
memory system.

**Last updated:** 2026-09-29.

## Overview

Standalone project for electrodialysis (ED) applied to metal recovery / acid
regeneration, split out of the ARWA AMD treatment-train repo because it's a
distinct effort (see ARWA repo's own memory for why). GitHub:
`SANTHAKUMAR09/ED-SIM`, private.

The project went through two design generations:
1. A first simple CEM/cation ED model (steady-state, Ni/Co only) — **parked**
   by explicit user decision, tagged `v1-cem-cation-model`.
2. The **current design**: a full electrodialysis metathesis (EDM) "quad
   stack" model based on real lab hardware topology the user provided, with
   a corrected membrane termination (CEM at both end positions) to avoid
   metal plating at the cathode and Cl₂ evolution at the anode. This is what
   all active work now targets.

## What exists (current EDM quad-stack design)

- **`MATLAB/EDM_QuadStack_Simulator.m`** — the reference model, ~2440 lines.
  27 equations, 15 tracked species (10 metals: Co/Ni/Fe/Al/Cr/Mn/Mg/Cu/Zn/Ca,
  plus Na/SO4/Cl/H/OH), 6 compartments (cathode rinse, D2, C2, D1/feed, C1,
  anode rinse), 5 membranes (CEM-AEM-CEM-AEM-CEM). Covers ion-competition
  transport (Faraday's law + Kohlrausch transport-number weighting),
  Nernst-Einstein-derived diffusivities, limiting current + over-limiting
  water-splitting, ionic-strength conductivity correction, multi-metal
  precipitation/scaling chemistry (Ksp, nucleation induction time,
  dissolution), suspended-solids/membrane fouling, Lévêque-correlation
  flow→boundary-layer coupling, electro-osmotic drag + osmosis (water
  transport), N-repeating-unit stack tiling, monovalent-selective and
  bipolar membrane branches, permselectivity/co-ion leakage,
  concentration-polarization film resistance, and a configurable initial
  pH per compartment (`cfg.pH0`/`H0_M`/`OH0_M`) with automatic acid/base
  counter-ion dosing (`PH_ACID_ANION`, default `'SO4'` for H2SO4 — any
  species in `sp_id` works, it's not tied to one specific acid).
  Source of truth for all physics — the web app is a verified 1:1 port
  of this file.
- **`MATLAB/EDM_QuadStack_Walkthrough.m`** — beginner-friendly flat script
  (no functions) version of the model, for understanding each module
  step-by-step. **Now in sync** with the extended 27-equation model (as of
  2026-09-29): all 15 species/10 metals, precipitation/scaling, water
  transport, fouling, monovalent/bipolar membranes, and the configurable
  Initial pH feature. Deliberately narrower than the reference file in two
  places (flagged in its own header): always N=1 (no repeating-unit
  tiling) and one fixed default configuration (no generic config-override
  mechanism) — neither is physics, both are reference-model engineering
  conveniences. Verified via Octave to produce numerically identical
  results to `EDM_QuadStack_Simulator.m` (see below).
- **`WebUI/edm_stack_simulator.html`** — self-contained interactive web app
  (vanilla JS, embedded fonts, canvas charts, no external libraries/CDN).
  Full port of `EDM_QuadStack_Simulator.m` — numerically verified identical
  across every KPI (mass balance, recovery, purity, current efficiency,
  voltage, full 9-solid scaling table to 3 decimal places). Features: live
  animated stack schematic (per-compartment species concentrations), an
  equations documentation section (all 27 equations as label + formula
  cards — no per-card descriptions, by explicit user request; a symbol
  legend below the cards is the only explanatory text left), a KPI
  dashboard (mass balance, Co²⁺ and Ni²⁺ recovery, avg. metal recovery,
  C1 purity, current efficiency, fouling/scale voltage, etc.),
  current-efficiency / mass-balance / voltage / scaling charts, topology
  safety badges (flags trace metal reaching cathode rinse or Cl⁻ reaching
  anode rinse — a real, small-magnitude consequence of co-ion leakage
  physics, not a bug), a per-compartment **Initial pH** input (with
  a selectable acid anion for the counter-ion it doses — Cl⁻/HCl or
  SO4²⁻/H2SO4, default SO4²⁻) mirroring the MATLAB model's `cfg.pH0`,
  and a full 6-selector membrane panel for N>1 repeating-unit tiling
  (not just the 5 N=1 selectors — see the corrected D2-start topology's
  cathode|D2 / D2|C2 / C2|D1 / D1|C1 / C1|D2-next-unit / C1|anode rows).
- **Live deployment**: **https://arwa-edm-sim.web.app** (Firebase Hosting,
  multi-site target `arwa` under GCP project `edm-sim-hosting`, free
  `.web.app` subdomain — no custom domain purchased). Deploy source is
  `firebase-deploy/public/index.html`, kept in sync with
  `WebUI/edm_stack_simulator.html` (copy before `firebase deploy --only
  hosting:arwa`). The old default Hosting site (`edm-sim-hosting.web.app`)
  is disabled (`firebase hosting:disable`), not deletable outright.
- **`assets/edm_sim_qr.png`** — QR code pointing to the live URL above,
  generated via Python `qrcode`, verified via `opencv-python-headless`.
- **`AMD_Electrodialysis_MetalRecovery.m`** + `Presets/` + `scripts/` +
  `References/electrodialysis_metal_recovery_papers.xlsx` — the **parked**
  first-generation CEM/cation model and its literature-review artifacts.
  Left in place for reference; not part of the active design. See
  `v1-cem-cation-model` git tag for the state at which this was parked.

## Latest session's changes (2026-09-29)

- Merged `EDM_QuadStack_Simulator13.m` (updated on another machine) into
  `MATLAB/EDM_QuadStack_Simulator.m`. The real change: a configurable
  **initial pH per compartment** (`cfg.pH0`, or `cfg.H0_M`/`cfg.OH0_M`
  directly) instead of a hard-wired pH 7 everywhere, with automatic
  acid/base counter-ion dosing (`PH_DOSE_COUNTERION`, default on) so a
  non-neutral start stays electroneutral, plus new dashboard reporting
  of initial pH/[H+]/[OH-]/dosing per compartment. Also restored this
  repo's accurate 27-equation header comment — the incoming file's
  header had reverted to a stale "same 8 equations as the web tool"
  template that no longer matched the actual (untouched) 27-equation
  code body.
- Ported the same initial-pH feature to the web app: each compartment's
  Solutions tab gets an **Initial pH** input with a live hint showing
  the resulting [H+]/[OH-] and any counter-ion dosed.
- Fixed a modeling gap the user caught: the counter-ion dosing was
  hardcoded to Cl⁻ (i.e. assumed HCl), but the user's actual acid is
  H2SO4. Made the acid anion **user-selectable** (Cl⁻/HCl or SO4²⁻/H2SO4)
  in both the web app (dropdown, default SO4²⁻) and MATLAB
  (`PH_ACID_ANION`, default changed from `'Cl'` to `'SO4'`) — and fixed
  the stoichiometry to divide by the anion's charge, since a divalent
  anion like SO4²⁻ needs half the moles Cl⁻ would for the same eq/L of
  acid. Removed the "dose counter-ion" on/off checkbox from the web UI
  (it was always checked, no real use case for turning it off) — dosing
  stays permanently on in both files.
- Removed more explanatory hint text per explicit user request: under
  Repeating units (N) and under the Membranes card.
- Reorganized the Solutions panel's per-compartment rows (Initial pH,
  Loop volume, Flow, Suspended solids): label + hint now group together
  on the left, with the input consistently pinned to the end of the row
  (previously the input sat sandwiched between the label and a
  variable-length hint with uneven spacing).
- Deployed all of the above to the live Firebase site
  (https://arwa-edm-sim.web.app) — confirmed live: acid-anion selector
  present and defaulting to SO4, old dosing checkbox gone, new row
  layout in place, Ni²⁺ KPI still present.
- Removed the remaining equation-number references from the left-panel
  labels ("(eq. 12–16)", "(eq. 19–20)", "(eq. 18)", "(eq. 9-11)" under
  Chemistry options and the Flow/Suspended solids rows).
- **Full rewrite of `MATLAB/EDM_QuadStack_Walkthrough.m`** to match the
  27-equation model (previously flagged as out of sync for several
  sessions — it still reflected the original 8-equation/2-metal
  version). Now covers all 15 species/10 metals, the Initial-pH feature,
  monovalent/bipolar membranes, permselectivity/co-ion leakage,
  concentration polarisation, precipitation/nucleation/dissolution for
  all 11 solids, suspended-solids fouling, and electro-osmotic-drag +
  osmosis water transport — one flat, heavily-commented script, same
  style as before. Deliberately narrower than the reference file in two
  places (flagged in its own header): always N=1, and one fixed default
  configuration rather than a generic config-override API — neither is
  physics. **Installed Octave (`brew install octave`) specifically to
  verify this** — no MATLAB available in this environment otherwise.
  Confirmed byte-for-byte-matching results against
  `EDM_QuadStack_Simulator.m` run the same way (mass balance 100.00%,
  Co recovered 47.0%, C1 purity 65.2%, current efficiency 26.2%, voltage
  2.83 V, and all 10 metals' recovery to 3+ significant figures), plus a
  smoke test of the new Initial-pH/TSS knobs (D1 to pH 2 + 200 mg/L TSS
  runs cleanly, mass/solids balance still 100.00%). One portability fix
  found by this verification: MATLAB allows local functions at the end
  of a script, Octave does not, so the three small helper functions
  (bisection, areal-deposit-resistance) are inlined in the loop instead.
- **Fixed N>1 membrane selection** (user-reported via a live test with
  screenshots): setting M1 to Mono-CEM at N=1 then raising Repeating
  units to N=2 silently reverted it to plain CEM. Root cause: N>1 tiles
  into a 6-row membrane table (cathode|D2, D2|C2, C2|D1, D1|C1, C1|D2 of
  the *next* unit, C1|anode — the "next unit" handoff is why it can't be
  a simple 5-slot linear chain), and the web app only ever exposed the 5
  N=1 selectors, always building the N>1 tile from a hardcoded default
  array. Added a second 6-selector panel (T1–T6) shown only when N>1 so
  every tiling row is actually choosable, matching what
  `EDM_QuadStack_Simulator.m` already supports via a 6-entry
  `membrane_types`. This also surfaced a second, independent display bug:
  the schematic's membrane badges always showed hardcoded defaults for
  N>1 regardless of selection (pre-existing, not caused by this feature) —
  fixed by mapping the 6-row tile table onto the 5 visible gaps (the
  internal C1|D2-next-unit row has no single gap in the lumped
  6-compartment picture, so it's the one row not shown visually).
  Verified in-browser: N=2 + T1=Mono-CEM now shows "Mono-CEM" in the
  schematic; switching back to N=1 restores the 5-selector panel; a run
  completes with mass balance at 100.00%.
- Deployed all of the above (walkthrough rewrite is MATLAB-only, not
  deployed) to the live Firebase site.

## Previous session's changes (2026-09-15)

- Stripped all explanatory text out of the Model Equations section: the
  intro paragraph under the "Model equations" heading, all 22 per-card
  `eq-desc` blocks (the 8 original + 14 added last session), and the
  trailing "Model scope" disclaimer note. Each of the 22 cards now shows
  only its label and formula; the symbol legend was kept since it's now
  the sole explanation of what each symbol means. Also removed the hint
  line under the Chemistry Options toggles ("Both default on, matching
  the MATLAB reference model...").
- Added an **"Ni²⁺ recovered → C2"** KPI tile to the dashboard, next to
  the existing Co²⁺ one — same `metalRecoveryAt()` helper already used
  for the 10-metal average, just applied to `'Ni'`.
- Deployed both changes to the live Firebase site
  (https://arwa-edm-sim.web.app) — confirmed live via a fresh page load
  (0 `eq-desc`/`eq-intro`/`footer-note` nodes, 22 eq-cards intact, Ni KPI
  present in the strip).

## Earlier session's changes (2026-09-09)

- Ported the full extended MATLAB model (27 equations / 15 species / 10
  metals, up from the original 8 equations / 2 metals) into the web app —
  a complete rewrite of the JS physics engine and UI, explicitly scoped by
  the user as "full port, everything." Verified line-for-line against
  MATLAB output.
- Fixed the live stack schematic: compartment box height increased
  120px → 210px (SVG height 225 → 320) so a compartment's full species list
  (up to 13 lines for D1/Feed) fits inside its box instead of overflowing
  past the bottom edge. Confirmed via direct browser testing at both
  desktop and mobile (375px) widths — the SVG uses `viewBox` scaling so the
  fix holds proportionally at any screen size.
- Expanded the equations documentation section: eq. 5b and 9–27 were
  previously shown only as a compact one-line summary table; replaced with
  14 full equation cards (formula + explanatory description each) matching
  the depth/style of the original cards 1–8, plus 20 new legend entries for
  the symbols they introduce (μ, ε, ℓ, Ksp/IAP, Sh/Re/Sc, n_drag, α, etc.).
- Deployed both fixes to the live Firebase site
  (https://arwa-edm-sim.web.app) — confirmed live via a fresh page load
  showing 22 eq-cards and no leftover summary table.
- Mobile-responsiveness check on the live site: single-column layout below
  980px (existing media query) and the `eq-grid`'s `auto-fit, minmax(340px,
  1fr)` naturally collapses to one column on phone widths — both changes
  above render correctly with no changes needed for mobile.

## Key literature reviewed in depth

1. **Isaksson et al. 2025** (Paper #1) — Ni/Co recovery via electrodialysis
   metathesis (EDM), EDTA-chelated leachate, Membranes journal.
   DOI: 10.3390/membranes15040097. Full methods extracted (membrane
   sequence, all 4 solution streams, flow/current/voltage) via PMC. Key
   result: 97.9% Ni / 96.6% Co separated at 0.10 M. This paper's real quad
   topology (with the termination correction) is what the current design
   generation is modeled on.
2. **Xing & Srinivasan 2023** — Li recovery via chelating-agent-facilitated
   bipolar-membrane ED (BMED) from real industrial LIB leachate. Chemical
   Engineering Journal, DOI: 10.1016/j.cej.2023.145306. Tested EDTA/HEDTA/
   GLDA/DTPA — DTPA best (63.9% Li recovery / 99.4% purity at optimum
   dosage). Reviewed for the parked first-generation model; not yet
   revisited against the current quad-stack design.

## Validated facts / model behavior (current quad-stack design)

- Reference-case run (0.30 A, 10 cm² area, N=1, default balanced feed):
  mass balance 100.00% closure, 47.0% Co recovered, 65.2% C1 purity, 26.2%
  current efficiency, 0/5 membranes over-limiting, 2.83 V stack voltage —
  matched exactly between MATLAB and the web app, including the full
  9-solid scaling table to 3 decimal places.
- Real lab CC test (0.5 A, up to 10 V) showed voltage rising from low to
  ~10 V and plateauing, not spiking — this is why the over-limiting-current
  + water-splitting mechanism (eq. 8) was implemented as real physics
  rather than a hardcoded cap, per explicit user request.
- Trace metal reaching the cathode rinse / trace Cl⁻ reaching the anode
  rinse over a long run is a genuine, small-magnitude (~1e-5–1e-3 M)
  consequence of the permselectivity/co-ion-leakage mechanism (eq. 23),
  confirmed by direct numeric inspection — not a bug.

## Open items / not yet done

- No acid-regeneration (bipolar-membrane ED) module in the current design.
- No train-level integration (this is intentionally independent of ARWA).
- (Parked v1 model only) Chemistry proxy factors for DTPA/HEDTA/GLDA/none
  and the pH curve shape were illustrative assumptions, not fit to a
  Ni/Co-specific dataset.

## How to run

```matlab
% Full extended quad-stack reference model
R = EDM_QuadStack_Simulator();   % see file header for options/signature

% Beginner walkthrough (now in sync with the reference model, see above)
run('MATLAB/EDM_QuadStack_Walkthrough.m')

% Parked first-generation model
R = AMD_Electrodialysis_MetalRecovery();
run('scripts/run_paper1_scale_match.m')
run('scripts/run_chelator_pH_secondaryfeed_sensitivity.m')
```

Octave (`brew install octave`) works for both `.m` files above and is
useful as a free sanity-check runner when MATLAB itself isn't at hand —
just note `plotDashboard_EDM`'s use of `tiledlayout` isn't supported on
older Octave, so the numeric report/warnings print fine but the dashboard
figure itself will error; the walkthrough's plain `subplot`-based plots
don't hit that issue.

```bash
# Web app: open locally
open WebUI/edm_stack_simulator.html

# Deploy to the live Firebase site (after copying the latest WebUI file
# into firebase-deploy/public/index.html)
cd firebase-deploy && firebase deploy --only hosting:arwa
```

## Workflow

Direct commits to `main` so far (no feature-branch+PR requirement set for
this repo, unlike ARWA — revisit if the user asks for that workflow here
too). Dashboard/plot PNGs and the web app's deploy copy are regenerated in
place and re-committed alongside the source that produced them. Firebase
CLI credentials expire periodically — `firebase login --reauth` needs to be
run interactively by the user in their own terminal when this happens.
