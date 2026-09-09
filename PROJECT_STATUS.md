# ED-SIM — Project Status

Portable session handoff for this repo. Update and commit this file at every
session checkpoint ("save this session") so any machine can `git pull` and
resume with full context — this repo has no dependency on ARWA or its
memory system.

**Last updated:** 2026-09-09.

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
  bipolar membrane branches, permselectivity/co-ion leakage, and
  concentration-polarization film resistance. Source of truth for all
  physics — the web app is a verified 1:1 port of this file.
- **`MATLAB/EDM_QuadStack_Walkthrough.m`** — beginner-friendly flat script
  (no functions) version of the model, for understanding each module
  step-by-step. **Out of sync**: still reflects the old 8-equation/2-metal
  version, not yet updated to match the extended 27-equation model. Known
  gap, not yet requested to be fixed.
- **`WebUI/edm_stack_simulator.html`** — self-contained interactive web app
  (vanilla JS, embedded fonts, canvas charts, no external libraries/CDN).
  Full port of `EDM_QuadStack_Simulator.m` — numerically verified identical
  across every KPI (mass balance, recovery, purity, current efficiency,
  voltage, full 9-solid scaling table to 3 decimal places). Features: live
  animated stack schematic (per-compartment species concentrations), an
  equations documentation section (all 27 equations as full formula+
  description cards), current-efficiency / mass-balance / voltage /
  scaling charts, topology safety badges (flags trace metal reaching
  cathode rinse or Cl⁻ reaching anode rinse — a real, small-magnitude
  consequence of co-ion leakage physics, not a bug).
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

## Latest session's changes (2026-09-09)

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

- `MATLAB/EDM_QuadStack_Walkthrough.m` (beginner flat-script version) is
  out of sync with the extended 27-equation model — still reflects the old
  8-equation/2-metal version.
- No acid-regeneration (bipolar-membrane ED) module in the current design.
- No train-level integration (this is intentionally independent of ARWA).
- (Parked v1 model only) Chemistry proxy factors for DTPA/HEDTA/GLDA/none
  and the pH curve shape were illustrative assumptions, not fit to a
  Ni/Co-specific dataset.

## How to run

```matlab
% Full extended quad-stack reference model
R = EDM_QuadStack_Simulator();   % see file header for options/signature

% Beginner walkthrough (currently out of sync — old 8-eq model)
run('MATLAB/EDM_QuadStack_Walkthrough.m')

% Parked first-generation model
R = AMD_Electrodialysis_MetalRecovery();
run('scripts/run_paper1_scale_match.m')
run('scripts/run_chelator_pH_secondaryfeed_sensitivity.m')
```

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
