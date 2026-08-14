# ED-SIM — Project Status

Portable session handoff for this repo. Update and commit this file at every
session checkpoint ("save this session") so any machine can `git pull` and
resume with full context — this repo has no dependency on ARWA or its
memory system.

## Overview

Standalone MATLAB simulation project for electrodialysis (ED) applied to
metal recovery / acid regeneration, split out of the ARWA AMD treatment-train
repo because it's a distinct effort (see ARWA repo's own memory for why).
GitHub: `SANTHAKUMAR09/ED-SIM`, private.

## What exists

- **`AMD_Electrodialysis_MetalRecovery.m`** — steady-state ED stack model,
  standalone (no upstream module dependency), synthetic feed. Tracks Ni/Co
  (recovered cations, CEM side) and SO4 (counter-anion, AEM side). Faraday's
  law transfer, limiting-current-density design check, stack voltage/power/
  specific-energy-consumption. Dashboard PNGs to
  `ElectrodialysisMetalRecovery_Plots/`.
  - Extended with a **chemistry proxy**: `p.chelator` / `p.secondary_feed` /
    `p.pH` apply assumption-flagged multiplier factors onto `eta_Ni`/`eta_Co`
    (see `chelatorFactor_ED`/`secondaryFeedFactor_ED`/`pHFactor_ED` in the
    file). Calibrated to factor=1.0 at Isaksson et al. 2025's exact recipe
    (EDTA, Na2SO4, pH 7). **Not a mechanistic chelation model** — only the
    EDTA/Na2SO4/pH-7 baseline is literature-anchored; other values are
    illustrative proxies. Do not present sensitivity results from this as
    validated quantitative correlations.
- **`Presets/AMD_ED_Preset_Isaksson2025_Paper1.m`** — reusable feed/stack
  config matching Paper #1's hardware (5 cell pairs, 64 cm² membrane area,
  200 A/m² initial current density, 0.10 M NiSO4+CoSO4 feed, 20 L/h flow)
  plus their batch-run metadata. Pattern for future paper presets:
  `Presets/AMD_ED_Preset_<Author><Year>_Paper<N>.m`.
- **`scripts/run_paper1_scale_match.m`** — runs the core model against the
  Paper #1 preset, then translates single-pass removal into an equivalent
  recirculating-batch pass count/time (~147 passes ≈ 110 min) to sanity-check
  against Paper #1's reported batch endpoint.
- **`scripts/run_chelator_pH_secondaryfeed_sensitivity.m`** — sweeps
  chelator × secondary-feed and pH × secondary-feed on the Paper #1 preset;
  writes 2 CSVs + a 2-panel dashboard to `ChemistrySensitivity_Plots/`.
- **`References/electrodialysis_metal_recovery_papers.xlsx`** — literature
  tracking workbook, 3 sheets:
  - **"Nickel and Cobalt Recovery"** (first tab) — isolated deep-dive on
    Paper #1 (full profile: membrane sequence, 4-stream solution sequence,
    flow, voltage/current, run endpoint, findings) + a comparison table
    against Xing & Srinivasan 2023 (BMED, Li recovery).
  - **"Papers Matrix"** — general literature tracking table (title, metals,
    method, membrane, performance, scale, link), 18 rows as of last count.
  - **"Summary"** — pre-existing summary sheet.

## Key literature reviewed in depth

1. **Isaksson et al. 2025** (Paper #1) — Ni/Co recovery via electrodialysis
   metathesis (EDM), EDTA-chelated leachate, Membranes journal.
   DOI: 10.3390/membranes15040097. Full methods extracted (membrane
   sequence, all 4 solution streams, flow/current/voltage) via PMC.
   Key result: 97.9% Ni / 96.6% Co separated at 0.10 M.
2. **Xing & Srinivasan 2023** — Li recovery via chelating-agent-facilitated
   bipolar-membrane ED (BMED) from real industrial LIB leachate. Chemical
   Engineering Journal, DOI: 10.1016/j.cej.2023.145306. Full text pulled
   from NTU's open-access institutional repository. Tested EDTA/HEDTA/GLDA/
   DTPA — DTPA best (63.9% Li recovery / 99.4% purity at optimum dosage).

## Validated facts / model behavior

- Default synthetic feed (Ni 250 mg/L, Co 60 mg/L, SO4 5000 mg/L) gives
  58.7% Ni / 67.0% Co single-pass removal, i_lim check passes at 63% margin.
- At tens-of-mg/L Ni (dilute AMD-effluent-level feeds), the limiting current
  density collapses below any practical i_app — single-stage ED is
  physically impractical at that concentration in this model. This is a
  real constraint (concentration-polarization physics), not a tunable bug.
- Paper #1 hardware-matched run (via the preset): 3.3% Ni / 0.9% Co removed
  per single pass, i_lim OK at 32% margin — consistent with Paper #1's
  cumulative batch-recirculation result once translated to ~147 equivalent
  passes (~110 min at their 20 L/h / 250 cm³ donor volume).
- `eta_Ni + eta_Co` must stay ≤ 1 (shared CEM-side current budget) — this
  was an actual bug caught and fixed early in the session (defaults were
  0.85/0.85, summing to 1.7).

## How to run

```matlab
% Default synthetic feed + stack, with dashboard
R = AMD_Electrodialysis_MetalRecovery();

% Paper #1 hardware-matched run
run('scripts/run_paper1_scale_match.m')

% Chemistry (chelator/pH/secondary-feed) sensitivity sweep
run('scripts/run_chelator_pH_secondaryfeed_sensitivity.m')
```

## Open items / not yet done

- No mechanistic EDM/BMED transport model (chelated anionic complex via
  AEM) — deferred in favor of the CEM/cation model + chemistry-proxy
  approach (explicit user decision this session).
- No acid-regeneration (bipolar-membrane ED) module.
- No train-level integration (this is intentionally independent of ARWA).
- Chemistry proxy factors for DTPA/HEDTA/GLDA/none and the pH curve shape
  are illustrative assumptions, not fit to a Ni/Co-specific dataset — flag
  this if the sensitivity results are ever presented outside this repo.
- Papers Matrix has more rows (~18) than have been deep-dived — only
  Paper #1 and the Xing & Srinivasan paper have full-text-level detail so
  far.

## Workflow

Direct commits to `main` so far (no feature-branch+PR requirement set for
this repo, unlike ARWA — revisit if the user asks for that workflow here
too). Dashboard/plot PNGs are regenerated in place by scripts and re-committed
alongside the scripts that produced them.
