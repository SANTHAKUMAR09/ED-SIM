# ED-SIM

MATLAB simulation work for electrodialysis (ED) applied to acid mine
drainage (AMD) metal recovery / acid regeneration. This is an independent
project — not part of, and not wired into, the ARWA AMD treatment-train
repository.

## Modules

### `AMD_Electrodialysis_MetalRecovery.m`

Standalone, steady-state model of an ED stack recovering Ni and Co from a
synthetic dilute feed into a recirculating concentrate loop, with SO4
tracked as the balancing counter-ion (AEM side).

**Model** (algebraic, steady-state — no ODE integration):

```
I            = i_app * A_mem                     [A, per cell pair]
Ndot_i       = N_cp * eta_i * I / (z_i * F) * 3600   [mol/h transferred]
C_i,dil_out  = C_i,in - Ndot_i / Q_d              [diluate outlet]
Q_c          sized so Ni hits CF_target           [concentrate loop flow]
C_i,conc_out = Ndot_i / Q_c                       [concentrate, fresh-water makeup]
i_lim_i      = z_i * F * D_i * C_i,dil_out / delta   [limiting current, worst-case outlet]
V_stack      = N_cp * (V_junction + i_app * R_area)  [stack voltage]
```

Key outputs: diluate in/out concentrations and removal %, concentrate loop
concentrations and achieved concentration factor, a limiting-current-density
design check (concentration polarization), stack voltage/power, and
specific energy consumption (kWh/kg metal, kWh/m³ treated).

**Usage:**

```matlab
% Default synthetic feed + stack design, with dashboard
R = AMD_Electrodialysis_MetalRecovery();

% Override feed and/or stack parameters
feed_ov.Ni_mgL  = 300;
feed_ov.Q_in    = 100;
stack_ov.N_cp   = 20;
stack_ov.i_app  = 2;
R = AMD_Electrodialysis_MetalRecovery(feed_ov, stack_ov, true);

% Skip dashboard rendering (faster, for sweeps)
R = AMD_Electrodialysis_MetalRecovery([], [], false);
```

Dashboard PNGs are written to `ElectrodialysisMetalRecovery_Plots/`
(full 2×2 dashboard + one PNG per panel).

## Design notes

- Default feed is a moderately concentrated stream (Ni 250 mg/L, Co 60 mg/L,
  SO4 5000 mg/L) — representative of something like an ion-exchange
  eluate/bleed, not raw dilute AMD effluent. At tens-of-mg/L Ni, the
  limiting current density collapses to a fraction of an A/m², making
  single-stage ED impractical; that's a physical constraint of the model
  (see the `i_lim` check / EP3 panel), not a tunable.
- `eta_Ni + eta_Co` must stay ≤ 1: both share the same cation-side current
  budget through the CEM. `eta_SO4` is an independent anion-side budget
  through the AEM. Any shortfall from either budget is implicitly carried
  by untracked ions (Ca, Mg, Na, H⁺ / Cl⁻, HCO₃⁻).
- The concentrate loop flow `Q_c` is sized to hit the target concentration
  factor for **Ni only**. Other tracked species (Co, SO4) land wherever
  their own transfer rate puts them at that flow — SO4 in particular can
  come out *under* 1× feed concentration in the same loop, since its molar
  transfer relative to its large feed concentration is small. Concentrating
  SO4/acid-side products (e.g. for acid regeneration) would need separate
  loop sizing — not yet modeled here.

## Status

Single steady-state metal-recovery module, verified end-to-end in MATLAB
(no ODE, single function call). No sensitivity analysis, no train-level
integration, no acid-regeneration (bipolar-membrane ED) module yet.
