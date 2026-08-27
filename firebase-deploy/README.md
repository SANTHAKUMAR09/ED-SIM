# Firebase Hosting deploy

Deploys `WebUI/edm_stack_simulator.html` to Firebase Hosting, project `edm-sim-hosting`
(dedicated GCP project, separate from any Skawix production projects).

Live site: **https://arwa-edm-sim.web.app** (target name `arwa`).

The project's default site (`edm-sim-hosting.web.app`, target `main`) was disabled via
`firebase hosting:disable --site edm-sim-hosting` -- it now 404s. Firebase doesn't allow
deleting a project's default site outright (only non-default sites can be deleted), so
"disabled" is the closest equivalent to removal; the site record still exists but serves
nothing. `firebase.json` only lists the `arwa` target, so a normal `firebase deploy` won't
accidentally re-enable it.

## To redeploy after updating the simulator

```bash
cp ../WebUI/edm_stack_simulator.html public/index.html
firebase deploy --only hosting --project edm-sim-hosting
```

## Adding another site later

```bash
firebase hosting:sites:create <new-site-id> --project edm-sim-hosting
firebase target:apply hosting <target-name> <new-site-id> --project edm-sim-hosting
# then add a matching entry to firebase.json's "hosting" array
```
