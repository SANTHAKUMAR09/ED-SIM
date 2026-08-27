# Firebase Hosting deploy

Deploys `WebUI/edm_stack_simulator.html` to Firebase Hosting, project `edm-sim-hosting`
(dedicated GCP project, separate from any Skawix production projects).

Two Hosting sites serve the same content:

| Site | URL | Target name |
|---|---|---|
| `edm-sim-hosting` (default) | https://edm-sim-hosting.web.app | `main` |
| `arwa-edm-sim` | **https://arwa-edm-sim.web.app** (preferred, shorter link) | `arwa` |

## To redeploy after updating the simulator

```bash
cp ../WebUI/edm_stack_simulator.html public/index.html
firebase deploy --only hosting --project edm-sim-hosting   # deploys BOTH sites
```

To deploy to just one site: `firebase deploy --only hosting:arwa` or `--only hosting:main`.

## Adding another site later

```bash
firebase hosting:sites:create <new-site-id> --project edm-sim-hosting
firebase target:apply hosting <target-name> <new-site-id> --project edm-sim-hosting
# then add a matching entry to firebase.json's "hosting" array
```
