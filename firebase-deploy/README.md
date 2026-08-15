# Firebase Hosting deploy

Deploys `WebUI/edm_stack_simulator.html` to Firebase Hosting, project `edm-sim-hosting`
(dedicated GCP project, separate from any Skawix production projects).

Live URL: https://edm-sim-hosting.web.app

## To redeploy after updating the simulator

```bash
cp ../WebUI/edm_stack_simulator.html public/index.html
firebase deploy --only hosting --project edm-sim-hosting
```
