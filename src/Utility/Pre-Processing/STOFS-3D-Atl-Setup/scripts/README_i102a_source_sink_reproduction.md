# I102a source/sink milestone reproduction

Run `reproduce_i102a_from_initial.py` on a compute node after the initial NWM
generation and USGS flow replacement have completed. It starts with the results
in `initial_nwm_usgs/results` and stops after the artificial-island post-patch.
It does not run constant-sink generation.

The runner stages the initial source/sink files in a new output directory. It
copies the USGS-adjusted `vsource.th` because pre-relocation zeroing rewrites
that file. The completed initial results and HJ's reference are not modified.

Use the same Python environment as `initial_nwm_usgs/run.sh`. From the
STOFS-3D-Atl-Setup repository directory:

```bash
export PYTHONPATH="$PWD/src"
/sciclone/home/feiye/mambaforge/envs/stofs/bin/python -u scripts/reproduce_i102a_from_initial.py \
  --preflight \
  --initial-results /sciclone/schism10/feiye/STOFS3D-v7.4/I102a_test/initial_nwm_usgs/results \
  --reference-input /sciclone/schism10/hjyoo/task/task10_Atlantic/stofs3d-setup_v7p4/I102a \
  --output /sciclone/schism10/feiye/STOFS3D-v7.4/I102a_test/I102a_full_chain_milestone
```

If preflight passes, run the same command without `--preflight`:

```bash
set -o pipefail
/sciclone/home/feiye/mambaforge/envs/stofs/bin/python -u scripts/reproduce_i102a_from_initial.py \
  --initial-results /sciclone/schism10/feiye/STOFS3D-v7.4/I102a_test/initial_nwm_usgs/results \
  --reference-input /sciclone/schism10/hjyoo/task/task10_Atlantic/stofs3d-setup_v7p4/I102a \
  --output /sciclone/schism10/feiye/STOFS3D-v7.4/I102a_test/I102a_full_chain_milestone \
  2>&1 | tee /sciclone/schism10/feiye/STOFS3D-v7.4/I102a_test/I102a_full_chain_milestone.log
```

The script refuses to overwrite an existing output directory. Budget several gigabytes
for the new output; the copied adjusted flow, zeroed flow, relocated forcing,
and post-patch forcing total roughly 2.5 GB before diagnostics and filesystem
overhead.

The final report is `I102a_full_chain_milestone/reproduction_result.json`.
It checks byte identity for the initial NWM/static files, pre-relocation
zeroed flow, relocation mapping, relocated forcing, and final
`source_sink.in`, `vsource.th`, `msource.th`, and
`vsink.th`. It compares every `source.nc` dimension and variable, allowing the
previously observed `1e-10 m3/s` absolute tolerance only for `vsource`.
Exit status 0 means all checks passed; status 2 means at least one comparison
failed. The report also confirms that no constant-sink directory was created.
