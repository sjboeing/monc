### Relevant additions

# General plumbing changes needed to add new components and their inputs in
- global_config: also fully disabled the IO server
- components/componentheaders.static
- components/componentregistations.static

# New components
- components/simple_tke: makes TKE and dissipation available for CASIM
- components/simple_2d_diagnostics: per-rank 2D diagnostic output to bypass IO server. Needs stitching afterward.
- components/simple_netcdf_dump: per-rank 3D diagnostic output to bypass IO server. Needs stitching afterward.
- model_core/src/utils/simple_io_utils.F90: underlying shared functionality
- model_core/makefile: changes needed to actually add the diagnostics.

# Changes related to allocation of new "state" fields for TKE and dissipation
- components/simplesetup/src/simplesetup.F90 
- model_core/src/state.F90
- components/checkpointer/src/readcheckpoint.F90: adding the fields, but also additional error catching 

# Other
- model_core/src/components/timestepper.F90: commented out calls to handle_sampling, as they rely on the io server. I realise this is a bit of a hack.
Note: we check the CFL every time step in order to make the CFL work without the sampling.

# Job submission
- job_example.sh: example continuation job (I had trouble modifying/reading the usual MONC continuation job so set this up)
  Invokes tiny scripts
  - check_continuation.py: to check if continuation is needed
  - submit_next.sh: submits the next job at the end of a run
  - get_checkpoint_from_mcf.py: to find the relevant checkpoint name

# Test case
- testcases/arm_diurnal_cycle/arm_sb.mcf: ARM case including CASIM and setup for new diagnostics

# Stitching
- stitch_v2.py: script to stitch python files. This script needs a python environment with dask. It seems to work better on ARCHER2 (in terms of dask parallelism) but for now it is OK on MONSOON3.
- stitch_chain.pbs: submits this as a job in continuation-job style 
