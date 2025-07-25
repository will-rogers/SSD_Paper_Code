# SSD_Paper_Code
 
Simultations:

    Simulate movements:

        run_steps.slurm:
            requires: "base_fxns.R", "sur.val.packed.RDS"
            returns: "output_steps/used_X.rds"


    Fit RSF and SSF models:

        run_rsf.slurm: 
            requires: "base_fxns.R", "sims.RDS", "output_steps/used_X.rds"
            returns: "output_rsfs/rsfs_X.rds"

        run_ssf.slurm: 
            requires: "base_fxns.R", "sims.RDS", "output_steps/used_X.rds"
            returns: "output_ssfs/ssfs_X.rds"


    Run simulations for validation:

        run_simulations.slurm:
            requires: "base_fxns.R", "sur.val.packed.RDS"
            returns: "output_rasters/used_X.rds"


    Run ABMs: 

        run_abm.slurm:
            requires: "base_fxns.R", "sims.RDS", "index.sheet.csv", "output_ssfs/ssfs_X.rds", "output_steps/used_x.rds"
            returns: "output_abms/abms_1_X.rds"

        run_abm_size.slurm:
            requires: "base_fxns.R", "abm.example.RDS"
            returns: "output_abm_size/"

        run_abm_example.slurm:
            requires: "base_fxns.R", "abm.example.RDS"
            returns: "output_abm_example/abm_steps_X.rds"


    Map RSF, SSF, and SSD predictions: 

        run_rsf_map.slurm: 
            requires: "base_fxns.R", "sims.RDS", "output_ssfs/ssfs_X.rds", "output_rsfs/rsf_X.rds"
            returns: "output_rsf_maps/rsfm_X.rds"

        run_ssf_map.slurm: 
            requires: "base_fxns.R", "sims.RDS", "output_ssfs/ssfs_X.rds", "output_steps/used_X.rds"
            returns: "output_ssf_maps/ssfm_X.rds"

        run_ssd.slurm:
            requires: "base_fxns.R", "sims.RDS", "output_ssfs/ssfs_X.rds", "output_steps/used_X.rds"
            returns: "output_ssds/ssds_X.rds"


    Merging Output: 

        tidy_output.R:
            requires: "base_fxns.R", "output_abm/", "output_ssfs/", "output_rsfs/"
            returns: "output_rsf_maps/rsfm_X.rds"


Case studies ABMs:

    run_abm_case.slurm:
        requires: "base_fxns.R", "sims.RDS", "index.sheet.csv", "output_ssfs/ssfs_X.rds", "output_steps/used_x.rds"
        returns: "output_abm_case/"

    merge_case_abms.R:
        requires: "base_fxns.R", "output_abm_case/", "output_abm/", "abm_rasters/"
        returns: "abm_compiled", "store.convergence.case.RDS", "store.convergence.simulations.RDS"
