This is code to support: *Choices to landscapes: Mechanisms of animal movement scale to landscape patterns*

Authors:

Will Rogers will.rogers@yale.edu https://orcid.org/0000-0001-9899-8064

Scott Yanco yancos@si.edu https://orcid.org/0000-0003-4717-9370

Walter Jetz walter.jetz@yale.edu https://orcid.org/0000-0002-1971-7277

WR and SY conceptualized the paper. WR wrote the code and carried out the analyses supporting the manuscript. WR wrote the initial draft with contributions from SY. All authors provided feedback on the concept and contributed to the final draft.

Please write with any questions, concerns, or inquiries to will.rogers@yale.edu, rogerswill47@gmail.com

All movement data are publically available and spatial data was collected from public sources, or are otherwise noted (see S2 of paper). All data is provided as R data objects (or CSV for density data) for ease, but the raw movement data (besides that contained in amt) are available on MoveBank as described in S2. A substantial amount of computational effort was performed on Yale HPC, and the analyses for ABMs would not be possible using local computational resources (e.g. 16GB RAM computer with 8 processors). We provide the scripts necessary to implement this pipeline on the HPC, but please contact if there are issues in implementation. Note: simulation SSF and RSF models and raw ABM simulations are not possible to provide here because of their size and Github limits. The precursors to and products from these objects are provided (we document below where these objects are missing in the analysis pipeline). Also, note that the African buffalo census data is sourced from SANParks Kruger National Parks Scientific Services. Reuse of this data for another purpose is prohibited, and inqueries about further use of this data should be directed to SANParks Kruger National Parks Scientific Services. In accordance with data sharing polices, we randomized the year associated with the spatial density data, as our results do not rely on the temporal structure of density but this allows us to acocunt for between-year density variation. 


To reproduce the results of the paper, you should only have to run the code within "Paper_Compilation.RMD" - this file will generate precursors for HPC simulations, and it pulls in HPC outputs as specified below. To run this code in RStudio, make sure the file can properly access the "Simulation", "Casestudies", and "Code" subfolders. The code within this file should run with 10-15 minutes on a local PC with 16GB of RAM and 8 processors. We also provide the compiled document as a HTML file (S3 in the text) for clarity.  

Plots will be generated in the "Plots" subfolder - note: the papers appearing in the paper were compiled from these figure components. So, you won't see Figure 3 as it appears in the text for example, but you will see all the subcomponents that were compiled to make Figure 3. 

If you wish to recreate the cluster-based simulations, you need to load `base_fxns.R`, `sur.val.packed.RDS`, `sims.RDS`, `index.sheet.csv`, and `tidy_output.R` into cluster storage. Then, we recommend working through the slurm calls in this order - they should generate their own outdirectories locally, but this can be very tricky.

    1. `run_steps.slurm` - generates the 10,000 steps used for modeling based on the 500 simulation contexts
    2. `run_rsf.slurm` - generates RSF models based on run_steps.slurm output - too large provide
    3. `run_ssf.slurm` - generates SSF models based on run_steps.slurm output - too large provide
    4. `run_simulations.slurm` - generates population sampling 5 million steps per simulation context
    5. `run_abm.slurm` - generates ABM simulation based on run_steps.slurm and run_ssf.slurm output - too large provide
    6. `run_abm_size.slurm` - generates ABM simulations for differently sized landscapes based on output from Paper_Compilation.RMD "abm.example.RDS"
    7. `run_abm_example.slurm` - generates ABM simulations for differently sized landscapes based on output from Paper_Compilation.RMD "abm.example.RDS"
    8. `run_rsf_map.slurm` - generates RSF maps based on run_rsf.slurm and run_ssf.slurm output
    9. `run_ssf_map.slurm` - generates naive SSF maps based on run_steps.slurm and run_ssf.slurm output
    10. `run_ssd.slurm` - generates SSD maps based on run_steps.slurm and run_ssf.slurm output
    11. `tidy_output` - ingests outputs of run_abm.slurm, run_rsf.slurm, and run_ssf.slurm to make smaller (storage) outputs
        
    
If you wish to recreate the cluster-based ABM maps for the case studies, you need to load "base_fxns.R", "reddeer.abm.RDS", "roedeer.abm.RDS", "fisher.focal.abm.RDS", "fisher.nonfocal.abm.RDS", "buffalo.abm.RDS", and "merge_case_abms.R" into cluster storage. To make the best use of the cluster structure at Yale, we only run a single individual in each call of "run_abm_case.slurm". This is tedious, and means that you must manually change the `OUTDIR="output_abm_case/X/Y"`, `total <- readRDS("X.abm.RDS")`, `model <- total[[2]][[Y]]` for each simulation case (X) and for each individual (Y) in each simulation case. 

    1. run_abm_case.slurm - generates the ABM simulation steps for each case study - too large provide
    2. merge_case_abms.slurm - ingests outputs of run_abm_case.slurm to make smaller (storage) outputs
    


The buffalo density data "buffalo_locations.csv" has four columns: 

    TOTAL: total number of mixed-herd buffalo detected at a given X,Y location during dry-season helicopter surveys during a given year, t.rand
    X: the Easting (in m) of the detection
    Y: the Northing (in m) of the detection
    t.rand: the year of the detection - year order is randomized for data sharing

`base_fxns.R` is a critical file for any kind of cluster work, and is the library of all functions and helper functions used in the paper. Most of these are copied directly from "Paper_Compilation.RMD", but some are specific to running the ABM. Some of the functions used in ABM simulations are adaptations from code included in `amt`. 

Below we outline the function of the files contained in the Code subfolder:

`Paper_Compilation.RMD`:
    
        Simulations: 
            requires: "output_abm_example", "output_ssds", "output_ssf_maps", "output_rsf_maps", "output_abm_sim", "output_rasters", "store.convergence.simulations.RDS", "ssf.coefs.rds", "rsf.coefs.rds" 
            requires, but generates: "pop.locations.demo.RDS", "ssd.by.neighborhood.RDS", "ram.RDS", "nnzero.elements.RDS", "locs.and.rasters.RDS", "abm_size_bhs.RDS"
            generates: "sims.RDS", "abm.example.RDS", "sur.val.packed.RDS"          
        
        Red Deer case study:
            requires: amt package data, "sim.pred/final_red_1.tif"
            requires, but generates: "deerrasters.RDS"
            generates: "reddeer.abm.RDS" ***Missing from files given size***

        Roe Deer case study:
            requires: "EuroDeer_ Roe deer in Italy 2005-2008.csv" (MoveBank), Landcover Data - See S2, "sim.pred/final_roe_X.tif"
            requires, but generates: "italy.deer.RDS", "italyrasters.RDS"
            generates: "roedeer.abm.RDS" ***Missing from files given size***
            
        Fisher case study:
            requires: "Martes pennanti LaPoint New York.csv" (MoveBank), Spatial Data - See S2, "sim.pred/final_ff_X.tif", "sim.pred/final_fnf_X.tif"
            requires, but generates: "fisher.data.RDS", "Focal.Landscape.RDS", "Nonfocal.Landscape.RDS"
            generates: "fisher.focal.abm.RDS", "fisher.nonfocal.abm.RDS" ***Missing from files given size***
        
        Buffalo case study:
            requires: "Kruger African Buffalo, GPS tracking, South Africa.csv" (MoveBank), "2017_kruger_national_park.shp", Modis NDVI Range, Streams, Dams and Water Holes, "sim.pred/final_sum_bX.tif", "buffalocensusR.RDS"
            requires, but generates: "buffalo.data.RDS", "buffalo.landscape.RDS", buffalo.sparse.neighbors.RDS
            generates: "buffalo.abm.RDS" ***Missing from files given size***


Simultations:

    Simulate movements:

        run_steps.slurm:
            requires: "base_fxns.R", "sur.val.packed.RDS"
            returns: "output_steps/used_X.rds"


    Fit RSF and SSF models:

        run_rsf.slurm: 
            requires: "base_fxns.R", "sims.RDS", "output_steps/used_X.rds"
            returns: "output_rsfs/rsfs_X.rds" ***Missing from files given size***

        run_ssf.slurm: 
            requires: "base_fxns.R", "sims.RDS", "output_steps/used_X.rds"
            returns: "output_ssfs/ssfs_X.rds" ***Missing from files given size***


    Run simulations for validation:

        run_simulations.slurm:
            requires: "base_fxns.R", "sur.val.packed.RDS"
            returns: "output_rasters/used_X.rds"


    Run ABMs: 

        run_abm.slurm:
            requires: "base_fxns.R", "sims.RDS", "index.sheet.csv", "output_ssfs/ssfs_X.rds", "output_steps/used_x.rds"
            returns: "output_abms/abms_1_X.rds" ***Missing from files given size***

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
            returns: "output_abm_sim", "ssf.coefs.RDS", "rsf.coefs.RDS"


Case studies ABMs:

    run_abm_case.slurm:
        requires: "base_fxns.R", "reddeer.abm.RDS", "roedeer.abm.RDS", "fisher.focal.abm.RDS", "fisher.nonfocal.abm.RDS", "buffalo.abm.RDS" ***Case Study RDS files are missing from files given size***
        returns: "output_abm_case/" *Missing from files given size

    merge_case_abms.R:
        requires: "base_fxns.R", "output_abm_case/"
        returns: "abm_compiled", "store.convergence.case.RDS", "store.convergence.simulations.RDS"
        
    
