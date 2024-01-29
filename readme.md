# Intro
The scripts here have been used for the following analysis:

1. Perform cell composition analysis. In the current implementation, a cell composition was manually created by Jana based on the annotated dataset. This was then analyused using the script <cell_comp_analysis_subcluster.R> and <cell_comp_analysis_main_cluster.R>. The analysis contains both statistics test based on beta-binomial distribution implemented either by the <countdata> R-package or as implemented in the MMCA paper (by CX) using <VGAM> R-package using vector generalized linear models. Just run the script like on the local computer: 

        Rscript cell_comp_analysis_xxx.R

2. scVelo analysis. In the current implementation, since Jana did not use the RDS object from count_to_seurat pipeline that contained spliced and unpsliced counts, the script <sc_velo_analysis/add_spliced_count_to_rds_for_velocity.R> has to be run to manaually add this info. After this, run the custom <sc_velo.nf> script to make the scVelo-based RNA velocity plots. Do this in the cluster (Note, the internal script was changed on 19th September and the nextflow automation may not work - needs fixing.):

        srun -p shortterm -c 2 --mem 100GB debug --pty bash
        conda activate scVelocity
        module load nextflow/v22.04.1

        # Move to SCRATCH all the relevant scripts.
        cp * $SCRATCH/
        cp -r ../src $SCRATCH/
        cd $SCRATCH

        nextflow run sc_velo.nf -params-file sc_velo_params.yaml --id ${SCRATCH/"/scratch/"/}



# Installations
    conda env create --file henck_lmnb1.yaml 
    # install CRAN package <countdata>
    conda activate henck_lmnb1
    R
    install.packages("countdata")
