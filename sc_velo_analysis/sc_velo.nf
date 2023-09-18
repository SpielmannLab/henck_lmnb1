// Launching an srun (slurm job), and go into $SCRATCH
// Activate conda environment scVelocity
// Submit example: nextflow run sc_velo.nf -params-file sc_velo_params.yaml --id ${SCRATCH/"/scratch/"/}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

/*--------------DEFINE GLOBAL PARAMETERS----------------*/
code_rdsToh5ad = "$SCRATCH/src/rds_h5ad.R"
code_sc_velo = "$SCRATCH/src/sc_velo.py"

/*--------------MAIN WORKFLOW----------------*/

workflow {
    write_params()
    
    scobj_channel=Channel.fromList(params.in_seurat_rds)

    rdsToh5ad(code_rdsToh5ad,
        scobj_channel)

    sc_velo(code_sc_velo,
        rdsToh5ad.out)

}

/*--------------PROCESSES----------------*/

process write_params{
    publishDir params.outfolder, mode: 'copy', overwrite: true
    output:
        path "*.txt"
    """
    echo \$(date) > parameters_sc_velo_${params.id}.txt
    echo "Parameters used to perform scVelo (RNA Velocity):" >> parameters_velo_${params.id}.txt
    echo ${params} | tr , '\n' >> parameters_velo_${params.id}.txt
    """
}

process rdsToh5ad{ // conversion of seurat object to h5ad format required for the scVelo python script
    publishDir params.outfolder+'scVelo', mode: 'copy', overwrite: true
    input:
        path "code.R"
        path input
    output:
        path "*.h5ad"
    """
    Rscript code.R --scobject=${input}
    """
}

process sc_velo{
    errorStrategy 'ignore'
    publishDir params.outfolder+'scVelo', mode: 'copy', overwrite: true
    input:
        path "code.py"
        file scobj
    output:
        file "figures/*${params.plot_format}"
    """
    python code.py ${scobj} ${params.plot_format} ${params.resolution}
    """
}
