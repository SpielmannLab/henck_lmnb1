// Launching an srun (slurm job), and go into $SCRATCH
// Activate conda environment scVelocity
// Submit example: nextflow run sc_velo.nf -params-file sc_velo_params.yaml --id ${SCRATCH/"/scratch/"/}

// Use dsl 2 to separate process declaration and workflow
nextflow.enable.dsl=2

/*--------------MAIN WORKFLOW----------------*/

workflow {
    
    scobj_channel=Channel.fromList(params.in_seurat_rds)

    rdsToh5ad(scobj_channel)

    sc_velo(rdsToh5ad.out, params.subset_key, params.subset_val)

}

/*--------------PROCESSES----------------*/
process rdsToh5ad{ // conversion of seurat object to h5ad format required for the scVelo python script
    publishDir params.outfolder+'scVelo', mode: 'copy', overwrite: true
    input:
        path input
    output:
        path "*.h5ad"
    """
    rds_h5ad.R --scobject=${input}
    """
}

process sc_velo{
    errorStrategy 'ignore'
    publishDir params.outfolder+'scVelo', mode: 'copy', overwrite: true
    input:
        file scobj
        val subset_key
        each subset_val
    output:
        file "figures/*${params.plot_format}"
    """
    sc_velo.py ${scobj} ${params.plot_format} ${params.subset_key} ${params.subset_val}
    """
}
