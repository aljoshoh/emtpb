#! /usr/bin/env nextflow

println "\n-->EXPERIMENT DESCR $params.desc\n-->Executing $params.file in conda environment $params.conda with in total $params.amount jobs.\n"

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --file /path/to/python/file --conda name/of/built/conda/environment

        Mandatory arguments:
         --conda                       name of virtual conda environment
         --file                        script to exectute (full path required)
	 --desc			       description of experiment, also makes this directory       

	Optional arguments:
	--amount                       amount of submitted jobs
        --outdir                       Output directory to place final output
        --options                      Additional options [-evalue 1e-3]
        --path 			       Conda path 
	--node 			       Node for slurm execution [ibis216-010-071]
	--help                         This usage statement
        """
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

a = params.amount
counter_channel = Channel.of(1..a)



dir = file("./.runName_$workflow.runName")
dir.mkdir()   
file = file(params.file)
file.copyTo(dir)


// Convert jupyter notebook
process convert {
        input:
        file notebook from file(params.file)

        script:
        """
        jupyter nbconvert --to script "${notebook}"
        """
}


// Executing jobs
process myjob {

	conda params.condapath
	
	publishDir "$dir"
	
	input:
    	val counter from counter_channel
        file notebook from file(params.file)

	// output:
    	// file "test.txt" into
	
	script:
    	"""
	echo '\nExecuting process with file $params.file and conda environment /home/icb/$USER/miniconda3/envs/$params.conda ...\n'
	echo '\n Running script $params.file'
	echo '\n RUN #$counter'
	touch $counter
	python "${notebook/%.ipynb/.py}" $counter

	"""
}

//report = file('report.html')
//report.moveTo(dir)
//timeline = file('timeline.html')
//timeline.moveTo(dir)

