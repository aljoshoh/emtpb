#! /usr/bin/env nextflow

println "\n-->EXPERIMENT DESCR $params.desc\n-->Executing $params.file in conda environment $params.conda with in total $params.amount jobs.\n"

nextflow.enable.dsl=2

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
process nf_emtpb_convert_and_unpack {
        
	input:
	val notebook
        // file notebook from file(params.file)
	
	output:
	val ""

        script:
        """
        jupyter nbconvert --to script '${notebook}'
	
	mkdir -p /localscratch/$USER/container;
 	cd /localscratch/$USER/container;

 	#Check if directory exists
 	if [[ ! -d "/localscratch/$USER/containers/emtpb" ]]
 	then
 		ch-tar2dir /lustre/groups/cbm01/code/alexander.ohnmacht/emtpb/metadata/emtpb.tar.gz /localscratch/$USER/container
 		mkdir -p /localscratch/$USER/container/emtpb/lustre/groups
 	fi
	mkdir -p /localscratch/$USER/container/emtpb/ch
	touch /localscratch/$USER/container/emtpb/ch/environment
        """
}


// Executing jobs
process nf_emtpb_benchmark {

	conda params.condapath
	
	publishDir "$dir"
	
	input:
    	val counter
	val notebook
	val next
	// val counter from counter_channel
        //file notebook from file(params.file)

	// output:
    	// file "test.txt" into
	
	script:
    	"""
	echo '\nExecuting process with file $params.file and conda environment /home/icb/$USER/miniconda3/envs/$params.conda ...\n'
	echo '\n Running script $params.file'
	echo '\n RUN #$counter'
	touch $counter
	python3 '/lustre/groups/cbm01/code/alexander.ohnmacht/emtpb/scripts/${notebook.baseName}.py' $counter

	"""
}


workflow {
	next = ""	

	notebook = file(params.file)
	next = nf_emtpb_convert_and_unpack(notebook)
 
	counter = counter_channel
	notebook = file(params.file)
	next = nf_emtpb_benchmark(counter,notebook,next)
}



//report = file('report.html')
//report.moveTo(dir)
//timeline = file('timeline.html')
//timeline.moveTo(dir)

