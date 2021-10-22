// uc davis farm config

singularity {
	enabled = true
}

docker {
	enabled = false
}

process {
    executor = 'slurm'
    // beforeScript = 'conda activate nextflow; module load singularity'
    clusterOptions = '--partition bmm'
}

params {
    max_memory = 500.GB
    max_cpus = 50
}


//    withLabel:process_low {
//      cpus = { check_max( 2 * task.attempt, 'cpus' ) }
//      memory = { check_max( 14.GB * task.attempt, 'memory' ) }
//      time = { check_max( 6.h * task.attempt, 'time' ) }
//    }
//    withLabel:process_medium {
//      cpus = { check_max( 6 * task.attempt, 'cpus' ) }
//      memory = { check_max( 42.GB * task.attempt, 'memory' ) }
//      time = { check_max( 8.h * task.attempt, 'time' ) }
//    }
//    withLabel:process_high {
//      cpus = { check_max( 12 * task.attempt, 'cpus' ) }
//      memory = { check_max( 84.GB * task.attempt, 'memory' ) }
//     time = { check_max( 10.h * task.attempt, 'time' ) }
//    }
//    withLabel:process_long {
//      cpus = { check_max(6 * task.attempt, 'cpus') }
//      memory = { check_max( 5.GB * task.attempt, 'memory' ) }
//     time = { check_max( 20.h * task.attempt, 'time' ) }
//  }

