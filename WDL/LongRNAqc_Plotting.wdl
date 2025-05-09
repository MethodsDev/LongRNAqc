version 1.0

task LongRNAqcPlottingTask {
    input {
        Array[String] sampleName
        Array[File] classificationFile
        String outputPrefix
        String type
        Int preemptible
        String docker
    }

    # Calculate total memory required
    Int total_file_size = ceil(size(classificationFile, "GiB") + 8)

    command {

           set -ex
        
           LongRNAqc_classification_plots.py \
            --sample_names '~{sep="," sampleName}' \
            --classification_files '~{sep="," classificationFile}' \
            --output ~{outputPrefix} \
            --type ~{type}

       gzip ~{outputPrefix}.read_lengths.tsv
    }

    output {
        File QC_categories_plots = "~{outputPrefix}.categories.pdf"
        File QC_read_lengths_plots = "~{outputPrefix}.read_length_distributions.pdf"
        File? QC_categories_tsv = "~{outputPrefix}.categories.tsv"
        File QC_read_lengths_tsv = "~{outputPrefix}.read_lengths.tsv.gz"
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 30
        disks: "local-disk " + total_file_size*2 + " HDD"
        cpu: 1
        memory: "4 GiB"
        preemptible: preemptible
    }
}


workflow LongRNAqcPlotting {
    meta {
        description: "Generate multi-sample QC plots using Sqanti3 outputs."
    }

    input {
        Array[String] sampleName
        Array[File] classificationFile
        String outputPrefix
        String type
        Int preemptible = 1
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-sqanti3/lrtools-sqanti3-plotting"
    }

    call LongRNAqcPlottingTask {
        input:
            sampleName = sampleName,
            classificationFile = classificationFile,
            outputPrefix = outputPrefix,
            type = type,
            preemptible = preemptible,
            docker=docker
    }

    output {
        File QC_categories_plots = LongRNAqcPlottingTask.QC_categories_plots
        File QC_read_lengths_plots = LongRNAqcPlottingTask.QC_read_lengths_plots
        File? QC_categories_tsv = LongRNAqcPlottingTask.QC_categories_tsv
        File QC_read_lengths_tsv = LongRNAqcPlottingTask.QC_read_lengths_tsv
        
    }
}

