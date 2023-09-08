# Create sample sheet, run aquamis and validate qc values


rule create_sample_sheet:
    output:
        outdir=directory("aquamis"),
        sample_sheet="aquamis/samples.tsv",
    params:
        fastq_folder=config["fastq_folder"],
        fastq_naming=config["fastq_naming"],
        aquamis_path=config["aquamis_path"],
    conda:
        "../envs/aquamis.yaml"
    log:
        "logs/create_sample_sheet.log",
    shell:
        """
        exec 2> {log}
        echo 
        bash {params.aquamis_path}/scripts/create_sampleSheet.sh \
            --mode {params.fastq_naming} \
            --fastxDir {params.fastq_folder} \
            --outDir {output.outdir} \
            --force
        """


# rule run_aquamis:


# rule autoqc:
