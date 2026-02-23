# COMP 483 Pipeline Project
This pipeline was developed to analyze transcriptomic data of Human cytomegalovirus from two participants at 2 and 6 days post-infection. It integrates output from a variety of tools, including kallisto, sleuth, Bowtie2, SPAdes, and NCBI BLAST/datasets. Links for information on how to download these dependencies are included under "Required Tools". 

The "sample_data" directory in this repo includes short sample data from the 4 transcriptomic datasets to test if the pipeline is working. Information on how this sample data was downloaded as well as how to run the pipline is included in the sections "Data" and "Running on Sample Data", respectively.

### Required Tools
* kallisto
* sleuth (R package)
* Bowtie2
* SPAdes
* NCBI BLAST/datasets 
* fasterq-dump or wget

### Data
Sample data is the first 10,000 reads from the following 4 SRA samples: SRR5660030, SRR5660033,SRR5660044, and SRR5660045. The full set of reads was downloaded from NCBI using the following command for each sample:

    fasterq-dump \<accession\>

These reads can also be downloaded using wget:

    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/<accession>/<accession>

These reads were used to generate the Finger_PipelineReport.txt file included in the repo.

To generate the shorter sample input data files, the first 10,000 reads from each sample were extracted using the following command:

    head -n 40000 \<sample\>.fastq \> shortened_\<sample\>.fastq

### Running on Sample Data
The file "Snakefile" and all other scripts are correctly formatted to run on the provided sample data with one command. First, clone the repo using:

    git clone https://github.com/fingerga/pipeline_project.git

Then, cd into the pipeline_project directory:

    cd pipeline_project

Then, run:

    snakemake -c 1

The pipeline should finish running on the sample data within 2 minutes. You can optionally increase the number of cores used to speed up the process.

