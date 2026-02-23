# COMP 483 Pipeline Project
This pipeline was developed to analyze transcriptomic data of Human cytomegalovirus from two participants at 2 and 6 days post-infection. It integrates output from a variety of tools, including kallisto, sleuth, Bowtie2, SPAdes, and NCBI BLAST/datasets. Links for information on how to download these dependencies are included under "Required Tools". 

The "sample_data" directory in this repo includes short sample data from the 4 transcriptomic datasets to test if the pipeline is working. Information on how this sample data was downloaded as well as how to run the pipeline is included in the sections "Data" and "Running on Sample Data", respectively.

### Required Tools
* Python version 3.12.3 
* R version 4.5.2 
* kallisto: https://pachterlab.github.io/kallisto/download.html
* sleuth (R package) https://pachterlab.github.io/sleuth/download
* Bowtie2 https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2
* SPAdes https://github.com/ablab/spades/releases/tag/v4.2.0
* NCBI BLAST+ https://www.ncbi.nlm.nih.gov/books/NBK52640/
* NCBI datasets https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
* fasterq-dump through SRA Toolkit https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
    * alternatively, wget https://www.gnu.org/software/wget/manual/

## Data
### Full Datasets
Sample data is the first 10,000 reads from the following 4 SRA samples: 
* SRR5660030 (Donor 1, 2dpi) 
* SRR5660033 (Donor 1, 6dpi)
* SRR5660044 (Donor 3, 2dpi)
* SRR5660045 (Donor 3, 6dpi)

The full set of reads was downloaded from NCBI using the following command for each sample:

    fasterq-dump <accession>

These reads can also be downloaded using wget:

    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/<accession>/<accession>
These commands were repeated for each sample, changing \<accession\> out for every SRA accession number.

The full set of these reads were used to generate the Finger_PipelineReport.txt file included in the repo.

### Generating Sample Data
To generate the shorter sample input data files, the first 10,000 reads from each sample were extracted using the following command:

    head -n 40000 \<sample\>.fastq \> shortened_\<sample\>.fastq
This command was repeated for each sample, changing \<sample\> out for every SRA accession number.
## Running on Sample Data
The file "Snakefile" and all other scripts are correctly formatted to run on the provided sample data with one command. 
### Clone repository
First, cd into your desired directory and clone the repo using:

    git clone https://github.com/fingerga/pipeline_project.git

### Running pipeline
Then, cd into the pipeline_project directory:

    cd pipeline_project

Then, run:

    snakemake -c 1

The pipeline should finish running on the sample data within 2 minutes. You can optionally increase the number of cores used to speed up the process.