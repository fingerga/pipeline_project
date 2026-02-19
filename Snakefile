index_accession = "GCF_000845245.1"
samples = ["SRR5660030","SRR5660033","SRR5660044","SRR5660045"]

rule all:
    input:
        outfile_report= "PipelineReport.txt",
        sleuth_done = "sleuth_check.txt",
        mapped1 = expand("mapped_reads/{sample}/{sample}.1.fastq", sample= samples),
        mapped2 = expand("mapped_reads/{sample}/{sample}.2.fastq", sample= samples)

rule fetch_index:
    output:
        outfile_zip= "ncbi_dataset.zip"
    shell:
        "datasets download genome accession {index_accession} --include cds --filename {output.outfile_zip}"

rule unzip_index:
    input:
        outfile_zip= "ncbi_dataset.zip"
    output:
        cds= "data/cds_from_genomic.fna"
    shell:
        '''
        unzip -p {input.outfile_zip} > {output.cds}
        '''

rule clean_cds:
    input:
        cds = "data/cds_from_genomic.fna"
    output:
        clean_cds= "data/cds_clean.fna",
        outfile_report= "PipelineReport.txt"
    shell:
        '''
        sed '/^>/s/ .*//' {input.cds} | sed 's/>.*cds_/>/g' > data/cds_clean.fna
        echo "The HCMV genome (GCF_000845245.1) has $(wc -l {output.clean_cds} | awk '{{print $1}}') CDs" >> {output.outfile_report}
        '''

rule build_index:
    input:
        clean_cds= "data/cds_clean.fna"
    output:
        ref_index= "index.idx"
    shell:
        '''
        kallisto index -i index.idx {input.clean_cds}
        '''

rule kallisto_on_reads:
    input: 
        r1="sample_data/{sample}/{sample}_1.fastq",
        r2="sample_data/{sample}/{sample}_2.fastq",
        ref_index= "index.idx"
    output:
        quant_reads = directory("quant_reads/{sample}")
    shell:
        '''
        kallisto quant -i {input.ref_index} -o {output.quant_reads} -b 10 -t 2 {input.r1} {input.r2}
        '''

rule sleuth:
    input:
        quant_reads = expand("quant_reads/{sample}",sample=samples),
    output:
        sleuth_done = "sleuth_check.txt"
    shell:
        '''
        Rscript sleuth_script.R 
        touch {output.sleuth_done}
        cat sleuth_out.txt >> PipelineReport.txt
        '''

rule cleanup:
    shell:
        '''
        rm -r quant_reads
        rm -r data
        rm index.idx
        rm ncbi_dataset.zip
        rm PipelineReport.txt
        rm sleuth_out.txt
        rm sleuth_check.txt
        rm -r mapped_reads
        rm ncbi_dataset_gen.zip
        rm bowtiebuild.done
        rm -r ref_gen
        rm dummy.sam
        '''

rule fetch_index_gen:
    output:
        outfile_zip_gen= "ncbi_dataset_gen.zip"
    shell:
        "datasets download genome accession {index_accession} --include genome --filename {output.outfile_zip_gen}"

rule unzip_index_gen:
    input:
        outfile_zip_gen= "ncbi_dataset_gen.zip"
    output:
        ref_genome= "data/HCMVgenome.fasta"
    shell:
        '''
        unzip -p {input.outfile_zip_gen}  > {output.ref_genome}
        sed -n '/^>/,$p' {output.ref_genome} | sed '/^{{/,$d' > genometemp.fasta
        cat genometemp.fasta > {output.ref_genome}
        rm genometemp.fasta
        '''

rule bowtie_build:
    input:
        ref_genome= "data/HCMVgenome.fasta"
    output: 
        ref_index_gen= touch("bowtiebuild.done")
    shell:
        '''
        mkdir ref_gen
        bowtie2-build {input.ref_genome} ref_gen/ref
        '''

rule bowtie_run:
    input:
        r1="sample_data/{sample}/{sample}_1.fastq",
        r2="sample_data/{sample}/{sample}_2.fastq",
        ref_index_gen= touch("bowtiebuild.done")
    output:
        mapped1 = "mapped_reads/{sample}/{sample}.1.fastq",
        mapped2 = "mapped_reads/{sample}/{sample}.2.fastq"
    shell:
        '''
        bowtie2 -x ref_gen/ref -1 {input.r1} -2 {input.r2} --al-conc mapped_reads/{wildcards.sample}/{wildcards.sample}.%.fastq -S dummy.sam
        '''
