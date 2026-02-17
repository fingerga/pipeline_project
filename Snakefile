index_accession = "GCF_000845245.1"
samples = ["SRR5660030","SRR5660033","SRR5660044","SRR5660045"]

rule all:
    input:
        quant_reads = expand("quant_reads/{sample}",sample=samples),
        outfile_report= "PipelineReport.txt"

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
        sed '/^>/s/ .*//' data/cds_from_genomic.fna | sed 's/>.*cds_/>/g' > data/cds_clean.fna
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
        r1="reads/{sample}/{sample}_1.fastq",
        r2="reads/{sample}/{sample}_2.fastq",
        ref_index= "index.idx"
    output:
        quant_reads = directory("quant_reads/{sample}")
    shell:
        '''
        kallisto quant -i index.idx -o {output.quant_reads} -b 10 -t 2 {input.r1} {input.r2}
        '''

rule cleanup:
    shell:
        '''
        rm -r quant_reads
        rm -r data
        rm index.idx
        rm ncbi_dataset.zip
        rm PipelineReport.txt
        '''
    
