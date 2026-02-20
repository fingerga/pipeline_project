index_accession = "GCF_000845245.1"
samples = ["SRR5660030","SRR5660033","SRR5660044","SRR5660045"]

rule all:
    input:
        outfile_report= "PipelineReport.txt",
        sleuth_done = "sleuth_check.txt",
        donecheck= expand("pipeline{sample}.done", sample= samples)

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
        rm ncbi_dataset.zip
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
        printf "\n" >> {output.outfile_report}
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
        printf "\n" >> PipelineReport.txt
        rm sleuth_out.txt
        '''

rule cleanup:
    shell:
        '''
        rm -r quant_reads
        rm -r data
        rm index.idx
        rm PipelineReport.txt
        rm -r mapped_reads
        rm -r ref_gen
        rm -r assemblies
        rm -r blast_db
        rm bowtiebuild.done
        rm sleuth_check.txt
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
        rm ncbi_dataset_gen.zip
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
        ref_index_gen= touch("bowtiebuild.done"),
        sleuth_done = "sleuth_check.txt"
    output:
        mapped1 = "mapped_reads/{sample}/{sample}.1.fastq",
        mapped2 = "mapped_reads/{sample}/{sample}.2.fastq",
        bowtie_map= "mapped_reads/{sample}/{sample}_out.txt",
    shell:
        '''
        bowtie2 --quiet -x ref_gen/ref -1 {input.r1} -2 {input.r2} --al-conc mapped_reads/{wildcards.sample}/{wildcards.sample}.%.fastq -S dummy.sam
        echo "Sample {wildcards.sample} had $(echo $(wc -l < {input.r1}) / 4 | bc) read pairs before and $(echo $(wc -l < {output.mapped1}) / 4 | bc) read pairs after Bowtie2 filtering." >> {output.bowtie_map}
        rm dummy.sam
        cat {output.bowtie_map} >> PipelineReport.txt
        '''

rule spades:
    input:
        mapped1 = "mapped_reads/{sample}/{sample}.1.fastq",
        mapped2 = "mapped_reads/{sample}/{sample}.2.fastq",
    output:
        assembled_gen= directory("assemblies/{sample}_assembly/"),
    shell:
        '''
        spades.py -k 127 -t 2 --only-assembler -1 {input.mapped1} -2 {input.mapped2} -o {output.assembled_gen}
        '''

rule fetch_blast:
    output:
        outfile_blast= "ncbi_dataset_blast.zip"
    shell:
        "datasets download virus genome taxon betaherpesvirinae --refseq --include genome --filename {output.outfile_blast}"

rule unzip_blast:
    input:
        outfile_blast= "ncbi_dataset_blast.zip"
    output:
        blast_genome= "data/blastdb.fasta"
    shell:
        '''
        unzip -p {input.outfile_blast}  > {output.blast_genome}
        sed -n '/^>/,$p' {output.blast_genome} | sed '/^The NCBI Datasets Project/,$d' > tempdb.fasta
        cat tempdb.fasta > {output.blast_genome}
        rm tempdb.fasta
        rm ncbi_dataset_blast.zip
        '''

rule build_blast:
    input:
        blast_genome= "data/blastdb.fasta"
    output:
        blast_db = directory("blast_db/")
    shell:
        '''
        makeblastdb -in {input.blast_genome} -out {output.blast_db}/betaherpesvirinae -title betaherpesvirinae -dbtype nucl
        '''
rule longest_contig:
    input:
        assembled_gen= directory("assemblies/{sample}_assembly")
    output:
        longest_contig = ("assemblies/{sample}_lc.fasta")
    shell:
        '''
        python blast_assembly.py -i {input.assembled_gen}/scaffolds.fasta -o {output.longest_contig}
        '''

rule run_blast:
    input:
        longest_contig = ("assemblies/{sample}_lc.fasta"),
        blast_db = directory("blast_db/")
    output:
        blast_out= "blast_out/{sample}"
    shell:
        '''
        blastn -query {input.longest_contig} -db {input.blast_db}/betaherpesvirinae -max_hsps 1 -out {output.blast_out} -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"
        '''


rule add_blast:
    input: 
        blast_out= "blast_out/{sample}"
    output:
        donecheck= touch('pipeline{sample}.done')
    shell:
        '''
        printf "\n" >> PipelineReport.txt
        echo "{wildcards.sample}:" >> PipelineReport.txt
        echo "sacc  pident  length  qstart  qend    sstart  send    bitscore    evalue  stitle" >> PipelineReport.txt
        head -n 5 {input.blast_out} >> PipelineReport.txt
        '''

#to do:
    #rule blast: call py file that extracts longest contig, using SeqIO probably, and uses a os.system call to blast, immediately >> to report file
    #rule remove_excess: remove extra files that aren't needed, keep: data, assemblies, mapped_reads, quant_reads
    #write comments
    #write README file/documentation
    #run with full sample files and SAVE REPORT