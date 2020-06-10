# _______
# This is a snakemake file/script for processing of Single End Illumina RNAseq generated from Mtb.

### Maximillian Marin (mgmarin@g.harvard.edu)




# Define PATH to the reference genome to be used: 
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
refGenome_GFF_PATH = config["RefGenome_GFF_PATH"]



# Define PATH of OUTPUT directory
MainDir = config["output_dir"]

#output_Dir = MainDir + "sample_Output/"
output_Dir = MainDir   # + "/"


import pandas as pd

# Read in data regarding input 
df = pd.read_csv( config["inputSampleData_TSV"], sep='\t')


# Create a python list of all input samples
input_SampleIDs = list( df["SampleID"].values )

# Save a list of SRA/ENA "Run" Accessions
input_SampleIDs_SRA_RunAcc = list( df["Run"].values )



rule all:
    input:
        expand(f"{output_Dir}/{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}.fastq.gz", sampleID=input_SampleIDs_SRA_RunAcc),



rule download_FQ-FromSRA_RunID:
    output: 
        f"{output_Dir}/{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}.fastq.gz",
    params:
        outdir = f"{output_Dir}/{sampleRunAcc}/IlluminaRNAseq/FASTQs/"
    shell:
        "fastq-dump --sra-id {wildcards.sampleRunAcc} --outdir {params.outdir} --gzip "




rule fastqc_BeforeTrim:
    input:
        f"{output_Dir}/{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}.fastq.gz",
    output:
        html = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/fastqc/{sampleRunAcc}_fq_fastqc.html",
        zip = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/fastqc/{sampleRunAcc}_fq_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        output_Dir+"logs/fastqc/{sampleRunAcc}_fq.log"
    wrapper:
        "0.38.0/bio/fastqc"




# adapter list from: https://github.com/stephenturner/adapters/blob/master/adapters_combined_256_unique.fasta
# can move adapter list to config file later
rule trimmomatic_Illumina_SE_Trimming:
    input:
         f"{output_Dir}/{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}.fastq.gz",
    output:
        output_Dir + "/{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}_trimmed.fastq.gz",
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:'references/CustomTrimmoatic_IlluminaWGS_AdapterList.fasta':2:30:10:2:true SLIDINGWINDOW:4:20 MINLEN:75"],
        compression_level="-9"
        # optional parameters
        # extra=" "
    threads: 8
    wrapper:
        "0.38.0/bio/trimmomatic/se"




rule fastqc_OfTrimmomaticReads:
    input:
        output_Dir + "/{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}_trimmed.fastq.gz",
    output:
        html = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/fastqc_Trimmomatic/{sampleRunAcc}_fq_fastqc.html",
        zip = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/fastqc_Trimmomatic/{sampleRunAcc}_fq_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    wrapper:
        "0.38.0/bio/fastqc"




rule bwa_map_IllPE_AlignTo_H37rv:
    input:
        H37rv_FA_PATH = refGenome_FA_PATH,
        fq_trimmed = output_Dir + "/{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}_trimmed.fastq.gz",
    output:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.sam"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    params:
        rg=r"@RG\tID:{sampleRunAcc}\tSM:{sampleRunAcc}"
    threads: 8
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input.H37rv_FA_PATH} {input.fq_trimmed} > {output}"






rule samtools_ViewAndSort_IllPE_AlignTo_H37rv:
    input:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.sam"
    output:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.bam"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools view -bS {input} | samtools sort - > {output}"



rule samtools_index_IllPE_AlignTo_H37rv:
    input:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.bam"
    output:
        output_Dir+ "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllPE.H37rv.bam.bai"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools index {input}"



rule featureCounts_H37rv_GFF:
    input:
        ""
    output:
        ""
    conda:
        ""
    shell:
        "featureCounts -t gene -g locus_tag -a {refGTF} -o {output.ReadCountsMatrix} {input.RNAseq_To_H37rv_BAM}"












