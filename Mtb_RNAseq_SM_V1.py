# _______
# This is a snakemake file/script for processing of Single End Illumina RNAseq generated from Mtb.

### Maximillian Marin (mgmarin@g.harvard.edu)




# Define PATH to the reference genome to be used: 
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
refGenome_GFF_PATH = config["RefGenome_GFF_PATH"]


Kraken2_DB_PATH = config["Kraken2_DB_PATH"]



# Define PATH of OUTPUT directory
MainDir = config["output_dir"]

#output_Dir = MainDir + "sample_Output/"
output_Dir = MainDir  + "/"


import pandas as pd

# Read in data regarding input 
df = pd.read_csv( config["inputSampleData_TSV"], sep=',')


# Save a list of SRA/ENA "Run" Accessions
input_SampleIDs_SRA_RunAcc = list( df["Run"].values )



rule all:
    input:
        expand(output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}.fastq.gz", sampleRunAcc=input_SampleIDs_SRA_RunAcc),
        expand(output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FeatureCounts/{sampleRunAcc}.RNAseq.H37rv.ReadCounts.tsv", sampleRunAcc=input_SampleIDs_SRA_RunAcc),
        expand(output_Dir + "{sampleRunAcc}/IlluminaRNAseq/fastqc_Trimmomatic/{sampleRunAcc}_fq_fastqc.html", sampleRunAcc=input_SampleIDs_SRA_RunAcc),
        expand(output_Dir + "{sampleRunAcc}/IlluminaRNAseq/fastqc/{sampleRunAcc}_fq_fastqc.html", sampleRunAcc=input_SampleIDs_SRA_RunAcc),


rule download_FQ_FromSRA_RunID:
    output: 
        FQ_GZ = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}.fastq.gz",
    params:
        target_DownloaDir = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FASTQs/"
    conda:
        "envs/sratools_2_10_7_Conda.yml"
    shell:
        "fastq-dump --stdout {wildcards.sampleRunAcc} | gzip > {output.FQ_GZ} \n"



rule fastqc_BeforeTrim:
    input:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}.fastq.gz",
    output:
        html = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/fastqc/{sampleRunAcc}_fq_fastqc.html",
        zip = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/fastqc/{sampleRunAcc}_fq_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    wrapper:
        "0.38.0/bio/fastqc"




# adapter list from: https://github.com/stephenturner/adapters/blob/master/adapters_combined_256_unique.fasta
# can move adapter list to config file later
rule trimmomatic_Illumina_SE_Trimming:
    input:
         output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}.fastq.gz",
    output:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}_trimmed.fastq.gz",
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:'references/CustomTrimmoatic_IlluminaWGS_AdapterList.fasta':2:30:10:2:true SLIDINGWINDOW:4:20 MINLEN:40"],
        compression_level="-9"
        # optional parameters
        # extra=" "
    threads: 8
    wrapper:
        "0.38.0/bio/trimmomatic/se"




rule fastqc_OfTrimmomaticReads:
    input:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}_trimmed.fastq.gz",
    output:
        html = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/fastqc_Trimmomatic/{sampleRunAcc}_fq_fastqc.html",
        zip = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/fastqc_Trimmomatic/{sampleRunAcc}_fq_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    wrapper:
        "0.38.0/bio/fastqc"




rule bwa_map_IllSE_AlignTo_H37rv:
    input:
        refGenome_FA = refGenome_FA_PATH,
        fq_trimmed = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}_trimmed.fastq.gz",
    output:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.sam"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    params:
        rg=r"@RG\tID:{sampleRunAcc}\tSM:{sampleRunAcc}"
    threads: 8
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input.refGenome_FA} {input.fq_trimmed} > {output}"


rule samtools_ViewAndSort_IllPE_AlignTo_H37rv:
    input:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.sam"
    output:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.bam"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools view -bS {input} | samtools sort - > {output}"



rule samtools_index_IllSE_AlignTo_H37rv:
    input:
        output_Dir + "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.bam"
    output:
        output_Dir+ "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.bam.bai"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools index {input}"



rule featureCounts_H37rv_GFF:
    input:
        RNAseq_ToH37rv_BAM = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.bam",
        RNAseq_ToH37rv_BAI = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/IlluminaSE_AlignedTo_H37rv/{sampleRunAcc}.IllSE.RNAseq.H37rv.bam.bai", 
        refGenome_GFF = refGenome_GFF_PATH,
    output:
        ReadCountsMatrix = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FeatureCounts/{sampleRunAcc}.RNAseq.H37rv.ReadCounts.tsv",
        FeatureCounts_Summary = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FeatureCounts/{sampleRunAcc}.RNAseq.H37rv.ReadCounts.tsv.summary"
    conda:
        "envs/subread_featurecounts_2_0_1_Conda.yml"
    shell:
        "featureCounts -t gene -g locus_tag -a {input.refGenome_GFF} -o {output.ReadCountsMatrix} {input.RNAseq_ToH37rv_BAM}"


rule Kraken2_Illumina_SE:
    input:
        Kraken2_DB = Kraken2_DB_PATH,
        fq1_trimmed = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FASTQs/{sampleRunAcc}.fastq.gz",
    output:
        ReadClassifications = output_Dir + "{sampleID_WiIll}/IlluminaRNAseq/Kraken2/{sampleID_WiIll}.Kraken2.Output.Reads.txt",
        Kraken2_Report = output_Dir + "{sampleID_WiIll}/IlluminaRNAseq/Kraken2/{sampleID_WiIll}.Kraken2.Output.Report.txt"
    conda:
        "envs/kraken2_2_0_8_Conda.yml"
    threads: 1
    shell:
        " kraken2 --use-names --threads {threads} --db {input.Kraken2_DB} --output {output.ReadClassifications} "
        " --report {output.Kraken2_Report} {input.fq1_trimmed}   "


"""
rule Filter_FeatureCountsOutput_ReadCountsOnly:
    input:
        FeatureCounts_OutputFile = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FeatureCounts/{sampleRunAcc}.RNAseq.H37rv.ReadCounts.tsv",
    output:
        ReadCounts_Col7_Only_TSV = output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FeatureCounts/{sampleRunAcc}.RNAseq.H37rv.ReadCounts.Col7_Only.tsv",
    shell:
        "cat {input} | sed '1d' | cut -f 7 > {output}"


rule get_GeneIDs_FromFeatureCounts_Output:
    input:
        expand(output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FeatureCounts/{sampleRunAcc}.RNAseq.H37rv.ReadCounts.Col7_Only.tsv", sampleRunAcc=input_SampleIDs_SRA_RunAcc),
    output: 
        geneIDs_TXT = output_Dir + "FeatureCounts_MergedReadCounts/geneids.txt",
    shell:
        "ls -1 {input} | head -1 | xargs cut -f 1 > {output} "


rule merge_All_FeatureCounts_ToMatrix:
    input:
        geneIDs_TXT = output_Dir + "FeatureCounts_MergedReadCounts/geneids.txt",
        setOfAll_CountsFiles = expand(output_Dir + "{sampleRunAcc}/IlluminaRNAseq/FeatureCounts/{sampleRunAcc}.RNAseq.H37rv.ReadCounts.Col7_Only.tsv", sampleRunAcc=input_SampleIDs_SRA_RunAcc),
    output: 
        FeatureCounts_ReadsCountMatrix_Merged = output_Dir + "FeatureCounts_MergedReadCounts/FeatureCounts.ReadsCountMatrix.Merged.tsv",
    shell:
        "paste {input.geneIDs_TXT} {input.setOfAll_CountsFiles} > {output.FeatureCounts_ReadsCountMatrix_Merged}"
"""






