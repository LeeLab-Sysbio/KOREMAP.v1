######################################################################
# Define paths to utility directories and tools
######################################################################
util=/data/utils
qc_dir=${util}/FastQC
tg_dir=${util}/TrimGalore-0.6.6
star_dir=${util}/STAR-2.7.3a/bin/Linux_x86_64
rsem_dir=${util}/RSEM-1.3.3

######################################################################
# STAR and RSEM INDEX generation
######################################################################
# Define reference directory and move to it
ref=/data/sequencing_data/ref/homo_sapiens
cd $ref

# Download the human genome FASTA file and decompress it
wget https://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download the human genome GTF file and decompress it
wget https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
gzip -d Homo_sapiens.GRCh38.84.gtf.gz &

# Define paths to genome FASTA and GTF files
gf=${ref}/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=${ref}/Homo_sapiens.GRCh38.84.gtf
rindex=${ref}/rsem
sindex=${ref}/star99

# Create directories for STAR and RSEM indexes
mkdir -p $sindex
mkdir -p $rindex

# Generate STAR genome index
${star_dir}/STAR --runThreadN 4 \
--runMode genomeGenerate  \
--genomeDir $sindex \
--genomeFastaFiles $gf \
--sjdbGTFfile $gtf \
--sjdbOverhang 99 

# Prepare RSEM reference
${rsem_dir}/rsem-prepare-reference -p 4 --gtf $gtf $gf ${rindex}/RSEM

######################################################################
# Define directories for main workflow
######################################################################
main=/data/RNAseq/
fq=${main}/fastq        # Directory for paired FASTQ files (R1.fq.gz and R2.fq.gz)
qout=${main}/QC         # Directory for quality control output files
tout=${main}/trimmed    # Directory for trimmed output files
bout=${main}/bam        # Directory for BAM files
eout=${main}/expression # Directory for expression data
mkdir -p ${qout}
mkdir -p ${bout}
mkdir -p ${eout}
mkdir -p ${tout}

######################################################################
# List unique file names in the fastq directory
######################################################################
ls -a $fq | cut -d '.' -f 1 | sort -u > $main/fname.txt

######################################################################
# Process each file: trimming, alignment, and expression quantification
######################################################################
while read name
do
    echo $name

    # Trim adapters and low-quality bases from reads
    ${tg_dir}/trim_galore \
    --adapter AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
    --adapter2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
    --phred33 --quality 20 --length 50 --paired \
    --cores 4 --gzip --output_dir ${tout} ${fq}/${name}.R1.fq.gz ${fq}/${name}.R2.fq.gz

    echo $name

    # Align trimmed reads to the reference genome using STAR
    ${star_dir}/STAR \
    --genomeDir ${sindex} \
    --readFilesIn ${tout}/${name}.R1_val_1.fq.gz ${tout}/${name}.R2_val_2.fq.gz \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM \
    --twopassMode Basic \
    --runThreadN 4 \
    --sjdbGTFfile ${gtf} \
    --genomeLoad NoSharedMemory \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --outSAMtype BAM Unsorted \
    --outSAMunmapped None \
    --outFileNamePrefix ${bout}/${name}.STAR.

    # Quantify gene expression using RSEM
    ${rsem_dir}/rsem-calculate-expression --time --num-threads 4 --alignments --paired-end --no-bam-output \
    ${bout}/${name}.STAR.Aligned.toTranscriptome.out.bam ${rindex}/index ${eout}/${name}.STAR.RSEM
    
done < $main/fname.txt
