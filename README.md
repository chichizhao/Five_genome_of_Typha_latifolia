

# Scripts were used to build the five genome of ***Typha latifolia***

# 1. Data source 
| Name     | TGS Accession | TGS Platform      | TGS Bases | NGS Accession | NGS Bases | Publication                  |
|----------|---------------|-------------------|-----------|---------------|-----------|------------------------------|
| CN_SC    | CRR784823     | OXFORD NANOPORE   | 36.6G     | CRR784824     | 15.99G    | GWHCBIL00000000                           |
| CA_ON    | SRR15383431   | PacBio SMRT       | 86.8G     | SRR15383432   | 138.6G    | (Widanagama et al., 2022)     |
| USA_PA1  | SRR15691249   | OXFORD NANOPORE   | 27.6G     | SRR15691252   | 18.6G     | (Aylward et al., 2023)        |
| USA_PA2  | SRR15691248   | OXFORD NANOPORE   | 8.4G      | SRR15680755   | 22.4G     | (Aylward et al., 2023)                           |
| USA_CA   | SRR15691261   | OXFORD NANOPORE   | 10.5G     | SRR15680754   | 21.7G     | (Aylward et al., 2023)                           |


# 2. Filt the third generation reads (TGS) and the next generation reads (NGS)
    Nanofilt: https://github.com/wdecoster/nanofilt

        gunzip -c TGS.fastq.gz | NanoFilt -q 10 | gzip > highQuality-reads.fastq.gz

    fastp: https://github.com/OpenGene/fastp

        fastp -i NGS.R1.fq.gz -I NGS.R2.fq.gz -o NGS_clean.R1.fq.gz -O NGS_clean.R2.fq.gz
    
# 3. Evaluate the genome size
    Jellyfish: https://github.com/gmarcais/Jellyfish
    GenomeScope: http://genomescope.org/genomescope2.0/

        jellyfish count -C -m 21 -s 1000000000 -t 10 *.fastq -o reads.jf
        jellyfish histo -t 10 reads.jf > reads.histo

# 4. Assemble the genome
    Canu: https://github.com/marbl/canu
        canu  -p CN_SC -d CN_SC useGrid=false genomeSize=250m -nanopore ./data/* 

    Minimap2: https://github.com/lh3/minimap2
    racon: https://github.com/isovic/racon
        minimap2 -t 16 -ax map-pb canu.draft.fasta highQuality-reads.fastq.gz > aln1.sam
        racon -t 16 highQuality-reads.fastq.gz aln1.sam canu.draft.fasta > racon1.fasta
        minimap2 -t 16 -ax map-pb racon1.fasta highQuality-reads.fastq.gz > aln2.sam
        racon -t 16 highQuality-reads.fastq.gz aln2.sam racon1.fasta > racon2.fasta
        rm aln2.sam
        minimap2 -t 16 -ax map-pb racon2.fasta highQuality-reads.fastq.gz > aln3.sam
        racon -t 16 highQuality-reads.fastq.gz aln3.sam racon2.fasta > racon3.fasta
        rm aln3.sam
    
    NextPolish: https://github.com/Nextomics/NextPolish
    SAMtools: https://github.com/lh3/samtools
    BWA: https://github.com/lh3/bwa
        round=2
        threads=32
        read1=NGS_clean.R1.fq.gz
        read2=NGS_clean.R2.fq.gz
        input=racon3.fasta
        for ((i=1; i<=${round};i++)); do
        #step 1:
        #index the genome file and do alignment
        bwa index ${input};
        bwa mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 32 -F 0x4 -b -|samtools fixmate -m --threads 32  - -|samtools sort -m 2g --threads 32 -|samtools markdup --threads 32 -r - sgs.sort.bam
        #index bam and genome files
        samtools index -@ ${threads} sgs.sort.bam;
        samtools faidx ${input};
        #polish genome file
        python /home/chichi/anaconda3/envs/nextpolish/share/nextpolish-1.4.1/lib/nextpolish1.py -g ${input} -t 1 -p ${threads} -s sgs.sort.bam > genome.polishtemp.fa;
        input=genome.polishtemp.fa;
        #step2:
        #index genome file and do alignment
        bwa index ${input};
        bwa mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 32 -F 0x4 -b -|samtools fixmate -m --threads 32  - -|samtools sort -m 2g --threads 32 -|samtools markdup --threads 32 -r - sgs.sort.bam
        #index bam and genome files
        samtools index -@ ${threads} sgs.sort.bam;
        samtools faidx ${input};
        #polish genome file
        python /home/chichi/anaconda3/envs/nextpolish/share/nextpolish-1.4.1/lib/nextpolish1.py -g ${input} -t 2 -p ${threads} -s sgs.sort.bam > genome.nextpolish.fa;
        input=genome.nextpolish.fa;
        done;
        #Finally polished genome file: genome.nextpolish.fa
    
    Ragtag: https://github.com/malonge/RagTag
        ragtag.py scaffold -o CN_SC --reference gapfree.fasta CN_SC.nextpolish.fasta
    
    TGS-GapClose: https://github.com/BGI-Qingdao/TGS-GapCloser
        tgsgapcloser --scaff ragtag.scaffold.fasta --reads TGS.fastq --output chinaclosegap --ne
        
    NextPolish: https://github.com/Nextomics/NextPolish
    SAMtools: https://github.com/lh3/samtools
    BWA: https://github.com/lh3/bwa
        round=2
        threads=32
        read1=NGS_clean.R1.fq.gz
        read2=NGS_clean.R2.fq.gz
        input=Tgs_gapclose.fasta
        for ((i=1; i<=${round};i++)); do
        #step 1:
        #index the genome file and do alignment
        bwa index ${input};
        bwa mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 32 -F 0x4 -b -|samtools fixmate -m --threads 32  - -|samtools sort -m 2g --threads 32 -|samtools markdup --threads 32 -r - sgs.sort.bam
        #index bam and genome files
        samtools index -@ ${threads} sgs.sort.bam;
        samtools faidx ${input};
        #polish genome file
        python /home/chichi/anaconda3/envs/nextpolish/share/nextpolish-1.4.1/lib/nextpolish1.py -g ${input} -t 1 -p ${threads} -s sgs.sort.bam > genome.polishtemp.fa;
        input=genome.polishtemp.fa;
        #step2:
        #index genome file and do alignment
        bwa index ${input};
        bwa mem -t ${threads} ${input} ${read1} ${read2}|samtools view --threads 32 -F 0x4 -b -|samtools fixmate -m --threads 32  - -|samtools sort -m 2g --threads 32 -|samtools markdup --threads 32 -r - sgs.sort.bam
        #index bam and genome files
        samtools index -@ ${threads} sgs.sort.bam;
        samtools faidx ${input};
        #polish genome file
        python /home/chichi/anaconda3/envs/nextpolish/share/nextpolish-1.4.1/lib/nextpolish1.py -g ${input} -t 2 -p ${threads} -s sgs.sort.bam > genome.nextpolish.fa;
        input=genome.nextpolish.fa;
        done;
        #Finally polished genome file: genome.nextpolish.fa

# 5. Evaluate the genome 
    compleasm: https://github.com/huangnengCSU/compleasm
        python compleasm.py run -a CN_SC.fasta -o CN_SC_compare -l embryophyta_odb10 -t 8
    
    tidk: https://github.com/tolkit/telomeric-identifier
        tidk search -s AAACCCT -out CN_SC CN_SC.fasta
    
# 6. Repeat and gene annotation
    EDTA: https://github.com/oushujun/EDTA
        EDTA.pl --genome CN_SC.fa --sensitive 1 --overwrite 1 --anno 1 --species others --threads 8
    
    RepeatMasker: https://github.com/Dfam-consortium/RepeatMasker
    RepeatModeler: https://github.com/Dfam-consortium/RepeatModeler
    
        BuildDatabase -name CN_SC CN_SC.fa
        RepeatModeler -database CN_SC -engine ncbi -threads 8
        RepeatMasker -pa 8 -s -lib /path/CN_SC-families.fa CN_SC.fa
    
    minimap2: https://github.com/lh3/minimap2
        minimap2 -t 8 -ax splice -uf CN_SC.fa.masked /path/.flRNA.fastq.gz | samtools view -@ 8 -bS - | samtools sort -@ 8 -o CN_SC.bam
    
    braker: https://github.com/Gaius-Augustus/BRAKER
        braker.pl --species=CN_SC --genome=/path/CN_SC.fa.masked --prot_seq=/path/proteins.fasta --softmasking --gff3 --threads=8 --bam=/path/CN_SC.bam \
        --GENEMARK_PATH /path/GeneMark-ETP/bin \
        --PROTHINT_PATH /patj/bin \
        --workingdir=/path/CN_SC_braker --min_contig=1000

# 7. Build the phylogentic tree
    ## the pep files of other species
    | Assembly  | Species   | Source |
    | :---: | :---: | :---: |
    |GCF_004353265.1    |*Vutis ripapria* |https://www.ncbi.nlm.nih.gov/assembly/GCF_004353265.1/|
    |GCA_902729315.2	|*Spirodela intermedia*   |https://www.ncbi.nlm.nih.gov/assembly/GCA_902729315.2/|
    |GCF_000001605.2	|*Sorghun bicolor*    |https://www.ncbi.nlm.nih.gov/assembly/GCF_000003195.3/|
    |GCF_001263595.1	|*Phalaenopsis equestris*   |https://www.ncbi.nlm.nih.gov/assembly/GCF_001263595.1/|
    |GCF_001433935.1	|*Oryza sativa*   |https://www.ncbi.nlm.nih.gov/assembly/GCF_001433935.1/|
    |GCF_000313855.2	|*Musa acuminata*  |https://www.ncbi.nlm.nih.gov/assembly/GCF_000313855.2/|
    |GCF_000442705.1	|*Elaeis guineensis*  |https://www.ncbi.nlm.nih.gov/assembly/GCF_000442705.1/|
    |GCF_001876935.1	|*Asparagus officinalis*  |https://www.ncbi.nlm.nih.gov/assembly/GCF_001876935.1/|
    |GCF_001540865.1	|*Ananas comosus* |https://www.ncbi.nlm.nih.gov/assembly/GCF_001540865.1/|
    |GCF_000130695.1	|*Amborella trichopoda*   |https://www.ncbi.nlm.nih.gov/assembly/GCF_000130695.1/|
    
    orthofinder: https://github.com/davidemms/OrthoFinder 
        orthofinder -f workdir -t 32 -S diamond

# 7. Align the single copy genes
    mafft
        mafft --auto single_copy_genes.fasta > aligned.fasta
    trimAl
        trimal -in aligned.fasta -out aligned.trim.fasta -automated1

    ## merge the aligned files with python script
    
# 8. Use BEAST estimate the divergence time
    Beast: https://beast.community/programs
        
# 9. Pan-genome construction with OrthoFinder
    orthofinder: https://github.com/davidemms/OrthoFinder
        orthofinder -f workdir2 -t 32 -S diamond

# 10. MCScanX
    MCScanX: https://github.com/wyp1125/MCScanX
        MCScanX Typha_latifolia 
    
# 11. Ka_Ks calculation
    Ka_Ks_Calculator: https://sourceforge.net/projects/kakscalculator2/
        KaKs_Calculator -i duplication.fasta -o duplication.kaks -m YN
    
# 12. SV, Indel, SNP detection
    minimap2: https://github.com/lh3/minimap2
    syri: https://github.com/schneebergerlab/syri

    minimap2 -ax asm5 --eqx CN_SC qrygenome.fa > out.sam
    syri -c out.sam -r CN_SC -q qrygenome.fa -o out
    ## snp and small indels
        BWA :
        GATK: 
        bwa index CN_SC.fasta
        bwa mem -R '@RG\tID:L001\tSM:CN_SC' -M CN_SC -t 8 read1.fastq read2.fastq > CN_SC.sam
        samtools view -bS CN_SC.sam | samtools sort -o CN_SC.bam
        gatk4 MarkDuplicates -I CN_SC.bam -O CN_SC.markdup.bam -M CN_SC.markdup.metrics.txt
        samtools index CN_SC.markdup.bam
        gatk4 HaplotypeCaller -R CN_SC.fasta -I CN_SC.markdup.bam -O CN_SC.g.vcf.gz -ERC GVCF
        
        gatk4 CombineGVCFs -R CN_SC.fasta --variant CN_SC.g.vcf.gz --variants CA_ON.g.vcf.gz --variant USA_PA1.g.vcf.gz --variant USA_PA2.g.vcf.gz --variant USA_CA.g.vcf.gz -O CN_SC_combined.g.vcf.gz
        gatk4 GenotypeGVCFs -R CN_SC.fasta -V CN_SC_combined.g.vcf.gz -O CN_SC_combined_raw.vcf.gz
        gatk4 variantFiltration -R CN_SC.fasta -V CN_SC_combined_raw.vcf.gz --filter-expression "DP < 30||DP>8000 || QD < 15.0 || MQ < 58.0 || FS > 2.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -4.0" --filter-name "my_snp_filter" -O CN_SC_combined_filtered.vcf.gz
        gatk4 SelectVariants -R CN_SC.fasta -V CN_SC_combined_filtered.vcf.gz -O raw.vcf.gz --exclude-filtered true
        gatk4 SelectVariants -R CN_SC.fasta -V CN_SC_combined_filtered.vcf.gz -O snp.vcf.gz --exclude-filtered true --select-type SNP
        gatk4 SelectVariants -R CN_SC.fasta -V CN_SC_combined_filtered.vcf.gz -O indel.vcf.gz --exclude-filtered true --select-type INDEL

# 13. Build the SNP based phylogenetic tree
    IQ-TREE: https://github.com/iqtree/iqtree3
        iqtree -s tree/raw_filter_variants.min4.fasta -T 32 -m TVM+F -bb 1000 -pre tree/boostrap_1000
        
# 14. PSMC
    bcftools: https://github.com/iqtree/iqtree3
    psmc: https://github.com/lh3/psmc

        bcftools mpileup -C50 -f CN_SC.fa CN_SC.rmdup.bam | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 200 | gzip >CN_SC_psmc.fq.gz
        fq2psmcfa -q20 CN_SC_psmc.fq.gz >CN_SC.psmcfa
        splitfa CN_SC.psmcfa > CN_SC_split.psmcfa
        seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o CN_SC_round-{}.psmc CN_SC_split.psmcfa | sh

# 15. Momi2 
    # all information record in momi2.ipynb