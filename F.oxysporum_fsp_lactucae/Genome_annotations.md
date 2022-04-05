# Analysis of F.oxysporum f.sp. lactucae genomes

### Data organization

```bash
# Raw data stored on NIAB HPC. Working directory on gruff
scp -r /main/projects/Fus_all_Ex/RNAseq/X204SC21111729-Z01-F001/raw_data/Folac_R1 agomez@gruffalo.cropdiversity.ac.uk:/home/agomez/scratch/fusarium_species/F.oxysporum_fsp_lactucae/raw_rna/X204SC21111729-Z01-F001/raw_data
scp -r /main/projects/Fus_all_Ex/RNAseq/X204SC21111729-Z01-F001/raw_data/Folac_R4 agomez@gruffalo.cropdiversity.ac.uk:/home/agomez/scratch/fusarium_species/F.oxysporum_fsp_lactucae/raw_rna/X204SC21111729-Z01-F001/raw_data
```

## RNA-Seq qc

```bash
# Run fastqc
    for Strain in Folac_R1_1.fq.gz Folac_R1_2.fq.gz; do
        for RawData in $(ls raw_rna/X204SC21111729-Z01-F001/raw_data/Folac_R1/$Strain); do
        echo $RawData
        OutDir=qc_rna/X204SC21111729-Z01-F001/race_1/AJ520/fastqc/raw
        ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p short $ProgDir/fastqc.sh $RawData $OutDir
        done
    done

    for Strain in Folac_R4_1.fq.gz Folac_R4_2.fq.gz; do
        for RawData in $(ls raw_rna/X204SC21111729-Z01-F001/raw_data/Folac_R4/$Strain); do
        echo $RawData
        OutDir=qc_rna/X204SC21111729-Z01-F001/race_4/AJ516/fastqc/raw
        ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p short $ProgDir/fastqc.sh $RawData $OutDir
    done
```

```bash
# Run fastq-mcf
    for Strain  in Folac_R1 Folac_R4; do
        for RNADir in $(ls -d raw_rna/X204SC21111729-Z01-F001/raw_data/$Strain); do
        FileNum=$(ls $RNADir/*_1.fq.gz | wc -l)
            for num in $(seq 1 $FileNum); do
            printf "\n"
            FileF=$(ls $RNADir/*1.fq.gz | head -n $num | tail -n1)
            FileR=$(ls $RNADir/*2.fq.gz | head -n $num | tail -n1)
            echo $FileF
            echo $FileR
            if [ $Strain == Folac_R1 ]; then
            OutDir=qc_rna/X204SC21111729-Z01-F001/race_1/AJ520
            else
            OutDir=qc_rna/X204SC21111729-Z01-F001/race_4/AJ516
            fi
        echo $OutDir
        IluminaAdapters=/home/agomez/scratch/apps/scripts/bioinformatics_tools/SEQdata_qc/illumina_full_adapters.fa
        ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p long $ProgDir/fastq-mcf_himem.sh $FileF $FileR $IluminaAdapters RNA $OutDir
        done
    done
```

```bash
# Decontamination of rRNA reads
    for Strain  in AJ520 AJ516; do
        for RNADir in $(ls -d qc_rna/X204SC21111729-Z01-F001/*/$Strain); do
        FileNum=$(ls $RNADir/F/*_1_trim.fq.gz | wc -l)
        for num in $(seq 1 $FileNum); do
        printf "\n"
        FileF=$(ls $RNADir/F/*trim.fq.gz | head -n $num | tail -n1)
        FileR=$(ls $RNADir/R/*trim.fq.gz | head -n $num | tail -n1)
        echo $FileF
        echo $FileR
        Ref=/home/agomez/scratch/apps/prog/bbmap/ribokmers.fa.gz
        ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/SEQdata_qc
        echo $RNADir
        sbatch -p himem $ProgDir/bbduk.sh $Ref "$RNADir"/cleaned $FileF $FileR $ProgDir $Strain
        done
    done
```
```bash
# Run fastqc
    for Strain  in AJ520 AJ516; do
        for RawData in $(ls qc_rna/X204SC21111729-Z01-F001/*/$Strain/cleaned/*/*cleaned.fq.gz); do
        echo $RawData
        if [ $Strain == AJ520 ]; then
        OutDir=qc_rna/X204SC21111729-Z01-F001/race_1/AJ520/fastqc/cleaned
        else
        OutDir=qc_rna/X204SC21111729-Z01-F001/race_4/AJ516/fastqc/cleaned
        fi
        ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/SEQdata_qc
        sbatch -p short $ProgDir/fastqc.sh $RawData $OutDir
        done
    done
```

## Star alignment

```bash
# race 1 and 4
    for Strain in AJ520 AJ516; do
        for Assembly in $(ls assembly/race_*/$Strain/"$Strain"_contigs_unmasked.fa); do
        echo "$Assembly"
        Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev) 
        Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
        echo "$Organism - $Strain"
            for FileF in $(ls qc_rna/X204SC21111729-Z01-F001/race_*/$Strain/cleaned/F/Folac_R*_1_cleaned.fq.gz)
            do
            FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_cleaned/_2_cleaned/g')
            echo $FileF
            echo $FileR
            Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_cleaned.fq.gz//g')
            OutDir=alignment/star/$Organism/$Strain/$Sample_Name
            ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/Genome_aligners
            sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir 11
            done
        done
    done
# race 4
    for Strain in AJ592 AJ705; do
        for Assembly in $(ls assembly/race_4/$Strain/"$Strain"*fasta); do
        echo "$Assembly"
        Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev) 
        Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
        echo "$Organism - $Strain"
            for FileF in $(ls qc_rna/X204SC21111729-Z01-F001/race_*/AJ516/cleaned/F/Folac_R*_1_cleaned.fq.gz)
            do
            FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1_cleaned/_2_cleaned/g')
            echo $FileF
            echo $FileR
            Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_cleaned.fq.gz//g')
            OutDir=alignment/star/$Organism/$Strain/$Sample_Name
            ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/Genome_aligners
            sbatch $ProgDir/star.sh $Assembly $FileF $FileR $OutDir 11
            done
        done
    done
```

## Gene prediction with Braker

```bash
    for Strain in AJ520 AJ516; do
        for Assembly in $(ls assembly/race_*/$Strain/"$Strain"_contigs_unmasked.fa); do
            Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
            Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
            OutDir=gene_pred/braker/$Organism/$Strain
            AcceptedHits=alignment/star/assembly/$Strain/*/*_aligmentAligned.sortedByCoord.out.bam 
            GeneModelName="$Organism"_"$Strain"_braker 
            ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/Gene_prediction
            sbatch -p himem $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
        done
    done

    for Strain in AJ592 AJ705; do
        for Assembly in $(ls assembly/race_4/$Strain/"$Strain"*fasta); do
            Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
            Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
            OutDir=gene_pred/braker/$Organism/$Strain
            AcceptedHits=alignment/star/race_4/$Strain/*/*_aligmentAligned.sortedByCoord.out.bam 
            GeneModelName="$Organism"_"$Strain"_braker 
            ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/Gene_prediction
            sbatch -p himem $ProgDir/braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
        done
    done
```


## Cufflinks

scp -r /projects/oldhome/armita/prog/cufflinks/ agomez@gruffalo.cropdiversity.ac.uk:/home/agomez/scratch

```bash
    for Strain in AJ520 AJ516; do
        for Assembly in $(ls assembly/race_*/$Strain/"$Strain"_contigs_unmasked.fa); do
            Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
            Organism=$(echo $Assembly| rev | cut -d '/' -f3 | rev) 
            echo "$Organism - $Strain"
            OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
            mkdir -p $OutDir
            AcceptedHits=alignment/star/race_*/$Strain/*/*_aligmentAligned.sortedByCoord.out.bam 
            ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/Gene_prediction
            sbatch $ProgDir/cufflinks.sh $AcceptedHits $OutDir
        done
    done

    for Strain in AJ592 AJ705; do
        for Assembly in $(ls assembly/race_4/$Strain/"$Strain"*fasta); do
            Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
            Organism=$(echo $Assembly| rev | cut -d '/' -f3 | rev) 
            echo "$Organism - $Strain"
            OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
            mkdir -p $OutDir
            AcceptedHits=alignment/star/race_4/$Strain/*/*_aligmentAligned.sortedByCoord.out.bam 
            ProgDir=/home/agomez/scratch/apps/scripts/bioinformatics_tools/Gene_prediction
            sbatch $ProgDir/cufflinks.sh $AcceptedHits $OutDir
        done
    done
```

```bash
for Strain in AJ520 AJ516; do
for Assembly in $(ls assembly/race_*/$Strain/"$Strain"_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
Organism=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/codingquary/$Organism/$Strain/
mkdir -p $OutDir
GTF=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim/transcripts.gtf
ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/codingquarry.sh $Assembly $GTF $OutDir
done
done
```