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
# Contigs sorted in reverse. Sort and rename with seqkit

    for Strain in AJ592 AJ705; do
        for Assembly in $(ls assembly/race_4/$Strain/"$Strain"_22022022.fasta); do
        OutDir=$(dirname $Assembly)
        seqkit sort --by-length --reverse $Assembly | seqkit replace --pattern '.+' --replacement 'contig_{nr}' > $OutDir/"$Strain"_renamed_unmasked.fa
        done
    done

    for Strain in AJ592 AJ705; do
        for Assembly in $(ls assembly/race_4/$Strain/"$Strain"_renamed_unmasked.fa); do
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
            ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Genome_aligners
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
OutDir=gene_pred/codingquarry/$Organism/$Strain/
mkdir -p $OutDir
GTF=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim/transcripts.gtf
ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/codingquarry.sh $Assembly $GTF $OutDir
done
done

for Strain in AJ592 AJ705; do
for Assembly in $(ls assembly/race_4/$Strain/"$Strain"*fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev) 
Organism=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/codingquarry/$Organism/$Strain/
mkdir -p $OutDir
GTF=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim/transcripts.gtf
ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Gene_prediction
sbatch $ProgDir/codingquarry.sh $Assembly $GTF $OutDir
done
done
```

## Add additional transcripts to Braker gene models.

Additional transcripts predicted by CodingQuarry are added to the final gene models.

```bash
# The following perl scripts requires the installation of some libraries. Run these commands in a perly environment.
# Install the required libraries (if any) using cpanm
# cpanm Bio::Perl
conda activate perly_env

# Perl imposes a maximum of 65,536 characters per line. If needed, fold your sequence

fold assembly/race_1/AJ520/AJ520_contigs_unmasked.fa > assembly/race_1/AJ520/AJ520_unmasked_folded.fasta
fold assembly/race_4/AJ516/AJ516_contigs_unmasked.fa > assembly/race_4/AJ516/AJ516_unmasked_folded.fasta
fold assembly/race_4/AJ516/AJ592_22022022.fasta > assembly/race_4/AJ592/AJ592_unmasked_folded.fasta
fold assembly/race_4/AJ705/AJ705_22022022.fasta > assembly/race_4/AJ705/AJ705_unmasked_folded.fasta
```

```bash
for Strain in AJ520 AJ516 AJ592 AJ705; do
for BrakerGff in $(ls gene_pred/braker/*/$Strain/augustus.hints.gff3); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f2 | rev)
Organism=$(echo $BrakerGff | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
Assembly=assembly/*/$Strain/"$Strain"_unmasked_folded.fasta
CodingQuarryGff=gene_pred/codingquarry/$Organism/$Strain/out/PredictedPass.gff3
PGNGff=gene_pred/codingquarry/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquarry/$Organism/$Strain/additional # Additional transcripts directory
FinalDir=gene_pred/codingquarry/$Organism/$Strain/final # Final directory
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

# Create a list with the additional transcripts in CondingQuarry gff (and CQPM) vs Braker gene models
bedtools intersect -v -a $CodingQuarryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList

# Creat Gff file with the additional transcripts
ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Gene_prediction
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuarryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff

# Create a final Gff file with gene features
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3

# Create fasta files from each gene feature in the CodingQuarry gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary

# Create fasta files from each gene feature in the Braker gff3
# If any "Possible precedence issue with control flow operator" error, change "or" for "||" operator
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker

# Combine both fasta files
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

# Combine both gff3 files
GffBraker=$FinalDir/final_genes_CodingQuary.gff3
GffQuary=$FinalDir/final_genes_Braker.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended

done
done
```
```bash
# Check the final number of genes
for Strain in AJ520 AJ516 AJ592 AJ705; do
for DirPath in $(ls -d gene_pred/codingquarry/*/$Strain/final); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
  ```

```
gene_pred/codingquarry/race_1/AJ520/final
18949
1384
20333

gene_pred/codingquarry/race_4/AJ516/final
21172
1890
23062

gene_pred/codingquarry/race_4/AJ592/final
20768
1903
22671

gene_pred/codingquarry/race_4/AJ705/final
21002
1732
22734
```


  ### Remove duplicate and rename genes.

  ```bash


for Strain in AJ520 AJ516 AJ592 AJ705; do
for GffAppended in $(ls gene_pred/codingquarry/*/$Strain/final/final_genes_appended.gff3); do
Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FinalDir=gene_pred/codingquarry/$Organism/$Strain/final
# Remove duplicated genes
GffFiltered=gene_pred/codingquarry/$Organism/$Strain/final/filtered_duplicates.gff
ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Gene_prediction
$ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
# Rename genes
GffRenamed=gene_pred/codingquarry/$Organism/$Strain/final/final_genes_appended_renamed.gff3
LogFile=gene_pred/codingquarry/$Organism/$Strain/final/final_genes_appended_renamed.log
$ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
rm $GffFiltered
# Create renamed fasta files from each gene feature   
Assembly=$(ls assembly/$Organism/$Strain/"$Strain"_unmasked_folded.fasta)
$ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed
# The proteins fasta file contains * instead of Xs for stop codons, these should be changed
sed -i 's/\*/X/g' $FinalDir/final_genes_appended_renamed.pep.fasta
done 
done
```

### SignalP and tmhmm 

```bash
conda activate annotation
```
```bash
# Version 6 is available and should be included here in the future.
for Strain in AJ520 AJ516 AJ592 AJ705; do
ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Feature_annotation
CurPath=$PWD
for Proteome in $(ls gene_pred/codingquarry/*/$Strain/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/final_genes_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_final_preds
$ProgDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName 
for File in $(ls $SplitDir/*_final_preds_*); do
sbatch $ProgDir/pred_signalP.sh $File signalp-4.1 
done
done
done
```





The batch files of predicted secreted proteins needed to be combined into a single file for each strain. This was done with the following commands:

 ```bash
for Strain in AJ520 AJ516 AJ592 AJ705; do
  for SplitDir in $(ls -d gene_pred/final_genes_split/*/$Strain); do
    Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
    Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
    InStringAA=''
    InStringNeg=''
    InStringTab=''
    InStringTxt=''
    SigpDir=final_genes_signalp-4.1
    for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do
      InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";
      InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";
      InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
      InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";
    done
    cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
    cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
    tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
    cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
  done
done
```

Proteins containing a transmembrane domain were identified:

```bash
for Strain in AJ520 AJ516 AJ592 AJ705; do
for Proteome in $(ls gene_pred/*/*/$Strain/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Feature_annotation
sbatch $ProgDir/TMHMM.sh $Proteome
done
done
 ```

 Those proteins with transmembrane domains were removed from lists of Signal peptide containing proteins

 ```bash
for Strain in AJ520 AJ516 AJ592 AJ705; do
for File in $(ls gene_pred/trans_mem/*/$Strain/*_TM_genes_neg.txt); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $TmHeaders
SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
OutDir=$(dirname $SigP)
ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Feature_annotation
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
done
```


### EffectorP

```bash
    for Strain in AJ520 AJ516 AJ592 AJ705; do
        for Proteome in $(ls gene_pred/codingquarry/*/$Strain/final/final_genes_appended_renamed.pep.fasta); do
        Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
        Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
        BaseName="$Organism"_"$Strain"_EffectorP
        Version=3.0 # Version 2.0 or 3.0
        OutDir=analysis/effectorP_"$Version"/$Organism/$Strain
        ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Feature_annotation
        sbatch $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir $Version
        done
    done
```

```bash
# List of effectors
for Strain in AJ520 AJ516 AJ592 AJ705; do
  for File in $(ls analysis/effectorP*/*/$Strain/*_EffectorP.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    Effversion=$(echo $File | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
      if [ $Effversion == "effectorP_2.0" ]; then 
      Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
      cat $File | grep 'Effector' | cut -f1 > $Headers
      else 
      CitHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_Cytoplasmic_Effectors_headers.txt/g')
      cat $File | grep 'Cytoplasmic effector\|/' | cut -f1 > $CitHeaders
      ApHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_Apoplastic_Effectors_headers.txt/g')
      cat $File | grep 'Apoplastic effector\|/' | cut -f1 > $ApHeaders
      fi
  done
done

# sed can be used to remove first line and any info after gene id
# cat *_EffectorP_headers.txt| sed -i s/\ //g | sed '1d' > Clean_EffectorP_headers.txt 

# Extract secreted effectors 
for Strain in AJ516 AJ592 AJ705; do
cat analysis/effectorP_3.0/race_4/$Strain/*_Cytoplasmic_Effectors_headers.txt | sed '1d' > analysis/effectorP_3.0/race_4/$Strain/cit.txt
cat analysis/effectorP_3.0/race_4/$Strain/*_Apoplastic_Effectors_headers.txt | sed '1d' > analysis/effectorP_3.0/race_4/$Strain/apo.txt
cat analysis/effectorP_3.0/race_4/$Strain/cit.txt analysis/effectorP_3.0/*/$Strain/apo.txt | sort | uniq > analysis/effectorP_3.0/race_4/$Strain/"$Strain"_Effector_header_combined.txt
done

for Strain in AJ520; do
cat analysis/effectorP_3.0/race_1/$Strain/*_Cytoplasmic_Effectors_headers.txt | sed '1d' > analysis/effectorP_3.0/race_1/$Strain/cit.txt
cat analysis/effectorP_3.0/race_1/$Strain/*_Apoplastic_Effectors_headers.txt | sed '1d' > analysis/effectorP_3.0/race_1/$Strain/apo.txt
cat analysis/effectorP_3.0/race_1/$Strain/cit.txt analysis/effectorP_3.0/*/$Strain/apo.txt | sort | uniq > analysis/effectorP_3.0/race_1/$Strain/"$Strain"_Effector_header_combined.txt
done
```

## Identification of MIMP-flanking genes

Miniature impala (mimp) sequeces are found in promotor regions of SIX genes in fusarium.

### Typical run

```bash
conda activate perly_env

for Strain in AJ520 AJ516 AJ592 AJ705; do
for Assembly in $(ls assembly/*/$Strain/"$Strain"_unmasked_folded.fasta); do
Organism=$(echo "$Assembly" | rev | cut -d '/' -f3 | rev)
Strain=$(echo "$Assembly" | rev | cut -d '/' -f2 | rev)
GeneGff=$(ls gene_pred/codingquarry/*/$Strain/final/final_genes_appended_renamed.gff3)
OutDir=analysis/mimps/$Organism/$Strain
mkdir -p "$OutDir"
echo "$Organism - $Strain"
ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Feature_annotation
$ProgDir/mimp_finder.pl $Assembly $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
$ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
echo "The number of mimps identified:"
cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
echo "The following transcripts intersect mimps:"
MimpProtsTxt=$OutDir/"$Strain"_prots_in_2kb_mimp.txt
MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq > $MimpProtsTxt
cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | cut -f1 -d '.'| sort | uniq > $MimpGenesTxt
cat $MimpProtsTxt | wc -l
cat $MimpGenesTxt | wc -l
echo ""
done
done
```








### CAZY proteins

Carbohydrte active enzymes were identified from the CAZy database

```bash
for Strain in AJ520 AJ516 AJ592 AJ705; do 
  for Proteome in $(ls gene_pred/codingquarry/*/$Strain/final/final_genes_appended_renamed.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/CAZY/$Organism/$Strain
    mkdir -p $OutDir
    Prefix="$Strain"_CAZY
    CazyHmm=../dbCAN/dbCAN-HMMdb-V8.txt # databases are in /projects
    ProgDir=/home/agomez/scratch/apps/git_repos/bioinformatics_tools/Feature_annotation
    sbatch $ProgDir/hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
  done
done
```
```bash
 for File in $(ls path/to/*CAZY.out.dm); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $File)
  echo "$Organism - $Strain"
  ProgDir=/projects/dbCAN # Script from dbCAN tools
  $ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps # Creates a file with CAZy module and gene

  CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
  cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders # Extract gene names
  echo "Number of CAZY genes identified:"
  cat $CazyHeaders | wc -l

  Gff=$(ls path/to/final/gff3/file/final_genes_appended_renamed.gff3)
  CazyGff=$OutDir/"$Strain"_CAZY.gff
  ProgDir=/home/gomeza/git_repos/scripts/bioinformatics_tools/Feature_annotation
  $ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff # Creates a gff for all CAZymes
  
  SecretedProts=$(ls path/to/secreted/proteins/*signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
  SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
  cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
  CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
  $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted # Creates a gff for secreted CAZymes
  echo "Number of Secreted CAZY genes identified:"
  cat $CazyGffSecreted | grep -w 'gene' | cut -f9 | tr -d 'ID=' | wc -l
  done
  ```
