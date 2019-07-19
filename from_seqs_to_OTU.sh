
#!/bin/bash

########################################################################################
########################################################################################

# This is a collection of bash command, used for the initial steps of data analysis;
# from raw sequences to OTUs. Those command were used for microbiota analysis of article

# "IL-13 mRNA tissue content identifies two subsets of adult ulcerative colitis patients 
# with different clinical and mucosa-associated microbiota profiles"

# By authors: Alessia Butera, Monica Di Paola, Francesco Vitali, Daniela De Nitto, 
# Francesco Covotta, Francesco Borrini, Roberta Pica, Carlotta De Filippo, Duccio Cavalieri, 
# Alessandro Giuliani, Annamaria Pronio, Monica Boirivan

# Script author: Francesco Vitali

# DISCLAIMER: this is the collection of bash command that were used for analysis. 
# However, for publication purpose, absolute filepath of my local PC were changed with 
# variables (declared at the beginning of the script). As such, they were not tested and 
# are not guaranteed to work

# DISCLAIMER: Use these command and the information contained here at your own risk!
# I'm not responsible for loss of data

########################################################################################
########################################################################################

# Raw reads files are assumed to be in $currpath/Raw_reads !!

#################################################
# SETTING UP VARIABLES AND CREATING DIRECTORIES #

currpath=$(pwd) # top project directory
core=6 # number of core to use

mkdir $currpath/cutadapt/
mkdir $currpath/quality_control/
mkdir $currpath/renamed/
mkdir $currpath/seq_fasta
mkdir $currpath/cutadapt/trimmed
mkdir $currpath/cutadapt/trimmed/assembled

cd $currpath

##########################################################################################
# Renaming files from the sequencing service provider are usually redundant and "messy". #

# Here samples are renamed. Code should be edited depending on the type of name
# name assumed here is like SampleID-V3-V4_S239_L001_R1_001.fastq.gz
# where SampleID and R1 are the interesting part. "cut" is used to divide
# the name string on a desired character

# output would be SampleID-V3-V4_R1.fastq.gz

name=0
este=".fastq.gz"
for i in *.fastq.gz
do
	sample=$(echo "$i" | cut -d "_" -f1)
	read=$(echo "$i" | cut -d "_" -f4)
	mv $i $currpath/renamed/$sample"_"$read$este
	echo -e "$i\t-->\t$sample"_"$read$este" >> log_renamer.txt
done

##################################
# Counting raw reads per sample #

for i in $currpath/renamed/*.fastq.gz
do
        echo -n $i ' \t ' >> seq_count_16S.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_16S.txt
done

#######################################################################
# Creating a text file with sample names, it will be used for looping #

for i in $currpath/renamed/*.fastq.gz
do 
echo "$i" | cut -d "_" -f1 | >> names.txt
sed 'n; d' names.txt > names_single.txt
done

# coping the file to other folders
cp names_single.txt $currpath/cutadapt/
cp names_single.txt $currpath/cutadapt/trimmed/


#############################################
# Removing sequencing primers with CUTADAPT # 

# Reference: http://journal.embnet.org/index.php/embnetjournal/article/view/200
# Guide: https://cutadapt.readthedocs.io/en/stable/guide.html

while read file
do
	echo "Running cutadapt on file "${file}""
	'#qui le seq per 16S della BMr, cambiare se altri primer'
	cutadapt -a CCTACGGGNBGCA...GGATTAGATACCCBNGTAGTC -A GACTACNVGGGTATCTAATCC...CTGSTGCVNCCCGTAGG \
	--trim-n -o $currpath/cutadapt/"${file}_R1_cutadapt.fastq.gz" -p $currpath/cutadapt/"${file}_R2_cutadapt.fastq.gz" \
	$currpath/renamed/"${file}_R1.fastq.gz" $currpath/renamed/"${file}_R2.fastq.gz" 1>> report_cutadapt_primer.txt
done < names_single.txt

cp  report_cutadapt_primer.txt  $currpath/quality_control

#################################################################
# Evaluating quality of the reads to choose trimming parameters #

fastqc $currpath/cutadapt/*.fastq.gz -o $currpath/quality_control/ -t $core
multiqc $currpath/quality_control/. -o $currpath/quality_control/


##################################################################
# TRIMMING: removing low quality region towards the end of reads #

# Reference: https://github.com/najoshi/sickle


while read file
do
	echo "Running sickle on file "${file}""
	echo "Running sickle on file "${file}"" >> stats_trim.txt
	sickle pe -f $currpath/cutadapt/"${file}_R1_cutadapt.fastq.gz" -r $currpath/cutadapt/"${file}_R2_cutadapt.fastq.gz" \
	-o $currpath/cutadapt/trimmed/"${file}_trimmed_R1.fastq" -p $currpath/cutadapt/trimmed/"${file}_trimmed_R2.fastq" -s $currpath/cutadapt/trimmed/"${file}_singles" \
	-t sanger -q 20 -l 200 -g 1>> stats_trim.txt
done < names_single.txt


#####################################
# Joining forward and reverse read  # 

# Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3933873/
# Guide: https://cme.h-its.org/exelixis/web/software/pear/doc.html

while read file
do
	echo "Running pear on file "${file}""
	pear -f $currpath/cutadapt/trimmed/"${file}_trimmed_R1.fastq" -r $currpath/cutadapt/trimmed/"${file}_trimmed_R2.fastq" \
	-o $currpath/cutadapt/trimmed/assembled/"${file}" -j $core -p 0.001 -u 0 -q 30 1>> stats_assembly.txt
done < names_single.txt


#################################################################
# Evaluating quality on the final, pre-treated and joined reads #

fastqc $currpath/cutadapt/trimmed/assembled/*.fastq.gz -o $currpath/quality_control/ -t $core
multiqc $currpath/quality_control/. -o $currpath/quality_control/


##########################
# Convert fastq in fasta #

for i in $currpath/cutadapt/trimmed/assembled/*.fastq
do
	echo "Converting "$i""
	name=$(echo "$i" | cut -d "." -f1,2)
	fastq_to_fasta -i ./$i -o $currpath/seq_fasta/$name.fasta
done

###################################
# Validate mapping file for QIIME #

validate_mapping_file.py -m $currpath/mapping_colitis.csv -p -b

####################
# Add QIIME header #

add_qiime_labels.py -i $currpath/seq_fasta -m $currpath/mapping_colitis.csv -c InputFileName -o $currpath/seqs_with_header

#####################
# Chimera filtering #

# Reference: https://www.ncbi.nlm.nih.gov/pubmed/27781170
# Guide: https://github.com/torognes/vsearch

vsearch -uchime_ref $currpath/seqs_with_header/combined_seqs.fna --db /usr/local/lib/python2.7/dist-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta \
--nonchimeras $currpath/combined_seqs_chimera_filtered.fna --threads $core

##############################
# Open reference otu picking #

cat > $currpath/otu_picking_param.txt
pick_otus:enable_rev_strand_match True
pick_otus:otu_picking_method uclust

pick_open_reference_otus.py -i $currpath/combined_seqs_chimera_filtered.fna -r /usr/local/lib/python2.7/dist-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta \
-o $currpath/open_ref_otu_picking_uclust -s 0.1 -m uclust -p $currpath/otu_picking_param.txt -a -O 4

########################################
# Discard low sequence coverage OTUs   #

#REFERENCE: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3531572/

filter_otus_from_otu_table.py -i $currpath/open_ref_otu_picking_uclust/otu_table_mc2_w_tax_no_pynast_failures.biom -o $currpath/otu_table_mc2_w_tax_no_pynast_failures_filtered.biom \
--min_count_fraction 0.00005
