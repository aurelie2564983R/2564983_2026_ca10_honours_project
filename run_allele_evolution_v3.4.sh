#!/bin/bash


############# SLURM SETTINGS #############
#SBATCH --account=project0033   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=allele_evolution      # some descriptive job name of your choice
#SBATCH --output=%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=0-08:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=16G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

module load apps/miniforge



###### Inputs ######

rgt_in_dir="<path/to/dir>"
#directory contaning the RGT fastq output files. Should not end in a /

full_read_dir="<path/to/dir>"
#directory containing the full length ONT reads, 1MB either side of CA10.

csv_in="<path/to/csv>"
#path to csv file containing unphased sample information

out_dir="<path/to/dir>"
#directory to save output csv file of phased alleles. Should not end in a /

ref_seq="<path/to/reference>"
#path to HG38 fasta reference genome

region_of_interest="chr17:51821668-51841732"
# region of genome in which SNPs should be analysed to build the tree



conda activate phase_ont


#setup output CSVs

echo "Sample,Haplotype 1 Label, Haplotype 2 Label,Superpopulation" > "$out_dir""/""allele_labels.csv"



# for every sample in the compared output csv, test whether the RGT file exists

IFS=','

tail -n +2 $csv_in | while read -r sample superpopulation allele_1_length allele_2_length het_hom identical_hom allele_1_struc  allele_1_abun rgt_allele_1_len  same_len_allele_1 allele_1_EH_RGT_Diff  allele_2_struc allele_2_abun rgt_allele_2_len same_len_allele_2 allele_2_EH_RGT_Diff EH_RGT_comparison_notes
do
        #if a rgt file exists and rgt lengths are not within one repeat of each other
        file=$rgt_in_dir"/""$sample""_spliced.fastq"
        if test -f $file; then
                diff_1_2=$((rgt_allele_2_len-rgt_allele_1_len))
                if [ $diff_1_2 -ge 1 ]; then
                        echo
                        echo INFO: Processing $sample

			#output the read names for the reads that are the correct repeat lengths (reads have ~12bp flank (potentially -2bp depending on the flanking variant))
                        allele_1_min=$(((rgt_allele_1_len*3)+12-2))
                        allele_1_max=$(((rgt_allele_1_len*3)+12))
                        allele_2_min=$(((rgt_allele_2_len*3)+12-2))
                        allele_2_max=$(((rgt_allele_2_len*3)+12))
                        seqkit seq $file -m $allele_1_min -M $allele_1_max -n -g -o "$out_dir""/""$sample""_allele1_reads.txt"
                        seqkit seq $file -m $allele_2_min -M $allele_2_max -n -g -o "$out_dir""/""$sample""_allele2_reads.txt"


			
			# phase the long reads
			conda deactivate
			conda activate longshot
			longshot --out_bam "$out_dir""/""$sample""_longshot.bam" -F --bam "$full_read_dir""/"*"$sample"*".bam" --ref $ref_seq -s $sample --out "$out_dir""/""$sample""_longshot.vcf"
			conda deactivate
			conda activate phase_ont






			# check which haplotype contains the allele 1 reads
			# Extract read IDs belonging to each allele from the longshot bam output, and count the number of reads that have been marked as belonging to haplotype 1 and 2

			allele1_hap1=$(samtools view -h "$out_dir""/""$sample""_longshot.bam" | fgrep -f "$out_dir""/""$sample""_allele1_reads.txt" | grep -c 'HP:i:1')
			allele1_hap2=$(samtools view -h "$out_dir""/""$sample""_longshot.bam" | fgrep -f "$out_dir""/""$sample""_allele1_reads.txt" | grep -c 'HP:i:2')
			allele2_hap1=$(samtools view -h "$out_dir""/""$sample""_longshot.bam" | fgrep -f "$out_dir""/""$sample""_allele2_reads.txt" | grep -c 'HP:i:1')
			allele2_hap2=$(samtools view -h "$out_dir""/""$sample""_longshot.bam" | fgrep -f "$out_dir""/""$sample""_allele2_reads.txt" | grep -c 'HP:i:2')
			

			
			#check which RGT allele corresponds to which longshot haplotype. If there is no output for a SNP, assume 0/0 genotype

			if [ $allele1_hap1 -gt $allele1_hap2 -a $allele2_hap1 -lt $allele2_hap2 ]; then
				echo Allele 1 is hap 1, allele 2 is hap 2
				
				hap1_id="$allele_1_struc"" ""$sample""_1_""$rgt_allele_1_len"
				hap2_id="$allele_2_struc""_""$sample""_2_""$rgt_allele_2_len"


			elif [ $allele1_hap1 -lt $allele1_hap2 -a $allele2_hap1 -gt $allele2_hap2 ]; then
				echo Allele 1 is hap 2, allele 2 is hap 1

				hap1_id="$allele_2_struc""_""$sample""_2_""$rgt_allele_2_len"
				hap2_id="$allele_1_struc""_""$sample""_1_""$rgt_allele_1_len"
			else
				hap1_id="ERROR"
				hap2_id="ERROR"
				echo Could not resolve haplotypes
			fi
			
			
			# write labels to csv file for later use
			echo "$sample"",""$hap1_id"",""$hap2_id"",""$superpopulation" >> "$out_dir""/""allele_labels.csv"

		
			


			# compress vcf to use bcftools
			bgzip "$out_dir""/""$sample""_longshot.vcf"
			bcftools index "$out_dir""/""$sample""_longshot.vcf.gz"

			rm "$out_dir""/""$sample""_allele1_reads.txt"
			rm "$out_dir""/""$sample""_allele2_reads.txt"
			rm "$out_dir""/""$sample""_longshot.bam"
			rm "$out_dir""/""$sample""_longshot.bam.bai"
			#rm "$out_dir""/""$sample""_longshot.vcf.gz"
			#rm "$out_dir""/""$sample""_longshot.vcf.gz.csi"
		
		else
                        echo $sample alleles are only $diff_1_2 repeats apart, not processing this sample
                fi
fi
done


# combine VCFs from all samples

find $out_dir -type f -name "*.vcf.gz" > "$out_dir""/vcf.list"
bcftools merge -O z -o "$out_dir""/all_merged.vcf.gz" --file-list "$out_dir""/vcf.list"
bcftools index "$out_dir""/all_merged.vcf.gz"
echo "Merged VCF files"





#for each variant within the region of interest in the VCF, add the snp information for each allele to the string


# get list of SNP positions in VCF to use as header for CSV file

snp_list=$(bcftools query -f '%CHROM:%POS' -r $region_of_interest "$out_dir""/all_merged.vcf.gz" | sed ':a;N;$!ba;s/\n/,/g')
echo $snp_list

echo "Allele Label,Superpopulation,Interruption,""$snp_list" > "$out_dir""/roi_all_snps_haploid.csv"



IFS=','

tail -n +2 "$out_dir""/""allele_labels.csv" | while read -r sample hap1_label hap2_label superpopulation
do
	# write the allele label and snp string to the csv file
	bcftools query -f '[%GT]\n' -r $region_of_interest -s $sample "$out_dir""/all_merged.vcf.gz" > "$out_dir""/genotype_temp.txt"
	hap1_string=""
	hap2_string=""
	while read l; do
		hap1_genotype=${l:0:1}
		hap2_genotype=${l:2:3}
	
		# if the line is empty or there is no variant called for a specific position, classify the snp as being 0


		if [ -z "${hap1_genotype}" ]; then
			echo no hap1 genotype found for sample $sample at line $l
			hap1_genotype=0
		fi

		if [ -z "${hap2_genotype}" ]; then
			echo no hap2 genotype found for sample $sample at line $l
			hap2_genotype=0
		fi
		
		if [ $hap1_genotype = "." ]; then
			echo hap1 genotype is classed as . for $sample at line $l
			hap1_genotype=0
		fi

		if [ $hap2_genotype = "." ]; then
			echo hap2 genotype is classed as . for $sample at line $l
			hap2_genotype=0
		fi


		hap1_string="$hap1_string"",""$hap1_genotype"
		hap2_string="$hap2_string"",""$hap2_genotype"


	# check whether or not there are interruptions present in each haplotype

	if [[ $hap1_label == *"CAC"* ]]; then
		hap1_interruption="CAT[CAG]1CAC"
	elif [[ $hap1_label == *"CAT"* ]]; then
		hap1_interruption="[CAG]nCAT[CAG]n"
	elif [[ $hap1_label == *"CCG"* ]]; then
		hap1_interruption="CCG"
	elif [[ $hap1_label == *"CGG"* ]]; then
		hap1_interruption="CGG"
	else
		hap1_interruption="No Interruption"
	fi


	if [[ $hap2_label == *"CAC"* ]]; then
		hap2_interruption="CAT[CAG]1CAC"
	elif [[ $hap2_label == *"CAT"* ]]; then
		hap2_interruption="[CAG]nCAT[CAG]n"
	elif [[ $hap2_label == *"CCG"* ]]; then
		hap2_interruption="CCG"
	elif [[ $hap2_label == *"CGG"* ]]; then
		hap2_interruption="CGG"
	else
		hap2_interruption="No Interruption"
	fi


	done < "$out_dir""/genotype_temp.txt"
	echo "$hap1_label"",""$superpopulation"",""$hap1_interruption""$hap1_string" >> "$out_dir""/roi_all_snps_haploid.csv"
	echo "$hap2_label"",""$superpopulation"",""$hap2_interruption""$hap2_string" >> "$out_dir""/roi_all_snps_haploid.csv"
done


# remove any rows containing errors
grep -v "ERROR" "$out_dir""/roi_all_snps_haploid.csv" > "$out_dir""/roi_all_snps_haploid_no_error.csv"
conda deactivate



# make tree using scikit
echo Making Tree

conda activate matplotlib
in_csv="$out_dir""/roi_all_snps_haploid_no_error.csv"

python <<EOF

#import packages
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.cluster.hierarchy as sch
from seaborn import load_dataset
import radialtree as rt
from omniplot import plot as op

#make dataframe from csv, and remove any samples with errors
df = pd.read_csv("$in_csv")


# for the allele labels only keep the allele structure information (remove sample ID and allele number)
#df.iloc[:, 0]= df.iloc[:, 0].str.replace('.*_(.*)", "$1"', '', regex=True)



X = df.iloc[:,3:]
Z = sch.linkage(X, method='ward')
length = df['Allele Label'].str.replace(r'^.*_', '', regex=True)
length = np.asarray(length)
length = pd.to_numeric(length)
superpop = df['Superpopulation'].tolist()
interruption = df['Interruption'].tolist()



# count how many times the interruptions occur

# find how often each interruption occurs

interruptions=("CAT","CAT[CAG]1CAC","CCG","CGG")

for i in interruptions:
	df[i] = list(map(lambda x: x.count(i), df['Allele Label']))

#find the interruption that occurs most frequently (all alleles contain only one interruption motif, so this should produce the true estimate and avoid CATCAGCAC alleles as being detected as 2 interruptions) and save as list


interruption_counts= df.iloc[:,-4:]

df['number_of_interruptions'] = interruption_counts.max(axis='columns')
interruption_num = df['number_of_interruptions'].tolist()
print(interruption_num)

labels=[''] * X.shape[0]

fig = plt.figure()
dend = sch.dendrogram(Z, labels=labels, no_plot=True)



# define colours for superpopulation and interruption labels


interruption_color_map = {
	'No Interruption': '#0000ff',
	'CAT[CAG]1CAC': '#008000',
	'[CAG]nCAT[CAG]n': '#00BEBE',
	'CGG': '#d9b3ff',
	'CCG': '#BC00BC'
}

superpopulation_color_map = {
	'AFR': 'g',
	'AMR': 'm',
	'EAS': 'b',
	'EUR': 'y',
	'SAS': 'c'
}



# set up colourmap for the number of interruptions

interruption_num_color_map = {
        0: 'w',
	1: '#cad6e8',
        2: '#b5caeb',
        3: '#a4c1ed',
        4: '#90b6f0',
        5: '#7daaf0',
        6: '#6da2f2',
        8: '#5795f2',
        9: '#478cf5',
        21: '#2a59a1',
        23: '#1d488a',
        31: '#183b70',
        36: '#0f3063',
        37: '#0a2959',
        42: '#051630',
        44: 'k'

}



# set up colourmap for total repeat length

length_color_map=matplotlib.cm.get_cmap("gist_ncar", 145)


# set up spacer colourmap

spacer_color_map = {
	'': 'w'
}


colors = {
	"Number of Interruptions": [interruption_num_color_map[x] for x in interruption_num],
	"Interruption": [interruption_color_map[x] for x in interruption]
}


colorlabels_legend = {

    "Number of Interruptions": {
        "colors": list(interruption_num_color_map.values()),
        "labels": list(interruption_num_color_map.keys())
        },
    "Interruption": {
        "colors": list(interruption_color_map.values()),
        "labels": list(interruption_color_map.keys())
    }
}


sample_classes = {
    "Number of Interruptions": interruption_num,
    "Interruption": interruption
}




# create plot

X = X.astype(float)



rt.plot(
    dend,
    sample_classes=sample_classes,
    colorlabels=colors,
    colorlabels_legend=colorlabels_legend,
    figsize=(20, 8)
)




plt.savefig('$out_dir/20kb_roi_allele_evolution_denodrogram_interruptions.png',dpi=300)
plt.close()


# create plot showing the superpopulation and the total repeat length


length_norm = (length - length.min()) / (length.max() - length.min())
repeat_length_colors = length_color_map(length_norm)


colors = {
	"Repeat Length": length_color_map((length - length.min()) / (length.max() - length.min())), 
	"Superpopulation": [superpopulation_color_map[x] for x in superpop]
}



#plot colourbar for number of interruptions
ax = plt.gca()
norm = matplotlib.colors.Normalize(vmin=0, vmax=145)
sm = matplotlib.cm.ScalarMappable(cmap=plt.cm.gist_ncar, norm=norm)
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label("Repeat Length")

plt.savefig('$out_dir/20kb_roi_allele_evolution_denodrogram_len_colourbar.png')
plt.close()


colorlabels_legend = {

    "Repeat Length": {
	"colors": list(length_color_map(np.linspace(0, 1, 145))),
        "labels": [str(i) for i in range(145)]
    },
    "Superpopulation": {
        "colors": list(superpopulation_color_map.values()),
        "labels": list(superpopulation_color_map.keys())
    }
}




sample_classes = {
    "Spacer": labels,
    "Number of Interruptions": interruption_num,
    "Interruption": interruption,
    "Superpopulation": superpop
}





# create plot

X = X.astype(float)



rt.plot(
    dend,
    sample_classes=sample_classes,
    colorlabels=colors,
    colorlabels_legend=colorlabels_legend,
    figsize=(20, 8)
)




plt.savefig('$out_dir/20kb_roi_allele_evolution_denodrogram_superpop_len.png',dpi=300)
plt.close()





EOF

conda deactivate

echo Finished processing all samples

echo Finished processing all samples
