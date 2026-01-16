#!/bin/bash


############# SLURM SETTINGS #############
#SBATCH --account=project0033   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=phase_ont      # some descriptive job name of your choice
#SBATCH --output=%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=0-08:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=64G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

module load apps/miniforge



###### Inputs ######

rgt_in_dir="<path/to/directory>"
#directory contaning the RGT output fastq files. Should not end in a /


full_read_dir="<path/to/directory>"
#directory containing the full length ONT reads 1MB either side of CA10. Should not end in a /

csv_in="<path/to/csv>"
#path to csv file containing unphased sample information

out_dir="<path/to/directory>"
#directory to save output csv file of phased alleles. Should not end in a /

ref_seq="<path/to/reference>"
#path to Hg38 fasta reference genome

conda activate phase_ont



#set up/overwrite a phased ont alleles csv

echo "Sample Name,Superpopulation,RGT Allele Structure,RGT Allele Abundance,RGT Allele Length,CF SNP Genotype,Insomnia SNP Genotype,Chronotype SNP Genotype,Morningness SNP Genotype" | cat > "$out_dir""/ca10_ont_phased.csv"
# for every sample in the compared output csv, test whether the RGT file exists

IFS=','

tail -n +2 $csv_in | while read -r sample superpopulation allele_1_length allele_2_length het_hom identical_hom allele_1_struc  allele_1_abun rgt_allele_1_len  same_len_allele_1 allele_1_EH_RGT_Diff  allele_2_struc allele_2_abun rgt_allele_2_len same_len_allele_2 allele_2_EH_RGT_Diff EH_RGT_comparison_notes
do
        #if a rgt file exists and rgt lengths are not within five repeats of each other
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
			longshot --out_bam "$out_dir""/""$sample""_longshot.bam" -F --bam "$full_read_dir""/"*"$sample"*".bam" --ref $ref_seq --out "$out_dir""/""$sample""_longshot.vcf"
			conda deactivate
			conda activate phase_ont



			# Get SNP information from the VCF file

                        # compress vcf to use bcftools
			bgzip "$out_dir""/""$sample""_longshot.vcf"
			bcftools index "$out_dir""/""$sample""_longshot.vcf.gz"
			
			
			# Extract SNP info

			cf_snp=$(bcftools query -r chr17:52183005-52183007 -f '[%GT]\n' "$out_dir""/""$sample""_longshot.vcf.gz")
			insomnia_snp=$(bcftools query -r chr17:52181781-52181783 -f '[%GT]\n' "$out_dir""/""$sample""_longshot.vcf.gz")
			chronotype_snp=$(bcftools query -r chr17:52014840-52014842 -f '[%GT]\n' "$out_dir""/""$sample""_longshot.vcf.gz")
			morningness_snp=$(bcftools query -r chr17:52129455-52129457 -f '[%GT]\n' "$out_dir""/""$sample""_longshot.vcf.gz")

			echo CF: $cf_snp
			echo Insomnia: $insomnia_snp
			echo Chronotype: $chronotype_snp
			echo Morningness: $morningness_snp

			# check which haplotype contains the allele 1 reads
			# Extract read IDs belonging to each allele from the longshot bam output, and count the number of reads that have been marked as belonging to haplotype 1 and 2

			allele1_hap1=$(samtools view -h "$out_dir""/""$sample""_longshot.bam" | fgrep -f "$out_dir""/""$sample""_allele1_reads.txt" | grep -c 'HP:i:1')
			allele1_hap2=$(samtools view -h "$out_dir""/""$sample""_longshot.bam" | fgrep -f "$out_dir""/""$sample""_allele1_reads.txt" | grep -c 'HP:i:2')
			allele2_hap1=$(samtools view -h "$out_dir""/""$sample""_longshot.bam" | fgrep -f "$out_dir""/""$sample""_allele2_reads.txt" | grep -c 'HP:i:1')
			allele2_hap2=$(samtools view -h "$out_dir""/""$sample""_longshot.bam" | fgrep -f "$out_dir""/""$sample""_allele2_reads.txt" | grep -c 'HP:i:2')
			

			
			#check which RGT allele corresponds to which longshot haplotype. If there is no output for a SNP, assume 0/0 genotype

			if [ $allele1_hap1 -gt $allele1_hap2 -a $allele2_hap1 -lt $allele2_hap2 ]; then
				echo Allele 1 is hap 1, allele 2 is hap 2
				
				if [ -z "${cf_snp}" ]; then
					cf_snp1=0
					cf_snp2=0
				else
					cf_snp1=${cf_snp:0:1}
					cf_snp2=${cf_snp:2:3}
				fi
				

				if [ -z "${insomnia_snp}" ]; then
                                        insomnia_snp1=0
                                        insomnia_snp2=0
                                else
                                        insomnia_snp1=${insomnia_snp:0:1}
                                        insomnia_snp2=${insomnia_snp:2:3}
                                fi
				
				if [ -z "${chronotype_snp}" ]; then
                                        chronotype_snp1=0
                                        chronotype_snp2=0
                                else
                                        chronotype_snp1=${chronotype_snp:0:1}
                                        chronotype_snp2=${chronotype_snp:2:3}
                                fi
				

				if [ -z "${morningness_snp}" ]; then
                                        morningness_snp1=0
                                        morningness_snp2=0
                                else
                                        morningness_snp1=${morningness_snp:0:1}
                                        morningness_snp2=${morningness_snp:2:3}
                                fi
			elif [ $allele1_hap1 -lt $allele1_hap2 -a $allele2_hap1 -gt $allele2_hap2 ]; then
				echo Allele 1 is hap 2, allele 2 is hap 1


                                if [ -z "${cf_snp}" ]; then
                                        cf_snp1=0
                                        cf_snp2=0
                                else
                                        cf_snp2=${cf_snp:0:1}
                                        cf_snp1=${cf_snp:2:3}
                                fi


                                if [ -z "${insomnia_snp}" ]; then
                                        insomnia_snp1=0
                                        insomnia_snp2=0
                                else
                                        insomnia_snp2=${insomnia_snp:0:1}
                                        insomnia_snp1=${insomnia_snp:2:3}
                                fi

                                if [ -z "${chronotype_snp}" ]; then
                                        chronotype_snp1=0
                                        chronotype_snp2=0
                                else
                                        chronotype_snp2=${chronotype_snp:0:1}
                                        chronotype_snp1=${chronotype_snp:2:3}
                                fi


                                if [ -z "${morningness_snp}" ]; then
                                        morningness_snp1=0
                                        morningness_snp2=0
                                else
                                        morningness_snp2=${morningness_snp:0:1}
                                        morningness_snp1=${morningness_snp:2:3}
                                fi
			else
				echo Could not resolve haplotypes
				cf_snp1="ERROR"
                                cf_snp2="ERROR"
				insomnia_snp1="ERROR"
                                insomnia_snp2="ERROR"
				chronotype_snp1="ERROR"
                                chronotype_snp2="ERROR"
				morningness_snp1="ERROR"
                                morningness_snp2="ERROR"
			fi
			
			
			# write outputs to file
				
			echo "$sample"",""$superpopulation"",""$allele_1_struc"",""$allele_1_abun"",""$rgt_allele_1_len"",""$cf_snp1"",""$insomnia_snp1"",""$chronotype_snp1"",""$morningness_snp1" | cat >> "$out_dir""/ca10_ont_phased.csv"
			echo "$sample"",""$superpopulation"",""$allele_2_struc"",""$allele_2_abun"",""$rgt_allele_2_len"",""$cf_snp2"",""$insomnia_snp2"",""$chronotype_snp2"",""$morningness_snp2" | cat >> "$out_dir""/ca10_ont_phased.csv"
			
			rm "$out_dir""/""$sample""_allele1_reads.txt"
			rm "$out_dir""/""$sample""_allele2_reads.txt"
			rm "$out_dir""/""$sample""_longshot.bam"
			rm "$out_dir""/""$sample""_longshot.bam.bai"
			rm "$out_dir""/""$sample""_longshot.vcf.gz"
			rm "$out_dir""/""$sample""_longshot.vcf.gz.csi"
		
		else
                        echo $sample alleles are only $diff_1_2 repeats apart, not processing this sample
                fi
fi
done


conda deactivate

echo Finished processing all samples
