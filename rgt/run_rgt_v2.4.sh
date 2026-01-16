#!/bin/bash -l


############# SLURM SETTINGS #############
#SBATCH --account=project0033   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=rgt      # some descriptive job name of your choice
#SBATCH --output=%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=0-07:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=8G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=1       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node


module load apps/miniforge 



#######Inputs########

input_dir=<path/to/directory>
#file path to where long-read bam files are stored. Should not end in a /

out_dir=<path/to/directory>/
#file path to suitable output directory location. Should end in a /

settings_file=<path/to/file>
#path to json settings file

####################

mkdir "$out_dir""fastq"

#run through each cram file in the input_directory

echo Starting to convert cram to fastq at `date`


for f in $input_dir/*.bam ; do

	sample_name=${f%-ONT-hg38*}
	sample_name="${sample_name##*/}"

	#reverse complement where needed, splice, and convert cram files to fastq files
	conda activate samtools

        samtools view -b -q 10 --threads 16 $f chr17:51831637-51831739 > "$out_dir""fastq/""$sample_name""_spliced.bam"

        samtools index "$out_dir""fastq/""$sample_name""_spliced.bam"
        conda deactivate

	conda activate jvarkit

	jvarkit pcrclipreads --bed <path/to/bed/containing/regions/to/keep> --samoutputformat BAM "$out_dir""fastq/""$sample_name""_spliced.bam" > "$out_dir""fastq/""$sample_name""_jvarkit_spliced.bam" 
	conda deactivate

	conda activate ngsutils
	bamutils removeclipping -f  "$out_dir""fastq/""$sample_name""_jvarkit_spliced.bam"  "$out_dir""fastq/""$sample_name""_spliced_no_clipped.bam"

	conda activate samtools 
	#add in another samtools view step to remove any reads that have low mapping quality?
	samtools view -b -q 10 --threads 16 "$out_dir""fastq/""$sample_name""_spliced_no_clipped.bam" > "$out_dir""fastq/""$sample_name""_spliced_no_clipped_high_mapq.bam"


	samtools sort "$out_dir""fastq/""$sample_name""_spliced_no_clipped_high_mapq.bam" > "$out_dir""fastq/""$sample_name""_spliced_no_clipped_sorted.bam"
	samtools index "$out_dir""fastq/""$sample_name""_spliced_no_clipped_sorted.bam"
	samtools fastq -f 16  "$out_dir""fastq/""$sample_name""_spliced_no_clipped_sorted.bam" > "$out_dir""fastq/""$sample_name""_spliced_reverse_strand.fastq"
	samtools fastq -F 16  "$out_dir""fastq/""$sample_name""_spliced_no_clipped_sorted.bam" > "$out_dir""fastq/""$sample_name""_spliced_forward_strand.fastq"

	conda deactivate
	conda activate seqtk
	seqtk seq -r "$out_dir""fastq/""$sample_name""_spliced_reverse_strand.fastq" > "$out_dir""fastq/""$sample_name""_spliced_reverse_strand_reverse_complement.fastq"
	conda deactivate
	cat "$out_dir""fastq/""$sample_name""_spliced_reverse_strand_reverse_complement.fastq" "$out_dir""fastq/""$sample_name""_spliced_forward_strand.fastq" > "$out_dir""fastq/""$sample_name""_spliced.fastq"

        rm "$out_dir""fastq/""$sample_name""_spliced_reverse_strand.fastq"
        rm "$out_dir""fastq/""$sample_name""_spliced_forward_strand.fastq"
        rm "$out_dir""fastq/""$sample_name""_spliced_reverse_strand_reverse_complement.fastq"
        rm "$out_dir""fastq/""$sample_name""_jvarkit_spliced.bam"
        rm "$out_dir""fastq/""$sample_name""_spliced.bam"
        rm "$out_dir""fastq/""$sample_name""_spliced.bam.bai"
	rm "$out_dir""fastq/""$sample_name""_spliced_no_clipped.bam"
	rm "$out_dir""fastq/""$sample_name""_spliced_no_clipped_high_mapq.bam"
done


echo Started running RGT at `date`

echo JSON settings file:
cat $settings_file

conda activate rgt

python /mnt/autofs/data/userdata/project0033/ERDA1/tools/RGT/main.py -i "$out_dir""fastq" -o $out_dir -s $settings_file

#convert output to csv file so I can read it in MARS
xlsx2csv "$out_dir""ResultsSummary.xlsx" "$out_dir""ResultsSummary.csv"

echo End of RGT at `date`
conda deactivate
