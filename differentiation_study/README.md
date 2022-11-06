# Differentiation analysis on single cell RNA-seq data  <br> 
[1] to prepare data using velocyto for running scvelo <br> 
[2] to analyse data using slingshot, monocle <br> 
In Nov 2022

<br>
<br>

##[1] to prepare data using velocyto for running scvelo <br>
from *.BAM to *.loom based on mice data 

### 1. To use velocyto  <br> 


	## To install 
	## -------------------------------- ##
	pip install velocyto -u
	conda install -c conda-forge loompy
	
		
	## In command line 
	## -------------------------------- ##
	#!/usr/../bin/zsh
	#SBATCH -J velocyto
	#SBATCH -t 100:00:00
	#SBATCH --mem-per-cpu=102400M
	#SBATCH --output=output.%J.txt
	
	/home/.../.pyenv/shims/velocyto run10x -m /../mm10_repeats_repeatMasker_allTracks.gtf /../SAMPLE_NAME /../refdata-cellranger-mm10-3.0.0/genes/genes.gtf --samtools-threads 20
	
	## "SAMPLE_NAME" is the folder which has the outputs from CellRanger, for example, _cmdline, _finalstate, outs, ...   
	## "*.loom" will be made in a new directory "velocyto" inside the folder "SAMPLE_NAME"	
	## For the mm10_repeats_repeatMasker_allTracks.gtf, I refer to this https://www.biostars.org/p/416125/


### 2. To combine loom  <br> 


	## In "combine.loom.py"
	## -------------------------------- ##
	import os,sys
	import loompy
	loompy.combine(sys.argv[1:-1], sys.argv[-1])
	
	
	## In "combine.loom.sh"
	## -------------------------------- ##
	#!/usr/../bin/zsh
	#SBATCH -J velocyto
	#SBATCH -t 100:00:00
	#SBATCH --mem-per-cpu=102400M
	#SBATCH --output=output.%J.txt
	
	file1=/.../SAMPLE_NAME1/velocyto/SAMPLE_NAME1.loom
	file2=/.../SAMPLE_NAME2/velocyto/SAMPLE_NAME2.loom
	file3=/.../SAMPLE_NAME3/velocyto/SAMPLE_NAME3.loom
	file4=/.../SAMPLE_NAME4/velocyto/SAMPLE_NAME4.loom
	
	python combine.loom.py $file1 $file2 $file3 $file4 ./ALL.SAMPLE.merged.loom
	


### 3. To run  <br>


	## In Slurm 
	## -------------------------------- ##
	sbatch combine.loom.sh



## [2] to analyse data using slingshot, monocle <br> 
./code/phate_slingshot_monocle.R <br>
https://github.com/genehaus/singlecellRNA-seq_2022/blob/main/differentiation_study/code/phate_slingshot_monocle.R




