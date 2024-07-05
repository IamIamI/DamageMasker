<p align="center">
  <img src="https://github.com/IamIamI/DamageMasker/blob/main/Images/Damage_Masker_logo.jpg" width="600"/>
</p>

***```DamageMasker``` is a lightweight python based tool usefull for masking damaged sites in ancient SAM/BAM files, to prevent it from being incorporated into genotyping.***

## Overview

```DamageMasker``` is a small tool written in Python3, using the Pysam and BioPython libraries.
It's designed to mask sequencing artefacts derived from deaminated cytosines, a common feature of DNA damage.
This problem is especially prevalent in ancient DNA, and mostly affects the 5' edges of DNA molecules.
Deaminated cytocines function similarly to a uracil, and during the library and amplification process,
these deaminated cytocines can be complimented with a adenine instead of a guanine. 


<p align="center">
  <img src="https://github.com/IamIamI/DamageMasker/blob/main/Images/Deamination.jpg" width="800"/>
</p>
Figure 1. Depiction of the cause and effect of deamination of a cytosine. A) A depiction of a small portion of a DNA molecule. B) Structural illustration of the two base pairs, C-G and T-A and their respective hydrogen bonds. C) A step-by-step illustration of what deamination entails and how it can result in mutations or in the case of aDNA, misincorporations by polymerases. D) The process by which an uracil-DNA glycosylase (UDG) enzyme can remove deaminated cytosines (uracil) from the DNA backbone, which in turn results in a cleavage of the strand when amplification or repair of the molecule occurs. Adapted from [3D DNA model by Holoxica hosted by sketchfab](https://skfb.ly/68N7T) under a CC license.

This process leads to the replacement of cytosines to thymines in a subset of reads through a adenine complementation. 
When mapping these C>T (or the G>A reverse compliment) reads to a reference genome, these can cause mismapping (loss of mappable reads), 
wrong genotype calls (Ts where there should have been Cs), or loss of stable genotypes due to hetrogeneity during SNP calling.
```DamageMasker``` is intended as a soft sollution to this problem if UDG treatment is not used during library creation. 


## Instalation

This tool is written for Python 3 and uses the biopython and pysam libraries.
biopython and pysam can be easily by running ```pip install biopython``` (or ```pip3 install biopython```) and ```pip install pysam``` (or```pip3 install pysam```). 

The ```DamageMasker``` can be simply cloned using git ```git clone https://github.com/IamIamI/DamageMasker.git```.
Or the python script can be downloaded directly using wget by typing ```wget https://raw.githubusercontent.com/IamIamI/DamageMasker/main/DamageMasker.py```

conda is WIP

## Features

```DamageMasker``` offers the following mapping strategies:

### Masking:
Hardmasking (```-m H```): 	Mask all Ts on forward mapped reads, all A's on reverse mapped reads

Edgemasking (```-m E```): 	Mask all Ts on the 5' edge of the forward reads, and A's on the 5' edge of the reverse reads. 
> [!TIP]
> The user can set how many nucleotides into the read will be masked by setting a value with options ```-e``` or ```--edge_count```

<br><br/>
### Reference guidance:
Allows reference guidance by supplying a path to the reference file used in the mapping by supplying the path to it with the ```-r``` or ```--ref_file``` option.\
For example: ```--ref_file genome.fasta```.

<br><br/>
### Library support:
Single stranded (```-s S```):		This will assume damage presents itself as Ts on forward mapped reads, and As on Reverse mapped reads. 		
> [!CAUTION]
> When combining this with Hardmasking, expect a 25% data loss (all Ts on Forward, all As on reverse).


Double stranded (```-s D```):		This will assume damage presents itself as both Ts on 5' side and As on 3' sides forward mapped reads, and As on the 5' side and Ts on the 3' side of reverse mapped reads. 
> [!CAUTION]
> When combining this with Hardmasking, expect a complete loss of As and Ts, and should likely not be used.

<br><br/>
### Read Filtering:
An additional feature is present which can filter reads that are either too short by setting a minimum length (in bp) cuttoff using the option ```-l``` or ```--len_cutoff``` followed by a value,
or remove reads from the output that have too low of a MapQ score using the option ```-q``` or ```--mapq_cutoff``` followed by a value.\
For example: ```--mapq_cutoff 20 --len_cutoff 35```

<br><br/>
The sofware has an overview of all options which can be called upon by typing 'python ```python DamageMasker.py -h``` or ```python DamageMasker.py --help```.  
  
An overview of the options are as followed:  
```
options:
  -m , --masking       Change masking behaviour. 'H' for HardMasking. 'E' for EdgeMasking. 'F' for only Filtering. (default: Hardmasking)
  -i , --input_file    Input BAM or SAM file (mandatory)
  -s', '--strandness   Determine strandness of dataset, 'S' for single stranded libraries, and 'D' for double stranded libraries (default: S, for sslib)
  -e , --edge_count    Number of 5' edges to be masked if --masking 'E' is turned on (default: 5)
  -r , --ref_file      Give the path to a reference genome file if you want to turn on reference guidance (default: turned off)
  -q , --mapq_cutoff   MAPQ cutoff value (default: 0)
  -l , --len_cutoff    Ignore reads below a certain length (default: 0)
  -o , --output_file   Output SAM file with modified reads (default: 'output_modified.sam/.bam')
  -h, --help           Show this help message and exit.
```

Most options can be combined in order to achieve the users goal. 

<p align="center">
  <img src="https://github.com/IamIamI/DamageMasker/blob/main/Images/DamageMasker_examples.jpg" width="600"/>
</p>
Figure 3. An overview of the impact different combination of settings can have when running ```DamageMasker```

## Expected results

An overview of the expected results of each method
<p align="center">
  <img src="https://github.com/IamIamI/DamageMasker/blob/main/Images/Examples_Masking_impact.jpg" width="600"/>
</p>
Figure 4. Example results that masking can yield depending on the situation of a given position. High coverage data will benefit more from harsher mapping strategies, while lower coverage data will benefit more from edge masking in combination with reference guidance. However, outlier will stil exist and legitimate SNPs can still be lost. 

The resulting Damage patterns that are to be expected from each masking type can be seen below
<p align="center">
  <img src="https://github.com/IamIamI/DamageMasker/blob/main/Images/DamagePlots.jpg" width="600"/>
</p>
Figure 5. Example damage patterns that can be expected when using ```DamageMasker```.

## Example parameters

- Use case: When setting to hardmasking, al T's on the forward strand and all A's on the reverse strand are masked regardless of anything else.

Example: ```python DamageMasker.py --input_file Sample.processed.bam --masking H --output_file Sample_Hardmasked.bam```  
  
- Use case: When adding reference guidance, all T's on the forward strand are masked if the reference has a C, and all A's on the reverse strand will be masked if the reference has a G on that position.
 
Example: ```python DamageMasker.py --input_file Sample.processed.bam --masking H --ref_file Reference.fasta --output_file Sample_RefGuidedmasked.bam```

- Use case: When using edge masking, only T's on the 5' and 3' of the forward read are masked, and only A's on 5' and 3' of the reverse strand will be masked, the user can specify how many bases into these edges the masking runs.

Example: ```python DamageMasker.py --input_file Sample.processed.bam --masking E --edge_count 3 --output_file Sample_Edgemasked.bam```  

- Use case: When using double stranded libraries, both Ts and As damage can occur as a result of complentation of damage of both strands of the molecule. When using Edge masking this means that not only the 5' side of the read is masked for Ts but also the 3' side is masked for As. This can be further reduced by applying reference guidandance.
 
Example: ```python DamageMasker.py --input_file Sample.processed.bam --strandness S --masking E --edge_count 3 --ref_file Reference.fasta --output_file Sample_RefGuidedmasked.bam```

- Use case: The user can also remove reads that are too short (default 0bp), or are not mapping with a high enough MapQ score (default 0).

Example: ```python DamageMasker.py --input_file Sample.processed.bam --output_file Sample_Hardmasked_Filtered.bam --masking F --mapq_cutoff 37 --len_cutoff 25```

- Use case: Based on testing these are some of the most clean vs conservative settings.

Example: ```python DamageMasker.py --input_file Sample.processed.bam --strandness S --masking E --edge_count 5 --ref_file Reference.fasta --output_file Sample_RefGuidedmasked.bam --mapq_cutoff 37 --len_cutoff 30```


