<p align="center">
  <img src="https://github.com/IamIamI/DamageMasker/blob/main/Damage_Masker_logo.jpg" width="600"/>
</p>

```DamageMasker``` is intended to be used to mask thymine and adenine sites in ancient SAM/BAM files, to prevent damage from being incorporated into genotyping.

Offers the following mapping strategies:
		Hardmasking (option -m 'H'): mask all Ts on forward mapped reads, all A's on reverse mapped reads
		Edgemasking (option -m 'E'): mask all Ts on the 5' edge of the forward reads, and A's on the 5' edge of the reverse reads. The user can set how many nucleotides into the read will be masked by setting a value with options -e or --edge_count
  
Allows reference guidance by supplying a path to a reference file using the -r or --ref_file option.

Allows masking of single stranded libraries (only mask on the 5') and double stranded libraries (mask both from the 3' and 5' of the read) using the option -s or --strandness. 
 	Single stranded libraries will only trim T's on the forward read, and A's on the reverse read. Hardmasking with double stranded libraries will result in ~25% data loss. Edgemasking mitigates this loss substantially.
 	Double stranded libraries will trim both T's and A's on both forward and reverse reads. Hardmasking with double stranded libraries will therefore result in ~50% data loss. Edgemasking mitigates this loss substantially.

This tool is written for Python 3 and uses the biopython and pysam libraries which can easily be installed by running ```pip install biopython``` (or ```pip3 install biopython```) and ```pip install pysam``` (or```pip3 install pysam```) or using a conda environment as below. 
  
There are three approaches, the most logical but also harshest is to mask all the 'T' on the forward strand and the 'A' on the reverse strand as we cannot determine which 'T' are biological and which artefactual, we call this 'Hard Masking'.  

The second approach is what we call 'Edge Masking', where only a fixed number of nucleotides from the edges of reads are hardmasked (default is 5bp). Since deamination apears to mostly happen on the first couple of bases on the 5' and 3', we can mask all the 'T' on the forwards strand if they are within the first and last couple of nucleotides, and 'A' if they are within the first and last couple of nucleotides of the reverse strand. The user can control this value, and this approach is maybe the most robust. One could use a damage profile as a guide to determine what value to set. 

The last approach is a reference guide approach where we also supply the reference genome to which the reads were mapped, and we only mask 'T' on the forward strand if the reference sequence has a 'C', and mask the 'A' on the reverse strand if the reference sequence has a 'G'. Although here we assume to know which nucleotides are damage, and will also delete biological C>T, the forward strand 'A' and reverse strand 'T' will still be able to corroborate the genotype if it's biological. Although we do expect that 'C'>'T'/'G'>'A' sites will drop substatially in coverage compared to all other sides, and with low coverage samples might result in a bias during SNP calling. When generating phylogenies we have seen some strange behaviour with this method, and it might be more suited for higher coverage genomes with substantial damage. 

Additionally, an option is added to only apply the length and mapq filtering parameters. Although this is likely not usefull for the audience interested in damage, it can be usefull to just reduce the overall sam/bam files by removing unmapped reads and those that are too short or have low quality mapping.

Masking just replaces the nucleotide with an 'N' which does not affect GATK's UnifiedGenotyper when calling genotypes, but might cause problems in other untested software packages, so use this software at your own discression. 
  
This script has error handling for most thinkable scenario's and should be relatively easy to run as followed:
- Use case: When setting to hardmasking, al T's on the forward strand and all A's on the reverse strand are masked regardless of anything else  
Example: ```python DamageMasker.py --input_file Sample.processed.bam --masking H --output_file Sample_Hardmasked.bam```  
  
- Use case: When setting the script to Reference guided masking, all T's on the forward strand are masked if the reference has a C, and all A's on the reverse strand will be masked if the reference has a G on that position  
Example: ```python DamageMasker.py --input_file Sample.processed.bam --masking R --ref_file Reference.fasta --output_file Sample_RefGuidedmasked.bam```

- Use case: When using edge masking, only T's on the 5' and 3' of the forward read are masked, and only A's on 5' and 3' of the reverse strand will be masked, the user can specify how many bases into these edges the masking runs.  
Example: ```python DamageMasker.py --input_file Sample.processed.bam --masking E --edge_count 3 --output_file Sample_Edgemasked.bam```  
  
- Use case: The user can also remove reads that are too short (default 0bp), or are not mapping with a high enough MapQ score (default 0)  
Example: ```python DamageMasker.py --input_file Sample.processed.bam --output_file Sample_Hardmasked_Filtered.bam --masking F --mapq_cutoff 37 --len_cutoff 25```

- Use case: Damage masking and short read removal etc can be easily combined 
Example: ```python DamageMasker.py --input_file Sample.processed.bam --output_file Sample_Hardmasked_Filtered.bam -masking E --edge_count 3 --mapq_cutoff 37 --len_cutoff 25```

Since the single stranded libraries are not synthetically amplified yet, they are assumed to not have artefactually complemented 'C'>'T'>'A' changes, and instead only have natural deamination artefacts. This means it's possible to mask elements that appear as damage based on the strand that exhibits it. Hardmasking using the single stranded library approach is estimated to result in ~25% data loss, this can be reduced by using the Edgemasking approach instead or taking advantage of the reference guided masking.

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

An overview of the expected results of each method
![.](https://github.com/IamIamI/DamageMasker/blob/main/DamageMasker.png)
