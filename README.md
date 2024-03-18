# DamageMasker
DamageMasker is written for Python 3 and uses the biopython and pysam libraries which can easily be installed by running ```pip install biopython``` and ```pip install pysam``` or using a conda environment. 
If Python 3 is installed under python 3 then use ```pip3 install biopython``` and ```pip3 install pysam```.

This tool is intended to be used on singel stranded (SS) libraries that have ancient "damage" (deaminated cytosines), and are not UDG treated.
This tool can mask forward strand 'T' and reverse strand 'A' in Sam/Bam files. The intended goal is to mask damage from the reads and prevent it from showing up as biological genotypes or prevent real biological genotypes from being called due to hetrozygocity.  
  
Since the single stranded libraries are not synthetically amplified yet, they are assumed to not have artefactually complemented 'C'>'T'>'A' changes, and instead only have natural deamination artefacts. This means it's possible to mask elements that appear as damage based on the strand that exhibits it.  
  
There are three approaches, the most logical but also harshest is to mask all the 'T' on the forward strand and the 'A' on the reverse strand as we cannot determine which 'T' are biological and which artefactual, we call this 'Hard Masking'.  

The second approach is what we call 'Edge Masking', where only a fixed number of nucleotides from the edges of reads are hardmasked (default is 5bp). Since deamination apears to mostly happen on the first couple of bases on the 5' and 3', we can mask all the 'T' on the forwards strand if they are within the first and last couple of nucleotides, and 'A' if they are within the first and last couple of nucleotides of the reverse strand. The user can control this value, and this approach is maybe the most robust. One could use a damage profile as a guide to determine what value to set. 

The last approach is a reference guide approach where we also supply the reference genome to which the reads were mapped, and we only mask 'T' on the forward strand if the reference sequence has a 'C', and mask the 'A' on the reverse strand if the reference sequence has a 'G'. Although here we assume to know which nucleotides are damage, and will also delete biological C>T, the forward strand 'A' and reverse strand 'T' will still be able to corroborate the genotype if it's biological. Although we do expect that 'C'>'T'/'G'>'A' sites will drop substatially in coverage compared to all other sides, and with low coverage samples might result in a bias during SNP calling. When generating phylogenies we have seen some strange behaviour with this method, and it might be more suited for higher coverage genomes with substantial damage. 

Masking just replaces the nucleotide with an 'N' which does not affect GATK's UnifiedGenotyper when calling genotypes, but might cause problems in other software packages, so use this software at your own discression. 
  
This script has error handling for most thinkable scenario's and should be relatively easy to run as followed:
- Use case: When setting to hardmasking, al T's on the forward strand and all A's on the reverse strand are masked regardless of anything else  
Example: ```python DamageMasker.py --input_file Sample.processed.bam --masking R --ref_file Reference.fasta --output_file Sample_RefGuidedmasked.bam```  
  
- Use case: When setting the script to Reference guided masking, all T's on the forward strand are masked if the reference has a C, and all A's on the reverse strand will be masked if the reference has a G on that position  
Example: ```python DamageMasker.py --input_file Sample.processed.bam --masking H --output_file Sample_Hardmasked.bam```  
  
- Use case: When using edge masking, only T's on the 5' and 3' of the forward read are masked, and only A's on 5' and 3' of the reverse strand will be masked, the user can specify how many bases into these edges the masking runs.  
Example: ```python DamageMasker.py --input_file Sample.processed.bam --masking E --edge_count 3 --output_file Sample_Edgemasked.bam```  
  
- Use case: The user can also remove reads that are too short (default 0bp), or are not mapping with a high enough MapQ score (default 0)  
Example: ```python DamageMasker.py --input_file Sample.processed.bam --output_file Sample_Hardmasked_Filtered.bam --mapq_cutoff 37 --len_cutoff 25```  
  
The sofware has an overview of all options which can be called upon by typing 'python SSLib_Masker.py -h' or 'python SSLib_Masker.py --help'  
  
An overview of the options are as followed:  
```
options:
  -m , --masking       Change masking behaviour. 'R' for Reference based Masking. 'H' for HardMasking. 'E' for EdgeMasking. (default: Hardmasking)
  -i , --input_file    Input BAM or SAM file (mandatory)
  -r , --ref_file      Input reference genome file in FASTA format (mandatory for --masking 'R')
  -e , --edge_count    Number of 5' edges to be masked if --masking 'E' is turned on (default: 5)
  -q , --mapq_cutoff   MAPQ cutoff value (default: 0)
  -l , --len_cutoff    Ignore reads below a certain length (default: 0)
  -o , --output_file   Output SAM file with modified reads (default: 'output_modified.sam')
  -h, --help           Show this help message and exit.
```

An overview of the expected results of each method
![alt text](https://github.com/IamIamI/Bioinformatics_scripts/blob/master/Python/ssLib_masker.jpg)
