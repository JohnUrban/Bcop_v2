# Bcop_v2
Repository for supplemental files associated with the Nov 2022 bioRxiv preprint, and subsequent paper, titled, "Chromosome-scale scaffolding of the fungus gnat genome (Diptera: Bradysia 2 coprophila)." Contains lifted over GFF files for genes and repeats, repeat libraries (FASTA), and files mapping Bcop_v1 to Bcop_v2 (BED and AGP).

# See preprint for details:
- https://doi.org/10.1101/2022.11.03.515061
	- version 1, Nov 2022
	- version 2, Dec 2023


# Version 2 of the Bradysia coprophila reference genome:
- Chromosome-scale scaffolds of all four somatic chromosomes (X, II, III, IV).
- Bcop_v2 was deposited to NCBI.
- The Whole Genome Shotgun project has been deposited at DDBJ/ENA/GenBank under the accession VSDI00000000: version VSDI02000000 (deposited 2022; released Jan 4, 2023). 
- See BioProjects PRJNA291918 and PRJNA672144; and BioSamples SAMN12533751 and SAMN20343824. 
- GenBank GCA_014529535.2; WGS VSDI02; release date 2023-01-04.
- Male Pupa Hi-C data was deposited to SRA: SRR23335771. 


# Current directory structure:
```
.
├── README.md
├── bcopv1_bcopv2_map
│   ├── BcopV1_corrected_contigs_on_BcopV2.bed.gz
│   └── Bcop_v2_from_Bcop_v1-ncbi_corrected.agp.gz
├── genes
│   ├── maker-Bcop_v1_to_v2-geneset-liftOff.gff.gz
│   └── ncbi-Bcop_v1_to_v2-geneset-liftOff.gff.gz
└── repeats
    ├── libraries
    │   ├── 01-Bradysia_coprophila_repeatmodeler_families_BcopV1-canu.fasta.gz
    │   └── 02-Bradysia_coprophila_repeatmodeler_families_AlternativeAssembly-falcon.fasta.gz
    ├── repeat-annotations-Bcop_v2.0.fasta.out.gff.gz
    └── repeat-annotations-Bcop_v2.0.fasta.tbl
```

# File explanation
- The "bcopv1_bcopv2_map" directory contains files that map Bcop_v1 to Bcop_v2.
- The "genes" directory contains GFF files of annotations produced on Bcop_v1 lifted over to Bcop_v2 chromosome-scale scaffolds.
- The "repeats" directory contains a GFF file of repeats across Bcop_v2 produced using the "comprehensive repeat library" described in Urban et al (2021)(Bcop_v1 paper).
	- The "libraries" sub-directory contains two B.coprophila-specific repeat libraries de novo modeled on two different assemblies.
		- 01 = Bcop_v1 (produced with Canu).
		- 02 = An assembly produced with the Falcon assembler, also described in Urban et al 2021.
		- Both are added together as parts of the comprehensive repeat library, which also has other components such as repeats from Dfam and RepBase, and previously-known repeat sequences from B. coprophila.


