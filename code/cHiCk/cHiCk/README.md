# cHiCk
- Various utilities while working with Hi-C data, bedpe, and bed files.
- What does "cHiCk" stand for, you ask? How dare you constrain me to acronyms. This is just a word that had "HiC" in it! And if you must force the "c" and "k" to mean something, then let it make the sound "seq" when brought together ("ck").
- If you really need meaning for "cHiCk", then chicks are baby chickens that need to figure out how to get out of eggs by pecking, and here you're trying to figure out the meaning of "cHiCk" by pecking at me. Get out of here you!

# Dependencies
- BEDtools
- SAMtools
- Ken utilities
- Awk
- Cpp / g++

# Applications
- Generate interaction counts in bins that emanate a way from a locus.
- Example loci might be:
	- a specific gene: e.g. comparing the interaction frequencies of a predicted HGT gene to a nearby highly conserved fly gene of similar length.
	- a specific contig end within a scaffold to visualize its Hi-C scaffold evidence in something like IGV.
	- a specific end of a specific contig to see the evidence for its attachment to the immediately neighboring contig end.


# Getting input files
- First need to run the "Phase Map" pipeline for mapping Hi-C reads, found in a parallel directory under ~/Bcop_v2/code/.
- Will be using the BEDPE output.


# Examples from Bcop_v2 paper
- Tested so-called "alien" p450 genes, compared to highly conserved fly genes.
- Tested contigs in scaffolds.
- Also tested TPS genes as part of another paper.
