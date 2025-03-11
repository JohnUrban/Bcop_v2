bedtools makewindows -g ~/base/data/sciara/phase/Bcop_v2.0/Bcop_v2.0_pkg/Bcop_v2.0.fasta.genome -w 1000 -s 1000 > ../windows-1kb.bed
bedtools makewindows -g ~/base/data/sciara/phase/Bcop_v2.0/Bcop_v2.0_pkg/Bcop_v2.0.fasta.genome -w 5000 -s 5000 > ../windows-5kb.bed
bedtools makewindows -g ~/base/data/sciara/phase/Bcop_v2.0/Bcop_v2.0_pkg/Bcop_v2.0.fasta.genome -w 10000 -s 10000 > ../windows-10kb.bed
bedtools makewindows -g ~/base/data/sciara/phase/Bcop_v2.0/Bcop_v2.0_pkg/Bcop_v2.0.fasta.genome -w 100000 -s 100000 > ../windows-100kb.bed
bedtools makewindows -g ~/base/data/sciara/phase/Bcop_v2.0/Bcop_v2.0_pkg/Bcop_v2.0.fasta.genome -w 500 -s 500 > ../windows-500bp.bed
