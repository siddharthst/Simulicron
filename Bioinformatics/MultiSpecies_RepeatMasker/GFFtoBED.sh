for f in ./RepeatMasker/*.fasta.out.gff; do gff2bed < "$f" > "${f%.*}".bed; done
