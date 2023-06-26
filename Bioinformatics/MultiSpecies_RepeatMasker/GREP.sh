for file in ./Shortlist/*.bed
do
  grep -f shortlisted_TEs.txt "$file" > ./"$file".bed2
done
