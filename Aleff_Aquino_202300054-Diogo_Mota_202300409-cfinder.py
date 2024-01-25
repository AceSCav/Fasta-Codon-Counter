# Usage: python3 Aleff_Aquino_202300054-Diogo_Mota_202300409-cfinder.py file.fasta
import sys
if len(sys.argv) != 2:
    sys.exit(1)
fasta_file = sys.argv[1]
with open(fasta_file, "rt") as fasta_reading:
    from collections import Counter
    dict_line = {}
    current_key = None
    sep_final = []
    codon_common = []
    n = 0
    for line in fasta_reading: #Convert fasta_reading in a dictionary
        line = line.strip()
        if line.startswith(">"):
            current_key = line
            dict_line[current_key] = ""
        else:
            dict_line[current_key] += line
    for separated in dict_line.values(): #Separe all codons in a list "sep_final"
        separated = [separated[i:i+3] for i in range(0, len(separated), 3)]
        sep_final.append(list(separated))
    for check_codon_incomplete in sep_final: #Remove the nucleotides that do not complete codons
        if len(check_codon_incomplete[-1]) !=3:
            check_codon_incomplete.pop()
        else:
            continue
    for codons_sep in sep_final: #Count codons for analysis
        counter_seq = Counter(codons_sep)
        codon_common.append(Counter(codons_sep).most_common(len(Counter(codons_sep))))
    while n < len(dict_line.keys()): #print results
        print(f"Sequence: {list(dict_line.keys())[n][1:]}", file=sys.stderr)
        print(f"Total codons: {len(sep_final[n])}", file=sys.stderr)
        print(f"Most Frequent Codon: {codon_common[n][0][0]} ({round(codon_common[n][0][1]/len(sep_final[n])*100 ,1)}%)", file=sys.stderr)
        print(f"Least Frequent Codon: {codon_common[n][-1][0]} ({round(codon_common[n][-1][1]/len(sep_final[n])*100 ,1)}%)", file=sys.stderr)
        print("")
        n = n + 1






