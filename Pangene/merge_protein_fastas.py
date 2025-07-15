def fasta_iter(handle):
    head = None
    seq = []
    for line in handle:
        line = line.rstrip()
        if line.startswith(">"):
            if head:
                yield head, seq
            head = line[1:]
            seq = []
        else:
            seq.append(line)
    if head:
        yield head, seq

# File paths
file1 = "/n/netscratch/edwards_lab/Lab/kelsielopez/miniprot/HemMar/redo/genes_protein.faa"
file2 = "/n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/hemMar.toga.merged_protein_output.aa"
outfile = "/n/netscratch/edwards_lab/Lab/kelsielopez/HemMar_annotation/toga/combined_headers.faa"

with open(file1) as f1, open(file2) as f2, open(outfile, "w") as out:
    iter1 = fasta_iter(f1)
    iter2 = fasta_iter(f2)
    while True:
        try:
            h1, s1 = next(iter1)
            h2, s2 = next(iter2)
        except StopIteration:
            break
        out.write(f">{h1}:{h2}\n")
        # Join lines to one single string and split to max 60 chars per line
        seq = ''.join(s2)
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")
