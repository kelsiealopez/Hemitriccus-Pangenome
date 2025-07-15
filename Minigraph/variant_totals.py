from collections import Counter

filename = "minigraph_sv_results_cleaned_HMRG_ingroup_only.tsv"
counts = Counter()

with open(filename) as f:
    header = f.readline()
    for line in f:
        if line.strip() == "" or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        subtype = fields[3]  # 0-based index: 0,1,2,3 = subtype
        counts[subtype] += 1

for subtype, count in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
    print(f"{subtype}\t{count}")

