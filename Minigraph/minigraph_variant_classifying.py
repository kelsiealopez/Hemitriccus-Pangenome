import re

# ======== CONFIGURE THESE ==========
VCF_FILE = "Minigraph_12haps_polarizing_edited_cleaned.vcf"
HMRG_samples = [
    "HMRG_6371",
    "HMRG_6386",
    "HMRG_6388",
    "HMRG_6431",
    "HMRG_6433"
]
OUTGROUP_SAMPLE = "VEFL_149044"
# =====================================

def parse_alen(info):
    """Parse ALEN values as list of ints"""
    m = re.search(r"ALEN=([^;\t]+)", info)
    if m:
        return [int(x) for x in m.group(1).split(",")]
    else:
        return []

def parse_end(info):
    m = re.search(r"END=([^;\t]+)", info)
    if m:
        return int(m.group(1))
    else:
        return None

def get_sample_indices(header):
    """Get column indices for samples of interest"""
    return {s: header.index(s) for s in HMRG_samples + [OUTGROUP_SAMPLE] if s in header}

def classify_sv(alen, allele_count, ref_index, base_allele_index, polarized):
    # Defensive index: if base_allele_index out of range, use 0 instead
    if base_allele_index >= len(alen) or base_allele_index < 0:
        base_len = alen[0]
    else:
        base_len = alen[base_allele_index]
    alt_lens = [alen[i] for i in range(len(alen)) if i != base_allele_index and i < len(alen)]
    alt_len_max = max(alt_lens) if alt_lens else base_len
    alt_len_min = min(alt_lens) if alt_lens else base_len

    max_bp = max(alt_len_min, alt_len_max, base_len)
    min_bp = min(alt_len_min, alt_len_max, base_len)

    if allele_count == 2:
        if polarized and max_bp >= 50:
            if base_len == 1 and alt_len_max >= 50:
                return "SVINS"
            elif base_len >= 50 and alt_len_max == 1:
                return "SVDEL"
            elif base_len > alt_len_max:
                return "SVDEL_Complex"
            elif base_len < alt_len_max:
                return "SVINS_Complex"
            else:
                return "SV_Complex"
        if polarized and max_bp < 50:
            if base_len == 1 and alt_len_max > 1:
                return "INS"
            elif base_len > 1 and alt_len_max == 1:
                return "DEL"
            elif base_len > alt_len_max:
                return "DEL_Complex"
            elif base_len < alt_len_max:
                return "INS_Complex"
            else:
                return "INDEL_Complex"
        if not polarized and max_bp >= 50:
            return "SV_Complex"
        if not polarized and max_bp < 50:
            return "INDEL_Complex"
    elif allele_count > 2:
        return "SV_Complex"
    return "Complex"

def main():
    with open(VCF_FILE) as in_vcf:
        header_seen = False
        for line in in_vcf:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                idx = get_sample_indices(header)
                # Print column titles as header
                print(
                    "chrom\tpos\tend\tsubtype\tbase_allele_index\tderived\tcalled\tmissing\t" +
                    "\t".join(HMRG_samples)
                )
                header_seen = True
                continue
            if not header_seen:
                continue  # skip until header found

            cols = line.strip().split('\t')
            chrom, pos, _id, ref, alt_field, qual, fltr, info = cols[:8]
            alen = parse_alen(info)
            end = parse_end(info)
            allele_count = len(alen)
            genotypes = cols[9:]  # per-sample columns, after FORMAT

            # Build sample->genotype dict, robust to out-of-bounds/missing
            sample_gts = {}
            for s in idx:
                ix = header.index(s) - 9
                if 0 <= ix < len(genotypes):
                    sample_gts[s] = genotypes[ix]
                else:
                    sample_gts[s] = "."

            # Handle missing/empty ALEN
            if not alen or allele_count == 0:
                print(
                    f"{chrom}\t{pos}\t{end if end is not None else ''}\tNA\tNA\tNA\tNA\tNA\t" +
                    "\t".join(sample_gts.get(s, ".").split(':')[0] for s in HMRG_samples)
                )
                continue

            # Outgroup determination
            out_gt = sample_gts.get(OUTGROUP_SAMPLE, ".").split(":")[0]
            base_allele_index = None
            polarized = False
            if out_gt and out_gt not in [".", "./.", ".|."]:
                gt_alleles = re.split("/|\|", out_gt)
                for x in gt_alleles:
                    if x != ".":
                        try:
                            idx_val = int(x)
                            if idx_val < len(alen) and idx_val >= 0:
                                base_allele_index = idx_val
                                polarized = True
                                break
                            else:
                                # Outgroup's allele not present among ALENs
                                # Optionally print warning:
                                # print(f"WARNING: base_allele_index {idx_val} out of range ({alen}) at {chrom}:{pos}; using 0")
                                base_allele_index = 0
                                polarized = False
                                break
                        except Exception:
                            # Genotype is not an integer (shouldn't happen, but be safe)
                            # print(f"WARNING: Could not parse outgroup genotype '{x}' at {chrom}:{pos}, using 0")
                            base_allele_index = 0
                            polarized = False
                            break
            if base_allele_index is None:
                base_allele_index = 0  # fallback to REF

            if base_allele_index >= len(alen) or base_allele_index < 0:
                # Defensive, just in case
                # print(f"WARNING: base_allele_index {base_allele_index} out of range ({alen}) at {chrom}:{pos}; using 0")
                base_allele_index = 0
                polarized = False

            # Tally derived/called/missing for HMRG
            derived = 0
            called = 0
            missing = 0
            for s in HMRG_samples:
                gt = sample_gts.get(s, ".").split(":")[0]
                if not gt or gt in [".", "./.", ".|."]:
                    missing += 1
                    continue
                for x in re.split("/|\|", gt):
                    if x == ".":
                        missing += 1
                        continue
                    try:
                        if int(x) != base_allele_index:
                            derived += 1
                        called += 1
                    except Exception:
                        missing += 1

            subtype = classify_sv(alen, allele_count, 0, base_allele_index, polarized)

            print(
                f"{chrom}\t{pos}\t{end if end is not None else ''}\t{subtype}\t{base_allele_index}\t{derived}\t{called}\t{missing}\t"
                + "\t".join(sample_gts.get(s, ".").split(':')[0] for s in HMRG_samples)
            )


if __name__ == "__main__":
    main()
