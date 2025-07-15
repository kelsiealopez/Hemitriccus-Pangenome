def process_vcf_line(line):
    parts = line.strip().split("\t")
    
    chrom, pos, _id, ref, alts, qual, filt, info, fmt, *samples = parts

    # Modify the samples to merge haplotypes into phased genotype format
    transformed_samples = []
    
    for idx in range(0, len(samples), 2):
        hap1 = samples[idx]
        hap2 = samples[idx + 1] if idx + 1 < len(samples) else "."
        
        # Convert the genotype from e.g., "0:0" to "0" and combine them with "|"
        gt1 = hap1.split(":")[0]
        gt2 = hap2.split(":")[0]
        
        phased_genotype = f"{gt1}|{gt2}"
        transformed_samples.append(phased_genotype)
    
    # Re-assemble the transformed line
    transformed_line = "\t".join([chrom, pos, _id, ref, alts, qual, filt, info, fmt] + transformed_samples)
    return transformed_line

def transform_vcf(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith("##"):
                # Keep metadata headers
                outfile.write(line)
            elif line.startswith("#CHROM"):
                # Transform the header line to consolidate haplotype columns
                header_parts = line.strip().split("\t")
                individuals = [header_parts[i] for i in range(9, len(header_parts), 2)]  # Every second sample header
                new_header = header_parts[:9] + individuals
                outfile.write("\t".join(new_header) + "\n")
            else:
                # Process and join haplotype data for each individual
                transformed_line = process_vcf_line(line)
                outfile.write(transformed_line + "\n")

# Specify your input and output files
input_vcf = "Minigraph_12haps_polarizing.vcf"
output_vcf = "Minigraph_12haps_polarizing_edited.vcf"
transform_vcf(input_vcf, output_vcf)
