# Pangenomes reveal extensive structural variation in a suboscine passerine bird, the Pearly-vented Tody-Tyrant (Hemitriccus margaritaceiventer)

## Overview

This repository contains data files related to the pangenome analysis of Hemitriccus margaritaceiventer genomes. The data includes:

* **PGGB Pangenome Graph Files (pggb_graphs.tar.gz)**: This tar file contains pangenome graphs generated per chromosome using the PanGenome Graph Builder (PGGB).
* **Final PGGB Variation File (`pggb_variation_overlaps_final.tab.gz`)**: A gzipped Variant Call Format (VCF)-like tab file containing variant information decomposed from the pangenome graphs (see supplemental methods in the paper).

### Processed Variant Information (`pggb_variation_overlaps_final.tab.gz`)

#### Column Descriptions
**chrom**: Chromosome identifier for the variant position (with PanSN-spec naming; see methods). 

**bedStart**: Start position of the variant on the chromosome (0-based).

**bedEnd**: End position of the variant on the chromosome (1-based).

**type**: bcftools type, one of snp,mnp,indel,other

**overlap**: Genomic region(s) overlapped by variant

**repeat**: Comma-separated list of repeatmasker annotations overlapped by variant (otherwise none)

**repeat_family**: Comma-separated list of the repeat family of the repeat overlapped by variant (otherwise none)

**subtype**: One of SNP, SV, SVINS, SVDEL, INDEL, DEL, INS, possibly with _Complex suffix. See definitions below.

**ref**: Reference allele nucleotide(s) at the variant position.

**alt**: Alternative allele nucleotide(s) at the variant position.

**aa**: Ancestral allele nucleotide(s) at the variant position.

**inv**: Inversion present? 

**polarized**: Is the variant polarized? False if aa == ".", True otherwise

**base_allele_len**: length of reference allele if aa == ".", otherwise length of aa allele

**alt_len_max**: Longest non-base allele length

**alt_len_min**: Shortest non-base allele length

**allele_count**: Total Allele Count for all samples. 

**HMRG_DAC**: Derived Allele Count for all samples. 

**HMRG_AN**: Number of alleles

**HMRG_MISS**: Missing genotype counts 

**HMRG_6371**: Genotype for HMRG_6371.

**HMRG_6386**: Genotype for HMRG_6386.

**HMRG_6388**: Genotype for HMRG_6388.

**HMRG_6431**: Genotype for HMRG_6431.

**HMRG_6433**: Genotype for HMRG_6433.
