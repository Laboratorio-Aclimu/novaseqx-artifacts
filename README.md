# Recurrent NovaSeq X sequencing artifacts at clinically interpreted positions

Laboratorio Aclimu · Genomics Unit · Buenos Aires, Argentina

Recurrent sequencing artifacts on the Illumina NovaSeq X can appear as variants at ClinVar pathogenic/likely
pathogenic positions and canonical splice sites, with apparent allele fractions of up to 0.42. They are supported by
reads on one strand only, mostly in the lowest base-quality bin, and present in only one of two overlapping mates.
They recur across patients, instruments and centers, and are essentially absent from NovaSeq 6000 data.

**Lookup page:** https://laboratorio-aclimu.github.io/novaseqx-artifacts/ (search by gene, HGVS or position, or
check a VCF in the browser; nothing is uploaded).

## Contents

| File | Description |
|---|---|
| `data/novaseqx_v8_artifacts_v1.1.tsv` | 1,633 alleles that met the artifact criterion in at least 3 of 17 clinical exomes (Agilent SureSelect Human All Exon V8, NovaSeq X), with ClinVar, gnomAD v4.1 and GIAB annotations |
| `data/novaseqx_exome_wide_catalog.tsv` | 77,356 recurrent artifact-like alleles across the capture region (72,026 specific to the NovaSeq X) |
| `nsx_artifact_check.py` | Read-level check of listed alleles in a BAM/CRAM file, with optional VCF annotation |
| `docs/index.html` | Source of the lookup page |

Coordinates are GRCh38 (no-alt analysis set). HGVS descriptions refer to MANE Select v1.5 transcripts.

## Read-level check

Requires Python 3 and pysam.

```bash
python nsx_artifact_check.py --list data/novaseqx_v8_artifacts_v1.1.tsv --bam sample.cram --ref GRCh38.fa \
    --vcf sample.vcf.gz --out sample.nsx.vcf --report sample.nsx.tsv
```

Each listed allele found in the sample is classified as `artifact` (strand-confined support mostly in the lowest
quality bin, or discordant overlapping mates), `likely_real` (support on both strands, at least 3 reads of base
quality 20 or higher and concordant mates) or `review`. Without `--vcf`, all listed positions are scanned, which is
useful when a pipeline genotypes or rescues known pathogenic positions.

## Interpretation

A listed allele is not automatically false: a patient may carry it. The list is meant to prompt review of the reads
and orthogonal confirmation before reporting, not automatic exclusion. Record the instrument control software version
and quality table of each run: between versions 1.2 and 1.3, the lowest informative quality bin changed from 12 to 9,
which places these bases on either side of the default minimum base quality of common variant callers.

## Citation

Laboratorio Aclimu. Recurrent sequencing artifacts at known variant positions: Illumina NovaSeq X + Agilent
SureSelect Human All Exon V8. Zenodo. https://doi.org/10.5281/zenodo.23956232

## License

Code: MIT License (`LICENSE`). Data: Creative Commons Attribution 4.0 International (`LICENSE-DATA.md`).

## Contact

genomica@aclimu.com.ar · https://laboratorioaclimu.com.ar
