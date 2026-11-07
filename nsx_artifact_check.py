#!/usr/bin/env python3
"""nsx-artifact-check: flags recurrent NovaSeq X sequencing artifacts at clinically interpreted positions.

Given the per-allele list (Zenodo 10.5281/zenodo.23956233) and the BAM/CRAM of one sample, it inspects the reads
supporting each listed allele and classifies it as:
  artifact     strand-confined support, mostly in the lowest informative quality bin of the run, and/or present in
               only one of two overlapping mates;
  likely_real  support on both strands, at least 3 reads of base quality >= 20 and concordant mates;
  review       anything else with at least 3 supporting reads.
With --vcf, listed alleles present in the VCF are annotated (INFO/NSX_STATUS, INFO/NSX_EVIDENCE) and every other
record is copied unchanged. Without --vcf, all listed positions are scanned (useful when a pipeline genotypes or
rescues known pathogenic positions). Requires pysam.

Usage:
  nsx_artifact_check.py --list list.tsv --bam sample.cram --ref GRCh38.fa [--vcf in.vcf.gz --out out.vcf] \
      --report report.tsv
"""
import argparse
import collections
import csv
import sys

import pysam


def lowest_bin(bam, n=200000):
    """Lowest informative quality value (>2) used by the run, from the first reads of the file."""
    q = collections.Counter()
    for i, a in enumerate(bam.fetch(until_eof=True)):
        if i >= n:
            break
        if a.query_qualities is not None:
            q.update(a.query_qualities)
    tot = sum(q.values())
    return min(k for k, v in q.items() if k > 2 and v / tot >= 0.001)


def evidence(bam, chrom, pos, alt, low):
    """Read-level evidence for one allele at a 1-based position."""
    alt_f = alt_r = alt_low = alt_hq = ref_n = 0
    frag = collections.defaultdict(list)
    for col in bam.pileup(chrom, pos - 1, pos, truncate=True, min_base_quality=0, min_mapping_quality=20,
                          stepper='nofilter', ignore_overlaps=False):
        for r in col.pileups:
            a = r.alignment
            if (r.is_del or r.is_refskip or r.query_position is None or a.is_duplicate or a.is_secondary
                    or a.is_supplementary or a.is_qcfail):
                continue
            b = a.query_sequence[r.query_position]
            frag[a.query_name].append(b)
            if b == alt:
                bq = a.query_qualities[r.query_position]
                alt_r += a.is_reverse; alt_f += not a.is_reverse
                alt_low += bq <= low; alt_hq += bq >= 20
            else:
                ref_n += 1
    con_alt = [l for l in frag.values() if alt in l]
    solapados = [l for l in con_alt if len(l) >= 2]
    discordantes = sum(1 for l in solapados if any(x != alt for x in l))
    n = alt_f + alt_r
    return {'alt_reads': n, 'fwd': alt_f, 'rev': alt_r, 'low_bin_frac': round(alt_low / n, 3) if n else 0,
            'hq_reads': alt_hq, 'depth': n + ref_n, 'overlap_fragments': len(solapados),
            'discordant_mates': discordantes}


def classify(e):
    n = e['alt_reads']
    if n < 3:
        return 'no_support'
    minority = min(e['fwd'], e['rev']) / n
    mate_discord = e['overlap_fragments'] >= 3 and e['discordant_mates'] / e['overlap_fragments'] >= 0.9
    if (minority <= 0.1 and e['low_bin_frac'] >= 0.8 and e['hq_reads'] <= 2) or (minority <= 0.1 and mate_discord):
        return 'artifact'
    if minority > 0.1 and e['hq_reads'] >= 3 and not mate_discord:
        return 'likely_real'
    return 'review'


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--list', required=True); ap.add_argument('--bam', required=True)
    ap.add_argument('--ref', required=True); ap.add_argument('--vcf'); ap.add_argument('--out')
    ap.add_argument('--report', required=True)
    a = ap.parse_args(argv)
    listed = {(r['chrom'], int(r['pos']), r['ref'], r['alt']): r for r in
              csv.DictReader(open(a.list), delimiter='\t')}
    bam = pysam.AlignmentFile(a.bam, reference_filename=a.ref)
    low = lowest_bin(pysam.AlignmentFile(a.bam, reference_filename=a.ref))
    rep = csv.writer(open(a.report, 'w'), delimiter='\t')
    rep.writerow(['chrom', 'pos', 'ref', 'alt', 'gene', 'hgvs_c', 'status', 'alt_reads', 'fwd', 'rev', 'low_bin_frac',
                  'hq_reads', 'depth', 'overlap_fragments', 'discordant_mates', 'lowest_bin', 'in_vcf'])

    def informe(k, en_vcf):
        e = evidence(bam, k[0], k[1], k[3], low)
        st = classify(e)
        f = listed[k]
        rep.writerow([*k, f.get('gene', ''), f.get('hgvs_c', ''), st, e['alt_reads'], e['fwd'], e['rev'],
                      e['low_bin_frac'], e['hq_reads'], e['depth'], e['overlap_fragments'], e['discordant_mates'],
                      low, int(en_vcf)])
        return st, e

    if a.vcf:
        vin = pysam.VariantFile(a.vcf)
        vin.header.add_line('##INFO=<ID=NSX_STATUS,Number=A,Type=String,Description="nsx-artifact-check status of '
                            'listed NovaSeq X artifact alleles (artifact, likely_real, review, no_support; . if not '
                            'listed)">')
        vin.header.add_line('##INFO=<ID=NSX_EVIDENCE,Number=A,Type=String,Description="alt reads|fwd|rev|low-bin '
                            'fraction|reads with BQ>=20|overlapping fragments|discordant mates">')
        vout = pysam.VariantFile(a.out or '-', 'w', header=vin.header)
        for r in vin:
            sts, evs = [], []
            for alt in r.alts or ():
                k = (r.chrom, r.pos, r.ref, alt)
                if k in listed:
                    st, e = informe(k, True)
                    sts.append(st)
                    evs.append('|'.join(str(e[x]) for x in ('alt_reads', 'fwd', 'rev', 'low_bin_frac', 'hq_reads',
                                                              'overlap_fragments', 'discordant_mates')))
                else:
                    sts.append('.'); evs.append('.')
            if any(s != '.' for s in sts):
                r.info['NSX_STATUS'] = sts; r.info['NSX_EVIDENCE'] = evs
            vout.write(r)
        vout.close()
    else:
        for k in sorted(listed):
            informe(k, False)
    return 0


if __name__ == '__main__':
    sys.exit(main())
