import os
import re
import sys
from collections import OrderedDict
from .RunCmdsMP import run_cmd, logger

BLASType = {
    'qseqid': str,
    'sseqid': str,
    'pident': float,
    'length': int,
    'mismatch': int,
    'gapopen': int,
    'qstart': int,
    'qend': int,
    'sstart': int,
    'send': int,
    'evalue': float,
    'bitscore': float,   # surrogate (matches) when using minimap2
    'qlen': int,
    'slen': int,
    'qcovs': float,
    'qcovhsp': float,
    'sstrand': str,
}

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

def _parse_cigar_stats(cigar):
    """Return (ins_bases, del_bases, gap_opens) from a CIGAR string."""
    if not cigar:
        return 0, 0, 0
    ins_bases = del_bases = gap_opens = 0
    prev_is_gap = False
    for length, op in CIGAR_RE.findall(cigar):
        l = int(length)
        if op == 'I':
            ins_bases += l
            gap_opens += 1
        elif op == 'D':
            del_bases += l
            gap_opens += 1
        # treat separate I/D runs as separate gap opens (BLAST-like)
    return ins_bases, del_bases, gap_opens

def minimap(db_seq, qry_seq, preset='asm20', idx_name=None, paf_out=None, mm_opts='', ncpu=64):
    """
    Map qry_seq (FASTA) to db_seq (FASTA) using minimap2 and write PAF.
    Returns path to PAF.
    """
    if idx_name is None:
        idx_name = os.path.splitext(os.path.basename(db_seq))[0] + '.mmi'
    if paf_out is None:
        paf_out = qry_seq + '.paf'

    # index targets (fast; safe to redo)
    cmd = f"minimap2 -x {preset} -d {idx_name} {db_seq}"
    run_cmd(cmd, logger=logger, fail_exit=True)

    # map: best hit only, include CIGAR and cs for downstream parsing
    # -c adds cg:Z: CIGAR in PAF; --cs=short adds cs:Z: diff string (optional)
    cmd = (
        f"minimap2 -x {preset} -t {ncpu} --secondary=no -N 1 "
        f"-c --cs=short -K 1G -2 {idx_name} {qry_seq} > {paf_out}"
    )
    if mm_opts:
        # allow extra user-provided options (placed before index/query)
        cmd = f"minimap2 -x {preset} -t {ncpu} --secondary=no -N 1 -c --cs=short -K 1G -2 {mm_opts} {idx_name} {qry_seq} > {paf_out}"
    run_cmd(cmd, logger=logger, fail_exit=True)
    return paf_out

class MinimapPAFtoBLAST6(object):
    """
    Iterate over PAF lines and yield BLAST6-formatted records with the requested fields.
    """
    def __init__(self, paf_path):
        self.paf_path = paf_path

    def __iter__(self):
        return self.parse()

    def parse(self):
        with open(self.paf_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                # PAF mandatory columns:
                # 0 qname, 1 qlen, 2 qstart, 3 qend, 4 strand, 5 tname, 6 tlen, 7 tstart, 8 tend, 9 nmatch, 10 alnlen, 11 mapq
                qname = fields[0]
                qlen = int(fields[1])
                qstart0 = int(fields[2])            # 0-based, inclusive
                qend0 = int(fields[3])              # 0-based, exclusive
                strand = fields[4]                  # '+' or '-'
                tname = fields[5]
                tlen = int(fields[6])
                tstart0 = int(fields[7])
                tend0 = int(fields[8])
                nmatch = int(fields[9])
                alnlen = int(fields[10])

                # optional tags
                tags = {}
                for t in fields[12:]:
                    if ':' in t:
                        k, typ, val = t.split(':', 2)
                        tags[k] = (typ, val)

                # cg (CIGAR), NM (edit distance) if present
                cigar = tags.get('cg', (None, None))[1]
                nm_val = tags.get('NM', (None, None))[1]
                NM = int(nm_val) if nm_val is not None else None

                ins_bases, del_bases, gap_opens = _parse_cigar_stats(cigar)

                # BLAST-ish stats:
                # pident
                pident = 100.0 * nmatch / alnlen if alnlen > 0 else 0.0

                # mismatches: BLAST counts mismatching bases (not gaps).
                # If NM available: NM = mismatches + insert_bases + del_bases
                if NM is not None:
                    mismatches = max(NM - (ins_bases + del_bases), 0)
                else:
                    # fallback: alnlen - matches - (ins+del)
                    mismatches = max(alnlen - nmatch - (ins_bases + del_bases), 0)

                # positions: convert to 1-based inclusive for BLAST format
                qstart = qstart0 + 1
                qend = qend0  # because qend0 is exclusive

                if strand == '+':
                    sstart = tstart0 + 1
                    send = tend0
                else:
                    # BLAST typically prints sstart > send for '-' strand
                    sstart = tend0
                    send = tstart0 + 1

                # evalue/bitscore: minimap2 doesn’t compute these.
                # Use 0 for evalue; surrogate bitscore = number of matches (monotonic with alignment quality).
                evalue = 0.0
                bitscore = float(nmatch)

                # coverages
                qcov = 100.0 * (qend0 - qstart0) / qlen if qlen > 0 else 0.0
                qcovs = qcov     # single-HSP → same
                qcovhsp = qcov   # ditto

                # assemble BLAST6 columns exactly in requested order
                out_values = [
                    qname,                 # qseqid
                    tname,                 # sseqid
                    f"{pident:.6f}",       # pident
                    str(alnlen),           # length
                    str(mismatches),       # mismatch
                    str(gap_opens),        # gapopen
                    str(qstart),           # qstart (1-based)
                    str(qend),             # qend   (1-based)
                    str(sstart),           # sstart
                    str(send),             # send
                    f"{evalue:.1e}",       # evalue
                    f"{bitscore:.1f}",     # bitscore (surrogate)
                    str(qlen),             # qlen
                    str(tlen),             # slen
                    f"{qcovs:.6f}",        # qcovs
                    f"{qcovhsp:.6f}",      # qcovhsp
                    strand,                # sstrand
                ]
                yield '\t'.join(out_values)

def minimap_to_blast6(db_seq, qry_seq, seqtype='nucl', db_name=None, blast_out=None,
                      preset='asm20', mm_opts='', ncpu=64):
    """
    Run minimap2 and write a BLAST6-formatted file with the desired columns.
    Returns path to the BLAST6-like output.
    """
    if seqtype not in ('nucl', 'prot'):
        # minimap2 here is used for nucleotide sequences (LTR-RTs). For protein, you'd need a different tool.
        raise ValueError('Only nucleotide sequences are supported with minimap2 in this wrapper.')
    if db_name is None:
        db_name = os.path.splitext(os.path.basename(db_seq))[0] + '.mmi'
    if blast_out is None:
        blast_out = qry_seq + '.blastout'

    paf_path = minimap(db_seq=db_seq, qry_seq=qry_seq, preset=preset,
                       idx_name=db_name, paf_out=None, mm_opts=mm_opts, ncpu=ncpu)

    # Convert PAF → BLAST6 (requested columns)
    with open(blast_out, 'w') as fout:
        for line in MinimapPAFtoBLAST6(paf_path):
            print(line, file=fout)
    return blast_out

# Back-compat wrapper so app.py can still `from .modules.Blast import blast, BlastOut`
def blast(db_seq, qry_seq, seqtype='nucl', blast_out=None, blast_outfmt=None, ncpu=4, **kwargs):
    """
    Backward-compatible wrapper that runs minimap2 and returns a BLAST6-like file path.
    Ignores blast_outfmt (we always write the requested fields).
    Additional kwargs (e.g. preset, mm_opts) get passed through.
    """
    return minimap_to_blast6(
        db_seq=db_seq,
        qry_seq=qry_seq,
        seqtype=seqtype,
        blast_out=blast_out,
        ncpu=ncpu,
        **kwargs
    )


# --- Backward-compatible “BlastOut” API using the same parser you had ---

class BlastOut(object):
    """
    Reads the BLAST6-like file produced by minimap_to_blast6 (or true BLAST outfmt 6).
    """
    def __init__(self, blast_out, outfmt=None):
        self.blast_out = blast_out
        if outfmt is not None:
            outfmt = outfmt.strip(''''"''')
        if outfmt is None:
            self.outfmt = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split()
        elif outfmt[0] == '6':
            self.outfmt = outfmt.split()[1:]
        else:
            raise ValueError('Only support for blast outfmt 6 = tabular, but {} input'.format(outfmt))
    def __iter__(self):
        return self.parse()
    def parse(self):
        for line in open(self.blast_out):
            values = line.strip().split('\t')
            yield BlastOutRecord(self.outfmt, values)
    def filter_besthit(self, fout=sys.stdout):
        d_best_hit = OrderedDict()
        for rc in self.parse():
            if rc.qseqid in d_best_hit:
                if rc.bitscore > d_best_hit[rc.qseqid].bitscore:
                    d_best_hit[rc.qseqid] = rc
            else:
                d_best_hit[rc.qseqid] = rc
        if fout is None:
            return d_best_hit
        for qseqid, rc in d_best_hit.items():
            rc.write(fout)
        return d_best_hit

class BlastOutRecord(object):
    def __init__(self, outfmt, values):
        self.values = values
        for key, value in zip(outfmt, values):
            setattr(self, key, BLASType[key](value))
    def write(self, fout=sys.stdout):
        print('\t'.join(self.values), file=fout)
    @property
    def scov(self):
        return 1.0 * (abs(self.send - self.sstart) + 1) / self.slen
