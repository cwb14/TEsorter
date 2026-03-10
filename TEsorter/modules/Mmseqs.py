import os
import sys
from collections import OrderedDict
from .RunCmdsMP import run_cmd, logger

def mmseqs_easy_search(db_seq, qry_seq, out_m8, tmpdir,
                       seqtype="nucl", ncpu=4,
                       min_seq_id=0.0, min_cov=0.0, cov_mode=2, min_aln_len=0,
                       sensitivity=None,
                       mmseqs_bin="mmseqs"):
    """
    Run mmseqs easy-search and write a BLAST-tab-like file.
    We request a custom format to make parsing easy.
    """

    os.makedirs(tmpdir, exist_ok=True)

    cleaned_qry = os.path.join(tmpdir, "query.cleaned.fa")
    clean_fasta_atcg_only(qry_seq, cleaned_qry)
    qry_seq = cleaned_qry

    # search-type: 3 nucleotide; 1 amino-acid
    if seqtype == "nucl":
        search_type = 3
    elif seqtype == "prot":
        search_type = 1
    else:
        raise ValueError(f"Unknown seqtype {seqtype}")

    # qstart/qend/tstart/tend are needed for the split-alignment merge step (post-processing)
    fmt = "query,target,fident,alnlen,qlen,qcov,bits,qstart,qend,tstart,tend"

    cmd = (
        f"{mmseqs_bin} easy-search "
        f"{qry_seq} {db_seq} {out_m8} {tmpdir} "
        f"--threads {ncpu} "
        f"--search-type {search_type} "
        f"--format-output \"{fmt}\" "
        f"--min-seq-id {min_seq_id} "
        f"-c {min_cov} --cov-mode {cov_mode} "
        f"--min-aln-len {min_aln_len} "
    )
    if sensitivity is not None:
        cmd += f"-s {float(sensitivity)} "

    run_cmd(cmd, logger=logger)
    return out_m8

def clean_fasta_atcg_only(in_fa, out_fa):
    """
    Remove all non-ATCG characters from sequences in a FASTA file.
    Keeps headers unchanged.
    """
    with open(in_fa) as fin, open(out_fa, "w") as fout:
        seq_buf = []
        header = None

        def flush():
            if header is None:
                return
            seq = "".join(seq_buf).upper()
            seq = "".join([b for b in seq if b in ("A", "T", "C", "G")])
            fout.write(header + "\n")
            fout.write(seq + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                flush()
                header = line
                seq_buf = []
            else:
                seq_buf.append(line)

        flush()

class MmseqsM8Record:
    __slots__ = ("qseqid", "sseqid", "fident", "alnlen", "qlen", "qcov", "bits", "qstart", "qend", "tstart", "tend")

    def __init__(self, line):
        # query,target,fident,alnlen,qlen,qcov,bits,qstart,qend,tstart,tend
        vals = line.rstrip("\n").split("\t")
        self.qseqid  = vals[0]
        self.sseqid  = vals[1]
        self.fident  = float(vals[2])      # 0..1
        self.alnlen  = int(vals[3])
        self.qlen    = int(vals[4])
        self.qcov    = float(vals[5])      # 0..1
        self.bits    = float(vals[6])
        self.qstart  = int(vals[7])        # 1-based
        self.qend    = int(vals[8])        # 1-based
        self.tstart  = int(vals[9])        # 1-based
        self.tend    = int(vals[10])       # 1-based


# Maximum gap (bp) on both query and target sides within which two blocks
# are considered part of the same MMseqs2-split alignment rather than two
# genuinely independent local alignments.
_MAX_SPLIT_GAP = 500


def _chain_blocks(records):
    """
    Partition records (all sharing the same query-target pair) into contiguous
    chains.  Two consecutive blocks (sorted by qstart) belong to the same chain
    when the gap between them is ≤ _MAX_SPLIT_GAP on *both* the query and the
    target.  A negative gap (overlap) always satisfies the criterion.

    Returns a list-of-lists, where each inner list is one contiguous chain.
    """
    blocks = sorted(records, key=lambda r: (r.qstart, r.qend))
    chains = [[blocks[0]]]
    for r in blocks[1:]:
        prev = chains[-1][-1]
        qgap = r.qstart - prev.qend - 1          # <0 = overlap
        tgap = abs(r.tstart - prev.tend) - 1      # strand-agnostic
        if qgap <= _MAX_SPLIT_GAP and tgap <= _MAX_SPLIT_GAP:
            chains[-1].append(r)
        else:
            chains.append([r])
    return chains


def _merge_chain(records):
    """
    Merge a list of MmseqsM8Record objects that form a single contiguous chain
    (same query-target, confirmed adjacent by _chain_blocks).

      - unions the aligned query intervals (handles overlaps correctly)
      - sums bit-scores
      - computes a length-weighted average fident
      - recomputes qcov = merged_alnlen / qlen
    """
    if len(records) == 1:
        return records[0]

    blocks = sorted(records, key=lambda r: (r.qstart, r.qend))

    # Union of query intervals → total non-overlapping aligned length
    merged_intervals = []
    for r in blocks:
        lo, hi = r.qstart, r.qend
        if merged_intervals and lo <= merged_intervals[-1][1]:
            prev_lo, prev_hi = merged_intervals[-1]
            merged_intervals[-1] = (prev_lo, max(prev_hi, hi))
        else:
            merged_intervals.append((lo, hi))

    merged_alnlen = sum(hi - lo + 1 for lo, hi in merged_intervals)

    # Weighted-average fident (weight by each block's alnlen)
    total_weight = sum(r.alnlen for r in blocks)
    weighted_fident = sum(r.fident * r.alnlen for r in blocks) / total_weight if total_weight else 0.0

    # Sum bit-scores (per user spec)
    merged_bits = sum(r.bits for r in blocks)

    ref = blocks[0]
    synthetic = MmseqsM8Record.__new__(MmseqsM8Record)
    synthetic.qseqid = ref.qseqid
    synthetic.sseqid = ref.sseqid
    synthetic.fident  = weighted_fident
    synthetic.alnlen  = merged_alnlen
    synthetic.qlen    = ref.qlen
    synthetic.qcov    = merged_alnlen / ref.qlen if ref.qlen else 0.0
    synthetic.bits    = merged_bits
    synthetic.qstart  = merged_intervals[0][0]
    synthetic.qend    = merged_intervals[-1][1]
    synthetic.tstart  = blocks[0].tstart
    synthetic.tend    = blocks[-1].tend
    return synthetic


def _merge_split_alignments(records):
    """
    Given a flat list of MmseqsM8Record objects (potentially many per query),
    group by (query, target), split each group into contiguous chains (so that
    genuinely independent local alignments to the same target are not merged),
    merge each chain, then return the results as a flat list.
    """
    from collections import defaultdict
    groups = defaultdict(list)
    for r in records:
        groups[(r.qseqid, r.sseqid)].append(r)

    merged = []
    for grp in groups.values():
        for chain in _chain_blocks(grp):
            merged.append(_merge_chain(chain))
    return merged


def parse_mmseqs_m8_besthit(m8_path):
    """
    Merge MMseqs2 split-alignment rows, then keep the best hit per query by
    bits (like BlastOut.filter_besthit()).
    """
    all_records = []
    with open(m8_path) as f:
        for line in f:
            if not line.strip():
                continue
            all_records.append(MmseqsM8Record(line))

    merged = _merge_split_alignments(all_records)

    best = OrderedDict()
    for rc in merged:
        if rc.qseqid not in best or rc.bits > best[rc.qseqid].bits:
            best[rc.qseqid] = rc
    return best
