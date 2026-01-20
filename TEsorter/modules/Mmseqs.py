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

    # search-type: 3 nucleotide; 1 amino-acid
    if seqtype == "nucl":
        search_type = 3
    elif seqtype == "prot":
        search_type = 1
    else:
        raise ValueError(f"Unknown seqtype {seqtype}")

    fmt = "query,target,fident,alnlen,qlen,qcov,bits"

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


class MmseqsM8Record:
    __slots__ = ("qseqid", "sseqid", "fident", "alnlen", "qlen", "qcov", "bits")

    def __init__(self, line):
        # query,target,fident,alnlen,qlen,qcov,bits
        vals = line.rstrip("\n").split("\t")
        self.qseqid = vals[0]
        self.sseqid = vals[1]
        self.fident = float(vals[2])      # 0..1
        self.alnlen = int(vals[3])
        self.qlen = int(vals[4])
        self.qcov = float(vals[5])        # usually 0..1
        self.bits = float(vals[6])


def parse_mmseqs_m8_besthit(m8_path):
    """
    Keep best hit per query by bits (like your BlastOut.filter_besthit()).
    """
    best = OrderedDict()
    with open(m8_path) as f:
        for line in f:
            if not line.strip():
                continue
            rc = MmseqsM8Record(line)
            if rc.qseqid not in best or rc.bits > best[rc.qseqid].bits:
                best[rc.qseqid] = rc
    return best
