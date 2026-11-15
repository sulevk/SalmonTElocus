#!/usr/bin/env python3
"""SalmonTElocus: per-locus transposable element quantification with Salmon.

A rewrite of SalmonTE (Jeong et al., PSB 2018) that quantifies every genomic
TE copy separately instead of one consensus sequence per subfamily, and
reports each copy with its genomic coordinates.

    SalmonTElocus.py index --genome GENOME.fa --te-gtf TE.gtf --out INDEX_DIR
    SalmonTElocus.py quant --index INDEX_DIR --out OUT_DIR FASTQ_or_DIR ...
    SalmonTElocus.py test  --in OUT_DIR --out DE_DIR --conditions control,treatment

Only TEs are in the index: reads are never assigned to genes.
"""
import argparse
import csv
import gzip
import json
import logging
import os
import re
import shutil
import subprocess
import sys

__version__ = "1.0.4"

LOCI_FILE = "loci.csv"
META_FILE = "salmontelocus_index.json"
LOCUS_COLS = ["locus_id", "chrom", "start", "end", "strand", "transcript_id",
              "subfamily", "family", "class", "length"]
SUBFAMILY_COLS = ["subfamily", "family", "class"]
FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
# <sample>[._-](R)1|2[_001].fastq[.gz]
MATE_RE = re.compile(r"^(?P<sample>.+?)[._-]R?(?P<mate>[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$",
                     re.IGNORECASE)
_REVCOMP = str.maketrans("ACGTRYKMBVDHN", "TGCAYRMKVBHDN")

log = logging.getLogger("SalmonTElocus")


def die(msg, *args):
    log.error(msg, *args)
    sys.exit(1)


def _open(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


# --------------------------------------------------------------------------
# index
# --------------------------------------------------------------------------

def _gtf_attr(attr, key):
    tag = key + ' "'
    i = attr.find(tag)
    if i < 0:
        return None
    i += len(tag)
    return attr[i:attr.find('"', i)]


def read_te_loci(te_gtf, feature="exon"):
    """Parse a TEtranscripts-format TE GTF (one line per TE copy, with gene_id,
    transcript_id, family_id and class_id attributes) into locus records.

    Locus IDs use the TElocal / TEfastR format
    transcript_id:gene_id:family_id:class_id, made unique with a _N suffix if a
    GTF repeats one (same rule as R's make.unique(sep = "_")).
    """
    loci, seen = [], {}
    with _open(te_gtf) as inp:
        for n, line in enumerate(inp, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                die("%s line %d: not a 9-column GTF line", te_gtf, n)
            if f[2] != feature:
                continue
            tid, sub, fam, cls = (_gtf_attr(f[8], k) for k in
                                  ("transcript_id", "gene_id", "family_id", "class_id"))
            if None in (tid, sub, fam, cls):
                die("%s line %d lacks gene_id/transcript_id/family_id/class_id "
                    "(expected a TEtranscripts-format TE GTF)", te_gtf, n)
            locus_id = ":".join((tid, sub, fam, cls))
            if locus_id in seen:
                seen[locus_id] += 1
                locus_id = "{}_{}".format(locus_id, seen[locus_id])
            else:
                seen[locus_id] = 0
            start, end = int(f[3]), int(f[4])
            loci.append({"locus_id": locus_id, "chrom": f[0], "start": start, "end": end,
                         "strand": f[6], "transcript_id": tid, "subfamily": sub,
                         "family": fam, "class": cls, "length": end - start + 1})
    if not loci:
        die("No '%s' lines found in %s", feature, te_gtf)
    return loci


def iter_fasta(path):
    """Yield (name, UPPERCASE sequence) one record at a time. Uppercasing
    matters: soft-masked genomes write repeats -- i.e. every TE -- in lowercase."""
    name, chunks = None, []
    with _open(path) as inp:
        for line in inp:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).upper()


def extract_loci(genome, loci, out_fa, min_length):
    """Write each locus' sequence (reverse-complemented on '-') to out_fa.
    Returns the loci actually written."""
    by_chrom = {}
    for rec in loci:
        by_chrom.setdefault(rec["chrom"], []).append(rec)
    written, short, found = [], 0, set()
    with open(out_fa, "w") as oup:
        for chrom, seq in iter_fasta(genome):
            if chrom not in by_chrom:
                continue
            found.add(chrom)
            log.info("  %s: %d loci", chrom, len(by_chrom[chrom]))
            for rec in by_chrom[chrom]:
                s = seq[rec["start"] - 1:rec["end"]]
                if len(s) < min_length:
                    short += 1
                    continue
                if rec["strand"] == "-":
                    s = s.translate(_REVCOMP)[::-1]
                oup.write(">{}\n{}\n".format(rec["locus_id"], s))
                written.append(rec)
    missing = sorted(set(by_chrom) - found)
    if missing:
        log.warning("%d GTF sequence(s) absent from the genome FASTA (chromosome naming "
                    "mismatch, e.g. chr1 vs 1?): %d loci skipped, e.g. %s", len(missing),
                    sum(len(by_chrom[c]) for c in missing), ", ".join(missing[:3]))
    if short:
        log.warning("%d loci shorter than %d bp skipped (too short for Salmon's k-mers)",
                    short, min_length)
    if not written:
        die("No TE loci extracted: do the GTF and genome use the same chromosome names?")
    return written


def write_csv(path, rows, cols):
    with open(path, "w", newline="") as oup:
        w = csv.DictWriter(oup, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def salmon_bin(args):
    exe = args.salmon or os.environ.get("SALMON_BIN") or shutil.which("salmon")
    if not exe:
        die("Salmon not found: put it on PATH, set SALMON_BIN, or pass --salmon")
    return exe


def run_cmd(cmd):
    log.info("Running: %s", " ".join(cmd))
    ret = subprocess.run(cmd).returncode
    if ret != 0:
        die("Command failed (exit %d): %s", ret, " ".join(cmd))


def cmd_index(args):
    for p in (args.genome, args.te_gtf):
        if not os.path.exists(p):
            die("Input file does not exist: %s", p)
    if os.path.exists(os.path.join(args.out, META_FILE)) and not args.force:
        die("%s already holds a SalmonTElocus index; pass --force to rebuild it", args.out)
    os.makedirs(args.out, exist_ok=True)
    exe = salmon_bin(args)

    log.info("Reading TE loci from %s", args.te_gtf)
    loci = read_te_loci(args.te_gtf, args.feature)
    log.info("%d TE loci", len(loci))

    ref_fa = os.path.join(args.out, "te_loci.fa")
    log.info("Extracting TE sequences from %s", args.genome)
    written = extract_loci(args.genome, loci, ref_fa, args.min_length)
    write_csv(os.path.join(args.out, LOCI_FILE), written, LOCUS_COLS)

    # --keepDuplicates: identical young copies stay separate targets; without it
    # Salmon collapses them and those loci disappear from the output.
    salmon_idx = os.path.join(args.out, "salmon_index")
    run_cmd([exe, "index", "-t", ref_fa, "-i", salmon_idx, "--keepDuplicates",
             "-k", str(args.kmer), "-p", str(args.threads)])
    if not args.keep_fasta:
        os.remove(ref_fa)

    meta = {"tool": "SalmonTElocus", "version": __version__, "genome": os.path.abspath(args.genome),
            "te_gtf": os.path.abspath(args.te_gtf), "n_loci_in_gtf": len(loci),
            "n_loci_indexed": len(written), "min_length": args.min_length, "kmer": args.kmer}
    with open(os.path.join(args.out, META_FILE), "w") as oup:
        json.dump(meta, oup, indent=2)
    log.info("Index ready: %d of %d TE loci in %s", len(written), len(loci), args.out)


# --------------------------------------------------------------------------
# quant
# --------------------------------------------------------------------------

def is_fastq(path):
    return os.path.isfile(path) and path.lower().endswith(FASTQ_EXTS)


def _strip_fastq_ext(name):
    low = name.lower()
    for ext in sorted(FASTQ_EXTS, key=len, reverse=True):
        if low.endswith(ext):
            return name[:-len(ext)]
    return name


def collect_fastqs(inputs):
    files = []
    for item in inputs:
        if os.path.isdir(item):
            found = sorted(os.path.join(item, f) for f in os.listdir(item))
            found = [f for f in found if is_fastq(f)]
            if not found:
                die("No FASTQ files (%s) in directory %s", "/".join(FASTQ_EXTS), item)
            files.extend(found)
        elif is_fastq(item):
            files.append(item)
        else:
            die("Not a FASTQ file or directory: %s", item)
    return sorted(set(os.path.abspath(f) for f in files))


def pair_samples(files):
    """Group FASTQs into samples by file name: <sample>_R1/_R2 (also _1/_2,
    .1/.2, -1/-2, optional Illumina _001). A lone mate file is single-end."""
    mates = {}
    single = []
    for f in files:
        m = MATE_RE.match(os.path.basename(f))
        if m:
            mates.setdefault(m.group("sample"), {}).setdefault(m.group("mate"), []).append(f)
        else:
            single.append(f)
    samples = {}

    def add(name, entry):
        if name in samples:
            die("Two inputs map to the same sample name '%s'; use --samples", name)
        samples[name] = entry

    for name, m in mates.items():
        if any(len(v) > 1 for v in m.values()):
            die("Sample '%s' has several files for one mate (lanes?); use --samples", name)
        if "1" in m and "2" in m:
            add(name, {"r1": m["1"][0], "r2": m["2"][0]})
        else:  # e.g. sample_1.fastq without sample_2: single-end
            f = next(iter(m.values()))[0]
            single.append(f)
    for f in single:
        add(_strip_fastq_ext(os.path.basename(f)), {"r": f})
    return samples


def read_sample_sheet(path):
    """CSV with columns sample,r1[,r2] (r2 empty/absent = single-end)."""
    samples = {}
    base = os.path.dirname(os.path.abspath(path))
    with open(path, newline="") as inp:
        for row in csv.DictReader(inp):
            name = (row.get("sample") or "").strip()
            r1 = (row.get("r1") or "").strip()
            r2 = (row.get("r2") or "").strip()
            if not name or not r1:
                die("Sample sheet %s needs columns sample,r1[,r2]", path)
            if name in samples:
                die("Sample '%s' listed twice in %s", name, path)
            fix = lambda p: p if os.path.isabs(p) else os.path.join(base, p)
            for p in filter(None, (r1, r2)):
                if not os.path.exists(fix(p)):
                    die("Sample '%s': file not found: %s", name, fix(p))
            samples[name] = {"r1": fix(r1), "r2": fix(r2)} if r2 else {"r": fix(r1)}
    if not samples:
        die("Sample sheet %s lists no samples", path)
    return samples


def run_salmon_quant(exe, index_dir, name, entry, out_dir, args):
    sample_dir = os.path.join(out_dir, "salmon", name)
    if os.path.exists(os.path.join(sample_dir, "quant.sf")) and not args.force:
        log.info("%s: quant.sf exists, skipping (use --force to redo)", name)
        return sample_dir
    cmd = [exe, "quant", "-i", os.path.join(index_dir, "salmon_index"), "-l", args.libtype,
           "-o", sample_dir, "-p", str(args.threads)]
    cmd += ["-1", entry["r1"], "-2", entry["r2"]] if "r1" in entry else ["-r", entry["r"]]
    cmd += args.salmon_args.split() if args.salmon_args else []
    run_cmd(cmd)
    return sample_dir


def mapping_info(sample_dir):
    path = os.path.join(sample_dir, "aux_info", "meta_info.json")
    try:
        with open(path) as inp:
            m = json.load(inp)
        return {k: m.get(k) for k in ("num_processed", "num_mapped", "percent_mapped",
                                      "library_types")}
    except (OSError, ValueError):
        return {}


GENE_COLS = ["gene_context", "overlapping_genes", "gene_types", "gene_strand"]


class _IntervalIndex:
    """Per-chromosome interval lookup: features sorted by start plus a running
    maximum of ends, so the features that can overlap [s, e] form one
    contiguous index range found by two binary searches."""

    def __init__(self, recs_by_chrom):
        import numpy as np
        self.np = np
        self.data = {}
        for chrom, recs in recs_by_chrom.items():
            recs = sorted(recs)
            starts = np.array([r[0] for r in recs], dtype=np.int64)
            ends = np.array([r[1] for r in recs], dtype=np.int64)
            self.data[chrom] = (starts, ends, np.maximum.accumulate(ends), recs)

    def overlaps(self, chrom, s, e):
        d = self.data.get(chrom)
        if d is None:
            return []
        starts, ends, cummax, recs = d
        lo = int(self.np.searchsorted(cummax, s, "left"))
        hi = int(self.np.searchsorted(starts, e, "right"))
        return [recs[i] for i in range(lo, hi) if ends[i] >= s]


def read_gene_models(gene_gtf):
    """Gene bodies and exons from a gene GTF (GENCODE/Ensembl). Records are
    (start, end, gene_name, gene_type, strand). Gene bodies are derived from
    exons when the GTF has no 'gene' lines."""
    genes, exons = {}, {}
    with _open(gene_gtf) as inp:
        for line in inp:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            name = _gtf_attr(f[8], "gene_name") or _gtf_attr(f[8], "gene_id")
            gtype = _gtf_attr(f[8], "gene_type") or _gtf_attr(f[8], "gene_biotype") or "NA"
            rec = (int(f[3]), int(f[4]), name, gtype, f[6])
            (genes if f[2] == "gene" else exons).setdefault(f[0], set()).add(rec)
    if not exons:
        die("No exon lines found in gene GTF %s", gene_gtf)
    if not genes:
        span = {}
        for chrom, recs in exons.items():
            for s, e, name, gtype, strand in recs:
                k = (chrom, name, gtype, strand)
                lo, hi = span.get(k, (s, e))
                span[k] = (min(lo, s), max(hi, e))
        for (chrom, name, gtype, strand), (s, e) in span.items():
            genes.setdefault(chrom, set()).add((s, e, name, gtype, strand))
    return _IntervalIndex(genes), _IntervalIndex(exons)


def annotate_gene_context(loci, gene_gtf):
    """Label each TE locus by the genes it overlaps (no reads are assigned to
    genes; this is annotation only). Returns a DataFrame with GENE_COLS,
    aligned to `loci`."""
    import pandas as pd
    log.info("Annotating %d loci with gene context from %s", len(loci), gene_gtf)
    gene_idx, exon_idx = read_gene_models(gene_gtf)
    shared = set(loci["chrom"]) & set(exon_idx.data)
    if not shared:
        log.warning("No chromosome names shared between the TE loci and %s "
                    "(chr1 vs 1?); every locus will be 'intergenic'", gene_gtf)
    rows = []
    for chrom, s, e, strand in zip(loci["chrom"], loci["start"], loci["end"], loci["strand"]):
        hits = exon_idx.overlaps(chrom, s, e)
        context = "exon"
        if not hits:
            hits = gene_idx.overlaps(chrom, s, e)
            context = "intron" if hits else "intergenic"
        if not hits:
            rows.append(("intergenic", "", "", ""))
            continue
        by_gene = {}
        for _, _, name, gtype, gstrand in hits:
            by_gene.setdefault(name, (gtype, set()))[1].add(gstrand)
        names = sorted(by_gene)
        rel = set()
        for name in names:
            for gs in by_gene[name][1]:
                if strand in "+-" and gs in "+-":
                    rel.add("sense" if gs == strand else "antisense")
        rows.append((context, ";".join(names), ";".join(by_gene[n][0] for n in names),
                     "both" if len(rel) == 2 else (rel.pop() if rel else "")))
    return pd.DataFrame(rows, columns=GENE_COLS, index=loci.index)


def collect_tables(index_dir, sample_dirs, out_dir, keep_zero, gene_gtf=None):
    import pandas as pd
    loci = pd.read_csv(os.path.join(index_dir, LOCI_FILE), dtype={"chrom": str})
    counts, tpm = {}, {}
    for name, d in sample_dirs.items():
        q = pd.read_csv(os.path.join(d, "quant.sf"), sep="\t", index_col="Name")
        counts[name] = q["NumReads"]
        tpm[name] = q["TPM"]
    samples = sorted(counts)

    def per_locus(values):
        tab = pd.DataFrame(values)[samples].reindex(loci["locus_id"]).fillna(0.0)
        return pd.concat([loci.reset_index(drop=True), tab.reset_index(drop=True)], axis=1)

    lc, lt = per_locus(counts), per_locus(tpm)
    for kind, tab in (("counts", lc), ("tpm", lt)):
        sub = tab.groupby(SUBFAMILY_COLS, sort=True)[samples].sum().reset_index()
        sub.to_csv(os.path.join(out_dir, "subfamily_{}.csv".format(kind)), index=False)
    expressed = lc[samples].sum(axis=1) > 0
    if not keep_zero:
        lc, lt = lc[expressed], lt[expressed]
    if gene_gtf:
        lc, lt = lc.copy(), lt.copy()
        ann = annotate_gene_context(lc, gene_gtf)
        pos = lc.columns.get_loc("length") + 1
        for j, col in enumerate(GENE_COLS):
            lc.insert(pos + j, col, ann[col])
            lt.insert(pos + j, col, ann[col])
        expr = lc[lc[samples].sum(axis=1) > 0]
        reads = expr[samples].sum(axis=1)
        for ctx in ("exon", "intron", "intergenic"):
            m = expr["gene_context"] == ctx
            log.info("  %-10s %8d expressed loci, %5.1f%% of TE reads", ctx, int(m.sum()),
                     100 * reads[m].sum() / max(reads.sum(), 1e-9))
    lc.to_csv(os.path.join(out_dir, "locus_counts.csv"), index=False)
    lt.to_csv(os.path.join(out_dir, "locus_tpm.csv"), index=False)

    # BED6 of loci with reads, for genome browsers (score = mean count, capped at 1000)
    bed = lc[lc[samples].sum(axis=1) > 0]
    mean = bed[samples].mean(axis=1)
    with open(os.path.join(out_dir, "expressed_loci.bed"), "w") as oup:
        for (_, r), m in zip(bed.iterrows(), mean):
            oup.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(r["chrom"], r["start"] - 1, r["end"],
                                                        r["locus_id"], min(1000, int(round(m))),
                                                        r["strand"] if r["strand"] in "+-" else "."))
    pd.DataFrame([dict(sample=s, **mapping_info(sample_dirs[s])) for s in samples]) \
        .to_csv(os.path.join(out_dir, "mapping_info.csv"), index=False)

    cond = os.path.join(out_dir, "conditions.csv")
    if not os.path.exists(cond):
        pd.DataFrame({"sample": samples, "condition": "NA"}).to_csv(cond, index=False)
    log.info("%d samples; %d of %d loci have reads", len(samples), int(expressed.sum()),
             len(loci))


def cmd_quant(args):
    if not os.path.exists(os.path.join(args.index, META_FILE)):
        die("%s is not a SalmonTElocus index (build one with 'SalmonTElocus.py index')", args.index)
    if args.gene_gtf and not os.path.exists(args.gene_gtf):
        die("Gene GTF not found: %s", args.gene_gtf)
    if args.samples:
        if args.inputs:
            die("Give either --samples or FASTQ inputs, not both")
        samples = read_sample_sheet(args.samples)
    else:
        if not args.inputs:
            die("No FASTQ inputs given")
        samples = pair_samples(collect_fastqs(args.inputs))
    n_pe = sum("r1" in e for e in samples.values())
    log.info("%d samples (%d paired-end, %d single-end): %s", len(samples), n_pe,
             len(samples) - n_pe, ", ".join(sorted(samples)))
    os.makedirs(args.out, exist_ok=True)
    exe = salmon_bin(args)
    sample_dirs = {name: run_salmon_quant(exe, args.index, name, entry, args.out, args)
                   for name, entry in sorted(samples.items())}
    collect_tables(args.index, sample_dirs, args.out, args.keep_zero, args.gene_gtf)
    log.info("Done. Tables are in %s", args.out)


# --------------------------------------------------------------------------
# test
# --------------------------------------------------------------------------

def cmd_test(args):
    table = os.path.join(args.inp, "{}_counts.csv".format(args.level))
    if not os.path.exists(table):
        die("%s not found: run 'quant' first", table)
    cond = args.condition_file or os.path.join(args.inp, "conditions.csv")
    if not os.path.exists(cond):
        die("Condition file %s not found", cond)
    rscript = shutil.which("Rscript")
    if not rscript:
        die("Rscript not found on PATH")
    if args.level == "locus" and not args.gene_counts:
        log.warning("Locus-level test without --gene-counts: size factors come from TE copies "
                    "expressed in every sample, a subset that can be strongly skewed and bias "
                    "all fold changes. Gene-based size factors (--gene-counts, see 'gene-table') "
                    "are recommended at locus level.")
    os.makedirs(args.out, exist_ok=True)
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "SalmonTElocus_DE.R")
    if args.gene_counts and not os.path.exists(args.gene_counts):
        die("Gene count table %s not found", args.gene_counts)
    cmd = [rscript, script, table, cond, args.out, args.level, args.conditions or "",
           str(args.min_count), args.gene_counts or ""]
    run_cmd(cmd)


def cmd_gene_table(args):
    """Sum Salmon NumReads per gene across sample folders -> gene x sample CSV,
    for `test --gene-counts`. Genes are never part of the TE quantification."""
    import pandas as pd
    tx2gene = None
    if args.tx2gene:
        t = pd.read_csv(args.tx2gene, sep=None, engine="python", header=None, dtype=str)
        tx2gene = dict(zip(t.iloc[:, 0], t.iloc[:, 1]))
    cols = {}
    for d in sorted(args.quant_dirs):
        qf = os.path.join(d, "quant.sf")
        if not os.path.exists(qf):
            log.warning("No quant.sf in %s, skipped", d)
            continue
        q = pd.read_csv(qf, sep="\t", usecols=["Name", "NumReads"])
        if tx2gene is not None:
            gene = q["Name"].map(tx2gene)
            if gene.isna().all():
                die("No transcript in %s matched --tx2gene", qf)
            gene = gene.fillna(q["Name"])
        elif q["Name"].str.contains("|", regex=False).any():
            gene = q["Name"].str.split("|").str[1]  # GENCODE: ENST|ENSG|...
        else:
            die("%s: transcript names carry no gene ID; pass --tx2gene", qf)
        cols[os.path.basename(os.path.normpath(d))] = q.groupby(gene)["NumReads"].sum()
    if not cols:
        die("No Salmon results found")
    tab = pd.DataFrame(cols).fillna(0.0)
    tab.index.name = "gene"
    tab.to_csv(args.out)
    log.info("Wrote %d genes x %d samples to %s", tab.shape[0], tab.shape[1], args.out)


# --------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(prog="SalmonTElocus.py", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version="SalmonTElocus " + __version__)
    sub = p.add_subparsers(dest="command", required=True)

    i = sub.add_parser("index", help="build a per-locus TE index from a genome + TE GTF")
    i.add_argument("--genome", required=True, help="genome FASTA (.fa or .fa.gz)")
    i.add_argument("--te-gtf", required=True,
                   help="TEtranscripts-format TE GTF, e.g. hg38_rmsk_TE.gtf (.gtf or .gtf.gz)")
    i.add_argument("--out", required=True, help="output index directory")
    i.add_argument("--feature", default="exon", help="GTF feature type to use [exon]")
    i.add_argument("--min-length", type=int, default=31,
                   help="skip loci shorter than this (bp) [31]")
    i.add_argument("--kmer", type=int, default=31, help="Salmon k-mer size [31]")
    i.add_argument("--threads", type=int, default=4)
    i.add_argument("--keep-fasta", action="store_true", help="keep the extracted te_loci.fa")
    i.add_argument("--force", action="store_true", help="overwrite an existing index")
    i.add_argument("--salmon", help="path to the salmon binary")
    i.set_defaults(func=cmd_index)

    q = sub.add_parser("quant", help="quantify FASTQ samples against a locus index")
    q.add_argument("inputs", nargs="*", help="FASTQ files and/or directories of FASTQ files")
    q.add_argument("--index", required=True, help="index directory from 'index'")
    q.add_argument("--out", required=True, help="output directory")
    q.add_argument("--samples", help="CSV sample sheet: sample,r1[,r2] (instead of inputs)")
    q.add_argument("--threads", type=int, default=4)
    q.add_argument("--libtype", default="A", help="Salmon library type [A = auto-detect]")
    q.add_argument("--salmon-args", default="",
                   help="extra arguments passed to 'salmon quant', e.g. \"--useEM --seqBias\"")
    q.add_argument("--gene-gtf",
                   help="gene GTF (e.g. GENCODE) used ONLY to label each TE locus as exon/intron/"
                        "intergenic with the overlapping gene names; reads are never counted for genes")
    q.add_argument("--keep-zero", action="store_true",
                   help="keep loci with no reads in any sample in the locus tables")
    q.add_argument("--force", action="store_true", help="re-run samples that already have results")
    q.add_argument("--salmon", help="path to the salmon binary")
    q.set_defaults(func=cmd_quant)

    t = sub.add_parser("test", help="differential expression with DESeq2")
    t.add_argument("--in", dest="inp", required=True, help="output directory of 'quant'")
    t.add_argument("--out", required=True, help="directory for results")
    t.add_argument("--level", choices=["locus", "subfamily"], default="locus")
    t.add_argument("--conditions",
                   help="two condition labels, reference first, e.g. control,treatment. "
                        "Omit when the condition column is numeric (e.g. age).")
    t.add_argument("--condition-file", help="CSV sample,condition [<in>/conditions.csv]")
    t.add_argument("--min-count", type=int, default=10,
                   help="drop features with fewer total reads than this [10]")
    t.add_argument("--gene-counts",
                   help="gene x sample count table (CSV/TSV, first column gene ID, e.g. from "
                        "'gene-table'): DESeq2 size factors are computed from genes instead of TEs, "
                        "so global TE increases are not normalised away")
    t.set_defaults(func=cmd_test)

    g = sub.add_parser("gene-table",
                       help="gene x sample counts from Salmon transcriptome runs (for test --gene-counts)")
    g.add_argument("quant_dirs", nargs="+", help="Salmon output folders, one per sample (folder name = sample)")
    g.add_argument("--out", required=True, help="output CSV")
    g.add_argument("--tx2gene", help="2-column transcript,gene file (not needed for GENCODE FASTA headers)")
    g.set_defaults(func=cmd_gene_table)
    return p


def main(argv=None):
    logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
