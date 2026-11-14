"""Tests for SalmonTElocus.  Run from the repository root:

    python3 -m unittest discover -s tests -v

The end-to-end tests need `salmon` on PATH (and Rscript + DESeq2 for `test`);
they are skipped otherwise.
"""
import csv
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import SalmonTElocus as stl  # noqa: E402

SCRIPT = os.path.join(ROOT, "SalmonTElocus.py")
HAVE_SALMON = shutil.which("salmon") is not None
HAVE_DESEQ2 = shutil.which("Rscript") is not None and subprocess.run(
    ["Rscript", "-e", "quit(status = !requireNamespace('DESeq2', quietly = TRUE))"],
    capture_output=True).returncode == 0


def revcomp(s):
    return s.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def gtf_line(chrom, start, end, strand, tid, sub, fam, cls):
    return ('{}\tt\texon\t{}\t{}\t.\t{}\t.\tgene_id "{}"; transcript_id "{}"; '
            'family_id "{}"; class_id "{}";\n').format(chrom, start, end, strand, sub, tid, fam, cls)


def read_csv(path):
    with open(path, newline="") as inp:
        return list(csv.DictReader(inp))


class TestPairing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def touch(self, *names):
        for n in names:
            open(os.path.join(self.tmp, n), "w").close()
        return [os.path.join(self.tmp, n) for n in names]

    def test_pair_patterns(self):
        files = self.touch("A_R1.fastq.gz", "A_R2.fastq.gz", "B_1.fq", "B_2.fq",
                           "C_S1_L001_R1_001.fastq.gz", "C_S1_L001_R2_001.fastq.gz",
                           "single.fastq", "lonely_1.fastq")
        s = stl.pair_samples(files)
        self.assertEqual(set(s), {"A", "B", "C_S1_L001", "single", "lonely_1"})
        self.assertTrue(s["A"]["r1"].endswith("A_R1.fastq.gz"))
        self.assertTrue(s["B"]["r2"].endswith("B_2.fq"))
        self.assertEqual(set(s["single"]), {"r"})
        self.assertEqual(set(s["lonely_1"]), {"r"})

    def test_directory_named_fq_is_not_a_fastq(self):
        d = os.path.join(self.tmp, "fq")
        os.mkdir(d)
        open(os.path.join(d, "x.fastq"), "w").close()
        self.assertFalse(stl.is_fastq(d))
        self.assertEqual([os.path.basename(f) for f in stl.collect_fastqs([d])], ["x.fastq"])


class TestIndexParts(unittest.TestCase):
    def test_gtf_ids_and_duplicates(self):
        with tempfile.NamedTemporaryFile("w", suffix=".gtf", delete=False) as f:
            f.write(gtf_line("chr1", 1, 50, "+", "X_dup1", "X", "F", "C"))
            f.write(gtf_line("chr1", 60, 99, "-", "X_dup1", "X", "F", "C"))
        try:
            loci = stl.read_te_loci(f.name)
        finally:
            os.remove(f.name)
        self.assertEqual([l["locus_id"] for l in loci], ["X_dup1:X:F:C", "X_dup1:X:F:C_1"])
        self.assertEqual(loci[1]["length"], 40)

    def test_extract_softmasked_and_minus_strand(self):
        tmp = tempfile.mkdtemp()
        try:
            genome = os.path.join(tmp, "g.fa")
            with open(genome, "w") as f:
                f.write(">chr1 some description\nNNNNacgtacgtAAAACCCCGGGG\nTTTT\n")
            loci = [{"locus_id": "a", "chrom": "chr1", "start": 5, "end": 12, "strand": "+"},
                    {"locus_id": "b", "chrom": "chr1", "start": 13, "end": 28, "strand": "-"},
                    {"locus_id": "c", "chrom": "chrX", "start": 1, "end": 5, "strand": "+"}]
            out = os.path.join(tmp, "o.fa")
            written = stl.extract_loci(genome, loci, out, min_length=1)
            with open(out) as f:
                seqs = f.read().split()
            self.assertEqual(seqs, [">a", "ACGTACGT", ">b", "AAAACCCCGGGGTTTT"])
            self.assertEqual([w["locus_id"] for w in written], ["a", "b"])
        finally:
            shutil.rmtree(tmp)


@unittest.skipUnless(HAVE_SALMON, "salmon not on PATH")
class TestEndToEnd(unittest.TestCase):
    """Synthetic soft-masked genome:
      L1A   chr1 (+)  unique copy, 2x higher in 'treat'
      AluB  chr1 (-)  unique copy
      ERVC  chr2 (+)  two IDENTICAL copies (reads must split, both loci kept)
      X     chrUn     absent from the genome (skipped with a warning)
    """

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="stl_")
        rng = random.Random(7)
        rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        teA, teB, teC = rnd(1200), rnd(900), rnd(800)
        c1 = rnd(1000) + teA.lower() + rnd(1000) + revcomp(teB).lower() + rnd(1000)
        c2 = rnd(500) + teC.lower() + rnd(700) + teC.lower() + rnd(500)
        with open(os.path.join(cls.tmp, "genome.fa"), "w") as f:
            f.write(">chr1\n" + "\n".join(c1[i:i + 60] for i in range(0, len(c1), 60)) + "\n")
            f.write(">chr2\n" + "\n".join(c2[i:i + 60] for i in range(0, len(c2), 60)) + "\n")
        with open(os.path.join(cls.tmp, "te.gtf"), "w") as f:
            f.write(gtf_line("chr1", 1001, 2200, "+", "L1A_dup1", "L1A", "L1", "LINE"))
            f.write(gtf_line("chr1", 3201, 4100, "-", "AluB_dup1", "AluB", "Alu", "SINE"))
            f.write(gtf_line("chr2", 501, 1300, "+", "ERVC_dup1", "ERVC", "ERV1", "LTR"))
            f.write(gtf_line("chr2", 2001, 2800, "+", "ERVC_dup2", "ERVC", "ERV1", "LTR"))
            f.write(gtf_line("chrUn", 1, 500, "+", "X_dup1", "X", "X", "DNA"))

        cls.truth = {}
        reads_dir = os.path.join(cls.tmp, "fq")  # deliberately named like a FASTQ extension
        os.mkdir(reads_dir)
        design = {"ctrl1": 100, "ctrl2": 110, "treat1": 200, "treat2": 220}
        for sample, nA in design.items():
            frags = [(teA, nA), (teB, 60), (teC, 80)]
            with open(os.path.join(reads_dir, sample + "_R1.fastq"), "w") as r1, \
                    open(os.path.join(reads_dir, sample + "_R2.fastq"), "w") as r2:
                k = 0
                for seq, n in frags:
                    for _ in range(n):
                        s = rng.randrange(0, len(seq) - 250)
                        frag = seq[s:s + 250]
                        r1.write("@{}_{}/1\n{}\n+\n{}\n".format(sample, k, frag[:75], "I" * 75))
                        r2.write("@{}_{}/2\n{}\n+\n{}\n".format(sample, k, revcomp(frag[-75:]), "I" * 75))
                        k += 1
            cls.truth[sample] = {"L1A": nA, "AluB": 60, "ERVC": 80}

        cls.index = os.path.join(cls.tmp, "index")
        cls.out = os.path.join(cls.tmp, "out")
        cls.sh([SCRIPT, "index", "--genome", os.path.join(cls.tmp, "genome.fa"),
                 "--te-gtf", os.path.join(cls.tmp, "te.gtf"), "--out", cls.index, "--threads", "2"])
        cls.sh([SCRIPT, "quant", "--index", cls.index, "--out", cls.out, "--threads", "2", reads_dir])
        # Gene annotation (exon lines only -> gene bodies derived from exons):
        # GENE1 (+) has an exon overlapping L1A and an intron containing AluB.
        cls.gene_gtf = os.path.join(cls.tmp, "genes.gtf")
        with open(cls.gene_gtf, "w") as f:
            for s, e in ((900, 1100), (5000, 5200)):
                f.write('chr1\tt\texon\t{}\t{}\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; '
                        'gene_type "protein_coding";\n'.format(s, e))
        cls.out_genes = os.path.join(cls.tmp, "out_genes")
        cls.sh([SCRIPT, "quant", "--index", cls.index, "--out", cls.out_genes, "--threads", "2",
                "--gene-gtf", cls.gene_gtf, reads_dir])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    @classmethod
    def sh(cls, cmd):
        p = subprocess.run([sys.executable] + cmd, capture_output=True, text=True)
        if p.returncode != 0:
            raise AssertionError("command failed: {}\n{}".format(" ".join(cmd), p.stderr[-3000:]))
        return p

    def test_index_loci(self):
        loci = read_csv(os.path.join(self.index, "loci.csv"))
        self.assertEqual([l["locus_id"] for l in loci],
                         ["L1A_dup1:L1A:L1:LINE", "AluB_dup1:AluB:Alu:SINE",
                          "ERVC_dup1:ERVC:ERV1:LTR", "ERVC_dup2:ERVC:ERV1:LTR"])

    def test_locus_counts_and_coordinates(self):
        rows = {r["locus_id"]: r for r in read_csv(os.path.join(self.out, "locus_counts.csv"))}
        a = rows["L1A_dup1:L1A:L1:LINE"]
        self.assertEqual((a["chrom"], a["start"], a["end"], a["strand"]), ("chr1", "1001", "2200", "+"))
        for s, t in self.truth.items():
            self.assertAlmostEqual(float(rows["L1A_dup1:L1A:L1:LINE"][s]), t["L1A"], delta=1)
            self.assertAlmostEqual(float(rows["AluB_dup1:AluB:Alu:SINE"][s]), t["AluB"], delta=1)
            c1 = float(rows["ERVC_dup1:ERVC:ERV1:LTR"][s])
            c2 = float(rows["ERVC_dup2:ERVC:ERV1:LTR"][s])
            self.assertAlmostEqual(c1 + c2, t["ERVC"], delta=1)
            self.assertAlmostEqual(c1, c2, delta=1)  # identical copies share evenly

    def test_subfamily_sums_and_side_files(self):
        rows = {r["subfamily"]: r for r in read_csv(os.path.join(self.out, "subfamily_counts.csv"))}
        self.assertAlmostEqual(float(rows["ERVC"]["ctrl1"]), 80, delta=1)
        self.assertNotIn("X", rows)
        with open(os.path.join(self.out, "expressed_loci.bed")) as f:
            bed = [l.split("\t") for l in f]
        self.assertIn(["chr1", "1000", "2200", "L1A_dup1:L1A:L1:LINE"], [b[:4] for b in bed])
        info = read_csv(os.path.join(self.out, "mapping_info.csv"))
        self.assertEqual([r["sample"] for r in info], ["ctrl1", "ctrl2", "treat1", "treat2"])

    def test_gene_context_annotation(self):
        rows = {r["locus_id"]: r for r in read_csv(os.path.join(self.out_genes, "locus_counts.csv"))}
        a, b = rows["L1A_dup1:L1A:L1:LINE"], rows["AluB_dup1:AluB:Alu:SINE"]
        self.assertEqual((a["gene_context"], a["overlapping_genes"], a["gene_types"], a["gene_strand"]),
                         ("exon", "GENE1", "protein_coding", "sense"))
        self.assertEqual((b["gene_context"], b["gene_strand"]), ("intron", "antisense"))
        self.assertEqual(rows["ERVC_dup1:ERVC:ERV1:LTR"]["gene_context"], "intergenic")
        # annotation only: counts are identical with and without --gene-gtf
        plain = {r["locus_id"]: r for r in read_csv(os.path.join(self.out, "locus_counts.csv"))}
        for lid in plain:
            self.assertEqual(plain[lid]["ctrl1"], rows[lid]["ctrl1"])
        header = list(rows[next(iter(rows))])
        self.assertEqual(header[10:14], ["gene_context", "overlapping_genes", "gene_types", "gene_strand"])

    @unittest.skipUnless(HAVE_DESEQ2, "Rscript/DESeq2 not available")
    def test_differential_expression_with_gene_context(self):
        cond = os.path.join(self.tmp, "cond2.csv")
        with open(cond, "w") as f:
            f.write("sample,condition\nctrl1,ctrl\nctrl2,ctrl\ntreat1,treat\ntreat2,treat\n")
        de = os.path.join(self.tmp, "de_genes")
        self.sh([SCRIPT, "test", "--in", self.out_genes, "--out", de, "--conditions", "ctrl,treat",
                 "--condition-file", cond])
        res = {r["locus_id"]: r for r in read_csv(os.path.join(de, "results_locus.csv"))}
        self.assertEqual(res["L1A_dup1:L1A:L1:LINE"]["overlapping_genes"], "GENE1")
        self.assertAlmostEqual(float(res["L1A_dup1:L1A:L1:LINE"]["log2FoldChange"]), 1.0, delta=0.35)
        self.assertTrue(os.path.exists(os.path.join(de, "summary_by_class_and_context_locus.csv")))

    @unittest.skipUnless(HAVE_DESEQ2, "Rscript/DESeq2 not available")
    def test_differential_expression(self):
        cond = os.path.join(self.tmp, "cond.csv")
        with open(cond, "w") as f:  # rows deliberately in a different order than the table
            f.write("sample,condition\ntreat2,treat\nctrl1,ctrl\ntreat1,treat\nctrl2,ctrl\n")
        de = os.path.join(self.tmp, "de")
        self.sh([SCRIPT, "test", "--in", self.out, "--out", de, "--level", "locus",
                  "--conditions", "ctrl,treat", "--condition-file", cond, "--min-count", "10"])
        res = {r["locus_id"]: r for r in read_csv(os.path.join(de, "results_locus.csv"))}
        lfc = float(res["L1A_dup1:L1A:L1:LINE"]["log2FoldChange"])
        self.assertAlmostEqual(lfc, 1.0, delta=0.35)  # 2x up in treat
        self.assertTrue(os.path.exists(os.path.join(de, "ma_plot_locus.pdf")))


if __name__ == "__main__":
    unittest.main()
