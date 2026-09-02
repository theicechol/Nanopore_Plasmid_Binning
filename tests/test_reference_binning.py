import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from nanopore_plasmid_bin_assemble_polish import classify_read, sequence_kmers
from plasmid_qc_report import HOST_CATALOG, classify_multimer


class ReferenceBinningTests(unittest.TestCase):
    def setUp(self):
        shared = "ACGT" * 30
        self.a_marker = "GATTACAAGGCTTACCGATCGTACGATGCTAGCTAGGCTA"
        self.b_marker = "TTGCAACCTTGGTACCAAGTTCGGAATCCGTTACCGGATC"
        k = 9
        all_kmers = {
            "a": sequence_kmers(shared + self.a_marker, k),
            "b": sequence_kmers(shared + self.b_marker, k),
        }
        self.unique = {
            "a": all_kmers["a"] - all_kmers["b"],
            "b": all_kmers["b"] - all_kmers["a"],
        }
        self.markers = {"a": {"AmpR": self.unique["a"]}, "b": {"origin": self.unique["b"]}}
        self.k = k

    def classify(self, sequence):
        record = SeqRecord(Seq(sequence), id="read")
        return classify_read(record, self.unique, self.markers, self.k, 3, 2.0)

    def test_assigns_unique_marker(self):
        assignment, _, _, markers = self.classify(self.a_marker)
        self.assertEqual(assignment, "a")
        self.assertEqual(markers, ["AmpR"])

    def test_reverse_complement_has_same_assignment(self):
        assignment, _, _, _ = self.classify(str(Seq(self.b_marker).reverse_complement()))
        self.assertEqual(assignment, "b")

    def test_shared_backbone_is_ambiguous(self):
        assignment, _, _, _ = self.classify("ACGT" * 20)
        self.assertEqual(assignment, "ambiguous")

    def test_multimer_classification(self):
        self.assertEqual(classify_multimer(7100, 7190, 0.15), "1-mer")
        self.assertEqual(classify_multimer(14362, 7190, 0.15), "2-mer")
        self.assertEqual(classify_multimer(3500, 7190, 0.15), "other")

    def test_expected_host_catalog_accessions(self):
        self.assertEqual(HOST_CATALOG["vibrio-natriegens"]["accessions"], ("CP009977.1", "CP009978.1"))
        self.assertEqual(HOST_CATALOG["bacillus-subtilis"]["accessions"], ("NC_000964.3",))


if __name__ == "__main__":
    unittest.main()
