import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from no_reference_plasmid_analysis import collapse_multimers


class NoReferenceAnalysisTests(unittest.TestCase):
    def test_dimer_peak_is_collapsed_into_monomer(self):
        monomer = "ACGTTGCAGATCCGATGCTAGCTAGGCTAACCGT" * 70
        dimer = monomer + monomer
        monomer_record = SeqRecord(Seq(monomer), id="m")
        dimer_record = SeqRecord(Seq(dimer), id="d")
        monomer_record.letter_annotations["phred_quality"] = [40] * len(monomer)
        dimer_record.letter_annotations["phred_quality"] = [40] * len(dimer)
        peaks = [
            {"center": len(monomer), "records": [monomer_record], "count": 20},
            {"center": len(dimer), "records": [dimer_record], "count": 10},
        ]
        clusters = collapse_multimers(peaks)
        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0]["multimer_peaks"][0]["copies"], 2)


if __name__ == "__main__":
    unittest.main()
