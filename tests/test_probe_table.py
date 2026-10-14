import csv
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class TestBridgeProbeTable(unittest.TestCase):
    def test_runtime_table_matches_update_csv(self):
        with (ROOT / "Bridge_Probes1-2578.csv").open(
            encoding="latin-1", newline=""
        ) as handle:
            expected = [
                (row["New ID"], row["Seq"])
                for row in csv.DictReader(handle)
            ]

        actual = [
            tuple(line.split())
            for line in (ROOT / "resources/probe_table.txt").read_text().splitlines()
        ]

        expected_ids = [f"BP_{number:04d}" for number in range(1, 1001)]
        expected_ids += [f"BP_{number:04d}" for number in range(1101, 2579)]

        self.assertEqual(actual, expected)
        self.assertEqual([bp_id for bp_id, _ in actual], expected_ids)
        self.assertTrue(all(re.fullmatch(r"[ACGT]{19}", seq) for _, seq in actual))
        self.assertEqual(len({seq for _, seq in actual}), len(actual))


if __name__ == "__main__":
    unittest.main()
