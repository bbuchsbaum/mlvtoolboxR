import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from PIL import Image

from score_curve_rank import gates, main, score, validate_roster


class CurveRankScoreTest(unittest.TestCase):
    def write_export(self, root, ident, method, drawing, sites=None):
        prefix = root / f"{ident}.{method}"
        Image.fromarray(np.asarray(drawing, dtype=np.uint8) * 255).save(
            str(prefix) + ".drawing.png"
        )
        (Path(str(prefix) + ".sites.csv")).write_text("x,y,type,directions\n")
        Path(str(prefix) + ".contours.json").write_text(
            json.dumps({"sites": sites or [], "contours": []})
        )

    def test_missing_candidate_is_explicit_and_fails_gate(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            staged, base, experiment = (
                root / "staged",
                root / "base",
                root / "experiment",
            )
            for folder in (staged, base, experiment):
                folder.mkdir()
            (root / "roster.csv").write_text("id,category,panel,fold\na,city,test,1\n")
            (staged / "panel.csv").write_text("id,category,panel\na,city,test\n")
            geometry = {
                "id": "a",
                "working": {"width": 4, "height": 4},
                "polylines": [],
                "junctions": [],
            }
            (staged / "a.json").write_text(json.dumps(geometry))
            picture = np.eye(4, dtype=bool)
            for method in ("artist", "fdog", "union_vector", "mlsd"):
                self.write_export(base, "a", method, picture)
            rows = score(validate_roster(root / "roster.csv"), staged, base, experiment)
            candidate = next(r for r in rows if r["method"] == "curve_rank")
            self.assertEqual(candidate["output_status"], "missing")
            self.assertFalse(gates(rows)["gate_pass"])
            report = root / "report"
            main(
                SimpleNamespace(
                    staged=staged,
                    base_outputs=base,
                    experiment_outputs=experiment,
                    roster=root / "roster.csv",
                    report=report,
                )
            )
            self.assertFalse(
                json.loads((report / "curve-rank-gate.json").read_text())["gate_pass"]
            )

    def test_roster_rejects_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "roster.csv"
            path.write_text("id,category,panel,fold\na,c,p,1\na,c,p,2\n")
            with self.assertRaisesRegex(ValueError, "distinct"):
                validate_roster(path)

    def test_reserve_test_allows_an_empty_fold_only_there(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "roster.csv"
            path.write_text("id,category,panel,fold\na,c,reserve-test,\n")
            self.assertEqual(validate_roster(path)[0]["fold"], "")
            path.write_text("id,category,panel,fold\na,c,development,\n")
            with self.assertRaisesRegex(ValueError, "empty fold"):
                validate_roster(path)


if __name__ == "__main__":
    unittest.main()
