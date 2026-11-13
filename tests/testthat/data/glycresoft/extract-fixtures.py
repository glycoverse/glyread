"""Reproduce unmodified fixture lines from the two cited Figshare ZIPs.

Usage: python3 extract-fixtures.py PXD009654.zip yeast.zip output-directory
Download URLs, attribution, checksums, and license are in SOURCE.md.
"""
import csv
import hashlib
import io
import pathlib
import sys
import zipfile

pxd, yeast, output = map(pathlib.Path, sys.argv[1:])
output.mkdir(parents=True, exist_ok=True)
archives = [
    (pxd, "8f0fca7aebdac30f63af12c205c3210c", {
        "20160301_Serum_SA_glycopeptide_HCC_2.glycopeptides.csv": ("HCC_2.csv", {"1", "81", "118", "371"}),
        "20160301_Serum_SA_glycopeptide_HCC_3.glycopeptides.csv": ("HCC_3.csv", {"1", "289"}),
        "20160301_Serum_SA_glycopeptide_Normal_3.glycopeptides.csv": ("Normal_3.csv", {"1", "331", "372"}),
    }),
    (yeast, "5b0790a2dc71eef22c22b649e4ff2731", {
        "cwq_mix2-1_726.glycopeptides.csv": ("cwq_mix2-1.csv", {"1", "481", "832"}),
    }),
]
for archive, checksum, selected in archives:
    assert hashlib.md5(archive.read_bytes()).hexdigest() == checksum
    with zipfile.ZipFile(archive) as z:
        for member, (name, groups) in selected.items():
            content = z.read(member)
            rows = list(csv.DictReader(io.StringIO(content.decode("utf-8"))))
            lines = content.splitlines(keepends=True)
            assert len(lines) == len(rows) + 1  # These source files have no multiline fields.
            chosen = [line for line, row in zip(lines[1:], rows) if row["group_id"] in groups]
            assert chosen
            result = lines[0] + b"".join(chosen)
            (output / name).write_bytes(result)
            print(name, hashlib.sha256(result).hexdigest())
