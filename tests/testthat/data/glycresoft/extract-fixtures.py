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
    (pxd, "8bf9b5e77116034398ff87d16e1c7992d7d2e3b41db24e3269db207ddb33b4d2", {
        "20160301_Serum_SA_glycopeptide_HCC_2.glycopeptides.csv": ("HCC_2.csv", {"1", "81", "118", "371"}),
        "20160301_Serum_SA_glycopeptide_HCC_3.glycopeptides.csv": ("HCC_3.csv", {"1", "289"}),
        "20160301_Serum_SA_glycopeptide_Normal_3.glycopeptides.csv": ("Normal_3.csv", {"1", "331", "372"}),
    }),
    (yeast, "c88f01f8fb78a9a95c9de67a7e54db09dcd95b7bf9a5939a788f5e66bcda6c35", {
        "cwq_mix2-1_726.glycopeptides.csv": ("cwq_mix2-1.csv", {"1", "481", "832"}),
    }),
]
for archive, checksum, selected in archives:
    if hashlib.sha256(archive.read_bytes()).hexdigest() != checksum:
        raise ValueError(f"SHA-256 mismatch for {archive}")
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
