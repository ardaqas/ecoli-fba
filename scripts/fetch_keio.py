"""Download the Keio collection supplementary tables into data/raw/keio/.

Usage:
    python scripts/fetch_keio.py

The tables are not committed to the repository. This script downloads the
supplementary archives of the two papers from Europe PMC, extracts only the
tables the analysis needs and checks each one against a known SHA-256, so
everybody works with exactly the same files.

The Europe PMC download can be slow (several minutes for the 2006 archive).
If it keeps failing, download the files listed in ``ecoli_fba.keio.RAW_FILES``
by hand from the articles' supplementary sections and put them in
data/raw/keio/; rerunning this script then only verifies the checksums.
"""

import hashlib
import shutil
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

from ecoli_fba.keio import RAW_DIR, RAW_FILES, SOURCE_ARCHIVES

CHUNK = 64 * 1024
TIMEOUT_S = 120
ATTEMPTS = 3


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def is_valid(filename: str) -> bool:
    path = RAW_DIR / filename
    return path.is_file() and sha256(path) == RAW_FILES[filename].sha256


def download(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
    with urllib.request.urlopen(request, timeout=TIMEOUT_S) as response, open(destination, "wb") as out:
        done = 0
        while block := response.read(CHUNK):
            out.write(block)
            done += len(block)
            print(f"\r  {done / 1e6:.1f} MB", end="", flush=True)
    print()


def fetch_archive(source: str, url: str, filenames: list[str]) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        archive = Path(tmp) / f"{source}.zip"
        for attempt in range(1, ATTEMPTS + 1):
            print(f"Downloading {source} supplementary files (attempt {attempt}/{ATTEMPTS})")
            try:
                download(url, archive)
                with zipfile.ZipFile(archive) as zf:
                    for filename in filenames:
                        with zf.open(filename) as member, open(RAW_DIR / filename, "wb") as out:
                            shutil.copyfileobj(member, out)
                return
            except (OSError, zipfile.BadZipFile, KeyError) as error:
                print(f"  failed: {error!r}")
    raise SystemExit(f"Could not download {source}. See the module docstring for the manual route.")


def main() -> int:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    for source, url in SOURCE_ARCHIVES.items():
        filenames = [name for name, info in RAW_FILES.items() if info.source == source]
        if all(is_valid(name) for name in filenames):
            print(f"{source}: already present, checksums OK")
            continue
        if not all((RAW_DIR / name).is_file() for name in filenames):
            fetch_archive(source, url, filenames)

    bad = [name for name in RAW_FILES if not is_valid(name)]
    for name in RAW_FILES:
        print(f"{'MISMATCH' if name in bad else 'ok':8s} {name}  ({RAW_FILES[name].description})")
    if bad:
        print("Checksum mismatch: these files differ from the ones the analysis was built on.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
