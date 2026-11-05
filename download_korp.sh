#!/bin/bash
# Downloads the KORP energy file (korp6Dv1.bin, ~316 MB) needed by --energy_file.
# Usage: ./download_korp.sh [destination_dir]   (default: src/libenergy)
set -e
URL="https://chaconlab.org/downloads/korp6Dv1.bin"   # redirects to /download/
EXPECTED_SIZE=331777205
DEST="${1:-$(dirname "$(readlink -f "$0")")/src/libenergy}"
OUT="$DEST/korp6Dv1.bin"

mkdir -p "$DEST"
if [ -f "$OUT" ] && [ "$(stat -c %s "$OUT")" = "$EXPECTED_SIZE" ]; then
  echo "$OUT already present"; exit 0
fi
echo "Downloading $URL -> $OUT"
if command -v curl >/dev/null; then
  curl -L --fail -o "$OUT.part" "$URL"
elif command -v wget >/dev/null; then
  wget -O "$OUT.part" "$URL"
else
  echo "Error: curl or wget is required" >&2; exit 1
fi
if [ "$(stat -c %s "$OUT.part")" != "$EXPECTED_SIZE" ]; then
  echo "Error: unexpected file size (expected $EXPECTED_SIZE bytes)" >&2
  rm -f "$OUT.part"; exit 1
fi
mv "$OUT.part" "$OUT"
echo "Done: $OUT"
