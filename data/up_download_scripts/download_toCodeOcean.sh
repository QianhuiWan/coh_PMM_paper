#!/usr/bin/env bash
# Purpose: Fetch input data from Zenodo into /tmp at runtime, extract, run analysis, then auto-clean.
# Usage:
#   Edit DATA_URL below (or export it as an env var), then run this script inside your capsule.
# Optional env vars:
#   ZENODO_TOKEN     : Bearer token for restricted records (omit if public)
#   SHA256SUMS_URL   : URL to checksums file (if you want integrity verification)
#   TOPDIR           : Top-level directory inside the archive (if known). If unset, auto-detect.
#   RSCRIPT          : Rscript binary (default: Rscript)
#   ANALYSIS_SCRIPT  : R script path (default: Figure_1/activation_score_S4_getL1Num.R)
#   OUTDIR           : Output directory (default: results)

set -euo pipefail

# ====== Configuration ======
DATA_URL="${DATA_URL:-https://zenodo.org/records/<RECORD_ID>/files/input.tar.gz?download=1}"
RSCRIPT="${RSCRIPT:-Rscript}"
ANALYSIS_SCRIPT="${ANALYSIS_SCRIPT:-Figure_1/activation_score_S4_getL1Num.R}"
OUTDIR="${OUTDIR:-results}"

# ====== Prepare temp workspace (auto-clean on exit/error) ======
TMPD="$(mktemp -d -p /tmp zen_in.XXXXXX)"
trap 'rm -rf "$TMPD"' EXIT

# Build curl auth header if token is provided
CURL_AUTH=()
if [[ -n "${ZENODO_TOKEN:-}" ]]; then
  CURL_AUTH=(-H "Authorization: Bearer ${ZENODO_TOKEN}")
fi

echo "[1/3] Fetching data into $TMPD ..."

if [[ -n "${SHA256SUMS_URL:-}" ]]; then
  # -- Integrity path: download archive to a file (needed for checksum)
  echo " - Downloading archive file (with retries) ..."
  curl -L --fail --retry 5 --retry-delay 3 "${CURL_AUTH[@]}" \
       -o "$TMPD/input.tar.gz" "$DATA_URL"

  echo " - Downloading SHA256SUMS ..."
  curl -L --fail --retry 5 --retry-delay 3 "${CURL_AUTH[@]}" \
       -o "$TMPD/SHA256SUMS" "$SHA256SUMS_URL"

  echo " - Verifying checksums ..."
  ( cd "$TMPD" && sha256sum -c SHA256SUMS )

  echo " - Extracting archive ..."
  tar -xzf "$TMPD/input.tar.gz" -C "$TMPD"
  # Remove big archive immediately to save /tmp space
  rm -f "$TMPD/input.tar.gz"
else
  # -- Fast path: stream extract (no intermediate 20GB file)
  echo " - Streaming download & extract (with retries) ..."
  # Note: --retry with a pipe requires curl 7.71+ to retry on transient errors; still helpful.
  curl -L --fail --retry 5 --retry-delay 3 "${CURL_AUTH[@]}" \
      "$DATA_URL" | tar -xz -C "$TMPD"
fi

# Determine the data root inside the archive.
# If TOPDIR is provided, use it; otherwise pick the single top-level directory if there is exactly one.
if [[ -n "${TOPDIR:-}" ]]; then
  DATA_PATH="$TMPD/$TOPDIR"
else
  # Auto-detect: if extraction produced exactly one top-level dir, use it; else default to TMPD.
  mapfile -t _top < <(find "$TMPD" -mindepth 1 -maxdepth 1 -type d | sort)
  if [[ ${#_top[@]} -eq 1 ]]; then
    DATA_PATH="${_top[0]}"
  else
    DATA_PATH="$TMPD"
  fi
fi

echo "[2/3] Running analysis ..."
mkdir -p "$OUTDIR"
"$RSCRIPT" "$ANALYSIS_SCRIPT" --input "$DATA_PATH" --out "$OUTDIR"

echo "[3/3] Done. Temporary files will be removed automatically."
