#!/usr/bin/env bash
#
# test_init.sh - Check scripts/init_config.py against a fixture manifest
#
# Description:
#   Builds a throwaway manifest whose references are empty placeholder files,
#   then checks that init_config.py (1) refuses to guess when answers are
#   missing or unknown, and (2) writes configs whose reference paths
#   run_pipeline.sh reads back (same grep/sed parser) as the manifest's, with
#   every other setting identical to the shipped species config.
#
# Usage:
#   pixi run test-init
#
# Author: Brent Chick
# Date: 2026-09-28

set -euo pipefail

readonly REPO="$(cd "$(dirname "$0")/.." && pwd)"
T="$(mktemp -d)"
trap 'rm -rf "${T}"' EXIT

fail() { echo "FAIL: $*" >&2; exit 1; }

# ==============================================================================
# FIXTURE MANIFEST
# ==============================================================================
R="${T}/resource"
mkdir -p "${R}/hs.STAR" "${R}/mm.STAR"
touch "${R}"/{hs,mm}.{fa,gtf,bed12} "${R}/mm.hk.bed12"
cat > "${T}/manifest.yaml" << EOF
manifest_version: 1
updated: test
genomes:
  hg38:
    species: human
    aliases: [GRCh38, human]
    rnaseq:
      description: "fixture human"
      genome_fasta: ${R}/hs.fa
      gtf: ${R}/hs.gtf
      star_index: ${R}/hs.STAR
      gene_model_bed: ${R}/hs.bed12
      ercc: false
  mm39:
    species: mouse
    aliases: [GRCm39, mouse]
    rnaseq:
      description: "fixture mouse"
      genome_fasta: ${R}/mm.fa
      gtf: ${R}/mm.gtf
      star_index: ${R}/mm.STAR
      gene_model_bed: ${R}/mm.bed12
      genebody_bed: ${R}/mm.hk.bed12
      ercc: true
  mm10:
    species: mouse
    status: ok
    fasta: ${R}/mm.fa
EOF

init() { python "${REPO}/scripts/init_config.py" --manifest "${T}/manifest.yaml" "$@" < /dev/null; }

# Same parser as run_pipeline.sh get_config()
get_config() {
    grep -E "^\s*${2##*.}:" "$1" | head -1 | sed 's/.*:\s*//' | sed 's/\s*#.*//' | tr -d '"' | tr -d "'"
}

# ==============================================================================
# REFUSALS
# ==============================================================================
set +e
init --dir "${T}/x" > /dev/null 2>&1;                 [[ $? == 2 ]] || fail "missing --genome should exit 2"
init --dir "${T}/x" --genome mm10 > /dev/null 2>&1;   [[ $? == 1 ]] || fail "genome without rnaseq should be refused"
set -e
init --list | grep -q '^mm39' || fail "--list should show mm39"

# ==============================================================================
# GENERATED CONFIGS
# ==============================================================================
for combo in hg38:human:hs mouse:mouse:mm; do
    IFS=: read -r genome species px <<< "${combo}"
    d="${T}/proj_${species}"
    init --dir "${d}" --genome "${genome}" > /dev/null
    cfg="${d}/config.yaml"
    [[ "$(get_config "${cfg}" genome_fasta)" == "${R}/${px}.fa" ]] || fail "${species}: genome_fasta"
    [[ "$(get_config "${cfg}" star_index)" == "${R}/${px}.STAR" ]] || fail "${species}: star_index"
    [[ "$(get_config "${cfg}" gene_model_bed)" == "${R}/${px}.bed12" ]] || fail "${species}: gene_model_bed"
    [[ "$(get_config "${cfg}" organism)" == "${species}" ]] || fail "${species}: organism"
    [[ -f "${d}/samples.tsv" ]] || fail "${species}: samples.tsv not copied"
    # Every non-reference line equals the shipped config (after init's 7-line header)
    shipped="${REPO}/config/config.${species}.yaml"
    diff <(grep -vE '^\s+(genome_fasta|gtf|star_index|gene_model_bed|genebody_bed):' "${shipped}") \
         <(tail -n +8 "${cfg}" | grep -vE '^\s+(genome_fasta|gtf|star_index|gene_model_bed|genebody_bed):') \
        > /dev/null || fail "${species}: non-reference settings differ from ${shipped}"
done
[[ "$(get_config "${T}/proj_mouse/config.yaml" genebody_bed)" == "${R}/mm.hk.bed12" ]] || fail "mouse: genebody_bed"

echo "test_init: refusals OK, human and mouse configs OK"
