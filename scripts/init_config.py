#!/usr/bin/env python3
"""
init_config.py - Write a project config from the lab's reference manifest

Asks which genome the samples are from, then writes <project>/config.yaml: a
copy of config/config.human.yaml or config/config.mouse.yaml whose reference
paths come from the resource manifest (/data/resource/manifest.yaml on the
Salk server, genome entry `rnaseq`). A copy of config/samples.tsv is placed
next to it to edit.

RNA-seq needs no blacklist, IgG control or spike-in choice. The ERCC92
spike-ins are built into the mouse reference (the Plasmidsaurus portal's) and
counted as features; the human reference has none.

Usage:
    pixi run init                                    # interactive
    pixi run init --dir /data/<user>/<project> --genome mm39
    pixi run init --list

The manifest is --manifest, else $MCBLA_RESOURCE_MANIFEST, else
/data/resource/manifest.yaml.
"""

import argparse
import datetime
import json
import os
import re
import shutil
import sys

import yaml

# ==============================================================================
# CONSTANTS
# ==============================================================================
DEFAULT_MANIFEST = "/data/resource/manifest.yaml"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
# pixi runs tasks from the repo root; INIT_CWD is where the user typed the command
CWD = os.environ.get("INIT_CWD") or os.getcwd()
# manifest `species` -> shipped config holding every non-reference setting
SPECIES_CONFIG = {"human": "config.human.yaml", "mouse": "config.mouse.yaml"}
REF_KEYS = ("genome_fasta", "gtf", "star_index", "gene_model_bed", "genebody_bed")


# ==============================================================================
# FUNCTIONS
# ==============================================================================
def die(msg, code=1):
    """Print an error and exit."""
    print(f"init: {msg}", file=sys.stderr)
    sys.exit(code)


def find_manifest(arg):
    """First existing manifest of --manifest, $MCBLA_RESOURCE_MANIFEST, default."""
    for p in (arg, os.environ.get("MCBLA_RESOURCE_MANIFEST"), DEFAULT_MANIFEST):
        if p and os.path.exists(p):
            return p
    die(
        "no resource manifest found (looked at --manifest, $MCBLA_RESOURCE_MANIFEST, "
        f"{DEFAULT_MANIFEST}). Outside the Salk server, copy "
        "config/config.template.yaml and fill in your reference paths by hand."
    )


def rnaseq_genomes(man, allow_unverified):
    """Manifest genomes that have an RNA-seq (STAR) reference for a known species."""
    return {
        k: g
        for k, g in man["genomes"].items()
        if g.get("rnaseq")
        and g.get("species") in SPECIES_CONFIG
        and (g.get("status", "ok") == "ok" or allow_unverified)
    }


def resolve_genome(genomes, answer):
    """Genome key for a key or alias (case-insensitive)."""
    a = answer.lower()
    for key, g in genomes.items():
        if a == key.lower() or a in (x.lower() for x in g.get("aliases", [])):
            return key
    die(f"unknown genome '{answer}'; choose one of: {', '.join(genomes)}")


def ask(question, options, default=None):
    """Numbered-choice prompt. options: list of (value, description)."""
    print(f"\n{question}")
    for i, (val, desc) in enumerate(options, 1):
        mark = "  (default)" if val == default else ""
        print(f"  {i}. {val:<8} {desc}{mark}")
    while True:
        raw = input(f"Choice [1-{len(options)}]: ").strip()
        if not raw and default:
            return default
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return options[int(raw) - 1][0]
        for val, _ in options:
            if raw.lower() == val.lower():
                return val
        print("  Please enter one of the numbers above.")


def ask_text(question, default):
    """Free-text prompt with a default."""
    raw = input(f"\n{question} [{default}]: ").strip()
    return raw or default


def describe(key, g):
    """One-line description of a genome's RNA-seq reference."""
    return f"{g.get('species')}: {g['rnaseq'].get('description', key)}"


def render_config(template, rs, header):
    """Shipped species config with its reference paths replaced by the manifest's."""
    out = []
    for line in template.splitlines():
        m = re.match(r"^(\s+)(" + "|".join(REF_KEYS) + r"):", line)
        if m:
            value = rs.get(m.group(2))
            if value:
                line = f'{m.group(1)}{m.group(2)}: "{value}"'
            else:
                line = f"{m.group(1)}# {m.group(2)}: not in the manifest for this genome"
        out.append(line)
    # genebody_bed may be absent from the template but present in the manifest
    if rs.get("genebody_bed") and not any(
        re.match(r"^\s+genebody_bed:", x) for x in out
    ):
        i = next(n for n, x in enumerate(out) if re.match(r"^\s+gene_model_bed:", x))
        out.insert(i + 1, f'  genebody_bed: "{rs["genebody_bed"]}"')
    return header + "\n".join(out) + "\n"


# ==============================================================================
# MAIN
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--dir", help="project directory (created if needed)")
    ap.add_argument("--genome", help="hg38 / human or mm39 / mouse (manifest key or alias)")
    ap.add_argument("--fastq-dir", help="FASTQ directory, used in the printed run command")
    ap.add_argument("--manifest", help=f"resource manifest (default {DEFAULT_MANIFEST})")
    ap.add_argument(
        "--allow-unverified",
        action="store_true",
        help="also offer entries not marked status: ok",
    )
    ap.add_argument("--force", action="store_true", help="overwrite an existing config.yaml")
    ap.add_argument("--list", action="store_true", help="print the choices and exit")
    ap.add_argument("--json", action="store_true", help="with --list: machine-readable")
    a = ap.parse_args()

    mpath = find_manifest(a.manifest)
    with open(mpath) as fh:
        man = yaml.safe_load(fh)
    genomes = rnaseq_genomes(man, a.allow_unverified)

    if a.list:
        if a.json:
            print(json.dumps({k: describe(k, g) for k, g in genomes.items()}, indent=2))
        else:
            for k, g in genomes.items():
                print(f"{k:<6} {describe(k, g)}")
        return

    interactive = sys.stdin.isatty()
    needed = [f for f, v in (("--dir", a.dir), ("--genome", a.genome)) if v is None]
    if needed and not interactive:
        die(
            "not a terminal, so every answer must be a flag. Missing: "
            + ", ".join(needed)
            + ".\nAsk the user which genome (human or mouse) -- do not guess. "
            "`pixi run init --list` shows the options.",
            2,
        )

    if a.dir is None:
        a.dir = ask_text("Project directory (created if needed)", ".")
    proj = os.path.abspath(os.path.join(CWD, a.dir))
    out = os.path.join(proj, "config.yaml")
    if os.path.exists(out) and not a.force:
        die(f"{out} exists; use --force to overwrite")

    if a.genome is None:
        a.genome = ask(
            "Which genome are the samples from?",
            [(k, describe(k, g)) for k, g in genomes.items()],
        )
    gkey = resolve_genome(genomes, a.genome)
    g = genomes[gkey]
    rs = g["rnaseq"]

    missing = [rs[k] for k in REF_KEYS if rs.get(k) and not os.path.exists(rs[k])]
    if missing:
        die(
            "these resources are missing on disk (tell the manifest maintainer):\n  "
            + "\n  ".join(missing)
        )

    with open(os.path.join(REPO, "config", SPECIES_CONFIG[g["species"]])) as fh:
        template = fh.read()
    header = (
        "# ==============================================================================\n"
        f"# Project config -- written by `pixi run init` on {datetime.date.today()}\n"
        f"# from {mpath} (manifest version {man['manifest_version']}, "
        f"updated {man.get('updated')})\n"
        f"# Genome: {gkey} -- {describe(gkey, g)}\n"
        f"# Based on config/{SPECIES_CONFIG[g['species']]}; edit edger.contrasts below.\n"
        "# ==============================================================================\n\n"
    )
    os.makedirs(proj, exist_ok=True)
    with open(out, "w") as fh:
        fh.write(render_config(template, rs, header))
    samples = os.path.join(proj, "samples.tsv")
    copied = not os.path.exists(samples)
    if copied:
        shutil.copy(os.path.join(REPO, "config", "samples.tsv"), samples)

    fastq = a.fastq_dir or "/path/to/fastq"
    print(f"\nWrote {out}")
    if copied:
        print(f"Copied the example samples.tsv to {samples} -- replace its rows.")
    print(
        "Set edger.contrasts in config.yaml (\"<condition>-<reference condition>\")."
    )
    print(
        f"\nNext, from {REPO}:\n"
        f"  pixi run bash run_pipeline.sh -i {fastq} -o {proj}/results "
        f"-c {out} -m {samples} --dry-run"
    )


if __name__ == "__main__":
    main()
