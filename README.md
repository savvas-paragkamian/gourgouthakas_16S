# Amplicon analysis of Gourgouthakas cave microbiome

Full-length 16S rRNA (PacBio HiFi / Kinnex) from Gourgouthakas cave, Crete —
a 1100 m deep vertical cave. 51 libraries: sediment and water down a depth
transect, plus mock and blank controls.

Upstream ASV inference is `HiFi-16S-workflow` (Nextflow + DADA2); downstream is
a tidyverse + vegan analysis. Everything runs in one podman container — the
host needs **only podman**.

## Setup

```bash
git clone --recurse-submodules https://github.com/savvas-paragkamian/gourgouthakas_16S.git
cd gourgouthakas_16S

podman build -t gourgouthakas-16s .        # ~45 min, mostly compiling R packages
./scripts/run_pipeline.sh download         # once: SILVA + GTDB + GG2 (11 GB)
```

`HiFi-16S-workflow` is a submodule pinned to the exact upstream revision this
analysis was run against — `--recurse-submodules` is required, or it clones as
an empty directory. If you already cloned without it:

```bash
git submodule update --init
```

Databases go to `/mnt/data/databases` by default; override with `DB_DIR=...`.

Raw reads are not in git (4.5 GB). Copy `data/PB482_SP/` across separately and
verify it:

```bash
cd data && md5sum -c md5sum-16s.txt        # expect 51/51 OK
```

## Run the pipeline

```bash
Rscript scripts/build_inputs.R             # regenerate samplesheet + metadata
DETACH=1 ./scripts/run_pipeline.sh run     # 51 libraries, several hours
podman logs -f gourgouthakas-run           # follow progress
```

`DETACH=1` runs it under podman rather than your shell, so it survives the
terminal closing — use it for any full run. Use `resume` instead of `run` to
continue an interrupted run; completed tasks are reused from cache.

```bash
./scripts/run_pipeline.sh resume
./scripts/run_pipeline.sh shell            # interactive shell in the container
```

Results land in `results/hifi/`.

## Analysis

R packages are pinned in `renv.lock`. Run R inside the container so renv
activates from the project root:

```bash
podman run --rm -it --userns=keep-id --security-opt label=disable \
  -e HOME=/tmp -v "$PWD:/work" -w /work localhost/gourgouthakas-16s R
```

Analysis scripts are `scripts/00_setup.R` … `scripts/10_figures.R`, each
runnable standalone and in order. They read `data/` and earlier `results/`,
and write only to `results/` and `plots/`.

## Layout

```
Containerfile          the image (layers on omics-16s.Containerfile)
scripts/               build_inputs.R, run_pipeline.sh, nextflow.config, 00–10
data/                  raw fastq (immutable), MIxS sheets, samplesheet + metadata
HiFi-16S-workflow/     upstream pipeline — separate repo, do not edit
results/ plots/        outputs
```

See `AGENTS.md` for dataset specifics, the departures from `PLAN.md`
(no phylogeny; depth, not geography, is the gradient), and the environment
traps already worked around (SELinux, wget2/IPv6).
