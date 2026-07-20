#!/usr/bin/env bash
#
# Run HiFi-16S-workflow on the Gourgouthakas cave dataset, inside the
# gourgouthakas-16s container.
#
#   ./scripts/run_pipeline.sh download    # fetch SILVA + GTDB + GG2 (once)
#   ./scripts/run_pipeline.sh run         # the actual 51-library run
#   ./scripts/run_pipeline.sh resume      # resume an interrupted run
#   ./scripts/run_pipeline.sh shell       # interactive shell in the image
#
# The container gets the repo at /work and the database directory at
# /databases. Nothing is installed on the host beyond podman.
#
# Runs from either side of the container boundary:
#   * on the host  — launches podman with the mounts and flags set up below
#   * inside the image (e.g. a gourgouthakas-dev shell) — detects that and execs
#     nextflow directly, instead of nesting podman inside podman
#
# Env overrides:
#   IMAGE=localhost/gourgouthakas-dev   use a different image
#   DB_DIR=/path/to/databases           databases live elsewhere
#   DETACH=1                            background the run (use for full runs)

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${IMAGE:-localhost/gourgouthakas-16s:latest}"
DB_DIR="${DB_DIR:-/mnt/data/databases}"
ACTION="${1:-run}"

# Are we already inside a container? If so this script must NOT call podman —
# that would mean nested rootless podman, which is painful and unnecessary,
# since every tool the pipeline needs is already on PATH here. Instead the
# nextflow command is exec'd directly, so the same script works unchanged from
# the host or from a shell inside gourgouthakas-16s / gourgouthakas-dev.
if [ -f /run/.containerenv ] || [ -f /.dockerenv ] || [ -n "${container:-}" ]; then
    IN_CONTAINER=1
else
    IN_CONTAINER=0
fi

mkdir -p "$ROOT/results" "$ROOT/work"
[ "$IN_CONTAINER" = "0" ] && mkdir -p "$DB_DIR"

# --userns=keep-id maps the caller to the same uid inside the container, so
# results and the work dir come back owned by the host user.
run_in_container() {
    # Already inside the image: run it here. The pipeline env is prepended to
    # PATH by scripts/nextflow.config, exactly as it would be under podman.
    if [ "$IN_CONTAINER" = "1" ]; then
        if [ ! -d /databases ]; then
            echo "warning: /databases is not mounted — start the container with -v \$DB_DIR:/databases" >&2
        fi
        cd "$ROOT"
        exec bash -lc "$1"
    fi

    # Allocate a TTY only when there is one; otherwise podman fails outright
    # when this is run from a script, a cron job or a CI log.
    local tty_flag=()
    [ -t 0 ] && tty_flag=(-it)

    # This host runs SELinux in enforcing mode, which denies the container all
    # access to bind mounts unless they are labelled. label=disable is used in
    # preference to the :z mount flag because :z relabels the mount recursively
    # — that would rewrite the SELinux context of the whole 4.5 GB raw dataset
    # on every run. Disabling the label for these two mounts touches nothing on
    # disk.
    # HOME is redirected off /work because --userns=keep-id otherwise sets it to
    # the mounted repo, so tools scatter dotfiles (.wget-hsts, .local/, .cache/)
    # into the project -- and renv refuses to initialise a project that looks
    # like a home directory.
    # DETACH=1 runs the container in the background under podman's own conmon,
    # so the pipeline is not a child of the calling shell. A full 51-library run
    # takes hours; attached, it dies with the terminal/ssh session that started
    # it. Logs go to `podman logs -f gourgouthakas-run`.
    if [ "${DETACH:-0}" = "1" ]; then
        podman run -d --name gourgouthakas-run \
            --userns=keep-id \
            --security-opt label=disable \
            -e HOME=/tmp \
            -v "$ROOT:/work" \
            -v "$DB_DIR:/databases" \
            -w /work \
            "$IMAGE" \
            bash -lc "$1"
        echo "detached as 'gourgouthakas-run' — follow with: podman logs -f gourgouthakas-run"
        return
    fi

    podman run --rm "${tty_flag[@]}" \
        --userns=keep-id \
        --security-opt label=disable \
        -e HOME=/tmp \
        -v "$ROOT:/work" \
        -v "$DB_DIR:/databases" \
        -w /work \
        "$IMAGE" \
        bash -lc "$1"
}

# HiFi-16S-workflow/nextflow.config interpolates a bare $HOME for the conda
# cacheDir, which Nextflow 26.04's strict config parser rejects outright. The
# legacy parser accepts it, and is preferred over patching the upstream repo —
# that is a separate checkout we do not own. Remove once upstream quotes it as
# env('HOME').
NF_COMMON="NXF_SYNTAX_PARSER=v1 nextflow run HiFi-16S-workflow/main.nf -c scripts/nextflow.config"

case "$ACTION" in
    download)
        # Populates /databases with the <db>/nb and <db>/vsearch layout that
        # main.nf expects. Note this is NOT what fetch-16s-databases.sh
        # produces — that script serves the container's own flat DB_DIR.
        run_in_container "$NF_COMMON \
            --download_db \
            --download_targets silva,gtdb,gg2 \
            --db_base_dir /databases \
            --outdir results/db_download"
        ;;
    run|resume)
        [ "$ACTION" = resume ] && RESUME="-resume" || RESUME=""
        run_in_container "$NF_COMMON $RESUME \
            --input data/samplesheet.tsv \
            --metadata data/metadata.tsv \
            --db_base_dir /databases \
            --outdir results/hifi"
        ;;
    shell)
        run_in_container "exec bash"
        ;;
    *)
        echo "usage: $0 {download|run|resume|shell}" >&2
        exit 1
        ;;
esac
