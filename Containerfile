# Gourgouthakas cave 16S — complete analysis environment, self-contained.
#
# Builds from a public base with no dependency on any locally-built image, so
# it reproduces on any machine with podman or docker:
#
#   podman build -t gourgouthakas-16s .          # ~45 min
#   ./scripts/run_pipeline.sh download
#   DETACH=1 ./scripts/run_pipeline.sh run
#
# This deliberately replaces the previous fedora:44 -> devbase -> r-tidyverse ->
# omics-16s -> here chain. Those images are personal dev environments: they
# carry neovim, tmux, w3m and an unpinned `npm install -g` of three AI CLIs,
# none of which this analysis uses, and the unpinned installs mean the image
# could not be rebuilt identically even on the original machine.
#
# The image holds TWO R environments on purpose:
#   * system R           — downstream analysis (dada2 1.40.0, vegan, ALDEx2, ...)
#   * /opt/conda/envs/hifi16s — HiFi-16S-workflow tasks, pinned to the versions in
#     that repo's env/*.yml (dada2 1.38.0). Kept OFF the default PATH;
#     scripts/nextflow.config prepends it for task execution only.
# Merging them silently changes the DADA2 version the ASVs are inferred with.

FROM fedora:44

LABEL org.opencontainers.image.title="gourgouthakas-16s" \
      org.opencontainers.image.description="PacBio HiFi 16S analysis environment for the Gourgouthakas cave dataset"

# ---------------------------------------------------------------------------
# System packages: R, a full compiler toolchain (gfortran is required by
# glmnet/pcaPP/Rfast and hence by the Maaslin2 dependency tree), and the -devel
# headers the tidyverse and Bioconductor stacks compile against.
#
# Java is pinned to 25 — the version Fedora 44 ships and the one Nextflow 26.04
# is known to run on here. Not java-latest-openjdk, which is already 26 and
# would move under you on every rebuild.
# ---------------------------------------------------------------------------
RUN dnf install -y \
        bash git gawk curl wget ca-certificates findutils procps-ng tar unzip bzip2 \
        R R-devel \
        gcc gcc-c++ gcc-gfortran make cmake \
        libcurl-devel openssl-devel libxml2-devel \
        fontconfig-devel freetype-devel harfbuzz-devel fribidi-devel \
        libjpeg-turbo-devel libpng-devel libtiff-devel libuv-devel mbedtls-devel \
        zlib-devel bzip2-devel xz-devel \
        java-25-openjdk-headless \
    && dnf clean all \
    && rm -rf /var/cache/dnf

# Pin a working CRAN mirror. The previous base image inherited an /root/.Rprofile
# pinning a University of Crete mirror that went dead; because a user .Rprofile
# is sourced AFTER Rprofile.site it silently overrode everything, and
# BiocManager then reported every CRAN dependency as "not available" — failing as
# limma -> edgeR -> DESeq2 -> ALDEx2 -> Maaslin2, which looks like a
# Bioconductor problem rather than a mirror problem. The assertion below turns a
# bad mirror into an immediate, obvious build failure.
RUN printf 'options(repos = c(CRAN = "https://cloud.r-project.org"))\n' \
        > /usr/lib64/R/etc/Rprofile.site \
    && Rscript -e "stopifnot(length(rownames(available.packages())) > 20000); \
                   cat('CRAN mirror reachable\n')"

# ---------------------------------------------------------------------------
# Analysis R packages (PLAN.md section 5).
#
# Not installed, deliberately (see CLAUDE.md):
#   phyloseq / qiime2R  — banned by PLAN.md; the data model is tidy tables
#   picante / GUniFrac  — no phylogeny is built, so Faith's PD and UniFrac
#                         are out of scope
#   sf / terra / usdm   — every sample shares one lat/lon (a single cave), so
#                         the gradient is depth, not geography; usdm would drag
#                         in GDAL/PROJ/GEOS for one function car::vif provides
# ---------------------------------------------------------------------------
RUN Rscript -e "\
    install.packages(c( \
        'tidyverse', 'here', 'renv', 'data.table', 'scales', \
        'vegan', 'ape', 'SRS', \
        'compositions', 'zCompositions', \
        'car', 'patchwork', 'ggpubr', 'ggrepel', 'pheatmap', \
        'BiocManager' \
    ), Ncpus=parallel::detectCores())"

RUN Rscript -e "\
    BiocManager::install(c( \
        'dada2', 'Biostrings', 'biomformat', \
        'decontam', 'ALDEx2', 'DESeq2', 'Maaslin2' \
    ), update=FALSE, ask=FALSE, Ncpus=parallel::detectCores())"

# install.packages() only warns on failure, so assert explicitly: a missing
# package must fail the build, not surface hours later mid-analysis.
RUN Rscript -e "\
    pkgs <- c('tidyverse','here','renv','vegan','ape','SRS','compositions', \
              'zCompositions','car','patchwork','ggpubr','ggrepel','pheatmap', \
              'dada2','Biostrings','biomformat','decontam','ALDEx2','DESeq2', \
              'Maaslin2'); \
    missing <- pkgs[!pkgs %in% rownames(installed.packages())]; \
    if (length(missing)) stop('Failed to install: ', paste(missing, collapse=', ')); \
    cat('all', length(pkgs), 'analysis packages present\n')"

# ---------------------------------------------------------------------------
# Pipeline tools, pinned to HiFi-16S-workflow/env/*.yml.
#
# csvtk is the one deliberate deviation: upstream uses 0.28.0 in
# inspect_metadata and 0.31.0 in collect_qc; 0.31.0 covers both.
#
# GNU wget is here because Fedora ships wget2 as /usr/bin/wget, and the database
# download processes call `wget -O`. wget2 hangs on zenodo.org outright, and
# neither wget does Happy Eyeballs — they try zenodo's IPv6 address first and
# block until timeout, since a rootless container typically has no IPv6 route.
# The symptom is a silent stall: retries forever, zero bytes, nothing on stdout.
# ---------------------------------------------------------------------------
ARG DADA2_PIPELINE_VERSION=1.38.0
ARG CUTADAPT_VERSION=5.2
ARG SEQKIT_VERSION=2.13.0
ARG CSVTK_VERSION=0.31.0
ARG VSEARCH_VERSION=2.31.0
ENV MAMBA_ROOT_PREFIX=/opt/conda

RUN set -eux; \
    case "$(uname -m)" in \
        x86_64)        mm_arch="linux-64" ;; \
        aarch64|arm64) mm_arch="linux-aarch64" ;; \
        *) echo "unsupported architecture: $(uname -m)" >&2; exit 1 ;; \
    esac; \
    curl -L "https://micro.mamba.pm/api/micromamba/${mm_arch}/latest" \
        | tar -xj -C /usr/local bin/micromamba; \
    micromamba create -y -n hifi16s -c conda-forge -c bioconda \
        "bioconductor-dada2=${DADA2_PIPELINE_VERSION}" \
        "cutadapt=${CUTADAPT_VERSION}" \
        "seqkit=${SEQKIT_VERSION}" \
        "csvtk=${CSVTK_VERSION}" \
        "vsearch=${VSEARCH_VERSION}" \
        r-data.table wget jq gawk; \
    micromamba clean -a -y; \
    /opt/conda/envs/hifi16s/bin/cutadapt --version; \
    /opt/conda/envs/hifi16s/bin/seqkit version; \
    /opt/conda/envs/hifi16s/bin/Rscript -e \
        "cat('pipeline dada2', as.character(packageVersion('dada2')), '\n')"

ENV HIFI16S_ENV=/opt/conda/envs/hifi16s

# Force IPv4 for wget rather than patching upstream's `wget -O` calls. This path
# is GNU wget's compiled-in *system* config for this env, so no WGETRC variable
# is needed — setting one too makes wget warn on every invocation.
RUN printf 'prefer-family = IPv4\ntries = 5\ntimeout = 60\n' \
        > /opt/conda/envs/hifi16s/etc/wgetrc

# ---------------------------------------------------------------------------
# Nextflow. NXF_HOME is under /opt so a bind-mounted /work is not polluted with
# plugin caches, and plugins are fetched at build time so a run needs no network.
# ---------------------------------------------------------------------------
ARG NEXTFLOW_VERSION=26.04.6
ENV NXF_HOME=/opt/nextflow
RUN set -eux; \
    mkdir -p /opt/nextflow; \
    curl -s https://get.nextflow.io | NXF_VER=${NEXTFLOW_VERSION} bash; \
    mv nextflow /usr/local/bin/nextflow; \
    chmod +x /usr/local/bin/nextflow; \
    nextflow -version; \
    chmod -R a+rwX /opt/nextflow

WORKDIR /work
CMD ["/bin/bash"]
