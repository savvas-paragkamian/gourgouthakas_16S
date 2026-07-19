# Downstream R analysis environment for the Gourgouthakas cave 16S dataset.
#
# Layered on omics-16s (which already carries tidyverse, dada2, vegan, SRS,
# Biostrings, vsearch, fastp/fastplong). This image adds only the packages
# PLAN.md needs for the ecology analysis, so rebuilds are cheap.
#
# Deliberately NOT installed:
#   - sf / terra / spatial stack — every sediment sample shares one lat/lon
#     (a single cave), so the gradient is depth, not geography. See CLAUDE.md.
#   - picante / GUniFrac — the HiFi workflow emits no phylogeny and we opted
#     out of building one, so Faith's PD and UniFrac are out of scope.
#   - phyloseq / qiime2R — PLAN.md is explicitly free of both.
#
# Build (from the repo root):
#   podman build -t gourgouthakas-16s .
# Run:
#   podman run -it --rm -v $(pwd):/work -w /work gourgouthakas-16s

FROM localhost/omics-16s:latest

# The base image inherits /root/.Rprofile from the host dotfiles, which pins
# CRAN to a University of Crete mirror that is dead — its index 404s and it
# advertises zero packages. A user .Rprofile is sourced AFTER Rprofile.site, so
# it wins over any site-level setting and must be replaced outright.
#
# Left in place it silently breaks BiocManager: the Bioconductor deps resolve
# fine while every CRAN dep ("statmod", "locfit", "Rfast", ...) comes back
# "not available", cascading into limma -> edgeR -> DESeq2 -> ALDEx2 -> Maaslin2.
RUN printf 'options(repos = c(CRAN = "https://cloud.r-project.org"))\n' \
        > /root/.Rprofile \
    && printf 'options(repos = c(CRAN = "https://cloud.r-project.org"))\n' \
        >> /usr/lib64/R/etc/Rprofile.site \
    && Rscript -e "stopifnot(getOption('repos')[['CRAN']] == 'https://cloud.r-project.org'); \
                   stopifnot(length(rownames(available.packages())) > 20000); \
                   cat('CRAN mirror reachable\n')"

# CRAN packages. Ncpus keeps the source builds parallel; the install is
# verified below rather than trusting install.packages' silent failures.
#
# usdm (VIF) is deliberately absent: it pulls terra, hence GDAL/PROJ/GEOS,
# for a single function that car::vif already provides.
RUN Rscript -e "\
    install.packages(c( \
        'ape', \
        'compositions', \
        'zCompositions', \
        'patchwork', \
        'ggpubr', \
        'ggrepel', \
        'car', \
        'pheatmap' \
    ), Ncpus=parallel::detectCores())"

# Bioconductor packages. All of these take plain matrices / data frames —
# no phyloseq object is ever constructed.
RUN Rscript -e "\
    BiocManager::install(c( \
        'biomformat', \
        'decontam', \
        'ALDEx2', \
        'DESeq2', \
        'Maaslin2' \
    ), update=FALSE, ask=FALSE, Ncpus=parallel::detectCores())"

# Fail the build loudly if anything above silently did not install.
RUN Rscript -e "\
    pkgs <- c('ape','compositions','zCompositions','patchwork','ggpubr', \
              'ggrepel','car','pheatmap','biomformat','decontam', \
              'ALDEx2','DESeq2','Maaslin2','vegan','SRS','dada2','tidyverse'); \
    missing <- pkgs[!pkgs %in% rownames(installed.packages())]; \
    if (length(missing)) stop('Failed to install: ', paste(missing, collapse=', ')); \
    cat('All', length(pkgs), 'packages present\n')"

# ---------------------------------------------------------------------------
# Upstream pipeline layer: Nextflow + the HiFi-16S-workflow process tools.
#
# The workflow is run from *inside* this image rather than from the host, so
# there is no nested containerisation and no conda-env creation at run time:
# nextflow executes locally with both `enable_conda` and `enable_container`
# false, and every process finds its tool on PATH.
#
# Versions are pinned to exactly what HiFi-16S-workflow/env/*.yml declares, so
# the run matches upstream despite bypassing the conda profile. The one
# deliberate deviation is csvtk: upstream uses 0.28.0 in inspect_metadata and
# 0.31.0 in collect_qc; 0.31.0 covers both (the subcommands used are stable).
# ---------------------------------------------------------------------------
ARG DADA2_PIPELINE_VERSION=1.38.0
ARG CUTADAPT_VERSION=5.2
ARG SEQKIT_VERSION=2.13.0
ARG CSVTK_VERSION=0.31.0

# The pipeline's dada2 lives in its own env, kept OFF the default PATH: the
# image's system R carries dada2 1.40.0 for the downstream analysis, and the
# pipeline must use the pinned 1.38.0. scripts/nextflow.config prepends
# this env to PATH for task execution only.
# GNU wget is in this env deliberately. Fedora 44 ships wget2 as /usr/bin/wget,
# and the database-download processes call `wget -O`. Two problems follow:
# wget2 hangs outright on zenodo.org, and neither wget implements Happy
# Eyeballs — they try zenodo's IPv6 address first and block until timeout,
# because the rootless container has no working IPv6 route. (curl succeeds and
# masks this, since it races v4/v6 and falls back in milliseconds.)
#
# The symptom is a silent stall: wget retries forever, no bytes are written,
# and the pipeline looks hung rather than failed.
RUN set -eux; \
    micromamba create -y -n hifi16s -c conda-forge -c bioconda \
        "bioconductor-dada2=${DADA2_PIPELINE_VERSION}" \
        "cutadapt=${CUTADAPT_VERSION}" \
        "seqkit=${SEQKIT_VERSION}" \
        "csvtk=${CSVTK_VERSION}" \
        r-data.table \
        vsearch \
        wget \
        jq \
        gawk; \
    micromamba clean -a -y; \
    /opt/conda/envs/hifi16s/bin/cutadapt --version; \
    /opt/conda/envs/hifi16s/bin/seqkit version; \
    /opt/conda/envs/hifi16s/bin/csvtk version; \
    /opt/conda/envs/hifi16s/bin/Rscript -e \
        "cat('pipeline dada2', as.character(packageVersion('dada2')), '\n')"

ENV HIFI16S_ENV=/opt/conda/envs/hifi16s

# Force IPv4 for wget rather than patching upstream's `wget -O` calls. This
# path is GNU wget's compiled-in *system* config for this env, so it applies
# without a WGETRC environment variable — setting one as well makes wget warn
# "Both system and user wgetrc point to ..." on every single invocation.
RUN printf 'prefer-family = IPv4\ntries = 5\ntimeout = 60\n' \
        > /opt/conda/envs/hifi16s/etc/wgetrc \
    && /opt/conda/envs/hifi16s/bin/wget --version | head -1

# Nextflow. Java is already present in the base; NXF_HOME is put under /opt so
# a bind-mounted /work is not polluted with plugin caches.
ARG NEXTFLOW_VERSION=26.04.6
ENV NXF_HOME=/opt/nextflow
RUN set -eux; \
    mkdir -p /opt/nextflow; \
    curl -s https://get.nextflow.io | NXF_VER=${NEXTFLOW_VERSION} bash; \
    mv nextflow /usr/local/bin/nextflow; \
    chmod +x /usr/local/bin/nextflow; \
    nextflow -version

# Nextflow bundles plugins on first run; do it at build time so the pipeline
# can run without reaching the network.
RUN nextflow -version && chmod -R a+rwX /opt/nextflow

# renv resolves the project library by copying from this system library
# (renv::hydrate) rather than recompiling everything — see scripts/00_setup.R.
ENV RENV_PATHS_CACHE=/work/.renv-cache

WORKDIR /work
CMD ["/bin/bash"]
