#!/usr/bin/env bash
# Shared helpers for the driver scripts. Source it, do not execute it.
#
# Locating this file is the caller's job and is not as simple as it looks:
# under sbatch the script runs from a copy in the scheduler's spool directory,
# so BASH_SOURCE does not point into the repository. Each driver therefore
# tries SV_REPO_ROOT, then its own location, then SLURM_SUBMIT_DIR, then $PWD,
# accepting the first that actually contains bin/common.sh.
#
# The drivers used to hardcode absolute paths to one filesystem, which made them
# unusable anywhere else. Everything machine-specific is now an environment
# variable with a portable default, and container images are given as registry
# URIs that are pulled on demand.

# Published, signed images. Pinned to explicit tags, never :latest, so a rebuild
# cannot silently change what a rerun produces.
TRUVARI_IMAGE_DEFAULT='library://blazv/benchmark-sv/truvari_modded:5.4.0-overlaps-numeric'
ANALYSIS_IMAGE_DEFAULT='library://blazv/benchmark-sv/python-r-analysis:py3.11-r4.4.1'

# Print the available container runtime. Apptainer is preferred; singularity is
# accepted because many sites still ship only that name.
container_engine() {
    if command -v apptainer >/dev/null 2>&1; then
        echo apptainer
    elif command -v singularity >/dev/null 2>&1; then
        echo singularity
    else
        echo "ERROR: neither apptainer nor singularity is on PATH" >&2
        return 1
    fi
}

# resolve_image <uri-or-path> <cache-dir>
# Print a usable local .sif path. A value containing "://" is treated as a
# registry URI and pulled into the cache once; anything else is treated as an
# existing local file and returned unchanged.
resolve_image() {
    local image=${1:?resolve_image needs an image}
    local cache=${2:?resolve_image needs a cache directory}
    local engine name target

    if [[ "$image" != *"://"* ]]; then
        if [[ ! -e "$image" ]]; then
            echo "ERROR: container image not found: $image" >&2
            return 1
        fi
        printf '%s\n' "$image"
        return 0
    fi

    engine=$(container_engine) || return 1
    # library://user/collection/name:tag -> name-tag.sif
    name=$(printf '%s' "$image" | sed 's|.*/||; s|:|-|g')
    target="$cache/$name.sif"

    mkdir -p "$cache"
    if [[ ! -f "$target" ]]; then
        echo "Pulling $image -> $target" >&2
        "$engine" pull "$target" "$image" >&2 || return 1
    fi
    printf '%s\n' "$target"
}

# Print the sha256 of an image, or a marker when it is not a local file, so the
# run manifests stay meaningful without assuming a path exists.
image_checksum() {
    local image=${1:?}
    if [[ -f "$image" ]]; then
        sha256sum "$image"
    else
        echo "unresolved  $image"
    fi
}
