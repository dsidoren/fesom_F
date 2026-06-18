#!/usr/bin/env bash
# Source this to load the FESOM3 build environment for the current host.
#   source env.sh <compiler>      # compiler in {intel, gnu}
# Sets FC/CC/CXX (mpi wrappers), FESOM_PLATFORM_STRATEGY, and runtime env.
#
# Anchor toolchain = intel (matches FESOM2 v2.7.3). gnu is the portability build.

_compiler="${1:-intel}"

# Ensure `module` is defined (it is a shell function, not inherited by child
# shells — so a `./configure.sh` subprocess would otherwise lack it).
if ! type module >/dev/null 2>&1; then
    for _init in /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash \
                 /etc/profile.d/modules.sh; do
        # shellcheck disable=SC1090
        [ -f "$_init" ] && source "$_init" && break
    done
fi

if [ -n "$BASH_VERSION" ]; then
    _SOURCE="${BASH_SOURCE[0]}"
elif [ -n "$ZSH_VERSION" ]; then
    _SOURCE=${(%):-%N}
fi
_DIR="$( cd "$( dirname "${_SOURCE}" )" && pwd )"

_host="$(hostname -f)"
if [[ "$_host" =~ ^levante ]] || [[ "$_host" =~ \.lvt\.dkrz\.de$ ]]; then
    export FESOM_PLATFORM_STRATEGY="levante.dkrz.de"
else
    echo "env.sh: unknown host '$_host'; only levante.dkrz.de is configured." >&2
    return 1 2>/dev/null || exit 1
fi

_shellfile="${_DIR}/env/${FESOM_PLATFORM_STRATEGY}/shell.${_compiler}"
if [[ ! -e "${_shellfile}" ]]; then
    echo "env.sh: no shell file for compiler '${_compiler}': ${_shellfile}" >&2
    return 1 2>/dev/null || exit 1
fi
echo "env.sh: sourcing ${_shellfile} (platform=${FESOM_PLATFORM_STRATEGY})"
source "${_shellfile}"
