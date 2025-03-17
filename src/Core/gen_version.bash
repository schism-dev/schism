#!/bin/bash

set -e

get_git_revision_hash() {
    git rev-parse --short HEAD 2>/dev/null || echo "none"
}

gen_version_git() {
    version_template="      character(LEN=32),parameter :: schism_version = '@{VERSION_SCHISM}', git_version = '@{VERSION_GIT}' "
    version_path="$(dirname "$0")"
    versionscratch_path="$version_path/_version"

    describe_re='^(v?[0-9]+\.[0-9]+\.[0-9]+)(-[0-9]+)?(-g[a-z0-9]+)?$'

    git describe --always --dirty > "$versionscratch_path" 2>/dev/null || echo "unknown" > "$versionscratch_path"
    version_raw=$(head -n 1 "$versionscratch_path" | tr -d '\n')
    rm -f "$versionscratch_path"

    is_dirty=false
    [[ "$version_raw" == *-dirty ]] && is_dirty=true && version_raw=${version_raw%-dirty}

    if [[ $version_raw =~ $describe_re ]]; then
        schism_version="${BASH_REMATCH[1]}"
        git_version="${BASH_REMATCH[3]#-g}"
        if [[ -n "${BASH_REMATCH[2]}" || "$is_dirty" == true ]]; then
            nmod="${BASH_REMATCH[2]#-}"
            schism_version+="mod"
            git_version+=" ($nmod commits since semantic tag, edits=$is_dirty)"
        fi
    else
        schism_version="semantic version not determined"
        git_version=$(get_git_revision_hash)
    fi
    echo "$git_version" "$schism_version"
}

gen_version_user() {
    schism_user_version_file="schism_version_user.txt"
    user_version_path="$(dirname "$0")/$schism_user_version_file"
    if [[ -f "$user_version_path" ]]; then
        user_version=$(head -n 1 "$user_version_path" | tr -d '\n')
        [[ ${#user_version} -ge 3 ]] || user_version="develop"
    else
        user_version="develop"
    fi
    echo "$user_version"
}

gen_version() {
    version_template="      character(LEN=32),parameter :: schism_version = '@{VERSION_SCHISM}', git_version = '@{VERSION_GIT}' "
    version_path="$(dirname "$0")"
    template_path="$version_path/schism_version.F90.template"
    versionfile_path=${1:-$version_path/schism_version.F90}

    read git_version schism_version < <(gen_version_git)
    [[ "$schism_version" == "semantic version not determined" ]] && schism_version=$(gen_version_user)
    [[ -z "$git_version" || "$git_version" == "none" ]] && git_version=$(get_git_revision_hash)

    echo " SCHISM version:  $schism_version"
    echo " GIT commit       $git_version"
    
    sed -e "s/@{VERSION_GIT}/$git_version/" -e "s/@{VERSION_SCHISM}/$schism_version/" "$template_path" > "$versionfile_path"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    gen_version "$1"
fi

