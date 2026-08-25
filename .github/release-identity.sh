#!/usr/bin/env bash
# Works out what this run publishes, if anything, and exports it for every
# later step. Both jobs source this, so the build and the publish cannot
# disagree about what is being produced.
#
#   push to main            release, version = next patch after the highest v* tag
#   push a v* tag           release, version = that tag
#   workflow_dispatch       artifacts only, no tag and no release
set -euo pipefail

manifest_version() {
  # The first version = line in Cargo.toml is the package's own.
  grep -m1 '^version = ' Cargo.toml | sed 's/.*"\(.*\)".*/\1/'
}

if [ "${GITHUB_EVENT_NAME}" = "workflow_dispatch" ]; then
  # A manual build of whatever branch was selected. Nothing is published, so
  # the version is only a label on the artifacts.
  channel=artifacts
  version="$(manifest_version)"
elif [ "${GITHUB_REF_TYPE}" = "tag" ]; then
  # An explicit tag is how a minor or major bump is made; the patch counter
  # below carries on from it.
  channel=release
  version="${GITHUB_REF_NAME#v}"
else
  channel=release
  latest="$(git tag -l 'v*' --sort=-v:refname | head -1)"
  if [ -z "$latest" ]; then
    # Nothing released yet: start from whatever the manifest declares.
    version="$(manifest_version)"
  else
    base="${latest#v}"
    major="${base%%.*}"
    rest="${base#*.}"
    minor="${rest%%.*}"
    patch="${rest#*.}"
    patch="${patch%%[-+]*}"
    version="${major}.${minor}.$((patch + 1))"
  fi
fi

{
  echo "CHANNEL=$channel"
  echo "VERSION=$version"
  echo "TAG=v$version"
} >> "$GITHUB_ENV"

echo "channel=$channel version=$version tag=v$version"
