#!/usr/bin/env bash
# Decides which of the two release channels this run belongs to, and exports
# the tag and version every later step uses. Shared by both jobs so they cannot
# disagree about what is being published.
set -euo pipefail

if [ "${GITHUB_EVENT_NAME}" = "workflow_dispatch" ]; then
  tag="${INPUT_TAG:-}"
  [ -n "$tag" ] || { echo "workflow_dispatch needs a tag input" >&2; exit 1; }
  channel=tag
elif [ "${GITHUB_REF_TYPE}" = "tag" ]; then
  tag="${GITHUB_REF_NAME}"
  channel=tag
else
  # A push to main: one rolling prerelease, always overwritten.
  tag="latest"
  channel=rolling
fi

if [ "$channel" = "tag" ]; then
  version="${tag#v}"
else
  version="main"
fi

{
  echo "CHANNEL=$channel"
  echo "TAG=$tag"
  echo "VERSION=$version"
} >> "$GITHUB_ENV"

echo "channel=$channel tag=$tag version=$version"
