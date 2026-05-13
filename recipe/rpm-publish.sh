#!/usr/bin/env bash
# DNF channel publisher for rpm.ikafssn.org.
#
# Regenerates the channel from a single release: drops every elN/ tree
# (preserving CNAME, README.md, and the public keyring), downloads the
# four .rpm release assets, re-signs each .rpm in place with
# rpm --addsign, runs createrepo_c, and writes a detached GPG signature of
# repomd.xml so that dnf can run with repo_gpgcheck=1.
#
# Required environment:
#   IKAFSSN_VERSION   release version (no leading "v")
#   GH_TOKEN          token with read access to the release assets
#   GITHUB_REPOSITORY owner/repo of the upstream project ("astanabe/ikafssn")
#   IKAFSSN_GPG_KEYID long key id of the imported signing key
#   GPG_PASSPHRASE    passphrase for the signing key (may be empty)
#
# Usage: recipe/rpm-publish.sh <channel-repo-dir>

set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: $0 <channel-repo-dir>" >&2
  exit 2
fi

CHANNEL_DIR=$1
: "${IKAFSSN_VERSION:?IKAFSSN_VERSION is required}"
: "${GH_TOKEN:?GH_TOKEN is required}"
: "${GITHUB_REPOSITORY:?GITHUB_REPOSITORY is required}"
: "${IKAFSSN_GPG_KEYID:?IKAFSSN_GPG_KEYID is required}"
GPG_PASSPHRASE=${GPG_PASSPHRASE-}

tag="v${IKAFSSN_VERSION}"

# Configure rpm so --addsign uses our imported key non-interactively.
cat > "${HOME}/.rpmmacros" <<MACRO
%_signature gpg
%_gpg_name ${IKAFSSN_GPG_KEYID}
%__gpg_sign_cmd %{__gpg} gpg --batch --pinentry-mode loopback --passphrase '${GPG_PASSPHRASE}' --no-armor --detach-sign -u "%{_gpg_name}" -o %{__signature_filename} %{__plaintext_filename}
MACRO

cd "${CHANNEL_DIR}"

# Wipe every elN/ tree so the channel only ever holds the latest release.
shopt -s nullglob
for d in el*/; do
  rm -rf "${d}"
done
shopt -u nullglob

for el_ver in 9 10; do
  for arch in x86_64 aarch64; do
    arch_dir="el${el_ver}/${arch}"
    mkdir -p "${arch_dir}"
    asset="ikafssn-${IKAFSSN_VERSION}.el${el_ver}.${arch}.rpm"
    gh release download "${tag}" \
      --repo "${GITHUB_REPOSITORY}" \
      --pattern "${asset}" \
      --dir "${arch_dir}"
    if [ ! -f "${arch_dir}/${asset}" ]; then
      echo "::error::Failed to download ${asset}" >&2
      exit 1
    fi

    rpm --addsign "${arch_dir}/${asset}"

    createrepo_c --general-compress-type zstd "${arch_dir}"

    rm -f "${arch_dir}/repodata/repomd.xml.asc"
    gpg --batch --yes --pinentry-mode loopback \
        --passphrase "${GPG_PASSPHRASE}" \
        --default-key "${IKAFSSN_GPG_KEYID}" \
        --detach-sign --armor \
        -o "${arch_dir}/repodata/repomd.xml.asc" \
        "${arch_dir}/repodata/repomd.xml"
  done
done
