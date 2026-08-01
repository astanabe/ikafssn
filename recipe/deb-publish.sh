#!/usr/bin/env bash
# APT channel publisher for deb.ikafssn.org.
#
# Regenerates the channel from a single release: drops the existing pool/
# and dists/ trees (preserving CNAME, README.md, and the public keyring),
# downloads the six .deb release assets, then rebuilds the Packages
# indices and a clearsigned InRelease + detached Release.gpg per suite.
#
# Required environment:
#   IKAFSSN_VERSION   release version (no leading "v")
#   GH_TOKEN          token with read access to the release assets
#   GITHUB_REPOSITORY owner/repo of the upstream project ("astanabe/ikafssn")
#   IKAFSSN_GPG_KEYID long key id of the imported signing key
#   GPG_PASSPHRASE    passphrase for the signing key (may be empty)
#
# Usage: recipe/deb-publish.sh <channel-repo-dir>

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

cd "${CHANNEL_DIR}"

# Wipe pool/ and dists/ so the channel only ever holds the latest release.
rm -rf pool dists

pool_dir="pool/main/i/ikafssn"
mkdir -p "${pool_dir}"

declare -A suite_for_ubuntu=(
  ["22.04"]="jammy"
  ["24.04"]="noble"
  ["26.04"]="resolute"
)

# One .deb per directory: dpkg-scanpackages writes Filename relative to the
# working directory when its argument is a directory, but absolutises it when
# the argument is a file.  APT requires Filename to be relative to the channel
# root, and each suite must advertise only its own build.
for ubuntu_ver in 22.04 24.04 26.04; do
  for arch in amd64 arm64; do
    asset="ikafssn_${IKAFSSN_VERSION}_ubuntu-${ubuntu_ver}_${arch}.deb"
    deb_dir="${pool_dir}/ubuntu-${ubuntu_ver}/${arch}"
    mkdir -p "${deb_dir}"
    gh release download "${tag}" \
      --repo "${GITHUB_REPOSITORY}" \
      --pattern "${asset}" \
      --dir "${deb_dir}"
    if [ ! -f "${deb_dir}/${asset}" ]; then
      echo "::error::Failed to download ${asset}" >&2
      exit 1
    fi
  done
done

for ubuntu_ver in 22.04 24.04 26.04; do
  suite=${suite_for_ubuntu[${ubuntu_ver}]}
  suite_dir="dists/${suite}"
  for arch in amd64 arm64; do
    bin_dir="${suite_dir}/main/binary-${arch}"
    mkdir -p "${bin_dir}"
    dpkg-scanpackages -m --multiversion \
      "${pool_dir}/ubuntu-${ubuntu_ver}/${arch}" > "${bin_dir}/Packages"
    if grep -q '^Filename: /' "${bin_dir}/Packages"; then
      echo "::error::${bin_dir}/Packages carries an absolute Filename; APT needs it relative to the channel root" >&2
      exit 1
    fi
    gzip -kf9 "${bin_dir}/Packages"
  done

  apt-ftparchive \
    -o "APT::FTPArchive::Release::Origin=ikafssn" \
    -o "APT::FTPArchive::Release::Label=ikafssn" \
    -o "APT::FTPArchive::Release::Suite=${suite}" \
    -o "APT::FTPArchive::Release::Codename=${suite}" \
    -o "APT::FTPArchive::Release::Components=main" \
    -o "APT::FTPArchive::Release::Architectures=amd64 arm64" \
    -o "APT::FTPArchive::Release::Description=ikafssn APT channel (${suite})" \
    release "${suite_dir}" > "${suite_dir}/Release"

  rm -f "${suite_dir}/InRelease" "${suite_dir}/Release.gpg"
  gpg --batch --yes --pinentry-mode loopback \
      --passphrase "${GPG_PASSPHRASE}" \
      --default-key "${IKAFSSN_GPG_KEYID}" \
      --clearsign \
      -o "${suite_dir}/InRelease" \
      "${suite_dir}/Release"
  gpg --batch --yes --pinentry-mode loopback \
      --passphrase "${GPG_PASSPHRASE}" \
      --default-key "${IKAFSSN_GPG_KEYID}" \
      --detach-sign --armor \
      -o "${suite_dir}/Release.gpg" \
      "${suite_dir}/Release"
done
