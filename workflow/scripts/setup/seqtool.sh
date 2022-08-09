#!/bin/bash

set -xeuo pipefail


if [ $# -eq 0 ]; then
  echo "Usage: $0 <outdir>" >&2
  exit 1
fi

outdir="$1"

tag=v0.3.0
prefix=https://github.com/markschl/seqtool/releases/download/$tag/seqtool-$tag

u="$(uname -s)"
case "${u}" in
    Linux*)     url=$prefix-$(arch)-unknown-linux-gnu.tar.gz;;
    Darwin*)    url=$prefix-$(arch)-apple-darwin.tar.gz;;
    *)          echo "Linux or Mac expected" >&2 && exit 1
esac

tmp=_seqtool_tmp
f=$tmp/seqtool.tar.gz
mkdir -p $tmp

wget -O $f "$url"
tar -xzf $f -C $tmp
mv $tmp/st "$outdir"

rm -R $tmp
