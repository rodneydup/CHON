#!/bin/bash
result=${PWD##*/}
if [ $result == "deployment" ]; then
  cd ..
fi

mkdir -p "deployment/Linux/"

if [ $# -eq 0 ]; then
  echo "Error: No version number provided. Version number required in format <MajorVersion>.<MinorVersion>"
  exit 1
fi

VERSION="$1""-1"
RELEASENAME="CHON-""$VERSION""-amd64"
BUILDLOCATION=$(cd deployment && pwd)

echo "Packaging $RELEASENAME..."

# make directory structure

BUILDDIR="deployment/Linux/$RELEASENAME"

mkdir -p "$BUILDDIR/DEBIAN"
mkdir -p "$BUILDDIR/usr/bin"
mkdir -p "$BUILDDIR/usr/share/applications/"
mkdir -p "$BUILDDIR/usr/share/doc/CHON/"
mkdir -p "$BUILDDIR/usr/share/pixmaps/"

# copy necessary files over

cd "$Dir"

objcopy --strip-debug --strip-unneeded bin/CHON "$BUILDDIR/usr/bin/CHON"

cp "deployment/icons/CHON.png" "$BUILDDIR/usr/share/pixmaps/CHON.png"

# Make .desktop file
echo "[Desktop Entry]
Name=CHON
Comment=Launch CHON
Exec=CHON
Icon=/usr/share/pixmaps/CHON.png
Terminal=false
Type=Application
Categories=Audio;Music;Science;
Name[en_US]=CHON" >>"$BUILDDIR/usr/share/applications/CHON.desktop"

# make Debian control file
echo "Package: CHON
Architecture: amd64
Section: sound
Priority: optional
Version:$VERSION
Maintainer:Rodney DuPlessis <rodney@rodneyduplessis.com>
Depends:libsndfile1, libc6
Homepage: https://github.com/rodneydup/CHON
Description: This package provides CHON, an application for generating sonic gestures and control signals based on the physical simulation of a coupled harmonic oscillator network." >>"$BUILDDIR/DEBIAN/control"

# make copyright file
echo "Format: https://www.debian.org/doc/packaging-manuals/copyright-format/1.0/
Upstream-Name: CHON
Source: https://github.com/rodneydup/CHON

Files: *
Copyright: 2020 Rodney DuPlessis <rodney@rodneyduplessis.com>
License: GPL-3+
 This program is free software; you can redistribute it
 and/or modify it under the terms of the GNU General Public
 License as published by the Free Software Foundation; either
 version 3 of the License, or (at your option) any later
 version.
 .
 This program is distributed in the hope that it will be
 useful, but WITHOUT ANY WARRANTY; without even the implied
 warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the GNU General Public License for more
 details.
 .
 You should have received a copy of the GNU General Public
 License along with this package; if not, write to the Free
 Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 Boston, MA  02110-1301 USA
 .
 On Debian systems, the full text of the GNU General Public
 License version 3 can be found in the file
 '/usr/share/common-licenses/GPL-3'" >>"$BUILDDIR/usr/share/doc/CHON/copyright"
DATE="$(date +'%a, %d %b %Y %H:%M:%S %Z')"
echo "CHON ($VERSION) stable; urgency=high
  * Initial Release
 -- Rodney DuPlessis <rodney@rodneyduplessis.com>  $DATE" >>"$BUILDDIR/usr/share/doc/CHON/changelog.Debian"
gzip -9 -n "$BUILDDIR/usr/share/doc/CHON/changelog.Debian"
echo "Packaging .deb at $BUILDLOCATION..."

# package .deb
cd deployment/Linux
fakeroot dpkg -b "$RELEASENAME" "$RELEASENAME.deb"

echo "Packaging Complete!"
