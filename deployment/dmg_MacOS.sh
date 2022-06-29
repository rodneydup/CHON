#!/bin/bash

rm -rf "deployment/CHON_MacOS/"
mkdir -p "deployment/CHON_MacOS/"

cp -r bin/CHON.app deployment/CHON_MacOS
cd deployment
hdiutil create -volname CHON.app -srcfolder ../deployment/CHON_MacOS -ov CHON.dmg
cd CHON_MacOS/

echo "CHON.dmg created in: $(pwd)"
