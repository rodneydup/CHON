#!/bin/bash

rm -rf "deployment/OSX/"
mkdir -p "deployment/OSX/"

cp -r bin/CHON.app deployment/OSX
cd deployment
hdiutil create -volname CHON.app -srcfolder ../deployment/OSX -ov CHON.dmg
mv CHON.dmg OSX/
cd OSX/

echo "CHON.dmg created in: $(pwd)"
