cd build
sudo make package -j8

cd _CPack_Packages/Darwin/PackageMaker/

# sudo /Applications/PackageMaker.app/Contents/MacOS/PackageMaker -build -p "radiance-Darwin.pkg" -f "radiance-5.2.e4d2f765dc-Darwin" -r "Resources" -i "Info.plist" -d "Description.plist" -v
pkgbuild --identifier radiance.pkg --root "radiance-5.2.e4d2f765dc-Darwin/usr/local/radiance" --install-location "/usr/local/radiance" radiance-5.2.e4d2f765dc-Darwin.pkg

zip -r Radiance.zip *.pkg 
cp -r Radiance.zip ../../../Radiance.zip
du -h -d=0 *
set -e