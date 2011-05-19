python2.6 compileRedHat.py build
pushd build
rm -rvf ./SupportMixRedHat
mv exe.linux-x86_64-2.6 SupportMixRedHat
tar -cvzf ../SupportMixRedHat.tar.gz ./SupportMixRedHat
popd
