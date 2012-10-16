VERSION=rsl-v1.43
rm -rf $VERSION
tar -xvzf $VERSION.tar.gz
cp patch/* $VERSION/
cd $VERSION
autoreconf --force --install
./configure --prefix=/home/sc8/python/pyart/io
make
make install
