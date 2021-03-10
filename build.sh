set -ex

# if [ "$1" = "--rm" ]; then
#     make clean
# fi

./configure --enable-shared --disable-demos --enable-double-precision --disable-asserts --enable-malloc --prefix=`pwd`/ode-cuda-release

BUILD="make -j$(nproc)"
#make -j$(nproc)
$BUILD
cd tests/UnitTest++/src/Posix && $BUILD libhelper.la && cd -
cd tests/UnitTest++/src && $BUILD libunittestpp.la && cd -
cd tests && $BUILD tests && cd -
#make install
