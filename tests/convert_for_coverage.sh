set -e
for filename in ../coverage_testing/source/*.h
do
    echo ./force_cover $filename -- $@ > $filename.temp
    ./force_cover $filename -- --language c++ -I/usr/include/c++/7 -I/usr/include/clang/7/include -I/usr/local/include -I/usr/bin/../lib/gcc/x86_64-linux-gnu/5.5.0/include -I/usr/include/x86_64-linux-gnu -I/usr/include $@ > $filename.temp
    mv $filename.temp $filename
done