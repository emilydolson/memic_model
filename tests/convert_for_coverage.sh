set -e
for filename in ../coverage_testing/source/*.h
do
    echo ./force_cover $filename -- $@ > $filename.temp
    ./force_cover $filename -- --language c++ $@ > $filename.temp
    mv $filename.temp $filename
done