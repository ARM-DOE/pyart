# compress the file
rm KLOT20030101_000921.bz2
bzip2 -z --best -k KLOT20030101_000921

# rename the file
mv KLOT20030101_000921.bz2 example_nexrad_archive_msg1.bz2
