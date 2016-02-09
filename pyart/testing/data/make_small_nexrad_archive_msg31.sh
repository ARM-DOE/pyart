# create the DUMMY file
python dummify_nexrad_file.py

# test the dummy file
nosetests -v check_nexrad_dummy.py

# compress the file
rm KATX20130717_195021_V06_DUMMY.bz2
bzip2 -z --best KATX20130717_195021_V06_DUMMY

# rename the file
mv KATX20130717_195021_V06_DUMMY.bz2 example_nexrad_archive_msg31.bz2
