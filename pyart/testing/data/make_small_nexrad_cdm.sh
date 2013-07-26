# create the DUMMY file
python dummify_nexrad_file.py

# test the dummy file
nosetests -v check_nexrad_dummy.py

# convert to CDM file using toolsUI
java -classpath toolsUI-4.3.jar ucar.nc2.FileWriter \
     -in KATX20130717_195021_V06_DUMMY -out example_nexrad_cdm
rm KATX20130717_195021_V06_DUMMY

# compress 
rm example_nexrad_cdm.bz2
bzip2 -z --best example_nexrad_cdm
