#!/bin/sh

for nex in ~/Dropbox/Crandall_tobo/Nexus_files/*.nex; do
dirn=${nex/.nex/}
echo "making directory" $dirn
mkdir $dirn
cd $dirn
echo "running modeltest on" $nex
java -jar /Applications/jModelTest-2.1.4/jModelTest.jar -d $nex -t BIONJ -s 5 -g 4 -f -BIC -v -tr 4 -o modeltest_output.txt
cd ..
echo
done