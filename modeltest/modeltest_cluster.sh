#!/bin/sh
cd modeltest
for nex in ../Nexus_files*.nex; do
dirn=${nex/.nex/}
echo "making directory" $dirn
mkdir $dirn
cd $dirn
echo "running modeltest on" $nex
~mpj/bin/mpjrun.sh -dev niodev -wdir $PWD/ -np 10 -jar ~/jModelTest-2.1.4/jModelTest.jar -d $nex -t BIONJ -s 11 -g 4 -f -BIC -a -v
cd ..
echo
done

