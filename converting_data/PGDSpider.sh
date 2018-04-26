#!/bin/sh
for arp in *.arp; do
nex=${arp/.arp/.nex} #make a .nex filename
mig=${arp/.arp/.mig} #make a .mig filename
fas=${arp/.arp/.fasta} #make a .mig filename
echo $arp
echo
echo "converting to nexus"
java -Xmx2048m -Xms512m -jar /Applications/PGDSpider_2.0.5.1/PGDSpider2-cli.jar -inputfile $arp -outputfile ~/Dropbox/Crandall_tobo/Nexus_files/$nex -spid ~/Dropbox/Crandall_tobo/scripts/arp_to_nex.spid 
echo
echo "converting to migrate"
java -Xmx2048m -Xms512m -jar /Applications/PGDSpider_2.0.5.1/PGDSpider2-cli.jar -inputfile $arp -outputfile ~/Dropbox/Crandall_tobo/Migrate_datafiles/$mig -spid ~/Dropbox/Crandall_tobo/scripts/arp_to_mig.spid 
echo
echo "converting to fasta"
java -Xmx2048m -Xms512m -jar /Applications/PGDSpider_2.0.5.1/PGDSpider2-cli.jar -inputfile $arp -outputfile ~/Dropbox/Crandall_tobo/FASTA_datafiles/$fas -spid ~/Dropbox/Crandall_tobo/scripts/arp_to_mig.spid 
echo
echo
echo
done



 