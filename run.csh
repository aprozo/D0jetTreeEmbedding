#!/bin/csh
echo "cleaning that directory"
rm *dataset
rm *session.xml
rm -r LocalLibraries.*
rm temp_gccflags.c
rm test*.root
rm lists/HIOverlay*list

echo "cleaning previous output"
# find /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/production/ -mtime +100 -name \*.log -ls -delete
# find /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/production/ -mtime +100 -name \*.root -ls -delete

# echo "additional out/ report/ csh/"
# find /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/out/ -mtime +100 -name \* -ls -delete
# find /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/report/ -mtime +100 -name \* -ls -delete
# find /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/csh/ -mtime +100 -name \* -ls -delete

rm -r /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/production
rm -r /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/out
rm -r /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/report
rm -r /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/csh
echo "production"
mkdir -p /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/production/
mkdir -p /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/out/
mkdir -p /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/report/
mkdir -p /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/csh/

echo "sending jobs"

star-submit send2014.xml

# wait for jobs to finish using condor_q
# then merge the output files
./wait_for_job.sh

setup 64b
setup root 6.20.08

hadd -f -k -j output_jets_part1.root /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/production/*_[0-4]*_jets.root
hadd -f -k -j output_jets_part2.root /gpfs01/star/pwg/prozorov/output/HeavyIonOverlayForD0Pythia/2014/production/*_[5-9]*_jets.root
hadd -f -k output_jets_full.root output_jets_part*.root
