#!/bin/bash

# make all charges states from solution equilibrated PDB structure

system='NIST'
charges='26 25 24 23 22 21'
charges=0

for charge in $charges
do
cp orig-noh-equil7-${system}.pdb noh-equil7-${system}.pdb
python ../DistributeCharge/DistributeOnPDB.py -n ${charge} -p equil7-${system}.pdb -s ${system}_residue_sasa.dat
./charge.sh
mv noh-equil7-${system}.pdb noh-equil7-${system}-${charge}.pdb
cd ../
sed "s/XXX/${charge}/g" < leaprc-${system}-gas > tmprc
tleap -f tmprc > tmp.out
value=`grep "charge of the unit" tmp.out | awk '{print $8}'`
echo "-----------LEAP-----------"
echo "$system real: $value target: $charge"
echo "--------------------------"
cd charge-all-pdbs/
done





