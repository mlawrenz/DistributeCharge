#!/bin/bash

# make all charges states from solution equilibrated PDB structure


system='AMG102'
charges='25'

for charge in $charges
do
cp orig-noh-equil5-${system}.pdb noh-equil5-${system}.pdb
python ../DistributeCharge/DistributeOnPDB.py -n ${charge} -p equil5-${system}.pdb -s ${system}_residue_sasa.dat
./charge.sh
mv noh-equil5-${system}.pdb noh-equil5-${system}-${charge}.pdb
cd ../
sed "s/XXX/${charge}/g" < leaprc-${system}-gas > tmprc
tleap -f tmprc > tmp.out
value=`grep "charge of the unit" tmp.out | awk '{print $8}'`
echo "-----------LEAP-----------"
echo "$system real: $value target: $charge"
echo "--------------------------"
cd charge-all-pdbs/
done

