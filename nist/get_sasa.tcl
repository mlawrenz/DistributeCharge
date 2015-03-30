

mol load pdb nist.amber.pdb
set all [ atomselect top "all" ]
set residues [ lsort -unique [ $all get resid ] ]
set filename "residue_sasa.dat"
set file [open $filename "w"]

foreach res $residues {
    set sel [ atomselect top [ format "resid %s" $res ] ]
    set sasa [ measure sasa 1.4 $all -restrict $sel ]
    puts $file "$res $sasa"
}
close $file
exit
