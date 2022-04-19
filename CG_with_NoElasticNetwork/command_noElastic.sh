mkdir NoENWork

cd NoENWork

cp ../Leptin_complex.pdb .

martinize2 -f Leptin_complex.pdb -o system.top -x Leptin_complex_cg.pdb -dssp mkdssp -p backbone -ff martini3001




























