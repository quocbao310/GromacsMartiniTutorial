mkdir NoENWork

cd NoENWork

cp ../Leptin_complex.pdb .

martinize2 -f Leptin_complex.pdb -o system.top -x Leptin_complex_cg.pdb -dssp mkdssp -p backbone -ff martini3001

gmx_mpi editconf -f Leptin_complex_cg.pdb -d 1.5 -bt dodecahedron -o Leptin_complex_box.gro

# Update system.top: #include "../mdp/martini_v3.0.0.itp"

sed -i 's|martini.itp|../mdp/martini_v3.0.0.itp|g' system.top

gmx_mpi grompp -p system.top -f ../mdp/minimization.mdp -c Leptin_complex_box.gro -o minimization-vac.tpr -r Leptin_complex_box.gro

gmx_mpi mdrun -s minimization-vac.tpr -deffnm minimization-vac -ntomp 2 -v

cp ../mdp/water.gro .

# Update system.top: #include "../mdp/martini_v3.0.0_solvents_v1.itp" and "../mdp/martini_v3.0.0_ions_v1.itp"
sed -i '/#include "molecule_0.itp"/i\#include "../mdp/martini_v3.0.0_solvents_v1.itp"' system.top

sed -i '/#include "molecule_0.itp"/i\#include "../mdp/martini_v3.0.0_ions_v1.itp"' system.top

gmx_mpi solvate -cp minimization-vac.gro -cs water.gro -radius 0.21 -o solvated.gro -p system.top

# include water number in system.top ex: W 7999

gmx_mpi grompp -f ../mdp/ions.mdp -c solvated.gro -r solvated.gro -p system.top -o ions.tpr

gmx_mpi genion -s ions.tpr -o Leptin_complex_solv_ions.gro -p system.top -pname NA -nname CL -neutral<<EOF
13
EOF

gmx_mpi grompp -p system.top -c Leptin_complex_solv_ions.gro -f ../mdp/minimization_W.mdp -o minimization.tpr -r Leptin_complex_solv_ions.gro

gmx_mpi mdrun -s minimization.tpr -deffnm minimization -ntomp 2 -v

gmx_mpi grompp -f ../mdp/npt_2fs.mdp -c minimization.gro -p system.top -o npt_2fs.tpr -r minimization.gro

gmx_mpi mdrun -s npt_2fs.tpr -deffnm npt_2fs -ntomp 8 -v

gmx_mpi grompp -f ../mdp/npt_10fs.mdp -c npt_2fs.gro -r npt_2fs.gro -p system.top -o npt_10fs.tpr

gmx_mpi mdrun -s npt_10fs.tpr -deffnm npt_10fs -ntomp 8 -v

gmx_mpi grompp -f ../mdp/npt_20fs.mdp -c npt_10fs.gro -r npt_10fs.gro -p system.top -o npt_20fs.tpr

gmx_mpi mdrun -s npt_20fs.tpr -deffnm npt_20fs -ntomp 8 -v

gmx_mpi grompp -f ../mdp/prod.mdp -c npt_20fs.gro -r npt_20fs.gro -p system.top -o prod.tpr

gmx_mpi mdrun -s prod.tpr -deffnm prod -ntomp 8 -v

# Fix pbc 
gmx_mpi trjconv -s prod.tpr -f prod.xtc -o prod_noPBC.xtc -pbc nojump -center -ur compact<<EOF
1
1
EOF

gmx_mpi trjconv -f prod_noPBC.xtc -s prod.tpr -dump 0 -o frame0.pdb<<EOF
1
EOF


























