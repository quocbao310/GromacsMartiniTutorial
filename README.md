# GromacsMartiniTutorial
This is a tutorial for using martini force field in Gromacs. Martini version is Martini3 and gromacs version is 2020.4. The biological system in this tutorial is a complex of Leptin- Leptin binding domain CRH2. 

# To create environment for MDAnalysis
conda env create -f env.yml

#It is important to note that the script in the MDAnalysis directory can only be used for CGMD simulation, and it will not work for ATMD simulation. 
#The script in MDAnalysis directory requires software and python libraries specified in env.yml script.
#The backmapping approach requires python2 environment. The python, bash script and mapping files required for backmapping were downloaded from http://cgmartini.nl/index.php/tutorials-general-introduction-gmx5/others-gmx5. 

