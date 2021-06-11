# LJOptimization

Contains scripts to optimize the LJ parameters

before running you need to:

git clone https://github.com/chetanrrk/LJOptimization/
the working directory must contain input data (check 3-6) and master script 'lj-opt.py'
define training molecules in the exp_set.txt (see 'exp_set.txt' for example) with molecule id same as defined in the experimental data file (see 'allparameters.dat'
define you initial emin and rmin in the 'vdw-tofit.txt' file (see 'vdw-to-fit.txt' for example)
place your current force field file in the 'ff_to_opt' dir
place your psf files in the 'psfs' directory
place your solvent boxes in the 'train_all_pdbs' dir
to run: submit run.sh to the compute queue (might need to change the queue specific settings)
