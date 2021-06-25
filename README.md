# LJOptimization

Contains scripts to optimize the LJ parameters

##Content of the repo:

1) halo_code: code to compute gradients for the halogens containing molecules

2) liq_code: code to compute gradients for the liquid phase

3) liquid_boxes: contains all the liquid boxes used in optimization

4) non_halo_code: code to compute gradients for non-halogens atoms containing molecules

5) solvated_boxes: contains all the solvated box (small molecules in a SWM4 box) used to compute FEH

6) allparameters.dat: contains all the experimental info about the molecules used in optimization

7) exp_set.txt: example input used in run (see below)

8) lj-opt.py: main python script to run the optimization (see below)

9) run.sh: example submission script (see below)

10) vdw-to-fit.txt: example input to define bounds for the parameters (see below)

11) SI: folders containing supplementary info about

	a) feh.xlsx: optimized and scaled FEH

	b) init_final_params.xlsx: initial and final parameters from optimization

	c) testing.xlsx: spread sheet containing molecules and computed properties for testing

	d) training.xlsx: spread sheet containing molecules and computed properties for training

##before running you need to:

1) git clone https://github.com/chetanrrk/LJOptimization/

2) the working directory must contain input data (check 3-6) and master script 'lj-opt.py'
define training molecules in the exp_set.txt (see 'exp_set.txt' for example) with molecule id same as defined in the experimental data file (see 'allparameters.dat'

3) define you initial emin and rmin in the 'vdw-tofit.txt' file (see 'vdw-to-fit.txt' for example)

4) place your current force field file in the 'ff_to_opt' dir

5) place your psf files in the 'psfs' directory

6) place your solvent boxes in the 'train_all_pdbs' dir

7) to run: submit run.sh to the compute queue (might need to change the queue specific settings)

