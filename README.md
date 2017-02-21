Full sky Power spectrum code (in development)

The code is in python with only the mode coupling matrix computed in fortran.

to compile the fortran bit run:

FFLAGS="-fopenmp -fPIC -Ofast -ffree-line-length-none" f2py -c -m mcm_code mcm_code.f90 wigner3j_sub.f -lgomp

Once you choose your parameters run:

python beam_and_mask.py global.dict

Then change the dictfile for incluing the path to the beam and the mask you created, then

python generate_sim.py global.dict

python mode_coupling.py global.dict

python runMC.py

python checkMC.py




global.dict is the parameter file, it contains:

nside (for healpix map)

beamsize  (arcminute)

mask_fsky (sky fraction in use)

nSims (number of sims)

ncomp (for the moment use 3 for T,Q,U)

beamfile (your beam file: l,bl)

maskfile (your mask)

plot (not in use)

nd (number of splits)

theoryPS (the theoretical power spectrum)

lmax (maximum multipole)

noise (add some white noise to the map)

binsize (size of the bin for your power spectrum)
