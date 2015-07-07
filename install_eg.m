fprintf('Compiling mex files...\n');
cd('manoptAuxiliary/pdManifold');
mex sqrtm_triu_real.c
mex sqrtm_triu_complex.c
cd('../../eg');
mex elsd.c
cd('..');
addpath(genpath(pwd))
fprintf('Done.\n');