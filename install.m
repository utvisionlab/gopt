% First build the mex files

% THIS WORKS ON A UNIX LIKE SYSTEM
% FOR WINDOWS, TRY USING INSTALL_EG

fprintf('Compiling mex files...\n');
cd('manoptAuxiliary/pdManifold');
mex -largeArrayDims -lmwblas -lmwlapack -DUNDERSCORE_LAPACK_CALL ...
    -DBLAS64 CFLAGS="-std=c99 -Wno-deprecated-declarations" sqrtm_triu_real.c
mex -largeArrayDims -lmwblas -lmwlapack -DUNDERSCORE_LAPACK_CALL ...
    -DBLAS64 CFLAGS="-std=c99 -Wno-deprecated-declarations" sqrtm_triu_complex.c
cd('../../eg');
mex -largeArrayDims -DBLAS64 CFLAGS="-std=c99 -Wno-deprecated-declarations" elsd.c
cd('..');
addpath(genpath(pwd))
fprintf('DONE\n');
