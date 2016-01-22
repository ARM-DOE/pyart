# This is a bit hacky and should most likely be handled by a setup file, but for quick
# compiling and re-compiling this is effective

# Fortran KDP brute force module
rm -f kdp_brute.so
f2py -m kdp_brute -h _kdp_brute.pyf src/kdp_brute.f90 --overwrite-signature
f2py --fcompiler=gfortran --f90flags='-O3 -Wall' -c _kdp_brute.pyf src/kdp_brute.f90