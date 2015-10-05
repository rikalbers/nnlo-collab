FC=gfortran
#
sed -i -e 's/^\s*FC\s*=.*/FC = '$FC'/' makefile-multi
make -f makefile-multi clean || true
$FC -c particles.f90
#
cd ToyAmp
sed -i -e 's/^\s*FC\s*=.*/FC = '$FC'/' makefile
rm libtoyamp.a || true
make -f makefile clean || true
make -f makefile
#
cd ../..
sed -i -e 's/^\s*FC\s*=.*/FC = '$FC'/' makefile_multi
rm libkaleu_multi.a || true
make -f makefile_multi clean || true
make -f makefile_multi
#
cd example
make -f makefile-multi
