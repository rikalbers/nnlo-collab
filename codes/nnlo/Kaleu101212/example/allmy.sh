FC=gfortran
#
sed -i -e 's/^\s*FC\s*=.*/FC = '$FC'/' makefile
make -f makefile clean || true
$FC -c particles.f90
#
cd ToyAmp
sed -i -e 's/^\s*FC\s*=.*/FC = '$FC'/' makefile
rm libtoyamp.a || true
make -f makefile clean || true
make -f makefile
#
cd ../..
sed -i -e 's/^\s*FC\s*=.*/FC = '$FC'/' makefile_example
rm libkaleu_example.a || true
make -f makefile_example clean || true
make -f makefile_example
#
cd example
make -f makefile
