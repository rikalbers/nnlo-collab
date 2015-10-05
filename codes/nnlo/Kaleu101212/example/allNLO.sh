FC=gfortran
#
sed -i -e 's/^\s*FC\s*=.*/FC = '$FC'/' makefileNLO
make -f makefileNLO clean || true
$FC -c particles.f90
#
cd ToyAmp
sed -i -e 's/^\s*FC\s*=.*/FC = '$FC'/' makefileNLO
rm libtoyamp_NLO.a || true
make -f makefileNLO clean || true
make -f makefileNLO
#
cd ../..
sed -i -e 's/^\s*FC\s*=.*/FC = '$FC'/' makefile_exampleNLO
rm libkaleu_exampleNLO.a || true
make -f makefile_exampleNLO clean || true
make -f makefile_exampleNLO
#
cd example
make -f makefileNLO
