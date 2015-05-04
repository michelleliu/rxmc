#gfortran -O3 -ffixed-line-length-132 *.f
gfortran -O3 -ffixed-line-length-132 -ffpe-trap=zero,invalid *.f
