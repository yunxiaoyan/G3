 ifort -ipo -O3 -mkl -o FCIone FCIone.f
 ifort -CU -CB -o gensps%tbme_simplepairing gensps%tbme_simplepairing.f

rm *.mod
