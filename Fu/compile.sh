clear

 ifort -CU -CB -o gensps%tbme_nuclearmatter gensps%tbme_nuclearmatter.f
 ifort -CU -CB -mkl -o CCone CCone.f

# ifort -ipo -O3 -parallel -o gensps%tbme_simplepairing gensps%tbme_simplepairing.f
# ifort -ipo -O3 -parallel -mkl -o FCIone FCIone.f
# ifort -ipo -O3 -parallel -mkl -o CCone CCone.f

rm *.mod
