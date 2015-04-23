#!/bin/bash

################################################################################
# function : do_dih
# description :
#   
#
################################################################################
istart=211
iend=221
dir=`pwd`

mkdir $dir/outs
 
for (( i=$istart; i<$iend; i++ ))
do
   for (( j=$i+1; j<=$iend; j++ ))
   do
      #cd $dir/md$mdn/no$i
      echo residue $i ---------------------------------------
      sed "s/__RES2__/${i}/g" prm.inp|sed "s/__RES3__/${j}/g" |sed "s/__OUTFILE__/outs\/${i}_${j}.dat/g" > ${i}_${j}.inp
      ./a.out ${i}_${j}.inp > outs/${i}_${j}.out
      sed "s/__STARTRES__/${i}_${j}/g" plot_pseudo.plt > plot_${i}_${j}.plt
      cd $dir/outs
      gnuplot ../plot_${i}_${j}.plt
      cd $dir
   done
done


