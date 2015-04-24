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
crystal=/work1/hnishi/3rdTrial/mcmd/md/md0_2/no3/out_rmsd.dat
gnuplot=/work1/hnishi/local/bin/gnuplot

mkdir $dir/outs
 
for (( i=$istart; i<$iend; i++ ))
do
   for (( j=$i+1; j<=$iend; j++ ))
   do
      #cd $dir/md$mdn/no$i
      echo residue $i ---------------------------------------
      sed "s/STARTRES/${i}_${j}/g" plot_dih_rmsd_pseudo.plt > plot_${i}_${j}_rmsd.plt
      cd $dir/outs
      bash ../do_paste_rmsd_pseudo.bash $dir/outs/${i}_${j}.dat $crystal ${i}_${j}_rmsd.dat
      #echo "-148.335 130.968 0" >> $dir/outs/${i}_rmsd.dat
      #gnuplot ../plot_${i}_rmsd.plt
      $gnuplot ../plot_${i}_${j}_rmsd.plt
      cd $dir
   done
done

