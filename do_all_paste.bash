#!/bin/bash

################################################################################
# function : do_all_paste
# description :
#   
#
################################################################################
istart=210
iend=222
dir=`pwd`
crystal=../../no3/out_rmsd.dat
gnuplot=/work1/hnishi/local/bin/gnuplot


mkdir $dir/outs
 
for (( i=$istart; i<=$iend; i++ ))
do
   #cd $dir/md$mdn/no$i
   echo residue $i ---------------------------------------
   #sed "s/__FILEDIH__/${i}.dat/g" do_paste_rmsd.bash > zzz.bash
   #chmod +x zzz.bash
   #./a.out ${i}.inp > outs/$i.out
   sed "s/STARTRES/$i/g" plot_dih_rmsd.plt > plot_${i}_rmsd.plt
   cd $dir/outs
   bash ../do_paste_rmsd.bash $dir/outs/$i.dat $crystal ${i}_rmsd.dat
   echo "-148.335 130.968 0" >> $dir/outs/${i}_rmsd.dat
   #gnuplot ../plot_${i}_rmsd.plt
   $gnuplot ../plot_${i}_rmsd.plt
   cd $dir
done



