set term png # size 1000, 1000  #canvas size
#set term pbm size 1000, 1000  #canvas size
set output "STARTRES_rmsd_pseudo.png"
#set size 0.95,0.95  #ratio of graph in canvas

set title "Pseudo-dihedral Plot of Residue 210_STARTRES_222"
set xlabel "Pseudo-dihedral Angle"
set ylabel "RMSD"

set xrange[-180:180]
#set yrange[-180:180]
#set zrange[0:3]
#set cbrange[0:3]
set xtics 45
#set ytics 45
#set grid

#set size square
#set size ratio 1 1
#set rmargin 2
#set lmargin 2
#set tmargin 2
#set bmargin 2

#plot "out_dih.dat"
plot "STARTRES_rmsd.dat" #u 1:2:3 lc palette

#plot "out_dih_rmsd.dat" u 1:2:3 lc palette

#pause -1
