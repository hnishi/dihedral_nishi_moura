set term png # size 1000, 1000  #canvas size
#set term pbm size 1000, 1000  #canvas size
set output "__STARTRES__.png"
#set size 0.95,0.95  #ratio of graph in canvas

set title "Pseudo-dihedral Plot of Residue 210___STARTRES___222"
set xlabel "Time"
set ylabel "Dihedral Angle"

#set xrange[-180:180]
set yrange[-180:180]
#set xtics 45
set ytics 45
#set grid

set size square
#set size ratio 1 1
#set rmargin 2
#set lmargin 2
#set tmargin 2
#set bmargin 2

#plot "out_dih.dat"
plot "__STARTRES__.dat" w l

#pause -1
