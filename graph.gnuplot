set terminal png
set xrange [0:3]
set yrange [0:3]
set ticslevel 0
# set dgrid3d 50,50
# set hidden3d
# set pm3d
splot "output.dat" u 1:2:3 with lines
