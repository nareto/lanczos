set title "Convergenza autovalori"
set key off
i = 50 #numero di iterazioni
eig = i #numero di autovalori da tracciare
set xlabel "Iterazione"
set xrange [1:i]
set style line 1 lt 1 lw 0.5 linecolor "black" pt 4 pointsize 1
set style arrow 1 nohead lw 0.1 linecolor rgb "red"
plot for [k=2:eig] "eigv.txt" using 1:(column(k)) with linespoints ls 1,\
     "correcteigv.txt" using (i):1 with points pt 7 ps 0.5 linecolor rgb "red",\
     "" using (0.0):1:(i):(0.0) with vectors arrowstyle 1
pause -1