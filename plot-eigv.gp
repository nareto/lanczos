set key off
i = `cat eigv.txt |grep "#iterations:" | awk '{print $2}'` #numero di iterazioni
eig = i #numero di autovalori da tracciare
set xlabel "Iterazione"
set xrange [1:i]
set style line 1 lt 1 lw 0.5 linecolor rgb "black" pt 4
set style arrow 1 nohead lw 0.1 linecolor rgb "red"
#set term png
#set output "eigv.png"
plot "correcteigv.txt" using (0.0):1:(i):(0.0) with vectors arrowstyle 1,\
     "correcteigv.txt" using (i):1 with points pt 7 ps 0.5 linecolor rgb "red",\
     for [k=2:eig] "eigv.txt" using 1:(column(k)) with linespoints ls 1
pause -1