set xlabel "Iterazione"
set ylabel "Errore assoluto"
i = `cat errors.txt |grep "#iterations:" | awk '{print $2}'` #numero di iterazioni
eig = i
#set logscale y
set key off
set style line 1 lt 1 lw 0.5 pt 3 pointsize 1
#set term png
#set output "errors.png"
plot for [k=2:eig] "errors.txt" using 1:(column(k)) with linespoints
pause -1
