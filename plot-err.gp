set xlabel "Iterazione"
set ylabel "Errore assoluto"
i = 50 #numero di iterazioni
#set logscale y
set key off
set style line 1 lt 1 lw 0.5 pt 3 pointsize 1
plot for [k=2:i] "errors.txt" using 1:(column(k)) with linespoints
pause -1
