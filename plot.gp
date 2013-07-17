set xlabel "iteration"
set ylabel "absolute error"
set logscale y
set key outside right bottom horizontal
plot "errors.txt" using 1:2 title "primo" , \
     "errors.txt" using 1:3 title "ultimo" ,\
     "errors.txt" using 1:4 title "secondo" ,\
     "errors.txt" using 1:5 title "penultimo" ,\
     "errors.txt" using 1:6 title "terzo" ,\
     "errors.txt" using 1:7 title "terzultimo" ,\
     "errors.txt" using 1:8 title "quarto" ,\
     "errors.txt" using 1:9 title "quartultimo" 
pause -1
