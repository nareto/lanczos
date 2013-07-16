set xlabel "iteration"
plot "eigv.txt" using 1:2 title "primo" pointtype 7 linetype 1, \
     "eigv.txt" using 1:3 title "secondo" pointtype 7 linetype 2,\
     "eigv.txt" using 1:4title "penultimo" pointtype 7 linetype 3,\
     "eigv.txt" using 1:5 title "ultimo" pointtype 7 linetype 4
pause -1
