binwidth = 0.01
binstart = 0

# set width of single bins in histogram
set boxwidth 0.5*binwidth
# set fill style of bins
set style fill solid 0.5
# define macro for plotting the histogram
hist1 = 'u (binwidth*(floor(($4-binstart)/binwidth)+0.5)+binstart):(1.0) smooth freq w boxes'
hist2 = 'u (binwidth*(floor(($4-binstart)/binwidth)+1.0)+binstart):(1.0) smooth freq w boxes'


set ylab "Counts"
set xlab "V_{ij}"
set tit 'Bare Coulomb'
plot 'results_dft/barecoulomb.txt' i 0 @hist1 ls 1 tit 'this run'
repl 'reference_dft/barecoulomb.txt' i 0 @hist2 ls 2 tit 'reference' 
pause -1 

set ylab "Counts"
set xlab "W_{ij}"
set tit 'Screened Coulomb'
set tit 'Bare Coulomb'
plot 'results_dft/screencoulomb.txt' i 0 @hist1 ls 1 tit 'this run'
repl 'reference_dft/screencoulomb.txt' i 0 @hist2 ls 2 tit 'reference' 
pause -1 


