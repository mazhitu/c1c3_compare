This folder has the code to compare the result from c1 and c3.

The c1_forssa.jl and c3_daily.jl are the ones that actually compute the c1 and prestack c3.

The dealwithdaily_c[13].jl are the ones to stack the daily file. The order of running the routines are at the end of dealwithdaily_c3.jl. It also includes a routine remove_large_std, which deletes the correlation that has a std larger than the others (by 10 times), usually this won't happen.

The compare_c1c3.jl makes the comparison. The main generates the plots of CC.pdf and sumcorr.pdf for GNW-STOR pair. The plot_allpairs computes the sumcorr for all pairs.
