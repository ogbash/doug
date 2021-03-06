# ctl.file.assembled
solver 2
solve_maxiters 300
method 1
levels  2
overlap -1 
smoothers 0
input_type 2
symmstruct T
symmnumeric T
# ###################
# aggregate level 1:
radius1 2
strong1 0.67e0
minasize1 2
#maxasize1 19
# aggregate level 2:
radius2 1
strong2 0.67e0
minasize2 2
#maxasize2 96
# ###################
matrix_type 1
number_of_blocks 1
initial_guess 2
start_vec_file ./NOT.DEFINED.start_vec_file
start_vec_type 2
solve_tolerance 1.0e-12
solution_format 2
solution_file solution.xdr
#debug -5
debug 0
verbose 3
plotting 0
assembled_mtx_file Lap16x16.txt
assembled_mtx_format 0
assembled_rhs_format 2
assembled_rhs_file rhs.xdr
