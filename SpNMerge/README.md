# Which commit is what:

"last-merged"                 : 16fe96e Fixed a bug in the previous submission
master                        : 0af1e75 Another fix for ranluxs
real merge base for master&spn: 9e944b8
spn                           : 3bbf0e7 

(real merge base is newer than "last-merged")


# Notes on the actual git merge (solving conflicts)
## `NodeNumber/mk_eigval.c`:

It is assumed that the last changes by Ed on the `spn` branch are correct.
The conflict was likely due to the fact that the code was reformatted 
for a small change on the master branch 
(only removed `+MPI_PID` in the argument to `rlxd_init` call,
but on the spn branch that call was removed altogether).

This might require some additional review.

## `NodeNumber/mk_modenumber.c`:

Same considerations apply to this file as well.
