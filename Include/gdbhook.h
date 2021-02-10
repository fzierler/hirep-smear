#ifndef __DEBUG_H_
#define __DEBUG_H_

// To be able to run GDB on MPI jobs, as desribed in
// https://www.open-mpi.org/faq/?category=debugging#serial-debuggers
void gdb_hook();

#endif // __DEBUG_H_
