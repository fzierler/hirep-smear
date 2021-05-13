#include "gdbhook.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/utsname.h>

// To be able to run GDB on MPI jobs, as desribed in
// https://www.open-mpi.org/faq/?category=debugging#serial-debuggers
// Using uname instead of gethostname 
// (gethostname is not compatible with c99)
void gdb_hook() {

  char *debugging_status = getenv("GDBMPI");
  if (0 == strcmp("1", debugging_status)) {

    volatile int i = 0;
    struct utsname unameData;
    uname(&unameData);
    printf("PID %d on %s ready for attach\n", getpid(), unameData.nodename);
    fflush(stdout);
    while (0 == i)
      sleep(5);
  }
}
