#include <time.h>
#include <sys/time.h>

#define INIT_TIMEIT(flags)                      \
  struct timespec __start, __end;               \
  double __sdiff = 0, __udiff = 0;              \
  int __clock_flags = flags

#define START_TIMEIT()  clock_gettime(__clock_flags, &__start)

#define END_TIMEIT()                                                    \
  clock_gettime(__clock_flags, &__end);                                 \
  __sdiff += (__end.tv_sec - __start.tv_sec);                           \
  __udiff += (__end.tv_nsec - __start.tv_nsec)

#define GET_TIMEIT()                            \
  __sdiff + __udiff * 1e-9

#define TIMEIT_FORMAT "%lf"
