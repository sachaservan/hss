bin_PROGRAMS = rms ddlog_bench exp_bench mul_bench
check_PROGRAMS = group_test

DDLOG = ddlog.c ddlog.h
ELGAMAL = elgamal.c elgamal.h
ENTROPY = entropy.c entropy.h
FBASE = fbase.c fbase.h
GROUP = group.c group.h
HSS = hss.c hss.h
TIMEIT = timeit.h
TESTS = $(check_PROGRAMS)

group_test_SOURCES = group_test.c $(GROUP) $(ENTROPY)

rms_SOURCES = rms.c $(DDLOG) $(ELGAMAL) $(ENTROPY) $(GROUP) $(FBASE) $(HSS) $(TIMEIT)
ddlog_bench_SOURCES = ddlog_bench.c $(DDLOG) $(ENTROPY) $(GROUP) $(TIMEIT)
exp_bench_SOURCES = exp_bench.c $(ENTROPY) $(GROUP) $(FBASE) $(TIMEIT)
mul_bench_SOURCES = mul_bench.c $(ENTROPY) $(GROUP)
