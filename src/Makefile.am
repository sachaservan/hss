BENCHMARKS = rms_bench ddlog_bench exp_bench group_bench

bin_PROGRAMS = $(BENCHMARKS) sanity_check
check_PROGRAMS = group_test ddlog_test

DDLOG = ddlog.c ddlog.h
ELGAMAL = elgamal.c elgamal.h
ENTROPY = entropy.c entropy.h
FBASE = fbase.c fbase.h
GROUP = group.c group.h
HSS = hss.c rms.c hss.h
TIMEIT = timeit.h
TESTS = $(check_PROGRAMS)

group_test_SOURCES = group_test.c \
	$(GROUP) $(ENTROPY)
ddlog_test_SOURCES = ddlog_test.c \
	$(GROUP) $(ENTROPY) $(DDLOG)

rms_bench_SOURCES = rms_bench.c \
	$(DDLOG) $(ELGAMAL) $(ENTROPY) $(GROUP) $(FBASE) $(HSS) $(TIMEIT)
ddlog_bench_SOURCES = ddlog_bench.c \
	$(DDLOG) $(ENTROPY) $(GROUP) $(TIMEIT)
exp_bench_SOURCES = exp_bench.c \
	$(ENTROPY) $(GROUP) $(FBASE) $(TIMEIT)
group_bench_SOURCES = group_bench.c \
	$(ENTROPY) $(GROUP)
group_bench_LDADD = -lcrypto -lssl

sanity_check_SOURCES = sanity_check.c \
	$(DDLOG) $(ELGAMAL) $(ENTROPY) $(GROUP) $(FBASE) $(HSS) $(TIMEIT)
