bin_PROGRAMS = rms ddlog_bench
#check_programs = test_ssl1

DDLOG = ddlog.c ddlog.h
ELGAMAL = elgamal.c elgamal.h
ENTROPY = entropy.c entropy.h
GROUP = group.c group.h
HSS = hss.c hss.h
#TESTS = $(check_programs)

#test_ssl1_SOURCES = test_ssl1.c
rms_SOURCES = rms.c $(DDLOG) $(ELGAMAL) $(ENTROPY) $(GROUP) $(HSS)
ddlog_bench_SOURCES = ddlog_bench.c $(DDLOG) $(ENTROPY) $(GROUP)
