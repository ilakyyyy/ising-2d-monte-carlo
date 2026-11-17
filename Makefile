CC = gcc
CFLAGS = -Wall -O3
LDFLAGS = -lm

ifeq ($(OS),Windows_NT)
	EXE = .exe
else
	EXE =
endif

COMMON_CORE = io.c memory.c random.c utils.c
SIM_CORE = $(COMMON_CORE) analysis.c metropolis.c observables.c

.PHONY: all clean debug

all: ising$(EXE) ising_parallel$(EXE)

ising$(EXE): $(SIM_CORE) ising.c hysteresis.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

ising_parallel$(EXE): $(SIM_CORE) ising_parallel_thread.c hysteresis_ilaky.c
	$(CC) $(CFLAGS) -fopenmp $^ $(LDFLAGS) -o $@

debug: analysis$(EXE) metropolis$(EXE) observables$(EXE) utils$(EXE)

analysis$(EXE): $(SIM_CORE)
	$(CC) -DDEBUG_ANALYSIS -Wall -g $^ $(LDFLAGS) -o $@

metropolis$(EXE): $(SIM_CORE)
	$(CC) -DDEBUG_METROPOLIS -Wall -g $^ $(LDFLAGS) -o $@

observables$(EXE): memory.c observables.c random.c utils.c
	$(CC) -DDEBUG_OBSERVABLES -Wall -g $^ $(LDFLAGS) -o $@

utils$(EXE): $(COMMON_CORE)
	$(CC) -DDEBUG_UTILS -Wall -g $^ $(LDFLAGS) -o $@

clean:
	rm -f *.o *$(EXE)
