RELEASE    = v0.3.2
FSEDIR     = ../lib

ifeq ($(shell $(CC) -v 2>&1 | grep -c "gcc version "), 1)
ALIGN_LOOP = # -falign-functions=16 # -falign-loops=32
else
ALIGN_LOOP =
endif

DESTDIR   ?=
CPPFLAGS  += -I$(FSEDIR)
CFLAGS    ?= -O3 $(ALIGN_LOOP)
CFLAGS    += -ansi -std=c99 -fno-exceptions -march=native -mtune=native\
	     -mavx512f -mavx512bw -mavx512dq -mavx512cd \
	     -Wall -Werror -Wno-long-long -Wno-unsed-function -Wno-unsed-parameter -Wshadow
FLAGS      = $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(MOREFLAGS)
FSETEST   ?=
FSEU16TEST?= $(FSETEST)

# Define *.exe as extension for Windows systems
ifneq (,$(filter Windows%,$(OS)))
EXT =.exe
else
EXT =
endif

.PHONY: fullbench fullbench32 fullbenchx32 clean

fullbench32: CFLAGS += -m32
fullbenchx32: CFLAGS += -mx32
fullbench fullbench32 fullbenchx32: benchmark.c entropy_common.c fse_compress.c fse_decompress.c hist.c huf_compress.c huf_decompress.c xxhash.c 
	$(CC) $(FLAGS) $^ -o $@$(EXT)

clean:
	rm -rf *.o
