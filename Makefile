CC = gcc
LD = gcc
CFLAGS = -Wall -Wextra -O2
LDFLAGS = -lpthread -lm
SRCS := $(wildcard *.c) # wildcard
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.dep)
EXEC = $(SRCS:.c=)
RM = rm -f


all: twin pred kmer_filter

pred.o: src/cmd_args.h src/fasta.h src/kmer.h src/pred.h

pred: pred.o
	$(LD) $(LDFLAGS) -o $@ $^

twin.o: src/cmd_args.h src/fasta.h src/hic.h src/kmer.h src/l2boost.h

twin: twin.o
	$(LD) $(LDFLAGS) -o $@ $^

kmer_filter.o: src/cmd_args.h src/fasta.h src/hic.h src/kmer.h src/l2boost.h

kmer_filter: kmer_filter.o
	$(LD) $(LDFLAGS) -o $@ $^

main: main.o
	$(LD) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(OBJS) $(EXEC) *~

.PHONY:
	all clean
