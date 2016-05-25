CC = gcc
LD = gcc
CFLAGS = -Wall -Wextra -O2
LDFLAGS = -lpthread -lm
SRCS := $(wildcard *.c) # wildcard
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.dep)
EXEC = $(SRCS:.c=)
RM = rm -f


all: twin

twin.o: src/cmd_args.h src/fasta.h src/hic.h src/kmer.h src/l2boost.h

twin: twin.o
	$(LD) $(LDFLAGS) -o $@ $^

main: main.o
	$(LD) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(OBJS) $(EXEC) *~

.PHONY:
	all clean
