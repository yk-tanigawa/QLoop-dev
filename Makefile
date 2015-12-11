CC = gcc
LD = gcc
CFLAGS = -Wall -O2
LDFLAGS = #-lpthread -lm
SRCS := $(wildcard *.c) # wildcard
OBJS = $(SRCS:.c=.o)
DEPS = $(SRCS:.c=.dep)
EXEC = $(SRCS:.c=)
RM = rm -f


all: prep mkLibsvmFile

prep: prep.o
	$(LD) $(LDFLAGS) -o $@ $^

mkLibsvmFile: mkLibsvmFile.o
	$(LD) $(LDFLAGS) -o $@ $^


chrom: chrom.o
	$(LD) $(LDFLAGS) -o $@ $^

clean:
	$(RM) $(OBJS) $(EXEC) *~

.PHONY:
	all clean
