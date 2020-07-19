
#CC = gcc
CC = g++

CFLAG = -g -O3 -Wall 

HDRS := $(wildcard *.h)

CSRC1 := main1.C func1.C
OBJS1 := $(CSRC1:%.C=%.o)
CSRC2 := main2.C degcore.C factorial.C
OBJS2 := $(CSRC2:%.C=%.o)
CSRC3 := main3.C degcore.C factorial.C
OBJS3 := $(CSRC3:%.C=%.o)

all : one two three

one: $(OBJS1)
	$(CC) $(CFLAG) -o $@ $^

two: $(OBJS2)
	$(CC) $(CFLAG) -o $@ $^

three: $(OBJS3)
	$(CC) $(CFLAG) -o $@ $^

%.o : %.C $(HDRS)
	$(CC) $(CFLAG) -c $<

clean: 
	rm -f one two three $(OBJS1) $(OBJS2) $(OBJS3)
