
#CC = gcc
CC = gcc

CFLAG = -g -Wall
		
h_addres :=header
HDRS := $(wildcard $(h_addres)/*.h)
CSRC := $(wildcard $(h_addres)/*.C)
OBJS := $(patsubst $(h_addres)/%.C,$(h_addres)/%.o,$(CSRC))
TARGET := ans1 ans2 ans3

all : ans1 ans2 ans3

ans1 : $(OBJS) main1.o
	$(CC) $(CFLAG) -o $@ $^ -lm
	gnuplot ans1_plot.p

plot :
	gnuplot observ.p

ans2 : $(OBJS) main2.o
	$(CC) $(CFLAG) -o $@ $^ -lm

ans3 : $(OBJS) main3.o
	$(CC) $(CFLAG) -o $@ $^ -lm

($(h_addres)/%.o) : $(h_addres)/%.C $(HDRS)
	$(CC) $(CFLAG) -c $< -lm

main%.o : main%.C $(HDRS)
	$(CC) $(CFLAG) -c $< -lm

test.o : test.C $(OBJS)
	$(CC) $(CFLGS) -c $< -lm

clean: 
	rm -f $(TARGET) $(OBJS) main1.o main2.o main3.o ans_1a.png obs.png
