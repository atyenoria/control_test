# Makefile for direct_flood

CC = gcc
CPP	= g++

#CFLAGS	= -O
CFLAGS	= -g -DDEBUG

LDFLAGS	=
INCLUDES = -I.
LIBS	= -lm
TARGET	= model_validation_previous

OBJS	= hrandom++.o $(TARGET).o

# 
all:	$(TARGET)

$(TARGET): $(OBJS)
	$(CPP) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	-rm -f $(TARGET) $(OBJS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<
.cpp.o:
	$(CPP) $(CFLAGS) $(INCLUDES) -c $<

hrandom++.o: hrandom++.h
$(TARGET).o:$(TARGET).cpp
