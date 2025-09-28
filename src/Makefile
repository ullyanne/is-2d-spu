CPP = g++
CFLAGS = -g -O3
OBJECTS = is.o

all: is

is: $(OBJECTS)
	$(CPP) $(CFLAGS) $(OBJECTS) -o is -lm

clean:
	rm -f is $(OBJECTS)
