CC=gcc
#CC_FLAGS= -O3 -Wextra -c -std=c++11 
#CC2_FLAGS= -O3 -Wextra -lm -lcrypto 
CC_FLAGS= -O3 -Wextra -lm 

EXE=main

all: $(EXE)


main: WireWorld_AN5D_256.c
	$(CC) -o WireWorld_AN5D_256  WireWorld_AN5D_256.c $(CC_FLAGS) 

clean:
	rm -rf WireWorld_AN5D_256   *.o

