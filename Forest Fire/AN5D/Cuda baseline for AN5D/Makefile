CC=gcc
#CC_FLAGS= -O3 -Wextra -c -std=c++11 
#CC2_FLAGS= -O3 -Wextra -lm -lcrypto 
CC_FLAGS= -O3 -Wextra -lm 

EXE=main

all: $(EXE)


main: ForestFire_AN5D_256.c
	$(CC) -o ForestFire_AN5D_256  ForestFire_AN5D_256.c $(CC_FLAGS) 

clean:
	rm -rf ForestFire_AN5D_256   *.o

