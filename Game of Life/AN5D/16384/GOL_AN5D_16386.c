#include <stdio.h>
#include <stdlib.h>

#define TIMESTEP 1024
#define BENCH_RAD 1
#define SIDE 16386
#define CELL_TYPE char

#define DEAD   0
#define ALIVE  1
#define CELL_NEIGHBOURS 8
#define SRAND_VALUE 1985

void main_computation (CELL_TYPE (*grid)[SIDE][SIDE])
{
	CELL_TYPE rule_table[2][CELL_NEIGHBOURS+1] = {
		    {DEAD,DEAD,DEAD,ALIVE, DEAD,DEAD,DEAD,DEAD,DEAD}, // DEAD is current state
		    {DEAD,DEAD,ALIVE,ALIVE, DEAD,DEAD,DEAD,DEAD,DEAD}  // ALIVE is current state
	};

#pragma scop
    for (int t = 0; t < TIMESTEP; t++)
      for (int i = 1; i < SIDE - 1; i++)
        for (int j = 1; j < SIDE - 1; j++) {
           grid[(t+1)%2][i][j] =
           rule_table [grid[t%2][i][j]] [ (grid[t%2][i-1][j] + grid[t%2][i+1][j] + grid[t%2][i][j-1] + grid[t%2][i][j+1] + 
		   grid[t%2][i-1][j-1] + grid[t%2][i-1][j+1] + grid[t%2][i+1][j-1] + grid[t%2][i+1][j+1]) ];
        }
#pragma endscop
}

/**/

void init_grid (CELL_TYPE (*grid)[SIDE][SIDE]) {
    srand(SRAND_VALUE);
    for(int i = 1; i<SIDE-1; i++) {
        for(int j = 1; j<SIDE-1; j++) {
            grid[0][i][j] = (CELL_TYPE) (rand() % 2);
        }
    }

}



long int print_total_alive (int ca, CELL_TYPE (*grid)[SIDE][SIDE]) {
	int i,j;
	long int total = 0;
	for (i = 1; i<SIDE-1; i++) {
		for (j = 1; j<SIDE-1; j++) {
			total += grid[ca][i][j];
		}
	}
	//printf("Total Alive: %d\n", total);
	return total;
}



int main(int argc, char* argv[])
{
    // Define variables
    int i,j;
    long int total = 0;
    CELL_TYPE (*grid)[SIDE][SIDE];
 

    // Allocate grid
    grid = (CELL_TYPE (*)[SIDE][SIDE]) malloc (2*sizeof(CELL_TYPE)*(SIDE)*(SIDE));  // grid[0][][] --> current CA state
                                                                                // grid[1][][] --> next CA state
  
    // Assign initial population randomly
    init_grid (grid);

 
    // Main GOL game loop
    main_computation (grid);


    // Sum up alive cells and print results
    printf("Total Alive CA 0: %ld\n", print_total_alive (0, grid));
    //printf("Total Alive CA 1: %ld\n", print_total_alive (1, grid));
 

    // Release memory
    free(grid);

}
/**/


// 256 Result in console: "Total Alive: 3281"
// 512 Result in console: "Total Alive: 11072"
// 1024 Result in console: "Total Alive: 45224"
// 2048 Result in console: "Total Alive: 182485"
// 4096 Result in console: "Total Alive: 724393"
// 8192 Result in console: "Total Alive: 2896683"

