#include <stdio.h>
#include <stdlib.h>

#define SIZE 1024

#define TIMESTEP 1024

#define CELL_NEIGHBOURS 8
#define SRAND_VALUE 1985
#define BLOCK_SIZE 32

#define DEAD 0
#define ALIVE 1
//#define CELL_TYPE char
typedef unsigned char CELL_TYPE;

 
__global__ void ghostRows(int dim, CELL_TYPE *grid)
{
    // We want id ∈ [1,dim]
    int id = blockDim.x * blockIdx.x + threadIdx.x + 1;
 
    if (id <= dim)
    {
        //Copy first real row to bottom ghost row
        grid[(dim+2)*(dim+1)+id] = grid[(dim+2)+id];
        //Copy last real row to top ghost row
        grid[id] = grid[(dim+2)*dim + id];
    }
}
 
__global__ void ghostCols(int dim, CELL_TYPE *grid)
{
    // We want id ∈ [0,dim+1]
    int id = blockDim.x * blockIdx.x + threadIdx.x;
 
    if (id <= dim+1)
    {
        //Copy first real column to right most ghost column
        grid[id*(dim+2)+dim+1] = grid[id*(dim+2)+1];
        //Copy last real column to left most ghost column 
        grid[id*(dim+2)] = grid[id*(dim+2) + dim];    
    }
}
 
__global__ void GOL(int dim, CELL_TYPE *grid, CELL_TYPE *newGrid, CELL_TYPE *lookup_table_)
{
    // We want id ∈ [1,dim]
    int iy = blockDim.y * blockIdx.y + threadIdx.y + 1;
    int ix = blockDim.x * blockIdx.x + threadIdx.x + 1;
    int id = iy * (dim+2) + ix;
    CELL_TYPE (*lookup_table)[CELL_NEIGHBOURS+1] = (CELL_TYPE (*)[CELL_NEIGHBOURS+1]) lookup_table_;

    if (iy>0 && iy <= dim && ix>0 && ix <= dim) {
           newGrid[id] =
           lookup_table [grid[id]] [  grid[id+(dim+2)] + grid[id-(dim+2)]   
                                    + grid[id+1] + grid[id-1]
                                    + grid[id+(dim+3)] + grid[id-(dim+3)]
                                    + grid[id-(dim+1)] + grid[id+(dim+1)]  ];

    }

    /*int numNeighbors;
 
    if (iy <= dim && ix <= dim) {
 
        // Get the number of neighbors for a given grid point
        numNeighbors = grid[id+(dim+2)] + grid[id-(dim+2)] //upper lower
                     + grid[id+1] + grid[id-1]             //right left
                     + grid[id+(dim+3)] + grid[id-(dim+3)] //diagonals
                     + grid[id-(dim+1)] + grid[id+(dim+1)];
 
        int cell = grid[id];
        // Here we have explicitly all of the game rules
        if (cell == 1 && numNeighbors < 2)
            newGrid[id] = 0;
        else if (cell == 1 && (numNeighbors == 2 || numNeighbors == 3))
            newGrid[id] = 1;
        else if (cell == 1 && numNeighbors > 3)
            newGrid[id] = 0;
        else if (cell == 0 && numNeighbors == 3)
            newGrid[id] = 1;
        else
            newGrid[id] = cell;
    }*/
}


void setup_lookup_table (CELL_TYPE *lookup_table_) {
    int j;
    CELL_TYPE (*lookup_table)[CELL_NEIGHBOURS+1] = (CELL_TYPE (*)[CELL_NEIGHBOURS+1]) lookup_table_;

	/*lookup_table[2][CELL_NEIGHBOURS+1] = {
		    {DEAD,DEAD,DEAD,ALIVE, DEAD,DEAD,DEAD,DEAD,DEAD}, // DEAD is current state
		    {DEAD,DEAD,ALIVE,ALIVE, DEAD,DEAD,DEAD,DEAD,DEAD}  // ALIVE is current state
	};*/

    for (j=0; j<CELL_NEIGHBOURS+1; j++) {
        lookup_table[0][j] = DEAD;
        lookup_table[1][j] = DEAD;
    }
    lookup_table[0][3] = ALIVE;
    lookup_table[1][2] = ALIVE;
    lookup_table[1][3] = ALIVE;
}


 
int main(int argc, char* argv[])
{
    int i,j,iter;
    CELL_TYPE* h_grid; //Grid on host
    CELL_TYPE* d_grid; //Grid on device
    CELL_TYPE* d_newGrid; //Second grid used on device only
    CELL_TYPE* d_tmpGrid; //tmp grid pointer used to switch between grid and newGrid
    CELL_TYPE* h_lookup_table;
    CELL_TYPE* d_lookup_table;
 
    //int dim = 16384; //Linear dimension of our grid - not counting ghost cells
    //int maxIter = 1<<10; //Number of game steps
 
    size_t bytes = sizeof(CELL_TYPE)*(SIZE+2)*(SIZE+2); //2 added for periodic boundary condition ghost cells
    size_t bytes2 = sizeof(CELL_TYPE)*2*(CELL_NEIGHBOURS+1);    // Lookup table size
    // Allocate host Grid used for initial setup and read back from device
    h_grid = (CELL_TYPE *)malloc(bytes);
    h_lookup_table = (CELL_TYPE *)malloc(bytes2);

    // Allocate device grids
    cudaMalloc(&d_grid, bytes);
    cudaMalloc(&d_newGrid, bytes);
    cudaMalloc(&d_lookup_table, bytes2);
 
    // Assign initial population randomly
    srand(SRAND_VALUE);
    for(i = 1; i<=SIZE; i++) {
        for(j = 1; j<=SIZE; j++) {
            h_grid[i*(SIZE+2)+j] = rand() % 2;
        }
    }
 
    // Init look-up table:
    setup_lookup_table (h_lookup_table);

    // Copy over initial game grid (Dim-1 threads)
    cudaMemcpy(d_grid, h_grid, bytes, cudaMemcpyHostToDevice);
    // Copy lookup table
    cudaMemcpy(d_lookup_table, h_lookup_table, bytes2, cudaMemcpyHostToDevice);
 
    dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE,1);
    int linGrid = (int)ceil(SIZE/(float)BLOCK_SIZE);
    dim3 gridSize(linGrid,linGrid,1);
 
    dim3 cpyBlockSize(BLOCK_SIZE,1,1);
    dim3 cpyGridRowsGridSize((int)ceil(SIZE/(float)cpyBlockSize.x),1,1);
    dim3 cpyGridColsGridSize((int)ceil((SIZE+2)/(float)cpyBlockSize.x),1,1);
 
    // Main game loop
    for (iter = 0; iter<TIMESTEP; iter++) {
 
        ghostRows<<<cpyGridRowsGridSize, cpyBlockSize>>>(SIZE, d_grid);
        ghostCols<<<cpyGridColsGridSize, cpyBlockSize>>>(SIZE, d_grid);
        GOL<<<gridSize, blockSize>>>(SIZE, d_grid, d_newGrid, d_lookup_table);
 
        // Swap our grids and iterate again
        d_tmpGrid = d_grid;
        d_grid = d_newGrid;
        d_newGrid = d_tmpGrid;
    }//iter loop
 
    // Copy back results and sum
    cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
 
    // Sum up alive cells and print results
    int total = 0;
    for (i = 1; i<=SIZE; i++) {
        for (j = 1; j<=SIZE; j++) {
            total += h_grid[i*(SIZE+2)+j];
        }
    }
    printf("Total Alive: %d\n", total);
 
    // Release memory
    cudaFree (d_grid);
    cudaFree (d_newGrid);
    cudaFree (d_lookup_table);
    free (h_grid);
    free (h_lookup_table);
 
    return 0;
}

// 256 Result in console: "Total Alive: 3281"
// 512 Result in console: "Total Alive: 11072"
// 1024 Result in console: "Total Alive: 45224"
// 2048 Result in console: "Total Alive: 182485"
// 4096 Result in console: "Total Alive: 724393"
// 8192 Result in console: "Total Alive: 2896683"
// 16384 Result in console: "Total Alive: 11547651"

