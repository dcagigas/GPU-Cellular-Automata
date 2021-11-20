#include <stdio.h>
#include <stdlib.h>
 
#define GRID_SIZE 256                       // Grid dimension: SIZExSIZE
#define MAXITER 1024                        // Number of game steps

typedef unsigned int CELL_TYPE;
#define ELEMENTS_PER_CELL 4                 // "int_32" data type per cell
#define ROW_SIZE GRID_SIZE/ELEMENTS_PER_CELL    // Real grid dimension

#define SRAND_VALUE 1985
#define BLOCK_SIZE 32

#define CELL_NEIGHBOURS 8
// Classic GOL:
#define MIN_NOF_NEIGH_FROM_ALIVE_TO_ALIVE 2
#define MAX_NOF_NEIGH_FROM_ALIVE_TO_ALIVE 3
#define MIN_NOF_NEIGH_FROM_DEAD_TO_ALIVE 3
#define MAX_NOF_NEIGH_FROM_DEAD_TO_ALIVE 3

#define ALIVE  1
#define DEAD   0


__global__ void kernel_init_rule_table (int *GPU_rule_table);
__forceinline__ unsigned char getSubCellH (CELL_TYPE cell, char pos);
__forceinline__ void setSubCellH (CELL_TYPE *cell, char pos, unsigned char subcell);
__device__ unsigned char getSubCellD (CELL_TYPE cell, char pos);
__device__ void setSubCellD (CELL_TYPE *cell, char pos, unsigned char subcell);
__global__ void ghostRows(CELL_TYPE *grid);
__global__ void ghostCols(CELL_TYPE *grid);
__global__ void GOL (CELL_TYPE *grid, CELL_TYPE *newGrid, int *GPU_rule_table);

void init_random_values (CELL_TYPE *h_grid);
long int print_total_alive (CELL_TYPE *h_grid);



int main(int argc, char* argv[])
{
    long int total = 0;
    int iter;
    CELL_TYPE* h_grid; //Grid on host
    CELL_TYPE* d_grid; //Grid on device
    CELL_TYPE* d_newGrid; //Second grid used on device only
    CELL_TYPE* d_tmpGrid; //tmp grid pointer used to switch between grid and newGrid
    int *GPU_rule_table;
  
    size_t bytes = sizeof(CELL_TYPE)*(GRID_SIZE+2)*(ROW_SIZE+2);//2 added for periodic boundary condition ghost cells
    // Allocate host Grid used for initial setup and read back from device
    h_grid = (CELL_TYPE*)malloc(bytes);
 
    // Allocate device grids
    cudaMalloc(&d_grid, bytes);
    cudaMalloc(&d_newGrid, bytes);
    cudaMalloc(&GPU_rule_table, sizeof(int)*2*(CELL_NEIGHBOURS+1));
 
    // Assign initial population randomly
    init_random_values (h_grid);

    //total = print_total_alive (h_grid);
    //printf("Intitial Alive: %ld\n", total);

    // Copy over initial game grid (Dim-1 threads)
    cudaMemcpy(d_grid, h_grid, bytes, cudaMemcpyHostToDevice);
 
    dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE,1);
    int linGridx = (int)ceil(ROW_SIZE/(float)BLOCK_SIZE);
    int linGridy = (int)ceil(GRID_SIZE/(float)BLOCK_SIZE);
    dim3 gridSize(linGridx,linGridy,1);
 
    dim3 cpyBlockSize(BLOCK_SIZE,1,1);
    dim3 cpyGridRowsGridSize((int)ceil(ROW_SIZE/(float)cpyBlockSize.x),1,1);
    dim3 cpyGridColsGridSize((int)ceil((GRID_SIZE+2)/(float)cpyBlockSize.x),1,1);
 
    // Init rule_table
    kernel_init_rule_table<<<1,blockSize>>>(GPU_rule_table);

    // Main game loop
    for (iter = 0; iter<MAXITER; iter++) {
        ghostRows<<<cpyGridRowsGridSize, cpyBlockSize>>>(d_grid);
        ghostCols<<<cpyGridColsGridSize, cpyBlockSize>>>(d_grid);
        GOL<<<gridSize, blockSize>>>(d_grid, d_newGrid, GPU_rule_table); 

        // Swap our grids and iterate again
        d_tmpGrid = d_grid;
        d_grid = d_newGrid;
        d_newGrid = d_tmpGrid;
    }//iter loop
 
    // Copy back results and sum
    cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
 
    // Sum up alive cells and print results

    total = print_total_alive (h_grid);
    printf("Total Alive: %ld\n", total);
 
    // Release memory
    cudaFree(d_grid);
    cudaFree(d_newGrid);
    cudaFree(GPU_rule_table);
    free(h_grid);
 
    return 0;
}



void init_random_values (CELL_TYPE *h_grid) {
	int i,j,k;
	CELL_TYPE aux;
    srand(SRAND_VALUE);
	for(i = 1; i<=GRID_SIZE; i++) {
		for(j = 1; j<=ROW_SIZE; j++) {
			aux = h_grid[i*(ROW_SIZE+2)+j];
			for (k=0; k<ELEMENTS_PER_CELL; k++) {
				setSubCellH (&aux, k, rand() % 2);
			}
			h_grid[i*(ROW_SIZE+2)+j] = aux;
		}
	}
}


__global__ void kernel_init_rule_table (int *GPU_rule_table) {
    int (*rule_table)[CELL_NEIGHBOURS+1] = (int (*)[CELL_NEIGHBOURS+1]) GPU_rule_table;

    if ( threadIdx.y < 2 && threadIdx.x < (CELL_NEIGHBOURS+1) ) {
        // Init rule_table for GOL
        // Classic B3S23 GOL:
	    //rule_table[cases] = {
		//    d,d,d,a, d,d,d,d,d // DEAD is current state
		//    d,d,a,a, d,d,d,d,d // ALIVE is current state
		    if (threadIdx.y==0) 
                if (threadIdx.x >= MIN_NOF_NEIGH_FROM_DEAD_TO_ALIVE && threadIdx.x <= MAX_NOF_NEIGH_FROM_DEAD_TO_ALIVE)
			        rule_table[threadIdx.y][threadIdx.x] = ALIVE;
		        else
			        rule_table[threadIdx.y][threadIdx.x] = DEAD; 
 
		    if (threadIdx.y==1) 
                if (threadIdx.x >= MIN_NOF_NEIGH_FROM_ALIVE_TO_ALIVE && threadIdx.x <= MAX_NOF_NEIGH_FROM_ALIVE_TO_ALIVE)
			        rule_table[threadIdx.y][threadIdx.x] =  ALIVE;
		        else
			        rule_table[threadIdx.y][threadIdx.x] = DEAD;  
    }
}



__forceinline__ unsigned char getSubCellH (CELL_TYPE cell, char pos)
{
	return (cell >> (ELEMENTS_PER_CELL - 1 - pos)*8);
}

__forceinline__ void setSubCellH (CELL_TYPE *cell, char pos, unsigned char subcell)
{
	CELL_TYPE mask = 0xFF;
	CELL_TYPE maskNewCell = subcell;
	
	// Erase pos content in cell:
	mask = mask << (ELEMENTS_PER_CELL - 1 - pos)*8;
	mask = ~mask;
	*cell = *cell & mask;
	
	// Add subcell content to cell in pos:
	*cell = *cell | (maskNewCell << (ELEMENTS_PER_CELL - 1 - pos)*8);
}

__device__ unsigned char getSubCellD (CELL_TYPE cell, char pos)
{
	return (cell >> (ELEMENTS_PER_CELL - 1 - pos)*8);
}

__device__ void setSubCellD (CELL_TYPE *cell, char pos, unsigned char subcell)
{
	CELL_TYPE mask = 0xFF;
	CELL_TYPE maskNewCell = subcell;
	
	// Erase pos content in cell:
	mask = mask << (ELEMENTS_PER_CELL - 1 - pos)*8;
	mask = ~mask;
	*cell = *cell & mask;
	
	// Add subcell content to cell in pos:
	*cell = *cell | (maskNewCell << (ELEMENTS_PER_CELL - 1 - pos)*8);
}



__global__ void ghostRows(CELL_TYPE *grid)
{
    // We want id ∈ [1,GRID_SIZE]
    int id = blockDim.x * blockIdx.x + threadIdx.x + 1;

    if (id <= ROW_SIZE)
    {
        //Copy first real row to bottom ghost row
        grid[(ROW_SIZE+2)*(GRID_SIZE+1)+id] = grid[(ROW_SIZE+2)+id];
        //Copy last real row to top ghost row
        grid[id] = grid[(ROW_SIZE+2)*GRID_SIZE + id];
    }
}
 
__global__ void ghostCols(CELL_TYPE *grid)
{
    // We want id ∈ [0,SIZE+1]
    int id = blockDim.x * blockIdx.x + threadIdx.x;
 
    if (id <= GRID_SIZE+1)
    {
        //Copy first real column to right most ghost column
        grid[id*(ROW_SIZE+2)+ROW_SIZE+1] = grid[id*(ROW_SIZE+2)+1];
        //Copy last real column to left most ghost column 
        grid[id*(ROW_SIZE+2)] = grid[id*(ROW_SIZE+2) + ROW_SIZE];    
    }
}


__global__ void GOL (CELL_TYPE *grid, CELL_TYPE *newGrid, int *GPU_rule_table)
{
    // We want id ∈ [1,SIZE]
    int iy = blockDim.y * blockIdx.y + threadIdx.y + 1;
    int ix = blockDim.x * blockIdx.x + threadIdx.x + 1;
    int id = iy * (ROW_SIZE+2) + ix;
    CELL_TYPE cell, new_cell=0; 
    CELL_TYPE up_cell, down_cell, right_cell, left_cell;                // Up,down,right,left cells
    CELL_TYPE upleft_cell, downleft_cell, upright_cell, downright_cell; // Diagonal cells
    unsigned char subcell;
 
    int k, numNeighbors;
    int (*rule_table)[CELL_NEIGHBOURS+1] = (int (*)[CELL_NEIGHBOURS+1]) GPU_rule_table;
 
    if (iy <= GRID_SIZE && ix <= ROW_SIZE) {
         cell = grid[id];

        // First (0) subcell:
        up_cell = grid[id-(ROW_SIZE+2)];
        down_cell = grid[id+(ROW_SIZE+2)];
        left_cell = grid[id-1];
        upleft_cell = grid[id-(ROW_SIZE+3)];
        downleft_cell = grid[id+(ROW_SIZE+1)];

        numNeighbors = getSubCellD (up_cell, 0) + getSubCellD (down_cell, 0); // upper lower
        numNeighbors += getSubCellD (left_cell, ELEMENTS_PER_CELL-1) + getSubCellD (cell, 1); // left right
        numNeighbors += getSubCellD (upleft_cell, ELEMENTS_PER_CELL-1) + getSubCellD (downleft_cell, ELEMENTS_PER_CELL-1); // diagonals left
        numNeighbors += getSubCellD (up_cell, 1) + getSubCellD (down_cell, 1); // diagonals right
        subcell = getSubCellD (cell, 0);
        setSubCellD (&new_cell, 0, rule_table[subcell][numNeighbors]);

        // Middle subcells:
        for (k=1; k<CELL_NEIGHBOURS-1; k++) {
            numNeighbors = getSubCellD (up_cell, k) + getSubCellD (down_cell, k); // upper lower
            numNeighbors += getSubCellD (cell, k-1) + getSubCellD (cell, k+1); // left right
            numNeighbors += getSubCellD (up_cell, k-1) + getSubCellD (down_cell, k-1); // diagonals left
            numNeighbors += getSubCellD (up_cell, k+1) + getSubCellD (down_cell, k+1); // diagonals right
            subcell = getSubCellD (cell, k);
            setSubCellD (&new_cell, k, rule_table[subcell][numNeighbors]);
        }

        // Last (CELL_NEIGHBOURS-1) subcell:
        right_cell = grid[id+1];
        upright_cell = grid[id-(ROW_SIZE+1)];
        downright_cell = grid[id+(ROW_SIZE+3)];

        numNeighbors = getSubCellD (up_cell, ELEMENTS_PER_CELL-1) + getSubCellD (down_cell, ELEMENTS_PER_CELL-1); // upper lower
        numNeighbors += getSubCellD (cell, ELEMENTS_PER_CELL-2) + getSubCellD (right_cell, 0); // left right
        numNeighbors += getSubCellD (up_cell, ELEMENTS_PER_CELL-2) + getSubCellD (down_cell, ELEMENTS_PER_CELL-2); // diagonals left
        numNeighbors += getSubCellD (upright_cell, 0) + getSubCellD (downright_cell, 0); // diagonals right
        subcell = getSubCellD (cell, ELEMENTS_PER_CELL-1);
        setSubCellD (&new_cell, ELEMENTS_PER_CELL-1, rule_table[subcell][numNeighbors]);


        // Copy new_cell to newGrid:
        newGrid[id] = new_cell;

/* 
        // Get the number of neighbors for a given grid point
        numNeighbors = grid[id+(SIZE+2)] + grid[id-(SIZE+2)] //upper lower
                     + grid[id+1] + grid[id-1]             //right left
                     + grid[id+(SIZE+3)] + grid[id-(SIZE+3)] //diagonals
                     + grid[id-(SIZE+1)] + grid[id+(SIZE+1)];
 
        CELL_TYPE cell = grid[id];
        newGrid[id] = rule_table[cell][numNeighbors];
*/

    }
}


long int print_total_alive (CELL_TYPE *h_grid) {
	int i,j,k;
	CELL_TYPE aux;
	long int total = 0;
	
	for (i = 1; i<=GRID_SIZE; i++) {
		for (j = 1; j<=ROW_SIZE; j++) {
			aux = h_grid[i*(ROW_SIZE+2)+j];
			for (k=0; k<ELEMENTS_PER_CELL; k++) {
				total += getSubCellH (aux, k);
			}
		}
	}
	return total;
}



// 256 Result in console: "Total Alive: 3281"  "Initial Alive: 32708"
// 512 Result in console: "Total Alive: 11072"
// 1024 Result in console: "Total Alive: 45224"
// 2048 Result in console: "Total Alive: 182485"
// 4096 Result in console: "Total Alive: 724393"
// 8192 Result in console: "Total Alive: 2896683"
// 16384 Result in console: "Total Alive: 11547651"


