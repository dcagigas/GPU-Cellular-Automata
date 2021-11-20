#include <stdio.h>
#include <stdlib.h>
 
#define SRAND_VALUE 1985
#define BLOCK_SIZE_x 32
#define BLOCK_SIZE_y 16
 
__global__ void ghostRows(int dim, int* grid)
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
__global__ void ghostCols(int dim, int* grid)
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
 
__global__ void GOL(int dim, int *grid, int *newGrid)
{
        int iy = (blockDim.y -2) * blockIdx.y + threadIdx.y;
        int ix = (blockDim.x -2) * blockIdx.x + threadIdx.x;
        int id = iy * (dim+2) + ix;
 
        int i = threadIdx.y;
        int j = threadIdx.x;
        int numNeighbors;
 
        // Declare the shared memory on a per block level
        __shared__ int s_grid[BLOCK_SIZE_y][BLOCK_SIZE_x];
 
       // Copy cells into shared memory
       if (ix <= dim+1 && iy <= dim+1)
           s_grid[i][j] = grid[id];
 
       //Sync all threads in block
        __syncthreads();
 
       if (iy <= dim && ix <= dim) {
           if(i != 0 && i !=blockDim.y-1 && j != 0 && j !=blockDim.x-1) {
 
               // Get the number of neighbors for a given grid point
               numNeighbors = s_grid[i+1][j] + s_grid[i-1][j] //upper lower
                            + s_grid[i][j+1] + s_grid[i][j-1] //right left
                            + s_grid[i+1][j+1] + s_grid[i-1][j-1] //diagonals
                            + s_grid[i-1][j+1] + s_grid[i+1][j-1];
 
                int cell = s_grid[i][j];
 
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
           }
       }
 
}
 
int main(int argc, char* argv[])
{
    int i,j,iter;
    int *h_grid; //Grid on host
    int *d_grid; //Grid on device
    int *d_newGrid; //Second grid used on device only
    int *d_tmpGrid; //tmp grid pointer used to switch between grid and newGrid
 
    int dim = 8192; //Linear dimension of our grid - not counting ghost cells
    int maxIter = 1<<10; //Number of game steps
 
    size_t bytes = sizeof(int)*(dim+2)*(dim+2);
    // Allocate host Grid used for initial setup and read back from device
    h_grid = (int*)malloc(bytes);
 
    // Allocate device grids
    cudaMalloc(&d_grid, bytes);
    cudaMalloc(&d_newGrid, bytes);
 
    // Assign initial population randomly
    srand(SRAND_VALUE);
    for(i = 1; i<=dim; i++) {
        for(j = 1; j<=dim; j++) {
            h_grid[i*(dim+2)+j] = rand() % 2;
        }
    }
 
    cudaFuncSetCacheConfig(GOL, cudaFuncCachePreferShared);
 
    // Copy over initial game grid (Dim-1 threads)
    cudaMemcpy(d_grid, h_grid, bytes, cudaMemcpyHostToDevice);
 
    dim3 blockSize(BLOCK_SIZE_x, BLOCK_SIZE_y,1);
    int linGrid_x = (int)ceil(dim/(float)(BLOCK_SIZE_x-2));
    int linGrid_y = (int)ceil(dim/(float)(BLOCK_SIZE_y-2));
    dim3 gridSize(linGrid_x,linGrid_y,1);
 
    dim3 cpyBlockSize(BLOCK_SIZE_x,1,1);
    dim3 cpyGridRowsGridSize((int)ceil(dim/(float)cpyBlockSize.x),1,1);
    dim3 cpyGridColsGridSize((int)ceil((dim+2)/(float)cpyBlockSize.x),1,1);
 
    // Main game loop
    for (iter = 0; iter<maxIter; iter++) {
 
        ghostRows<<<cpyGridRowsGridSize, cpyBlockSize>>>(dim, d_grid);
        ghostCols<<<cpyGridColsGridSize, cpyBlockSize>>>(dim, d_grid);
        GOL<<<gridSize, blockSize>>>(dim, d_grid, d_newGrid);
 
        // Swap our grids and iterate again
        d_tmpGrid = d_grid;
        d_grid = d_newGrid;
        d_newGrid = d_tmpGrid;
    }//iter loop
 
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
        printf("CUDA error %s\n",cudaGetErrorString(error));
 
    // Copy back results and sum
    cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
 
    // Sum up alive cells and print results
    int total = 0;
    for (i = 1; i<=dim; i++) {
        for (j = 1; j<=dim; j++) {
            total += h_grid[i*(dim+2)+j];
        }
    }
    printf("Total Alive: %d\n", total);
 
    cudaFree(d_grid);
    cudaFree(d_newGrid);
    free(h_grid);
 
    return 0;
}
