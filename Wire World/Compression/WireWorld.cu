#include <stdio.h>
#include <stdlib.h>
#include <openssl/rand.h> 			// For random numbers. OpenSSL library must be installed.

#define GRID_SIZE 256
#define TIMESTEP 1024

#define SEED 1
#define BLOCK_SIZE 32
#define GRID_ROWS GRID_SIZE/2

#define EMPTY 0
#define HEAD 3
#define TAIL 1
#define CONDUCTOR 2

#define STRIDE 8        // Used to draw the initial wire world (conductor concentric squares)
#define P_HEAD 0.001    // Probability of setting a conductor cell to head in the initial wire world

/*
#define _getSubCellHost_H(x)    (0x0F & x>>4)
#define _getSubCellHost_L(x)    (0x0F & cell)
#define _setSubCellHost_H(x, v) 
#define _setSubCellHost_L(x,v) 
*/

// Uncomment for grid snapshots (image) writing:
//#define __VIDEO

typedef unsigned char CELL_TYPE;

long int init_wire_world (CELL_TYPE *grid);
long int count_total_heads (CELL_TYPE *h_grid);
float r4_uniform_01 ( int *seed );
uint8_t* generate_rgb (int width, int height, CELL_TYPE *grid_, uint8_t *rgb);
void write_CA_screenshot_image (const char *filename, int width, int height, uint8_t *img);

__global__ void wire_world (CELL_TYPE *grid, CELL_TYPE *newGrid);

__forceinline__ CELL_TYPE getSubCellHost_H (CELL_TYPE cell);
__forceinline__ CELL_TYPE getSubCellHost_L (CELL_TYPE cell);
__forceinline__ void setSubCellHost_H (CELL_TYPE *cell, CELL_TYPE value);
__forceinline__ void setSubCellHost_L (CELL_TYPE *cell, CELL_TYPE value);


__device__ CELL_TYPE getSubCellDevice_H (CELL_TYPE cell);
__device__ CELL_TYPE getSubCellDevice_L (CELL_TYPE cell);
__device__ void setSubCellDevice_H (CELL_TYPE *cell, CELL_TYPE value);
__device__ void setSubCellDevice_L (CELL_TYPE *cell, CELL_TYPE value);


int seed = SEED; 	// Seed for random uniform number output.


__global__ void wire_world (CELL_TYPE *grid, CELL_TYPE *newGrid)
{
    /*
    Cells behave as follows:
    empty → empty,
    (electron) head → (electron) tail,
    (electron tail) → conductor,
    conductor → (electron) head if exactly one or two of the neighbouring cells are (electron heads), otherwise remains conductor. 
    */
    // We want id ∈ [1,dim]
    int iy = blockDim.y * blockIdx.y + threadIdx.y;
    int ix = blockDim.x * blockIdx.x + threadIdx.x;
    int id = iy * (GRID_ROWS) + ix;
    //CELL_TYPE (*grid)[GRID_ROWS] = (CELL_TYPE (*)[GRID_ROWS]) grid_;
 
    int numHeadNeighbors = 0;
 
    if ( (iy>0 && iy < GRID_SIZE-1) && (ix < GRID_SIZE) ) {
        CELL_TYPE cell = grid[id];
        CELL_TYPE subcell, new_cell=0;

        // HIGH SUBCELL PROCESS:
        subcell = getSubCellDevice_H (cell);
        if (ix > 0) {
            // Here we have explicitly all of the rules
            if (subcell == HEAD) {
                setSubCellDevice_H (&new_cell, TAIL);
            } else if (subcell == TAIL) {
               setSubCellDevice_H (&new_cell, CONDUCTOR);
            } else if (subcell == CONDUCTOR) {
                // Get the number of neighbors for a given grid point
                numHeadNeighbors = (getSubCellDevice_H(grid[id+(GRID_ROWS)])==HEAD) + (getSubCellDevice_H(grid[id-(GRID_ROWS)])==HEAD)  // upper lower
                     + (getSubCellDevice_L(cell)==HEAD) + (getSubCellDevice_L(grid[id-1])==HEAD)                                        // right left
                     + (getSubCellDevice_L(grid[id+(GRID_ROWS)])==HEAD) + (getSubCellDevice_L(grid[id-(GRID_ROWS)])==HEAD)              // diagonals
                     + (getSubCellDevice_L(grid[id-(GRID_ROWS)-1])==HEAD) + (getSubCellDevice_L(grid[id+(GRID_ROWS)-1])==HEAD);
                if (numHeadNeighbors == 1 || numHeadNeighbors == 2)
                    setSubCellDevice_H (&new_cell, HEAD);
                else 
                    setSubCellDevice_H (&new_cell, CONDUCTOR);
            } else {
                // subcell is EMPTY:
                setSubCellDevice_H (&new_cell, EMPTY);
            }
        } else {
            // Left Border of the grid. It must be copied but not computed as a regular subcell.
            setSubCellDevice_H (&new_cell, subcell);
        }

        // LOW SUBCELL PROCESS:
        subcell = getSubCellDevice_L (cell);
        if (ix < GRID_SIZE-1) {
            // Here we have explicitly all of the rules
            if (subcell == HEAD) {
                setSubCellDevice_L (&new_cell, TAIL);
            } else if (subcell == TAIL) {
               setSubCellDevice_L (&new_cell, CONDUCTOR);
            } else if (subcell == CONDUCTOR) {
                // Get the number of neighbors for a given grid point
                numHeadNeighbors = (getSubCellDevice_L(grid[id+(GRID_ROWS)])==HEAD) + (getSubCellDevice_L(grid[id-(GRID_ROWS)])==HEAD)      // upper lower
                     + (getSubCellDevice_H(grid[id+1])==HEAD) + (getSubCellDevice_H(cell)==HEAD)                                       // left right 
                     + (getSubCellDevice_H(grid[id+(GRID_ROWS)+1])==HEAD) + (getSubCellDevice_H(grid[id-(GRID_ROWS)+1])==HEAD)          // diagonals
                     + (getSubCellDevice_H(grid[id-(GRID_ROWS)])==HEAD) + (getSubCellDevice_H(grid[id+(GRID_ROWS)])==HEAD);
                if (numHeadNeighbors == 1 || numHeadNeighbors == 2)
                    setSubCellDevice_L (&new_cell, HEAD);
                else 
                    setSubCellDevice_L (&new_cell, CONDUCTOR);
            } else {
                setSubCellDevice_L (&new_cell, EMPTY);
            }
        } else {
            // Right Border of the grid. It must be copied but not computed as a regular subcell.
            setSubCellDevice_L (&new_cell, subcell);
        }

        newGrid[id] = new_cell;
    }
}



// Host auxiliary functions: 

__forceinline__ CELL_TYPE getSubCellHost_H (CELL_TYPE cell) {
    return (0x0F & cell>>4);
}

__forceinline__ CELL_TYPE getSubCellHost_L (CELL_TYPE cell) {
    return (0x0F & cell);
}

__forceinline__ void setSubCellHost_H (CELL_TYPE *cell, CELL_TYPE value) {
    *cell = (*cell & 0x0F) | (value <<4);
}

__forceinline__ void setSubCellHost_L (CELL_TYPE *cell, CELL_TYPE value) {
    *cell = (*cell & 0xF0) | (0x0F & value);
}


// Device auxiliary functions: 

__device__ CELL_TYPE getSubCellDevice_H (CELL_TYPE cell) {
    return (0x0F & cell>>4);
}

__device__ CELL_TYPE getSubCellDevice_L (CELL_TYPE cell) {
    return (0x0F & cell);
}

__device__ void setSubCellDevice_H (CELL_TYPE *cell, CELL_TYPE value) {
    *cell = (*cell & 0x0F) | (value <<4);
}

__device__ void setSubCellDevice_L (CELL_TYPE *cell, CELL_TYPE value) {
    *cell = (*cell & 0xF0) | (0x0F & value);
}



long int init_wire_world (CELL_TYPE *grid_) {
    int i, j, ii, jj, k;
    //long int conductors;
    long int heads = 0;

    CELL_TYPE (*grid)[GRID_ROWS] = (CELL_TYPE(*)[GRID_ROWS]) grid_;

    // Set grid to empty
    for(i = 0; i<GRID_SIZE; i++) {
        for(j = 0; j<GRID_ROWS; j++) {
            //grid[i][j] = EMPTY;
            setSubCellHost_H (&grid[i][j], EMPTY);
            setSubCellHost_L (&grid[i][j], EMPTY);
        }
    }

    // Set a "conductor cross" in the middle of the grid
    jj = GRID_ROWS/2;
    // Up-Down
	for (i = 0; i<GRID_SIZE; i++) {
        //aux = grid[i][jj];
        //setSubCellHost_H (&aux, CONDUCTOR);
        //grid[i][jj] = aux;
        setSubCellHost_H (&grid[i][jj], CONDUCTOR);
    }
    // Left-Right
    ii = GRID_SIZE/2;
	for (j = 0; j<GRID_ROWS; j++) {
        setSubCellHost_H (&grid[ii][j], CONDUCTOR);
        setSubCellHost_L (&grid[ii][j], CONDUCTOR);
    }


    // Set concentric squares in the grid with STRIDE padding.
    for (k = STRIDE; k < (GRID_SIZE)/2; k = k + STRIDE) {
        // Draw square:
        // Draw left side and right side:
        jj = k/2;
        for (i=k; i < GRID_SIZE-k; i++) {
            // Left side:
            //aux = grid[i][jj];
            //setSubCellH (&aux, kk, CONDUCTOR);
            //grid[i][jj] = aux;
            setSubCellHost_H (&grid[i][jj], CONDUCTOR);
            // Right side:
            setSubCellHost_L (&grid[i][GRID_ROWS-jj], CONDUCTOR);
        }
        // Draw up side and down side:
        for (j=jj; j < GRID_ROWS - jj+1; j++) {
            // Top side:
            setSubCellHost_H (&grid[k][j], CONDUCTOR);
            setSubCellHost_L (&grid[k][j], CONDUCTOR);
            // Down side:
            //aux = grid[GRID_SIZE-k-1][j];
            setSubCellHost_H (&grid[GRID_SIZE-k][j], CONDUCTOR);
            setSubCellHost_L (&grid[GRID_SIZE-k][j], CONDUCTOR);
        }
    }
    
    // Set initial random heads
    for(i = 0; i<GRID_SIZE; i++) {
        for(j = 0; j<GRID_ROWS; j++) {
            if ( getSubCellHost_H (grid[i][j]) == CONDUCTOR && r4_uniform_01(&seed) < P_HEAD) {
                setSubCellHost_H (&grid[i][j], HEAD);
//printf ("\t (%d, %d) %x\n", i, j, grid[i][j]);
                heads++;
            }
            if ( getSubCellHost_L (grid[i][j]) == CONDUCTOR && r4_uniform_01(&seed) < P_HEAD) {
                setSubCellHost_L (&grid[i][j], HEAD);
//printf ("\t (%d, %d) %x\n", i, j, grid[i][j]);
                heads++;
            }
        }
    }
    
    return heads;
}



 
int main(int argc, char* argv[])
{
    int iter;
    long unsigned int total_heads = 0;
    CELL_TYPE* h_grid; //Grid on host
    CELL_TYPE* d_grid; //Grid on device
    CELL_TYPE* d_newGrid; //Second grid used on device only
    CELL_TYPE* d_tmpGrid; //tmp grid pointer used to switch between grid and newGrid
 
    //int dim = GRID_SIZE; //Linear dimension of our grid - not counting ghost cells
 
    size_t bytes = sizeof(CELL_TYPE)*(GRID_SIZE)*(GRID_ROWS);
    // Allocate host Grid used for initial setup and read back from device
    h_grid = (CELL_TYPE *)malloc(bytes);

    #ifdef __VIDEO
        char image_name[80];
        uint8_t *rgb= (uint8_t *)malloc (3 * sizeof(uint8_t) *GRID_SIZE*GRID_SIZE);	
    #endif
 
    // Allocate device grids
    cudaMalloc(&d_grid, bytes);
    cudaMalloc(&d_newGrid, bytes);
 
    total_heads = init_wire_world (h_grid);

    printf("Initial Heads: %ld \n", total_heads);
 
    // Copy over initial game grid (GRID_SIZE threads)
    cudaMemcpy(d_grid, h_grid, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_newGrid, h_grid, bytes, cudaMemcpyHostToDevice);
 
    dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE,1);
    int linGridx = (int)ceil((GRID_ROWS)/(float)BLOCK_SIZE);
    int linGridy = (int)ceil(GRID_SIZE/(float)BLOCK_SIZE);
    dim3 gridSize(linGridx,linGridy,1);
 
    // Main game loop
    for (iter = 0; iter<TIMESTEP; iter++) {
    //for (iter = 0; iter<6; iter++) {
    #ifdef __VIDEO
		if (iter==0 || iter==50 || iter==100 || iter==600 || iter==300 || iter==1000) {
		//if (iter==0 || iter==1 || iter==2 || iter==3 || iter==4 || iter==5) {
            cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
			rgb = generate_rgb (GRID_SIZE, GRID_SIZE, h_grid, rgb);
			sprintf (image_name, "%s%d.ppm", "ca_step_", iter);
			write_CA_screenshot_image (image_name, GRID_SIZE, GRID_SIZE, rgb);
		}
    #endif
        wire_world<<<gridSize, blockSize>>>(d_grid, d_newGrid);
        cudaDeviceSynchronize();

        // Swap our grids and iterate again
        d_tmpGrid = d_grid;
        d_grid = d_newGrid;
        d_newGrid = d_tmpGrid;
    }//iter loop
 
    // Copy back results and sum
    cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
 
    // Sum up head cells and print results
    total_heads = count_total_heads (h_grid);
    printf("Final heads: %ld \n", total_heads);
 
    // Release memory
    cudaFree(d_grid);
    cudaFree(d_newGrid);
    free(h_grid);
 
    return 0;
}


long int count_total_heads (CELL_TYPE *h_grid) {
	int i,j;
	CELL_TYPE aux;
	long int total = 0;
	
	for (i = 1; i<GRID_SIZE-1; i++) {
		for (j = 0; j<GRID_ROWS; j++) {
			aux = h_grid[i*(GRID_ROWS)+j];
            if (getSubCellHost_H (aux) == HEAD) 
				    total++;
            if (getSubCellHost_L (aux) == HEAD) 
				    total++;
		}
	}
	return total;
}


/******************************************************************************/
// Used for random numbers:
float r4_uniform_01 ( int *seed )
	
	/******************************************************************************/
	/*
	https://people.sc.fsu.edu/~jburkardt/c_src/normal/normal.html
	Purpose:
	
	R4_UNIFORM_01 returns a unit pseudorandom R4.
	
	Discussion:
	
	This routine implements the recursion
	
	seed = 16807 * seed mod ( 2^31 - 1 )
	r4_uniform_01 = seed / ( 2^31 - 1 )
	
	The integer arithmetic never requires more than 32 bits,
	including a sign bit.
	
	If the initial seed is 12345, then the first three computations are
	
	Input     Output      R4_UNIFORM_01
	SEED      SEED
	
	12345   207482415  0.096616
	207482415  1790989824  0.833995
	1790989824  2035175616  0.947702
	
	Licensing:
	
	This code is distributed under the GNU LGPL license. 
	
	Modified:
	
	16 November 2004
	
	Author:
	
	John Burkardt
	
	Reference:
	
	Paul Bratley, Bennett Fox, Linus Schrage,
	A Guide to Simulation,
	Springer Verlag, pages 201-202, 1983.
	
	Pierre L'Ecuyer,
	Random Number Generation,
	in Handbook of Simulation
	edited by Jerry Banks,
	Wiley Interscience, page 95, 1998.
	
	Bennett Fox,
	Algorithm 647:
	Implementation and Relative Efficiency of Quasirandom
	Sequence Generators,
	ACM Transactions on Mathematical Software,
	Volume 12, Number 4, pages 362-376, 1986.
	
	Peter Lewis, Allen Goodman, James Miller,
	A Pseudo-Random Number Generator for the System/360,
	IBM Systems Journal,
	Volume 8, pages 136-143, 1969.
	
	Parameters:
	
	Input/output, int *SEED, the "seed" value.  Normally, this
	value should not be 0.  On output, SEED has been updated.
	
	Output, float R4_UNIFORM_01, a new pseudorandom variate, strictly between
	0 and 1.
	*/
{
	const int i4_huge = 2147483647;
	int k;
	float r;
	
	k = *seed / 127773;
	
	*seed = 16807 * ( *seed - k * 127773 ) - k * 2836;
	
	if ( *seed < 0 )
	{
		*seed = *seed + i4_huge;
	}
	/*
	Although SEED can be represented exactly as a 32 bit integer,
	it generally cannot be represented exactly as a 32 bit real number!
	*/
	r = ( float ) ( *seed ) * 4.656612875E-10;
	
	return r;
}




// Creates a grid (health states) rgb screenshot.
uint8_t* generate_rgb (int width, int height, CELL_TYPE *grid_, uint8_t *rgb) {
	int x, y, jj, k, cur;
    CELL_TYPE *cell;
    CELL_TYPE value;
	
	// Transform "grid" (CA states grid) into a similar rgb image format.
	CELL_TYPE (*grid)[GRID_ROWS] = (CELL_TYPE (*)[GRID_ROWS]) grid_;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			cur = 3 * (y * width + x);
            jj = x/2 ;
            k = x % 2;
            cell = (CELL_TYPE *)&grid[y][jj];
            if (k==0) 
                value = getSubCellHost_H (*cell);
            else
                value = getSubCellHost_L (*cell);
			if (value == EMPTY) {
				// Empty: Grey
				rgb[cur + 0] = 128;
				rgb[cur + 1] = 128;
				rgb[cur + 2] = 128;
			} else if (value == HEAD) {
				// Head: Blue
				rgb[cur + 0] = 0;
				rgb[cur + 1] = 0;
				rgb[cur + 2] = 255;
			} else if (value == TAIL) {
				// Tail: Red
				rgb[cur + 0] = 255;
				rgb[cur + 1] = 0;
				rgb[cur + 2] = 0;
			} else if (value == CONDUCTOR) {
				// Conductor: Yellow
				rgb[cur + 0] = 255;
				rgb[cur + 1] = 255;
				rgb[cur + 2] = 0;
			} else {
				// Any other state: White
				rgb[cur + 0] = 255;
				rgb[cur + 1] = 255;
				rgb[cur + 2] = 255;
			}
		}
	}
	return rgb;
}


// Writes to disk a CA state screenshot.
void write_CA_screenshot_image (const char *filename, int width, int height, uint8_t *img)
{
	FILE *fp;
	//open file for output
	fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "Unable to open file '%s'\n", filename);
		exit(1);
	}
	
	//write the header file
	//image format
	fprintf(fp, "P6\n");
	
	//comments
	fprintf(fp, "# Created by %s\n","WireWorld.cu");
	
	//image size
	fprintf(fp, "%d %d\n", width, height);
	
	// rgb component depth
	fprintf(fp, "%d\n",255);
	
	// pixel data
	fwrite(img, 3 * width, height, fp);
	fclose(fp);
}



// 256 Result in console: initial heads: 17, final heads: 2216
// 512 Result in console: initial heads: 42, final heads: 2508
// 1024 Result in console: initial heads: 139, final heads: 2230
// 2048 Result in console: initial heads: 538, final heads: 2647
// 4096 Result in console: initial heads: 2132, final heads: 4131
// 8192 Result in console: initial heads: 8442, final heads: 5751
// 16384 Result in console: initial heads: 33484, final heads: 15860

