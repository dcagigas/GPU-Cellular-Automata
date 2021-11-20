#include <stdio.h>
#include <stdlib.h>
#include <openssl/rand.h> 			// For random numbers. OpenSSL library must be installed.

#define GRID_SIZE 256
#define TIMESTEP 1024

#define SEED 1
#define BLOCK_SIZE 32

#define EMPTY 0
#define HEAD 17
#define TAIL 1
#define CONDUCTOR 2

#define LOOKUP_TABLE_LIMIT (HEAD*8+1)   // One cell can be surrounded by 8 head cells. That's the maxium sum.

#define STRIDE 8        // Used to draw the initial wire world (conductor concentric squares)
#define P_HEAD 0.001    // Probability of setting a conductor cell to head in the initial wire world

// Uncomment for grid snapshots (image) writing:
//#define __VIDEO

typedef char CELL_TYPE;

void init_lookup_table (CELL_TYPE *lookup_table);
float r4_uniform_01 ( int *seed );
uint8_t* generate_rgb (int width, int height, CELL_TYPE *grid_, uint8_t *rgb);
void write_CA_screenshot_image (const char *filename, int width, int height, uint8_t *img);


int seed = SEED; 	// Seed for random uniform number output.


void init_lookup_table (CELL_TYPE *lookup_table_) {
    int i, j;
    CELL_TYPE (*lookup_table)[LOOKUP_TABLE_LIMIT] = (CELL_TYPE(*)[LOOKUP_TABLE_LIMIT]) lookup_table_;

    // Init lookup_table:
    // 0: Empty. If a cell is EMPTY remains EMPTY 
    for (j=0; j<LOOKUP_TABLE_LIMIT; j++) {
        lookup_table[EMPTY][j] = EMPTY;
    }
    // 2: Head. A HEAD cell turns always into TAIL 
    for (j=0; j<LOOKUP_TABLE_LIMIT; j++) {
        lookup_table[HEAD][j] = TAIL;
    }
    // 3: Tail. A TAIL cell turns always into CONDUCTOR 
    for (j=0; j<LOOKUP_TABLE_LIMIT; j++) {
        lookup_table[TAIL][j] = CONDUCTOR;
    }
    // 4: Conductor. A CONDUCTOR cell turns into HEAD if exactly one or two of the neighbouring cells are HEADS), 
    //    otherwise remains CONDUCTOR.
    for (j=0; j<LOOKUP_TABLE_LIMIT; j++) {
        lookup_table[CONDUCTOR][j] = CONDUCTOR;
    }
    for (j=HEAD; j<(HEAD*2 + CONDUCTOR*6 + 1); j++) {
        lookup_table[CONDUCTOR][j] = HEAD;
    }

    // Fill with EMPTY not used lookup_table entries
    for (i=CONDUCTOR+1; i<HEAD-1; i++)
        for (j=0; j<LOOKUP_TABLE_LIMIT; j++) 
            lookup_table[i][j] = EMPTY;
}



__global__ void wire_world (CELL_TYPE *grid, CELL_TYPE *newGrid, CELL_TYPE *lookup_table_)
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
    int id = iy * (GRID_SIZE) + ix;
 
    CELL_TYPE (*lookup_table)[LOOKUP_TABLE_LIMIT] = (CELL_TYPE(*)[LOOKUP_TABLE_LIMIT]) lookup_table_;
    int numHeadNeighbors = 0;
 
    if ( (iy>0 && iy < GRID_SIZE-1) && (ix>0 && ix < GRID_SIZE-1) ) {
        // Get the number of neighbors for a given grid point
        numHeadNeighbors = (grid[id+(GRID_SIZE)]) + (grid[id-(GRID_SIZE)])              // upper lower
                     + (grid[id+1]) + (grid[id-1])                              // right left
                     + (grid[id+(GRID_SIZE)+1]) + (grid[id-(GRID_SIZE)+1])     // diagonals
                     + (grid[id-(GRID_SIZE)-1]) + (grid[id+(GRID_SIZE)-1]);
        CELL_TYPE cell = grid[id];
        newGrid[id] = lookup_table[cell][numHeadNeighbors];
    }
}



long int  init_wire_world (CELL_TYPE *grid) {
    int i, j, k;
    //long int conductors;
    long int heads = 0;

    // Set grid to empty
    for(i = 0; i<GRID_SIZE; i++) {
        for(j = 0; j<GRID_SIZE; j++) {
            grid[i*(GRID_SIZE)+j] = EMPTY;
        }
    }
    
    // Set a "conductor cross" in the middle of the grid
    for(i = 0; i<GRID_SIZE; i++) {
        grid[i*(GRID_SIZE)+(GRID_SIZE/2)] = CONDUCTOR;
    }
    for(j = 0; j<GRID_SIZE; j++) {
        grid[(GRID_SIZE/2)*(GRID_SIZE) + j] = CONDUCTOR;
    }

    // Set concentric squares in the grid with STRIDE padding.
    for (k = STRIDE; k < (GRID_SIZE)/2; k = k + STRIDE) {
        // Draw square:
        // Draw left side and right side:
        for (i=k*(GRID_SIZE)+k; i<  (GRID_SIZE)*(GRID_SIZE)-k*(GRID_SIZE); i+=GRID_SIZE) {
            grid[i] = CONDUCTOR;
            grid[i+(GRID_SIZE)-k-k] = CONDUCTOR;
         }
        // Draw up side and down side:
        for (j=k*(GRID_SIZE)+k; j< k*(GRID_SIZE)+(GRID_SIZE)-k; j++) {
            grid[j] = CONDUCTOR;
        }
        for (j=(GRID_SIZE)*(GRID_SIZE)-k*(GRID_SIZE)+(GRID_SIZE)-k; j >= (GRID_SIZE)*(GRID_SIZE)-k*(GRID_SIZE)+k ; j--) {
            grid[j] = CONDUCTOR;
        }        
    }
    
    // Set initial random heads
    for(i = 0; i<GRID_SIZE; i++) {
        for(j = 0; j<GRID_SIZE; j++) {
            if ( grid[i*(GRID_SIZE)+j] == CONDUCTOR && r4_uniform_01(&seed) < P_HEAD) {
                grid[i*(GRID_SIZE)+j] = HEAD;
                heads++;
            }
        }
    }
    
    return heads;
}



 
int main(int argc, char* argv[])
{
    int i,j,iter;
    long unsigned int total_heads = 0;
    CELL_TYPE* h_grid; //Grid on host
    CELL_TYPE* d_grid; //Grid on device
    CELL_TYPE* d_newGrid; //Second grid used on device only
    CELL_TYPE* d_tmpGrid; //tmp grid pointer used to switch between grid and newGrid
 
    CELL_TYPE* h_lookup_table; // Look-up table
    CELL_TYPE* d_lookup_table; 

    //int dim = GRID_SIZE; //Linear dimension of our grid - not counting ghost cells
 
    size_t bytes = sizeof(CELL_TYPE)*(GRID_SIZE)*(GRID_SIZE);//2 added for periodic boundary condition ghost cells
    // Allocate host Grid used for initial setup and read back from device
    h_grid = (CELL_TYPE *)malloc(bytes);

    #ifdef __VIDEO
        char image_name[80];
        uint8_t *rgb= (uint8_t *)malloc (3 * sizeof(uint8_t) *GRID_SIZE*GRID_SIZE);	
    #endif
 
    // Allocate device grids
    cudaMalloc(&d_grid, bytes);
    cudaMalloc(&d_newGrid, bytes);
 
    // Alocate look-up tables:
    size_t bytes2 = sizeof(CELL_TYPE)*(HEAD+1)*(LOOKUP_TABLE_LIMIT);
    h_lookup_table =  (CELL_TYPE *)malloc(bytes2);
    cudaMalloc(&d_lookup_table, bytes2);

    // Init look-up table:
    init_lookup_table (h_lookup_table);

    total_heads = init_wire_world (h_grid);

    printf("Initial Heads: %ld \n", total_heads);
 
    // Copy over initial grid (GRID_SIZE threads)
    cudaMemcpy(d_grid, h_grid, bytes, cudaMemcpyHostToDevice);
 
    // Copy lookup-table
    cudaMemcpy(d_lookup_table, h_lookup_table, bytes2, cudaMemcpyHostToDevice);

    dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE,1);
    int linGrid = (int)ceil(GRID_SIZE/(float)BLOCK_SIZE);
    dim3 gridSize(linGrid,linGrid,1);
  
    // Main loop
    for (iter = 0; iter<TIMESTEP; iter++) {
    #ifdef __VIDEO
		if (iter==0 || iter==50 || iter==100 || iter==600 || iter==300 || iter==1000) {
            cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
			rgb = generate_rgb (GRID_SIZE, GRID_SIZE, h_grid, rgb);
			sprintf (image_name, "%s%d.ppm", "ca_step_", iter);
			write_CA_screenshot_image (image_name, GRID_SIZE, GRID_SIZE, rgb);
		}
    #endif
        wire_world<<<gridSize, blockSize>>>(d_grid, d_newGrid, d_lookup_table);

        // Swap our grids and iterate again
        d_tmpGrid = d_grid;
        d_grid = d_newGrid;
        d_newGrid = d_tmpGrid;
    }//iter loop
 
    // Copy back results and sum
    cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
 
    total_heads = 0; 
    // Sum up head cells and print results
    for (i = 0; i<GRID_SIZE; i++) {
        for (j = 1; j<GRID_SIZE; j++) {
            if (h_grid[i*(GRID_SIZE)+j] == HEAD)
                total_heads++;
        }
    }
    printf("Final heads: %ld \n", total_heads);
 
    // Release memory
    cudaFree(d_grid);
    cudaFree(d_newGrid);
    free(h_grid);
    free (h_lookup_table);
    cudaFree (d_lookup_table);

    return 0;
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
	int x, y, cur;
	
	// Transform "grid" (CA states grid) into a similar rgb image format.
	CELL_TYPE (*grid)[height] = (CELL_TYPE (*)[height]) grid_;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			cur = 3 * (y * width + x);
			if (grid[y][x] == EMPTY) {
				// Empty: Black
				rgb[cur + 0] = 128;
				rgb[cur + 1] = 128;
				rgb[cur + 2] = 128;
			} else if (grid[y][x] == HEAD) {
				// Head: Blue
				rgb[cur + 0] = 0;
				rgb[cur + 1] = 0;
				rgb[cur + 2] = 255;
			} else if (grid[y][x] == TAIL) {
				// Tail: Red
				rgb[cur + 0] = 255;
				rgb[cur + 1] = 0;
				rgb[cur + 2] = 0;
			} else if (grid[y][x] == CONDUCTOR) {
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



// 256 Result in console: initial heads: 17, final heads: 1908
// 512 Result in console: initial heads: 42, final heads: 2092
// 1024 Result in console: initial heads: 139, final heads: 1911
// 2048 Result in console: initial heads: 538, final heads: 2612
// 4096 Result in console: initial heads: 2131, final heads: 3833
// 8192 Result in console: initial heads: 8442, final heads: 6512
// 16384 Result in console: initial heads: 33482, final heads: 15822

