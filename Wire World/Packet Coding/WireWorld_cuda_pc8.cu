#include <stdio.h>
#include <stdlib.h>
#include <openssl/rand.h> 			// For random numbers. OpenSSL library must be installed.

#define GRID_SIZE 256                       // Grid dimension: SIZExSIZE
#define MAXITER 1024                        // Number of game steps

// Uncomment for grid snapshots (image) writing:
//#define __VIDEO

typedef unsigned long int CELL_TYPE;
typedef unsigned char SUBCELL_TYPE;
#define ELEMENTS_PER_CELL 8                 // "int_64" data type per cell
#define ROW_SIZE GRID_SIZE/ELEMENTS_PER_CELL    // Real grid dimension

#define SEED 1
#define BLOCK_SIZE 32

#define CELL_NEIGHBOURS 8

// WireWorld cell states:
#define EMPTY 0
#define HEAD 17
#define TAIL 1
#define CONDUCTOR 2

#define LOOKUP_TABLE_LIMIT (HEAD*8+1)   // One cell can be surrounded by 8 head cells. That's the maxium sum.

#define STRIDE 8        // Used to draw the initial wire world (conductor concentric squares)
#define P_HEAD 0.001    // Probability of setting a conductor cell to head in the initial wire world

void init_lookup_table (CELL_TYPE *lookup_table);
long int  init_wire_world (CELL_TYPE *grid);
long int count_total_heads (CELL_TYPE *h_grid);
float r4_uniform_01 ( int *seed );
uint8_t* generate_rgb (int width, int height, CELL_TYPE *grid_, uint8_t *rgb);
void write_CA_screenshot_image (const char *filename, int width, int height, uint8_t *img);

__forceinline__ SUBCELL_TYPE getSubCellH (CELL_TYPE cell, char pos);
__forceinline__ void setSubCellH (CELL_TYPE *cell, char pos, SUBCELL_TYPE subcell);
__device__ unsigned char getSubCellD (CELL_TYPE cell, char pos);
__device__ void setSubCellD (CELL_TYPE *cell, char pos, SUBCELL_TYPE subcell);
__global__ void wire_world (CELL_TYPE *grid, CELL_TYPE *newGrid, SUBCELL_TYPE *lookup_table);

int seed = SEED; 	// Seed for random uniform number output.


void init_lookup_table (SUBCELL_TYPE *lookup_table_) {
    int i, j;
    SUBCELL_TYPE (*lookup_table)[LOOKUP_TABLE_LIMIT] = (SUBCELL_TYPE(*)[LOOKUP_TABLE_LIMIT]) lookup_table_;

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
    for (i=CONDUCTOR+1; i<HEAD; i++)
        for (j=0; j<LOOKUP_TABLE_LIMIT; j++) 
            lookup_table[i][j] = EMPTY;
}



long int  init_wire_world (CELL_TYPE *grid_) {
    int i, j, k, h;
    int ii,jj,kk;
    int dirty;
    long int heads = 0;
    CELL_TYPE aux;

    CELL_TYPE (*grid)[ROW_SIZE] = (CELL_TYPE(*)[ROW_SIZE]) grid_;

    // Set grid to empty
	for (i = 0; i<GRID_SIZE; i++) {
		for (j = 0; j<ROW_SIZE; j++) {
			aux = grid[i][ROW_SIZE+0];
			for (k=0; k<ELEMENTS_PER_CELL; k++) {
                setSubCellH (&aux, k, EMPTY);
			}
			grid[i][ROW_SIZE+0] = aux;
		}
	}

    // Set concentric squares in the grid with STRIDE padding.
    for (k = STRIDE; k < (GRID_SIZE)/2; k = k + STRIDE) {
        // Draw square:
        // Draw left side and right side:
        jj = k/ELEMENTS_PER_CELL;
        kk = k % ELEMENTS_PER_CELL;
        for (i=k; i < GRID_SIZE-k; i++) {
            // Left side:
            aux = grid[i][jj];
            setSubCellH (&aux, kk, CONDUCTOR);
            grid[i][jj] = aux;
            // Right side:
            aux = grid[i][ROW_SIZE - jj - 0];
            setSubCellH (&aux, ELEMENTS_PER_CELL-1, CONDUCTOR);
            grid[i][ROW_SIZE - jj - 1] = aux;
        }
        // Draw up side and down side:
        for (j=jj; j < ROW_SIZE - jj; j++) {
            // Top side:
            aux = grid[k][j];
            for (h=0; h < ELEMENTS_PER_CELL; h++)
                setSubCellH (&aux, h, CONDUCTOR);
            grid[k][j] = aux;
            // Down side:
            aux = grid[GRID_SIZE-k-1][j];
            for (h=0; h < ELEMENTS_PER_CELL; h++)
                setSubCellH (&aux, h, CONDUCTOR);
            grid[GRID_SIZE-k-1][j] = aux;
        }
    }


    // Set a conductor cross in the grid:
    jj = GRID_SIZE/2;
    k = jj % ELEMENTS_PER_CELL;
    jj = jj / ELEMENTS_PER_CELL;
    // Up-Down
	for (i = 0; i<GRID_SIZE; i++) {
        aux = grid[i][jj];
        setSubCellH (&aux, k, CONDUCTOR);
        grid[i][jj] = aux;
    }
    // Left-Right
    ii = GRID_SIZE/2;
	for (j = 0; j<ROW_SIZE; j++) {
        aux = grid[ii][j];
        for (k=0; k < ELEMENTS_PER_CELL; k++) {
            setSubCellH (&aux, k, CONDUCTOR);
        }
        grid[ii][j] = aux;
    }


    // Set radom initial heads:
	for (i = 0; i<GRID_SIZE; i++) {
		for (j = 0; j<ROW_SIZE; j++) {
			aux = grid[i][j+0];
            dirty = 0;
			for (k=0; k<ELEMENTS_PER_CELL; k++) {
                if ( (getSubCellH (aux, k) == CONDUCTOR) && (r4_uniform_01(&seed) < P_HEAD) ) {              
                    setSubCellH (&aux, k, HEAD);
                    heads++;
                    dirty = 1;
                }
			}
            if (dirty)
			    grid[i][j+0] = aux;
		}
	}

    return heads;
}



int main(int argc, char* argv[])
{
    long unsigned int total_heads = 0;
    int iter;
    CELL_TYPE* h_grid; //Grid on host
    CELL_TYPE* d_grid; //Grid on device
    CELL_TYPE* d_newGrid; //Second grid used on device only
    CELL_TYPE* d_tmpGrid; //tmp grid pointer used to switch between grid and newGrid

    SUBCELL_TYPE* h_lookup_table; // Look-up table
    SUBCELL_TYPE* d_lookup_table; 
  
    size_t bytes = sizeof(CELL_TYPE)*(GRID_SIZE+0)*(ROW_SIZE+0);//2 added for periodic boundary condition ghost cells
    // Allocate host Grid used for initial setup and read back from device
    h_grid = (CELL_TYPE*)malloc(bytes);
 
    #ifdef __VIDEO
        char image_name[80];
        uint8_t *rgb= (uint8_t *)malloc (3 * sizeof(uint8_t) *GRID_SIZE*GRID_SIZE);	
    #endif

    // Allocate device grids
    cudaMalloc(&d_grid, bytes);
    cudaMalloc(&d_newGrid, bytes); 

    // Alocate look-up tables:
    size_t bytes2 = sizeof(SUBCELL_TYPE)*(HEAD+1)*(LOOKUP_TABLE_LIMIT);
    h_lookup_table =  (SUBCELL_TYPE *)malloc(bytes2);
    cudaMalloc(&d_lookup_table, bytes2);

    dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE,1);
    int linGridx = (int)ceil((ROW_SIZE+0)/(float)BLOCK_SIZE);
    int linGridy = (int)ceil(GRID_SIZE/(float)BLOCK_SIZE);
    dim3 gridSize(linGridx,linGridy,1);
 
    //dim3 cpyBlockSize(BLOCK_SIZE,1,1);
    //dim3 cpyGridRowsGridSize((int)ceil(ROW_SIZE/(float)cpyBlockSize.x),1,1);
    //dim3 cpyGridColsGridSize((int)ceil((GRID_SIZE+0)/(float)cpyBlockSize.x),1,1);
 
    // Init look-up table:
    init_lookup_table (h_lookup_table);

    total_heads = init_wire_world (h_grid);

    printf("Initial Heads: %ld \n", total_heads);

    // Copy over initial grid (Dim-1 threads)
    cudaMemcpy(d_grid, h_grid, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_newGrid, h_grid, bytes, cudaMemcpyHostToDevice);

    // Copy lookup-table
    cudaMemcpy(d_lookup_table, h_lookup_table, bytes2, cudaMemcpyHostToDevice);

    // Main loop
    for (iter = 0; iter<MAXITER; iter++) {
        wire_world <<<gridSize, blockSize>>> (d_grid, d_newGrid, d_lookup_table); 
        cudaDeviceSynchronize();
    #ifdef __VIDEO
		if (iter==0 || iter==50 || iter==100 || iter==600 || iter==300 || iter==1000 || iter==1023) {
		//if (iter==0 || iter==1 || iter==2 || iter==3 || iter==4 || iter==5 || iter == 600) {
            cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
			rgb = generate_rgb (GRID_SIZE, GRID_SIZE, h_grid, rgb);
			sprintf (image_name, "%s%d.ppm", "ca_step_", iter);
			write_CA_screenshot_image (image_name, GRID_SIZE, GRID_SIZE, rgb);
		}
    #endif

        // Swap our grids and iterate again
        d_tmpGrid = d_grid;
        d_grid = d_newGrid;
        d_newGrid = d_tmpGrid;
    }//iter loop
 
    // Copy back results and sum
    cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
 
    total_heads = count_total_heads (h_grid);
    printf("Final heads: %ld \n", total_heads);
 
    // Release memory
    cudaFree(d_grid);
    cudaFree(d_newGrid);
    free (h_lookup_table);
    cudaFree(d_lookup_table);
    free(h_grid);
 
    return 0;
}



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




__forceinline__ SUBCELL_TYPE getSubCellH (CELL_TYPE cell, char pos)
{
	return (cell >> (ELEMENTS_PER_CELL - 1 - pos)*8);
}

__forceinline__ void setSubCellH (CELL_TYPE *cell, char pos, SUBCELL_TYPE subcell)
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

__device__ void setSubCellD (CELL_TYPE *cell, char pos, SUBCELL_TYPE subcell)
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



__global__ void wire_world (CELL_TYPE *grid, CELL_TYPE *newGrid, SUBCELL_TYPE *lookup_table_)
{
    // We want id ∈ [1,SIZE]
    int iy = blockDim.y * blockIdx.y + threadIdx.y + 0;
    int ix = blockDim.x * blockIdx.x + threadIdx.x + 0;
    int id = iy * (ROW_SIZE+0) + ix;
    CELL_TYPE cell, new_cell=0; 
    CELL_TYPE up_cell, down_cell, right_cell, left_cell;                // Up,down,right,left cells
    CELL_TYPE upleft_cell, downleft_cell, upright_cell, downright_cell; // Diagonal cells
    SUBCELL_TYPE subcell;
 
    int k, numNeighbors;
    SUBCELL_TYPE (*lookup_table)[LOOKUP_TABLE_LIMIT] = (SUBCELL_TYPE (*)[LOOKUP_TABLE_LIMIT]) lookup_table_;
 
    if (iy>0 && iy < GRID_SIZE-1 && ix < ROW_SIZE) {
    //if (iy>0 && iy < GRID_SIZE-1 && ix> 0 && ix < ROW_SIZE-1) {
         cell = grid[id];
         //new_cell = cell;

        // First (0) subcell:
        up_cell = grid[id-(ROW_SIZE+0)];
        down_cell = grid[id+(ROW_SIZE+0)];
        if (ix>0) {
            left_cell = grid[id-1];
            upleft_cell = grid[id-(ROW_SIZE)-1];
            downleft_cell = grid[id+(ROW_SIZE)-1];

            numNeighbors = getSubCellD (up_cell, 0) + getSubCellD (down_cell, 0);   // upper lower
            numNeighbors += getSubCellD (left_cell, ELEMENTS_PER_CELL-1);           // left
            numNeighbors += getSubCellD (cell, 1);                                  // right
            numNeighbors += getSubCellD (upleft_cell, ELEMENTS_PER_CELL-1) + getSubCellD (downleft_cell, ELEMENTS_PER_CELL-1); // diagonals left
            numNeighbors += getSubCellD (up_cell, 1) + getSubCellD (down_cell, 1); // diagonals right
            subcell = getSubCellD (cell, 0);
            setSubCellD (&new_cell, 0, lookup_table[subcell][numNeighbors]);
        } else {
            // Not process the grid border
            setSubCellD (&new_cell, 0, getSubCellD(cell, 0));
        }


        // Middle subcells:
        for (k=1; k<CELL_NEIGHBOURS-1; k++) {
            numNeighbors = getSubCellD (up_cell, k) + getSubCellD (down_cell, k); // upper lower
            numNeighbors += getSubCellD (cell, k-1) + getSubCellD (cell, k+1); // left right
            numNeighbors += getSubCellD (up_cell, k-1) + getSubCellD (down_cell, k-1); // diagonals left
            numNeighbors += getSubCellD (up_cell, k+1) + getSubCellD (down_cell, k+1); // diagonals right
            subcell = getSubCellD (cell, k);
            setSubCellD (&new_cell, k, lookup_table[subcell][numNeighbors]);
        }

        // Last (CELL_NEIGHBOURS-1) subcell:
        if (ix < ROW_SIZE-1) {
            right_cell = grid[id+1];
            upright_cell = grid[id-(ROW_SIZE)+1];
            downright_cell = grid[id+(ROW_SIZE)+1];

            numNeighbors = getSubCellD (up_cell, ELEMENTS_PER_CELL-1) + getSubCellD (down_cell, ELEMENTS_PER_CELL-1); // upper lower
            numNeighbors += getSubCellD (cell, ELEMENTS_PER_CELL-2) + getSubCellD (right_cell, 0); // left right
            numNeighbors += getSubCellD (up_cell, ELEMENTS_PER_CELL-2) + getSubCellD (down_cell, ELEMENTS_PER_CELL-2); // diagonals left
            numNeighbors += getSubCellD (upright_cell, 0) + getSubCellD (downright_cell, 0); // diagonals right
            subcell = getSubCellD (cell, ELEMENTS_PER_CELL-1);
            setSubCellD (&new_cell, ELEMENTS_PER_CELL-1, lookup_table[subcell][numNeighbors]);
        } else {
            // Not process the grid border
            setSubCellD (&new_cell, ELEMENTS_PER_CELL-1, getSubCellD(cell, ELEMENTS_PER_CELL-1));
        }

        // Copy new_cell to newGrid:
        newGrid[id] = new_cell;
    }
}


long int count_total_heads (CELL_TYPE *h_grid) {
	int i,j,k;
	CELL_TYPE aux;
	long int total = 0;
	
	for (i = 1; i<GRID_SIZE-1; i++) {
		for (j = 0; j<ROW_SIZE; j++) {
			aux = h_grid[i*(ROW_SIZE+0)+j];
			for (k=0; k<ELEMENTS_PER_CELL; k++) {
                if (getSubCellH (aux, k) == HEAD)
				    total++;
			}
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
	int x, y, jj, k,cur;
    SUBCELL_TYPE *cell;
	
	// Transform "grid" (CA states grid) into a similar rgb image format.
	CELL_TYPE (*grid)[ROW_SIZE] = (CELL_TYPE (*)[ROW_SIZE]) grid_;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			cur = 3 * (y * width + x);
            jj = x/ELEMENTS_PER_CELL ;
            k = x % ELEMENTS_PER_CELL;
            cell = (SUBCELL_TYPE *)&grid[y][jj];
            //cell = cell + k;
            cell = cell + (ELEMENTS_PER_CELL - k -1);
			if (*cell == EMPTY) {
				// Empty: Black
				rgb[cur + 0] = 128;
				rgb[cur + 1] = 128;
				rgb[cur + 2] = 128;
			} else if (*cell == HEAD) {
				// Head: Blue
				rgb[cur + 0] = 0;
				rgb[cur + 1] = 0;
				rgb[cur + 2] = 255;
			} else if (*cell == TAIL) {
				// Tail: Red
				rgb[cur + 0] = 255;
				rgb[cur + 1] = 0;
				rgb[cur + 2] = 0;
			} else if (*cell == CONDUCTOR) {
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



// 256 Result in console: initial heads: 17, final heads: 2254
// 512 Result in console: initial heads: 42, final heads: 2612
// 1024 Result in console: initial heads: 138, final heads: 2421
// 2048 Result in console: initial heads: 538, final heads: 2722
// 4096 Result in console: initial heads: 2127, final heads: 3682
// 8192 Result in console: initial heads: 8439, final heads: 6719
// 16384 Result in console: initial heads: 33476, final heads: 15925

