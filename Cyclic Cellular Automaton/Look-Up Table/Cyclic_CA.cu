#include <stdio.h>
#include <stdlib.h>
//#include <openssl/rand.h> 			// For random numbers. OpenSSL library must be installed.

#define GRID_SIZE 256
#define TIMESTEP 1024

#define SEED 1
#define BLOCK_SIZE 32
//#define CELL_TYPE char
typedef  unsigned char CELL_TYPE;
typedef  unsigned char LUT_TYPE;

#define STATE1 0
#define STATE2 1
#define STATE3 2
#define STATE4 3
#define STATE5 4
#define STATE6 5
#define STATE7 6
#define STATE8 7
#define STATE9 8
#define STATE10 9
#define STATE11 10
#define STATE12 11
#define STATE13 12
#define STATE14 13
#define STATE15 14

#define N 15 // Number os states

#define P_STATE (1.0/N)



//#define __VIDEO

float r4_uniform_01 ( int *seed );
char* generate_rgb (int width, int height, CELL_TYPE *grid_, char *rgb);
void write_CA_screenshot_image (const char *filename, int width, int height, char *img);


int seed = SEED; 	// Seed for random uniform number output.

/*
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
*/

__global__ void cyclic_ca (int dim, CELL_TYPE *grid, CELL_TYPE *newGrid, LUT_TYPE *lut_)
{
    /*
    This auto-reproductive cell automaton is two-dimensional, and its cells can take four states (red, brown, purple, purple). 
    The state of a cell at time t + 1 depends on its state at time t and the state of its 4 neighbors (Von Neumann's neighborhood). 
    A cell moves from a state i to a state i + 1 in the state cycle when the state i + 1 is present in at least two neighboring cells
    */
    // We want id ∈ [1,dim]
    int iy = blockDim.y * blockIdx.y + threadIdx.y;
    int ix = blockDim.x * blockIdx.x + threadIdx.x;
    int id = iy * (dim+2) + ix;
 
    LUT_TYPE (*lut)[2] = (LUT_TYPE (*)[2]) lut_;
 
    if (iy>0 && iy < dim+1 && ix>0 && ix < dim+1) {
        char cell = grid[id];

        int transition = (grid[id+(dim+2)] == (cell+1)%N || grid[id-(dim+2)] == (cell+1)%N || grid[id+1] == (cell+1)%N || grid[id-1] == (cell+1)%N);
        newGrid[id] = lut[cell][transition];
    }
}


void setup_lookup_table (LUT_TYPE *lut_) {
    int i;
    LUT_TYPE (*lut)[2] = (LUT_TYPE (*)[2]) lut_;

    for (i=0; i<N; i++) {
        lut[i][0] = i;
        lut[i][1] = (i+1)%N;
    }

}

 
int main(int argc, char* argv[])
{
    int i,j,iter;
    float p;
    long unsigned int total_state1 = 0, total_state2 = 0, total_state3 = 0, total_state4 = 0, 
                      total_state5 = 0, total_state6 = 0, total_state7 = 0, total_state8 = 0;
    CELL_TYPE* h_grid; //Grid on host
    CELL_TYPE* d_grid; //Grid on device
    LUT_TYPE* h_lut; // Lookup_table on host
    LUT_TYPE* d_lut; // Lookup_table on device
    CELL_TYPE* d_newGrid; //Second grid used on device only
    CELL_TYPE* d_tmpGrid; //tmp grid pointer used to switch between grid and newGrid
 
    int dim = GRID_SIZE; //Linear dimension of our grid - not counting ghost cells
    int maxIter = TIMESTEP; //Number of game steps
 
    size_t bytes = sizeof(CELL_TYPE)*(dim+2)*(dim+2);//2 added for periodic boundary condition ghost cells
    size_t bytes_lut = sizeof(LUT_TYPE)*N*2;
    // Allocate host Grid used for initial setup and read back from device
    h_grid = (CELL_TYPE *)malloc(bytes);
    h_lut = (LUT_TYPE *)malloc(bytes_lut);

    #ifdef __VIDEO
        char image_name[80];
        char *rgb= (char *)malloc (3 * sizeof(char) *(GRID_SIZE+2)*(GRID_SIZE+2));	
    #endif
 
    // Allocate device grids
    cudaMalloc(&d_grid, bytes);
    cudaMalloc(&d_newGrid, bytes);
    cudaMalloc(&d_lut, bytes_lut);
 

    // Assign an initial randomly
    for(i = 0; i<dim+2; i++) {
        for(j = 0; j<dim+2; j++) {
            p = r4_uniform_01(&seed);
            if (p < P_STATE) {
                h_grid[i*(dim+2)+j] = STATE1;
                total_state1++;
            } else if (p < P_STATE*2) {
                h_grid[i*(dim+2)+j] = STATE2;
                total_state2++;
            } else if (p <  P_STATE*3) {
                h_grid[i*(dim+2)+j] = STATE3;
                total_state3++;
            } else if (p <  P_STATE*4) {
                h_grid[i*(dim+2)+j] = STATE4;
                total_state4++;
            } else if (p <  P_STATE*5) {
                h_grid[i*(dim+2)+j] = STATE5;
                total_state5++;
            } else if (p <  P_STATE*6) {
                h_grid[i*(dim+2)+j] = STATE6;
                total_state6++;
            } else if (p <  P_STATE*7) {
                h_grid[i*(dim+2)+j] = STATE7;
                total_state7++;
            } else if (p <  P_STATE*8) {
                h_grid[i*(dim+2)+j] = STATE8;
                total_state8++;
            } else if (p <  P_STATE*9) {
                h_grid[i*(dim+2)+j] = STATE9;
            } else if (p <  P_STATE*10) {
                h_grid[i*(dim+2)+j] = STATE10;
            } else if (p <  P_STATE*11) {
                h_grid[i*(dim+2)+j] = STATE11;
            } else if (p <  P_STATE*12) {
                h_grid[i*(dim+2)+j] = STATE12;
            } else if (p <  P_STATE*13) {
                h_grid[i*(dim+2)+j] = STATE13;
            } else if (p <  P_STATE*14) {
                h_grid[i*(dim+2)+j] = STATE14;
            } else {
                h_grid[i*(dim+2)+j] = STATE15;
            }
        }
    }


    //printf("Initial state. State1: %ld State2: %ld State3: %ld State4: %ld State5: %ld State6: %ld State7: %ld State8: %ld\n", total_state1, total_state2, total_state3, total_state4, total_state5, total_state6, total_state7, total_state8);

    // Setup LookUp table:
    setup_lookup_table (h_lut);
    cudaMemcpy(d_lut, h_lut, bytes_lut, cudaMemcpyHostToDevice);
 
    // Copy over initial game grid (Dim-1 threads)
    cudaMemcpy(d_grid, h_grid, bytes, cudaMemcpyHostToDevice);
 
    dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE,1);
    int linGrid = (int)ceil(dim/(float)BLOCK_SIZE);
    dim3 gridSize(linGrid,linGrid,1);
 
    //dim3 cpyBlockSize(BLOCK_SIZE,1,1);
    //dim3 cpyGridRowsGridSize((int)ceil(dim/(float)cpyBlockSize.x),1,1);
    //dim3 cpyGridColsGridSize((int)ceil((dim+2)/(float)cpyBlockSize.x),1,1);
 
    // Main game loop
    for (iter = 0; iter<maxIter; iter++) {
    #ifdef __VIDEO
		if (iter==0 || iter==100 || iter==200 || iter==300 || iter==400 || iter==1000) {
            cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
			rgb = generate_rgb (GRID_SIZE+2, GRID_SIZE+2, h_grid, rgb);
			sprintf (image_name, "%s%d.ppm", "ca_step_", iter);
			write_CA_screenshot_image (image_name, GRID_SIZE+2, GRID_SIZE+2, rgb);
		}
    #endif
        //ghostRows<<<cpyGridRowsGridSize, cpyBlockSize>>>(dim, d_grid);
        //ghostCols<<<cpyGridColsGridSize, cpyBlockSize>>>(dim, d_grid);
        cyclic_ca<<<gridSize, blockSize>>>(dim, d_grid, d_newGrid, d_lut);

        // Swap our grids and iterate again
        d_tmpGrid = d_grid;
        d_grid = d_newGrid;
        d_newGrid = d_tmpGrid;
    }//iter loop
 
    // Copy back results and sum
    cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
 
    total_state1 = 0; 
    total_state2 = 0; 
    total_state3 = 0; 
    total_state4 = 0; 
    total_state5 = 0; 
    total_state6 = 0; 
    total_state7 = 0; 
    total_state8 = 0;    
    // Sum up alive cells and print results
    for (i = 0; i<dim+2; i++) {
        for (j = 0; j<dim+2; j++) {
            if (h_grid[i*(dim+2)+j] == STATE1)
                total_state1++;
            else if (h_grid[i*(dim+2)+j] == STATE2)
                total_state2++;
            else if (h_grid[i*(dim+2)+j] == STATE3)
                total_state3++;
            else if (h_grid[i*(dim+2)+j] == STATE4)
                total_state4++;
            else if (h_grid[i*(dim+2)+j] == STATE5)
                total_state5++;
            else if (h_grid[i*(dim+2)+j] == STATE6)
                total_state6++;
            else if (h_grid[i*(dim+2)+j] == STATE7)
                total_state7++;
            else if (h_grid[i*(dim+2)+j] == STATE8)
                total_state8++;
        }
    }
    printf("Final. State1: %ld State2: %ld State3: %ld State4: %ld State5: %ld State6: %ld State7: %ld State8: %ld\n", total_state1, total_state2, total_state3, total_state4, total_state5, total_state6, total_state7, total_state8);
 
    // Release memory
    cudaFree(d_grid);
    cudaFree(d_newGrid);
    cudaFree(d_lut);
    free(h_grid);
    free(h_lut);
 
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
char* generate_rgb (int width, int height, CELL_TYPE *grid_, char *rgb) {
	int x, y, cur;
	
	// Transform "grid" (CA states grid) into a similar rgb image format.
	CELL_TYPE (*grid)[height] = (CELL_TYPE (*)[height]) grid_;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			cur = 3 * (y * width + x);
			if (grid[y][x] == STATE1) {
				// Purple light
				rgb[cur + 0] = (char) 177;
				rgb[cur + 1] = (char) 25;
				rgb[cur + 2] = (char) 251;
			} else if (grid[y][x] == STATE2) {
				// Purple light
				rgb[cur + 0] = (char) 177;
				rgb[cur + 1] = (char) 25;
				rgb[cur + 2] = (char) 251;
			} else if (grid[y][x] == STATE3) {
				// Purple light
				rgb[cur + 0] = (char) 165;
				rgb[cur + 1] = (char) 30;
				rgb[cur + 2] = (char) 233;
			} else if (grid[y][x] == STATE4) {
				// Purple light
				rgb[cur + 0] = (char) 153;
				rgb[cur + 1] = (char) 38;
				rgb[cur + 2] = (char) 215;
			} else if (grid[y][x] == STATE5) {
				// Purple dark
				rgb[cur + 0] = (char) 139;
				rgb[cur + 1] = (char) 51;
				rgb[cur + 2] = (char) 196;
			} else if (grid[y][x] == STATE6) {
				// Purple dark
				rgb[cur + 0] = (char) 127;
				rgb[cur + 1] = (char) 64;
				rgb[cur + 2] = (char) 180;
			} else if (grid[y][x] == STATE7) {
				// Purple dark
				rgb[cur + 0] = (char) 115;
				rgb[cur + 1] = (char) 78;
				rgb[cur + 2] = (char) 162;
			} else if (grid[y][x] == STATE8) {
				// Green dark
				rgb[cur + 0] = (char) 103;
				rgb[cur + 1] = (char) 93;
				rgb[cur + 2] = (char) 144;
			} else if (grid[y][x] == STATE9) {
				// Green dark
				rgb[cur + 0] = (char) 91;
				rgb[cur + 1] = (char) 108;
				rgb[cur + 2] = (char) 127;
			} else if (grid[y][x] == STATE10) {
				// Green dark
				rgb[cur + 0] = (char) 79;
				rgb[cur + 1] = (char) 124;
				rgb[cur + 2] = (char) 109;
			} else if (grid[y][x] == STATE11) {
				// Green dark
				rgb[cur + 0] = (char) 68;
				rgb[cur + 1] = (char) 138;
				rgb[cur + 2] = (char) 92;
			} else if (grid[y][x] == STATE12) {
				// Green light
				rgb[cur + 0] = (char) 58;
				rgb[cur + 1] = (char) 153;
				rgb[cur + 2] = (char) 76;
			} else if (grid[y][x] == STATE13) {
				// Green light
				rgb[cur + 0] = (char) 50;
				rgb[cur + 1] = (char) 169;
				rgb[cur + 2] = (char) 63;
			} else if (grid[y][x] == STATE14) {
				// Green light
				rgb[cur + 0] = (char) 45;
				rgb[cur + 1] = (char) 185;
				rgb[cur + 2] = (char) 50;
			} else if (grid[y][x] == STATE15) {
				// Green light
				rgb[cur + 0] = (char) 42;
				rgb[cur + 1] = (char) 200;
				rgb[cur + 2] = (char) 42;
			} else {
				// Any other state: White
				rgb[cur + 0] = (char) 255;
				rgb[cur + 1] = (char) 255;
				rgb[cur + 2] = (char) 255;
			}
		}
	}
	return rgb;
}


// Writes to disk a CA state screenshot.
void write_CA_screenshot_image (const char *filename, int width, int height, char *img)
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
	fprintf(fp, "# Created by %s\n","Cyclic_CA.cu");
	
	//image size
	fprintf(fp, "%d %d\n", width, height);
	
	// rgb component depth
	fprintf(fp, "%d\n",255);
	
	// pixel data
	fwrite(img, 3 * width, height, fp);
	fclose(fp);
}



// 256 Result in console: 
// 512 Result in console: 
// 1024 Result in console: 
// 2048 Result in console: 
// 4096 Result in console: 
// 8192 Result in console: 
// 16384 Result in console: 

