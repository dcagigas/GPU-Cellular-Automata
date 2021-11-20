#include <stdio.h>
#include <stdlib.h>
 
#define GRID_SIZE 256                       // Grid dimension: SIZExSIZE
#define MAXITER 1024                        // Number of game steps

typedef unsigned long int CELL_TYPE;
typedef unsigned char SUBCELL_TYPE;
typedef unsigned char LUT_TYPE;
#define ELEMENTS_PER_CELL 8                 // "int_64" data type per cell
#define ROW_SIZE GRID_SIZE/ELEMENTS_PER_CELL    // Real grid dimension

#define BLOCK_SIZE 32

//#define CELL_NEIGHBOURS 4
// Clyclic Automaton states:
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

// Uncomment to write images (CA screenshots) to disk.
//#define __VIDEO


int seed = 1; 	// Seed for random uniform number output.

void setup_lookup_table (LUT_TYPE *lut_);
float r4_uniform_01 ( int *seed );
__forceinline__ SUBCELL_TYPE getSubCellH (CELL_TYPE cell, char pos);
__forceinline__ void setSubCellH (CELL_TYPE *cell, char pos, SUBCELL_TYPE subcell);
__device__ SUBCELL_TYPE getSubCellD (CELL_TYPE cell, char pos);
__device__ void setSubCellD (CELL_TYPE *cell, char pos, SUBCELL_TYPE subcell);
__global__ void ghostRows(CELL_TYPE *grid);
__global__ void ghostCols(CELL_TYPE *grid);
__global__ void cyclic_ca (CELL_TYPE *grid, CELL_TYPE *newGrid, LUT_TYPE *GPU_lookup_table);

void init_random_values (CELL_TYPE *h_grid);
char* generate_rgb (int width, int height, CELL_TYPE *grid_, char *rgb);
void write_CA_screenshot_image (const char *filename, int width, int height, char *img);



int main(int argc, char* argv[])
{
    //long int total = 0;
    int iter;
    CELL_TYPE* h_grid; //Grid on host
    CELL_TYPE* d_grid; //Grid on device
    CELL_TYPE* d_newGrid; //Second grid used on device only
    CELL_TYPE* d_tmpGrid; //tmp grid pointer used to switch between grid and newGrid
    LUT_TYPE* h_lut; // Lookup_table on host
    LUT_TYPE* d_lut; // Lookup_table on device
  
    size_t bytes = sizeof(CELL_TYPE)*(GRID_SIZE+2)*(ROW_SIZE+2);//2 added for periodic boundary condition ghost cells
    size_t bytes_lut = sizeof(LUT_TYPE)*N*2;
    // Allocate host Grid used for initial setup and read back from device
    h_grid = (CELL_TYPE*)malloc(bytes);
    h_lut = (LUT_TYPE *)malloc(bytes_lut);
 
    #ifdef __VIDEO
        char image_name[80];
        char *rgb= (char *)malloc (3 * sizeof(char) *(GRID_SIZE+2)*(GRID_SIZE+2));	
    #endif

    // Allocate device grids
    cudaMalloc(&d_grid, bytes);
    cudaMalloc(&d_newGrid, bytes);
    cudaMalloc(&d_lut, bytes_lut);
 
    // Assign initial population randomly
    init_random_values (h_grid);

    // Init lookup_table
    setup_lookup_table (h_lut);
    cudaMemcpy(d_lut, h_lut, bytes_lut, cudaMemcpyHostToDevice);

    // Copy over initial game grid (Dim-1 threads)
    cudaMemcpy(d_grid, h_grid, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_newGrid, h_grid, bytes, cudaMemcpyHostToDevice);
 
    dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE,1);
    int linGridx = (int)ceil((ROW_SIZE+2)/(float)BLOCK_SIZE);
    int linGridy = (int)ceil((GRID_SIZE+2)/(float)BLOCK_SIZE);
    dim3 gridSize(linGridx,linGridy,1);
 
    dim3 cpyBlockSize(BLOCK_SIZE,1,1);
    dim3 cpyGridRowsGridSize((int)ceil(ROW_SIZE/(float)cpyBlockSize.x),1,1);
    dim3 cpyGridColsGridSize((int)ceil((GRID_SIZE+2)/(float)cpyBlockSize.x),1,1);
 
    // Main game loop
    for (iter = 0; iter<MAXITER; iter++) {
    #ifdef __VIDEO
		if (iter==0 || iter==100 || iter==200 || iter==300 || iter==400 || iter==1000) {
		//if (iter==0 || iter==1 || iter==2) {
            cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
			rgb = generate_rgb (GRID_SIZE+2, GRID_SIZE+2, h_grid, rgb);
			sprintf (image_name, "%s%d.ppm", "ca_step_", iter);
			write_CA_screenshot_image (image_name, GRID_SIZE+2, GRID_SIZE+2, rgb);
		}
    #endif
        ghostRows<<<cpyGridRowsGridSize, cpyBlockSize>>>(d_grid);
        ghostCols<<<cpyGridColsGridSize, cpyBlockSize>>>(d_grid);
        cyclic_ca<<<gridSize, blockSize>>>(d_grid, d_newGrid, d_lut); 

        // Swap our grids and iterate again
        d_tmpGrid = d_grid;
        d_grid = d_newGrid;
        d_newGrid = d_tmpGrid;
    }//iter loop
 
    // Copy back results and sum
    cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
 
    // Release memory
    cudaFree(d_grid);
    cudaFree(d_newGrid);
    cudaFree(d_lut);
    free(h_grid);
    free(h_lut);
 
    return 0;
}



void init_random_values (CELL_TYPE *h_grid) {
	int i,j,k;
    float p;
    // Assign an initial randomly
    for(i = 0; i<GRID_SIZE+2; i++) {
        for(j = 0; j<ROW_SIZE+2; j++) {
            for(k = 0; k<ELEMENTS_PER_CELL; k++) {
                p = r4_uniform_01(&seed);
                if (p < P_STATE) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE1);
                } else if (p < P_STATE*2) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE2);
                } else if (p <  P_STATE*3) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE3);
                } else if (p <  P_STATE*4) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE4);
                } else if (p <  P_STATE*5) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE5);
                } else if (p <  P_STATE*6) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE6);
                } else if (p <  P_STATE*7) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE7);
                } else if (p <  P_STATE*8) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE8);
                } else if (p <  P_STATE*9) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE9);
                } else if (p <  P_STATE*10) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE10);
                } else if (p <  P_STATE*11) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE11);
                } else if (p <  P_STATE*12) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE12);
                } else if (p <  P_STATE*13) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE13);
                } else if (p <  P_STATE*14) {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE14);
                } else {
                    setSubCellH (&h_grid[i*(ROW_SIZE+2)+j], k, STATE15);
                }
            }
        }
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

__device__ SUBCELL_TYPE getSubCellD (CELL_TYPE cell, char pos)
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


__global__ void cyclic_ca (CELL_TYPE *grid, CELL_TYPE *newGrid, LUT_TYPE *lookup_table_)
{
    // We want id ∈ [1,SIZE]
    int iy = blockDim.y * blockIdx.y + threadIdx.y;
    int ix = blockDim.x * blockIdx.x + threadIdx.x;
    int id = iy * (ROW_SIZE+2) + ix;
    CELL_TYPE cell, new_cell=0; 
    CELL_TYPE up_cell, down_cell, right_cell, left_cell;                // Up,down,right,left cells
    SUBCELL_TYPE subcell;
 
    int k, transition;
    LUT_TYPE (*lookup_table)[2] = (LUT_TYPE (*)[2]) lookup_table_;
 
    if (iy>0 && iy <= GRID_SIZE && ix <= ROW_SIZE) {
         cell = grid[id];

        // First (0) subcell:
        subcell = getSubCellD (cell, 0);
        up_cell = grid[id-(ROW_SIZE+2)];
        down_cell = grid[id+(ROW_SIZE+2)];
        left_cell = grid[id-1];

        transition = (getSubCellD (up_cell, 0) == (subcell+1)%N) || (getSubCellD (down_cell, 0) == (subcell+1)%N); // upper lower
        transition = transition || (getSubCellD (left_cell, ELEMENTS_PER_CELL-1) == (subcell+1)%N) || (getSubCellD (cell, 1) == (subcell+1)%N); // left right

        setSubCellD (&new_cell, 0, lookup_table[subcell][transition]);
        // Middle subcells:
        for (k=1; k<ELEMENTS_PER_CELL-1; k++) {
            subcell = getSubCellD (cell, k);
            transition = (getSubCellD (up_cell, k) == (subcell+1)%N) || (getSubCellD (down_cell, k) == (subcell+1)%N); // upper lower
            transition = transition || (getSubCellD (cell, k-1) == (subcell+1)%N) || (getSubCellD (cell, k+1) == (subcell+1)%N); // left right
            setSubCellD (&new_cell, k, lookup_table[subcell][transition]);
        }

        // Last (CELL_NEIGHBOURS-1) subcell:
        right_cell = grid[id+1];

        subcell = getSubCellD (cell, ELEMENTS_PER_CELL-1);
        transition = (getSubCellD (up_cell, ELEMENTS_PER_CELL-1) == (subcell+1)%N) || (getSubCellD (down_cell, ELEMENTS_PER_CELL-1) == (subcell+1)%N); // upper lower
        transition = transition || (getSubCellD (cell, ELEMENTS_PER_CELL-2) == (subcell+1)%N) || (getSubCellD (right_cell, 0) == (subcell+1)%N); // left right
        setSubCellD (&new_cell, ELEMENTS_PER_CELL-1, lookup_table[subcell][transition]);
        // Copy new_cell to newGrid:
        newGrid[id] = new_cell;
    }
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
char* generate_rgb (int width, int height, CELL_TYPE *grid, char *rgb) {
	int x, xx, y, k, cur;
    SUBCELL_TYPE subcell;
    CELL_TYPE cell; 
	
	// Transform "grid" (CA states grid) into a similar rgb image format.
	//CELL_TYPE (*grid)[height] = (CELL_TYPE (*)[height]) grid_;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
            xx = x / ELEMENTS_PER_CELL;
            cell = grid[y*(ROW_SIZE+2)+xx];
            k = x % ELEMENTS_PER_CELL;
            subcell = getSubCellH (cell, k);
			cur = 3 * (y * width + x);
			    if (subcell == STATE1) {
				    // Purple light
				    rgb[cur + 0] = (char) 177;
				    rgb[cur + 1] = (char) 25;
				    rgb[cur + 2] = (char) 251;
			    } else if (subcell == STATE2) {
				    // Purple light
				    rgb[cur + 0] = (char) 177;
				    rgb[cur + 1] = (char) 25;
				    rgb[cur + 2] = (char) 251;
			    } else if (subcell == STATE3) {
				    // Purple light
				    rgb[cur + 0] = (char) 165;
				    rgb[cur + 1] = (char) 30;
				    rgb[cur + 2] = (char) 233;
			    } else if (subcell == STATE4) {
				    // Purple light
				    rgb[cur + 0] = (char) 153;
				    rgb[cur + 1] = (char) 38;
				    rgb[cur + 2] = (char) 215;
			    } else if (subcell == STATE5) {
				    // Purple dark
				    rgb[cur + 0] = (char) 139;
				    rgb[cur + 1] = (char) 51;
				    rgb[cur + 2] = (char) 196;
			    } else if (subcell == STATE6) {
				    // Purple dark
				    rgb[cur + 0] = (char) 127;
				    rgb[cur + 1] = (char) 64;
				    rgb[cur + 2] = (char) 180;
			    } else if (subcell == STATE7) {
				    // Purple dark
				    rgb[cur + 0] = (char) 115;
				    rgb[cur + 1] = (char) 78;
				    rgb[cur + 2] = (char) 162;
			    } else if (subcell == STATE8) {
				    // Green dark
				    rgb[cur + 0] = (char) 103;
				    rgb[cur + 1] = (char) 93;
				    rgb[cur + 2] = (char) 144;
			    } else if (subcell == STATE9) {
				    // Green dark
				    rgb[cur + 0] = (char) 91;
				    rgb[cur + 1] = (char) 108;
				    rgb[cur + 2] = (char) 127;
			    } else if (subcell == STATE10) {
				    // Green dark
				    rgb[cur + 0] = (char) 79;
				    rgb[cur + 1] = (char) 124;
				    rgb[cur + 2] = (char) 109;
			    } else if (subcell == STATE11) {
				    // Green dark
				    rgb[cur + 0] = (char) 68;
				    rgb[cur + 1] = (char) 138;
				    rgb[cur + 2] = (char) 92;
			    } else if (subcell == STATE12) {
				    // Green light
				    rgb[cur + 0] = (char) 58;
				    rgb[cur + 1] = (char) 153;
				    rgb[cur + 2] = (char) 76;
			    } else if (subcell == STATE13) {
				    // Green light
				    rgb[cur + 0] = (char) 50;
				    rgb[cur + 1] = (char) 169;
				    rgb[cur + 2] = (char) 63;
			    } else if (subcell == STATE14) {
				    // Green light
				    rgb[cur + 0] = (char) 45;
				    rgb[cur + 1] = (char) 185;
				    rgb[cur + 2] = (char) 50;
			    } else if (subcell == STATE15) {
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
			//}
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
	fprintf(fp, "# Created by %s\n","Cyclic_CA_int64.cu");
	
	//image size
	fprintf(fp, "%d %d\n", width, height);
	
	// rgb component depth
	fprintf(fp, "%d\n",255);
	
	// pixel data
	fwrite(img, 3 * width, height, fp);
	fclose(fp);
}



