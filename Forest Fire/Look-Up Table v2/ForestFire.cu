#include <stdio.h>
#include <stdlib.h>
#include <openssl/rand.h> 			// For random numbers. OpenSSL library must be installed.

#define GRID_SIZE 256
#define TIMESTEP 1024

#define SEED 1
#define BLOCK_SIZE 32
#define CELL_TYPE char

#define EMPTY 0
#define TREE 1
#define ASH 2
#define FIRE 3

#define P_EMPTY 0.20
#define P_TREE 0.80
//#define P_FIRE (1.0/(GRID_SIZE*GRID_SIZE))
#define P_FIRE (1.5/((GRID_SIZE)*GRID_SIZE))

#define LOOKUP_TABLE_LIMIT 2

//#define __VIDEO

void init_lookup_table (CELL_TYPE *lookup_table);
float r4_uniform_01 ( int *seed );
uint8_t* generate_rgb (int width, int height, CELL_TYPE *grid_, uint8_t *rgb);
void write_CA_screenshot_image (const char *filename, int width, int height, uint8_t *img);


int seed = SEED; 	// Seed for random uniform number output.


__global__ void forest_fire (int dim, CELL_TYPE *grid, CELL_TYPE *newGrid, CELL_TYPE *lookup_table_)
{
    /*
    A tree box (1) with at least one neighbor on fire (3) becomes on fire, 
    a burning box (3) becomes ash (2), 
    an ash box (2) becomes empty (0) if it does not have a neighbor on fire (3), 
    and an empty box (0) remains empty (0). 
    */
    // We want id ∈ [1,dim]
    int iy = blockDim.y * blockIdx.y + threadIdx.y;
    int ix = blockDim.x * blockIdx.x + threadIdx.x;
    int id = iy * (dim+2) + ix;
 
    CELL_TYPE (*lookup_table)[LOOKUP_TABLE_LIMIT] = (CELL_TYPE(*)[LOOKUP_TABLE_LIMIT]) lookup_table_;
 
    if ( (iy>0 && iy < dim+1) && (ix>0 && ix < dim+1) ) {
 
        CELL_TYPE cell = grid[id];
        int there_is_fire;

        //Check if there is a FIRE cell:
        there_is_fire = (grid[id+(dim+2)]==FIRE || grid[id-(dim+2)]==FIRE || grid[id+1]==FIRE || grid[id-1]==FIRE);

        newGrid[id] = lookup_table[cell][there_is_fire];
    }
}


void init_lookup_table (CELL_TYPE *lookup_table_) {
    CELL_TYPE (*lookup_table)[LOOKUP_TABLE_LIMIT] = (CELL_TYPE(*)[LOOKUP_TABLE_LIMIT]) lookup_table_;

    // Init lookup_table:
    // 0: Empty. If a cell is EMPTY remains EMPTY 
    lookup_table[EMPTY][0] = EMPTY;
    lookup_table[EMPTY][1] = EMPTY;

    // 1: Tree. A TREE cell with at least one neighbor on fire becomes on FIRE 
    lookup_table[TREE][0] = TREE;
    lookup_table[TREE][1] = FIRE;

    // 2: Ash. An ASH cell becomes EMPTY if it does not have a neighbor on FIRE 
    lookup_table[ASH][0] = EMPTY;
    lookup_table[ASH][1] = ASH;

    // 3: Fire. A FIRE cell becomes ASH 
    lookup_table[FIRE][0] = ASH;
    lookup_table[FIRE][1] = ASH;
}


 
int main(int argc, char* argv[])
{
    int i,j,iter;
    float p;
    long unsigned int total_tree = 0, total_fire = 0;
    CELL_TYPE* h_grid; //Grid on host
    CELL_TYPE* d_grid; //Grid on device
    CELL_TYPE* d_newGrid; //Second grid used on device only
    CELL_TYPE* d_tmpGrid; //tmp grid pointer used to switch between grid and newGrid

    CELL_TYPE* h_lookup_table; // Look-up table
    CELL_TYPE* d_lookup_table; 

    int dim = GRID_SIZE; //Linear dimension of our grid - not counting ghost cells
    int maxIter = TIMESTEP; //Number of game steps
 
    size_t bytes = sizeof(CELL_TYPE)*(dim+2)*(dim+2);//2 added for periodic boundary condition ghost cells
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
    size_t bytes2 = sizeof(CELL_TYPE)*(FIRE+1)*(LOOKUP_TABLE_LIMIT);
    h_lookup_table =  (CELL_TYPE *)malloc(bytes2);
    cudaMalloc(&d_lookup_table, bytes2);

    // Init look-up table:
    init_lookup_table (h_lookup_table);

    // Assign an initial forest randomly
    for(i = 1; i<dim+1; i++) {
        for(j = 1; j<dim+1; j++) {
            p = r4_uniform_01(&seed);
            if (p < P_FIRE) {
                h_grid[i*(dim+2)+j] = FIRE;
                total_fire++;
            } else if (p < P_EMPTY) {
                h_grid[i*(dim+2)+j] = EMPTY;
            } else {
                h_grid[i*(dim+2)+j] = TREE;
                total_tree++;
            }
        }
    }
    // Borders are set to empty
    for(i = 0; i<dim+2; i++) {
        h_grid[i*(dim+2)] = EMPTY;
        h_grid[i*(dim+2)+(dim+1)] = EMPTY;
    }
    for(j = 0; j<dim+2; j++) {
        h_grid[j] = EMPTY;
        h_grid[(dim+2)*(dim+1)+j] = EMPTY;
    }


    printf("Initial Tree: %ld Fire: %ld\n", total_tree, total_fire);
 
    // Copy over initial game grid (Dim-1 threads)
    cudaMemcpy(d_grid, h_grid, bytes, cudaMemcpyHostToDevice);
 
    // Copy lookup-table
    cudaMemcpy(d_lookup_table, h_lookup_table, bytes2, cudaMemcpyHostToDevice);

    dim3 blockSize(BLOCK_SIZE, BLOCK_SIZE,1);
    int linGrid = (int)ceil(dim/(float)BLOCK_SIZE);
    dim3 gridSize(linGrid,linGrid,1);
 
    dim3 cpyBlockSize(BLOCK_SIZE,1,1);
    dim3 cpyGridRowsGridSize((int)ceil(dim/(float)cpyBlockSize.x),1,1);
    dim3 cpyGridColsGridSize((int)ceil((dim+2)/(float)cpyBlockSize.x),1,1);
 
    // Main game loop
    for (iter = 0; iter<maxIter; iter++) {
    #ifdef __VIDEO
		if (iter==0 || iter==50 || iter==100 || iter==200 || iter==300) {
            cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
			rgb = generate_rgb (GRID_SIZE+2, GRID_SIZE+2, h_grid, rgb);
			sprintf (image_name, "%s%d.ppm", "ca_step_", iter);
			write_CA_screenshot_image (image_name, GRID_SIZE+2, GRID_SIZE+2, rgb);
		}
    #endif
        forest_fire<<<gridSize, blockSize>>>(dim, d_grid, d_newGrid, d_lookup_table);

        // Swap our grids and iterate again
        d_tmpGrid = d_grid;
        d_grid = d_newGrid;
        d_newGrid = d_tmpGrid;
    }//iter loop
 
    // Copy back results and sum
    cudaMemcpy(h_grid, d_grid, bytes, cudaMemcpyDeviceToHost);
 
    total_tree = 0; 
    total_fire = 0;
    // Sum up alive cells and print results
    for (i = 1; i<=dim; i++) {
        for (j = 1; j<=dim; j++) {
            if (h_grid[i*(dim+2)+j] == TREE)
                total_tree++;
            else if (h_grid[i*(dim+2)+j] == FIRE)
                total_fire++;
        }
    }
    printf("Final Tree: %ld Fire: %ld\n", total_tree, total_fire);
 
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
			if (grid[y][x] == ASH) {
				// Ash: Grey
				rgb[cur + 0] = 125;
				rgb[cur + 1] = 125;
				rgb[cur + 2] = 125;
			} else if (grid[y][x] == TREE) {
				// Tree: Green
				rgb[cur + 0] = 0;
				rgb[cur + 1] = 255;
				rgb[cur + 2] = 0;
			} else if (grid[y][x] == FIRE) {
				// Fire: Red
				rgb[cur + 0] = 255;
				rgb[cur + 1] = 0;
				rgb[cur + 2] = 0;
			} else if (grid[y][x] == EMPTY) {
				// Empty: Brown
				rgb[cur + 0] = 102;
				rgb[cur + 1] = 0;
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
	fprintf(fp, "# Created by %s\n","ForestFire.cu");
	
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

