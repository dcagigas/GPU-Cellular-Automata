#include <stdio.h>
#include <stdlib.h>
//#include <openssl/rand.h> 			// For random numbers. OpenSSL library must be installed.

#define SIDE 512
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

//#define __VIDEO

typedef char CELL_TYPE;

float r4_uniform_01 ( int *seed );
char* generate_rgb (int width, int height, CELL_TYPE *grid_, char *rgb);
void write_CA_screenshot_image (const char *filename, int width, int height, char *img);

int seed = SEED; 	// Seed for random uniform number output.


void main_computation (CELL_TYPE *grid_)
{
    int t,i,j;
	CELL_TYPE (*grid)[SIDE][SIDE] = (CELL_TYPE (*)[SIDE][SIDE]) grid_;
    CELL_TYPE lookup_table[HEAD+1][LOOKUP_TABLE_LIMIT];

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

    #ifdef __VIDEO
	    char *rgb = NULL;
	    rgb= (char *)malloc(3 * sizeof(char) *SIDE*SIDE);
        char image_name[80];
    #endif

    //setup_lookup_table (lookup_table);

#pragma scop
    for (t = 0; t < TIMESTEP; t++)
      for (i = 1; i < SIDE - 1; i++)
        for (j = 1; j < SIDE - 1; j++) {
           grid[(t+1)%2][i][j] = lookup_table [grid[t%2][i][j]] [ (grid[t%2][i-1][j] + grid[t%2][i+1][j] + grid[t%2][i][j-1] + grid[t%2][i][j+1] + 
		                         grid[t%2][i-1][j-1] + grid[t%2][i-1][j+1] + grid[t%2][i+1][j-1] + grid[t%2][i+1][j+1]) ];
            #ifdef __VIDEO
		        if (t==0 || t== 300 || t==600 || t==1000 ) {
			        rgb = generate_rgb (SIDE, SIDE, (CELL_TYPE *)grid[0], rgb);
			        sprintf (image_name, "%s%d.ppm", "ca_step_", t);
			        write_CA_screenshot_image (image_name, SIDE, SIDE, rgb);
		        }
            #endif
        }
#pragma endscop
}



long int  init_wire_world (int dim, CELL_TYPE *grid) {
    int i, j, k;
    //long int conductors;
    long int heads = 0;

    // Set grid to empty
    for(i = 0; i<dim+2; i++) {
        for(j = 0; j<dim+2; j++) {
            grid[i*(dim+2)+j] = EMPTY;
        }
    }
    
    // Set a "conductor cross" in the middle of the grid
    for(i = 0; i<dim+2; i++) {
        grid[i*(dim+2)+((dim+1)/2)] = CONDUCTOR;
    }
    for(j = 0; j<dim+2; j++) {
        grid[((dim+1)/2)*(dim+2) + j] = CONDUCTOR;
    }
    //conductors = (dim+2)*2;

    // Set concentric squares in the grid with STRIDE padding.
    for (k = STRIDE; k < (dim+2)/2; k = k + STRIDE) {
        // Draw square:
        // Draw left side and right side:
        for (i=k*(dim+2)+k; i<  (dim+2)*(dim+2)-k*(dim+2); i+=(dim+2)) {
            grid[i] = CONDUCTOR;
            grid[i+(dim+2)-k-k] = CONDUCTOR;
         }
        // Draw up side and down side:
        for (j=k*(dim+2)+k; j< k*(dim+2)+(dim+2)-k; j++) {
            grid[j] = CONDUCTOR;
        }
        for (j=(dim+2)*(dim+2)-k*(dim+2)+(dim+2)-k; j >= (dim+2)*(dim+2)-k*(dim+2)+k ; j--) {
            grid[j] = CONDUCTOR;
        }        
    }
    
    // Set initial random heads
    for(i = 0; i<dim+2; i++) {
        for(j = 0; j<dim+2; j++) {
            if ( grid[i*(dim+2)+j] == CONDUCTOR && r4_uniform_01(&seed) < P_HEAD) {
                grid[i*(dim+2)+j] = HEAD;
                heads++;
            }
        }
    }
    
    return heads;
}




//void print_final_statistics (int ca, CELL_TYPE (*grid)[SIDE][SIDE]) {
void print_final_statistics (int ca, CELL_TYPE *grid_) {
	int i,j;
    long unsigned int total_heads = 0;
	CELL_TYPE (*grid)[SIDE][SIDE] = (CELL_TYPE (*)[SIDE][SIDE]) grid_;

	for (i = 1; i<SIDE-1; i++) {
		for (j = 1; j<SIDE-1; j++) {
            if (grid[ca][i][j] == HEAD) {
			    total_heads++;
            }
		}
	}
    printf("Final Heads: %ld\n", total_heads);
}



int main(int argc, char* argv[])
{
    // Define variables
    CELL_TYPE *grid;
 

    // Allocate grid
    grid = (CELL_TYPE *) malloc (2*sizeof(CELL_TYPE)*(SIDE)*(SIDE));  // grid[0][][] --> current CA state
                                                                                // grid[1][][] --> next CA state
    // Assign initial population randomly
    long int heads = init_wire_world (SIDE-2, grid);
    printf("Initial Heads: %ld\n", heads);

    // Main FireForest loop
    main_computation (grid);

    // Print results
    print_final_statistics (0, grid); 

    // Release memory
    free(grid);

}
/**/



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
	fprintf(fp, "# Created by %s\n","WireWorld");
	
	//image size
	fprintf(fp, "%d %d\n", width, height);
	
	// rgb component depth
	fprintf(fp, "%d\n",255);
	
	// pixel data
	fwrite(img, 3 * width, height, fp);
	fclose(fp);
}




// 256 Result in console: 
// 512 Result in console: initial heads: 42, final heads: 493
// 1024 Result in console: 
// 2048 Result in console: 
// 4096 Result in console: 
// 8192 Result in console: 

