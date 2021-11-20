#include <stdio.h>
#include <stdlib.h>
#include <openssl/rand.h> 			// For random numbers. OpenSSL library must be installed.

typedef unsigned char CELL_TYPE;   // Used for CA grid states

#define TIMESTEP 1024
#define BENCH_RAD 1
#define SIDE 258
//#define CELL_TYPE char

#define CELL_NEIGHBOURS 4   // Von Neumann distance
#define SEED 1985

/*
const CELL_TYPE EMPTY = 0;
const CELL_TYPE TREE = 1;
const CELL_TYPE ASH = 2;
const CELL_TYPE FIRE = 9;
*/
#define EMPTY 0
#define TREE 1
#define ASH 2
#define FIRE 9


#define P_EMPTY 0.20
#define P_TREE 0.80
#define P_FIRE (1.0/(SIDE*SIDE))

#define LOOKUP_TABLE_LIMIT (4*FIRE+1)

#define __VIDEO



float r4_uniform_01 ( int *seed );
uint8_t* generate_rgb (int width, int height, CELL_TYPE *grid_, uint8_t *rgb);
void write_CA_screenshot_image (const char *filename, int width, int height, uint8_t *img);

int seed = SEED; 	// Seed for random uniform number output.


void setup_lookup_table (CELL_TYPE *lookup_table_) {
    int i,j;

    CELL_TYPE (*lookup_table)[LOOKUP_TABLE_LIMIT] = (CELL_TYPE (*)[LOOKUP_TABLE_LIMIT]) lookup_table_;

/**/
    // 0: Empty. If a cell is EMPTY remains EMPTY 
    for (j=0; j<LOOKUP_TABLE_LIMIT; j++) {
        lookup_table[EMPTY][j] = EMPTY;
    }
    // 1: Tree. A TREE cell with at least one neighbor on fire becomes on FIRE 
    for (j=FIRE; j<LOOKUP_TABLE_LIMIT; j++) {
        lookup_table[TREE][j] = FIRE;
    }
    for (j=0; j<FIRE; j++) {
        lookup_table[TREE][j] = TREE;
    }
    // 2: Ash. An ASH cell becomes EMPTY if it does not have a neighbor on FIRE 
    for (j=FIRE; j<LOOKUP_TABLE_LIMIT; j++) {
        lookup_table[ASH][j] = ASH;
    }
    for (j=0; j<FIRE; j++) {
        lookup_table[ASH][j] = EMPTY;
    }
    // 9: Fire. A FIRE cell becomes ASH 
    for (j=0; j<LOOKUP_TABLE_LIMIT; j++) {
        lookup_table[FIRE][j] = ASH;
    }
    // Fill with EMPTY not used lookup_table entries
    for (i=ASH+1; i<FIRE-1; i++)
        for (j=0; j<LOOKUP_TABLE_LIMIT; j++) 
            lookup_table[i][j] = EMPTY;
/**/
}



void main_computation (CELL_TYPE *grid_, CELL_TYPE *lookup_table_)
{
    int t,i,j;
    CELL_TYPE (*lookup_table)[LOOKUP_TABLE_LIMIT] = (CELL_TYPE (*)[LOOKUP_TABLE_LIMIT]) lookup_table_;
	CELL_TYPE (*grid)[SIDE][SIDE] = (CELL_TYPE (*)[SIDE][SIDE]) grid_;

#ifdef __VIDEO
	uint8_t *rgb = NULL;
	rgb= (uint8_t *)malloc(3 * sizeof(uint8_t) *SIDE*SIDE);
    char image_name[80];
#endif


#pragma scop
    for (t = 0; t < TIMESTEP; t++)
      for (i = 1; i < SIDE - 1; i++)
        for (j = 1; j < SIDE - 1; j++) {
           grid[(t+1)%2][i][j] =
           lookup_table [grid[t%2][i][j]] [ grid[t%2][i-1][j] + grid[t%2][i+1][j] + grid[t%2][i][j-1] + grid[t%2][i][j+1] ];
            #ifdef __VIDEO
		        if (t==300) {
			        rgb = generate_rgb (SIDE, SIDE, (CELL_TYPE *)grid[0], rgb);
			        sprintf (image_name, "%s%d.ppm", "ca_step_", t);
			        write_CA_screenshot_image (image_name, SIDE, SIDE, rgb);
		        }
            #endif
        }
#pragma endscop
}

/**/


//void init_grid (CELL_TYPE (*grid)[SIDE][SIDE]) {
void init_grid (CELL_TYPE *grid_) {
    int i,j;
    float p;
    long unsigned int total_tree = 0, total_fire = 0;
	CELL_TYPE (*grid)[SIDE][SIDE] = (CELL_TYPE (*)[SIDE][SIDE]) grid_;

    // Assign empty to borders
    for (j = 0; j<SIDE; j++) {
        grid[0][0][j] = EMPTY;
        grid[0][SIDE-1][j] = EMPTY;
    }
    for (i = 0; i<SIDE; i++) {
        grid[0][i][0] = EMPTY;
        grid[0][i][SIDE-1] = EMPTY;
    }

    // Assign an initial forest randomly
    for (i = 1; i<SIDE-1; i++) {
        for (j = 1; j<SIDE-1; j++) {
            p = r4_uniform_01(&seed);
            if (p < P_FIRE) {
                grid[0][i][j] = FIRE;
                total_fire++;
            } else if (p < P_EMPTY) {
                grid[0][i][j] = EMPTY;
            } else {
                grid[0][i][j] = TREE;
                total_tree++;
            }
        }
    }

    printf("Initial Tree: %ld Fire: %ld\n", total_tree, total_fire);
}



//void print_final_statistics (int ca, CELL_TYPE (*grid)[SIDE][SIDE]) {
void print_final_statistics (int ca, CELL_TYPE *grid_) {
	int i,j;
    long unsigned int total_tree = 0, total_fire = 0;
	CELL_TYPE (*grid)[SIDE][SIDE] = (CELL_TYPE (*)[SIDE][SIDE]) grid_;

	for (i = 1; i<SIDE-1; i++) {
		for (j = 1; j<SIDE-1; j++) {
            if (grid[ca][i][j] == FIRE) {
			    total_fire++;
            } else if (grid[ca][i][j] == TREE) {
			    total_tree++;
            }
		}
	}
    printf("Final Tree: %ld Fire: %ld\n", total_tree, total_fire);
}



int main(int argc, char* argv[])
{
    // Define variables
    //CELL_TYPE (*grid)[SIDE][SIDE];
    //CELL_TYPE (*lookup_table)[FIRE][LOOKUP_TABLE_LIMIT];
    CELL_TYPE *grid;
    CELL_TYPE *lookup_table;
 

    // Allocate grid
    //grid = (CELL_TYPE (*)[SIDE][SIDE]) malloc (2*sizeof(CELL_TYPE)*(SIDE)*(SIDE));  // grid[0][][] --> current CA state
    grid = (CELL_TYPE *) malloc (2*sizeof(CELL_TYPE)*(SIDE)*(SIDE));  // grid[0][][] --> current CA state
                                                                                // grid[1][][] --> next CA state
    // Allocate Lookup_table
    //lookup_table = (CELL_TYPE (*)[FIRE][LOOKUP_TABLE_LIMIT]) malloc ( sizeof(CELL_TYPE)*(FIRE)*(LOOKUP_TABLE_LIMIT) ); 
    lookup_table = (CELL_TYPE *) malloc (sizeof(CELL_TYPE)*(FIRE+1)*(LOOKUP_TABLE_LIMIT) ); 

  
    // Assign initial population randomly
    init_grid (grid);

    // Init lookup_table
    setup_lookup_table (lookup_table);

    // Main FireForest loop
    main_computation (grid, lookup_table);

    // Print results
    print_final_statistics (0, grid); 

    // Release memory
    free(lookup_table);
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
	fprintf(fp, "# Created by %s\n","ForestFire");
	
	//image size
	fprintf(fp, "%d %d\n", width, height);
	
	// rgb component depth
	fprintf(fp, "%d\n",255);
	
	// pixel data
	fwrite(img, 3 * width, height, fp);
	fclose(fp);
}




// 256 Result in console: "Total Alive: 3281"
// 512 Result in console: "Total Alive: 11072"
// 1024 Result in console: "Total Alive: 45224"
// 2048 Result in console: "Total Alive: 182485"
// 4096 Result in console: "Total Alive: 724393"
// 8192 Result in console: "Total Alive: 2896683"

