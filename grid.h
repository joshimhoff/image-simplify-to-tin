// Color type designed to store RGB values
typedef struct _color {
    int red;
    int green;
    int blue;
} Color;

// Grid type designed to store in memory what the fields in a .asc file
typedef struct _grid {
    int ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value;
    Color** color;
} Grid;

// Allocates a two dimensional array of ints to be used by the Grid struct
// @param grid, a pointer to a grid structure with the two dimensional array unallocated
// @paran nrows, the size of the array of int *
// @param ncols, the size of the arrays of ints
// @returns 0 if failure, 1 if successful
int allocateGrid(Grid* grid, int nrows, int ncols);

// Frees the two dimensional array used by the Grid struct
// @param grid, a pointer to a grid structure with the two dimensional array allocated
void freeGrid(Grid* grid);

// Read .asc file into a Grid structure
// @param grid, a pointer to the grid structure to fill with data
// @param filename, the name of the file on the filesystem
// @returns 0 if failure, 1 if successful
int readFileIntoGrid(Grid* grid, char* filename);

// Save Grid structure as .asc file
// @param grid, a pointer to the grid structure to be saved to the filesystem
// @param filename, the name of the file on the filesystem
// @returns 0 if failure, 1 if successful
int saveGridToFile(Grid* grid, char* filename);

// Copies all fields except values of the from grid to the to grid
// @param to, a pointer to the grid that we are copying the header to
// @param from, a pointer to the grid that we are copying the header from
void copyGridHeader(Grid* to, Grid* from);

// Print the grid
// @param grid, a pointer to the grid structure to be printed
void printGrid(Grid* grid);

// Check whether (i,j) is a valid index into the grid
// @param grid, a pointer to the grid to check bounds on
// @param i, is grid->values[i][j] within the bounds of the two dim array?
// @param j, is grid->values[i][j] within the bounds of the two dim array?
// @returns a boolean indicating whether (i,j) is a valid index
int inBoundsOnGrid(Grid* grid, int i, int j);

// Check whether there is actual data at point (i,j) on grid
// @param grid, a pointer to the grid to check
// @param i, is grid->values[i][j] equal to grid->NODATA_value
// @param j, is grid->values[i][j] equal to grid->NODATA_value
// @returns a boolean indicating whether there is actual data at point (i,j)
int noDataAtPoint(Grid* grid, int i, int j);
