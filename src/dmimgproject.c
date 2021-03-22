/*
**  Copyright (C) 2021  Smithsonian Astrophysical Observatory
*/

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */

#include <math.h>
#include <float.h>

#include <dslib.h>
#include <dsnan.h>
#include <histlib.h>
#include <cxcregion.h>

#include <dmfilters.h>
#include <dmimgio.h>


/* Hold info for an input image */
typedef struct {
    void *data;        // pixel values
    dmDataType dt;     // pixel datatype
    long *lAxes;       // axis lenghts
    short *mask;        // mask of valid pixels
    dmDescriptor *xdesc;  // X (or sky) coordinate descriptor
    dmDescriptor *ydesc;  // Y coordinate descriptor
    dmBlock *block; // The block image came from
} Image;

typedef struct {
    char name[30];     // Name of stat, also used for column name
    double *values;    // Output values per bin
} Stat;

typedef struct {
    double x_min;       // min x in rotated frame
    double x_max;       // max x in rotated frame
    double binsize;     // binsize

    double cos_angle;   // Save cos(angle) to only compute once
    double sin_angle;   // ditto for sin(angle)

    long num_bins;      // number of bins
    long *buflen;       // length of each buffer per bin
    double **buffers;   // buffer to hold pixel values in each bin

    double *xmid;       // x coordinate in rotated frame
    double *ymid;       // y coordinate in rotated frame    

    long *grid_image;   // image showing which grid each pixel belong to.

} Grid;

typedef struct {
    char infile[DS_SZ_PATHNAME]; // Input file name
    char outfile[DS_SZ_PATHNAME]; // ooutput file name
    double rotang;      // rotation angle
    double binsize;     // bin size
    short verbose;      // chatter level
    short clobber;      // rm outptu file if it exists?
} Parameters;

// ------------------------
// Function prototypes

Image* load_infile(char *infile);
int dmimgproject2();
int convert_coords( Image *image, double x_in, double y_in, 
        double *x_out, double *y_out);
Grid *get_output_grid(Image *image, double rotang, double binsize);
int setup_buffers(Image *image, Grid *grid);
int fill_buffers(Image *image, Grid *grid);
Stat *compute_stats(Grid *grid);
int write_output(char *outfile, Image *image, Grid *grid, Stat *stats);
int get_coordinates(Image *image, Grid *grid);


// ---------


#define NULL_GRID_VALUE -999

// ----------------------
// Functions

Image* load_infile(char *infile)
{
    // Load image

    Image *image;
    if (NULL == (image = calloc(1,sizeof(Image)))) {
        err_msg("ERROR: Cannot allocate memory for image\n");
        return(NULL);
    }
    
    if (NULL == (image->block = dmImageOpen(infile))) {
        err_msg("ERROR: Cannot load infile '%s'\n",infile);
        return(NULL);
    }

    // dmimgio
    regRegion *dss = NULL;
    long null_value;
    short has_null;
    image->dt = get_image_data(image->block, &(image->data), 
                    &(image->lAxes), &dss, &null_value, &has_null);
    get_image_wcs(image->block, &(image->xdesc), &(image->ydesc));
    image->mask = get_image_mask(image->block, image->data, 
                    image->dt, image->lAxes, dss, null_value,
                    has_null, image->xdesc, image->ydesc);

    if (dss != NULL){
        regFree(dss);
        dss=NULL;
    }
    
    return(image);
}


int convert_coords( Image *image, double x_in, double y_in, 
        double *x_out, double *y_out)
{
    // Convert from 0-based array index to physical coords
    if (image->xdesc == NULL) {
        *x_out = x_in;
        *y_out = y_in;
        return(0);        
    }

    double logical[2];
    double physical[2];
    
    logical[0] = x_in + 1;
    logical[1] = y_in + 1;
    dmCoordCalc_d(image->xdesc, logical, physical);
    if (image->ydesc){
        dmCoordCalc_d(image->ydesc, logical+1, physical+1);
    }
    *x_out = physical[0];
    *y_out = physical[1];

    return(0);
}


Grid *get_output_grid(Image *image, double rotang, double binsize)
{
    // Compute the grid in the rotated frame
    
    Grid *grid;
    if (NULL == (grid = calloc(1, sizeof(Grid)))) {
        err_msg("ERROR: Problem allocating memory");
        return(NULL);
    }

    double pi = acos(-1);

    grid->cos_angle = cos(rotang*pi/180.0);
    grid->sin_angle = sin(rotang*pi/180.0);
    grid->binsize = binsize;

    grid->x_min = DBL_MAX;
    grid->x_max = -DBL_MAX;
    
    double x_corners[4] = {0,0,image->lAxes[0]-1,image->lAxes[0]-1};
    double y_corners[4] = {0,image->lAxes[1]-1,image->lAxes[1]-1,0};
    
    int ii;
    for (ii=0;ii<4;ii++) {
        double xx,yy;

        xx = x_corners[ii];
        yy = y_corners[ii];
        
        // I only need to know the X coord in rotated coord system
        double rotx;
        rotx =  xx * grid->cos_angle + yy * grid->sin_angle;

        if (rotx < grid->x_min) {
            grid->x_min = rotx;
        }
        if (rotx > grid->x_max) {
            grid->x_max = rotx;
        }

    }   // end for ii

    grid->num_bins = floor((grid->x_max - grid->x_min)/binsize)+1;
        
    return(grid);
}


int setup_buffers(Image *image, Grid *grid)
{
    // Setup buffers for each of the grid bins
    
    if (NULL == (grid->buflen = calloc(grid->num_bins, sizeof(long)))) {
        err_msg("ERROR: problem allocating memory");
        return(1);
    }

    if (NULL == (grid->buffers = calloc(grid->num_bins, sizeof(double*)))) {
        err_msg("ERROR: problem allocating memory");
        return(1);        
    }
    

    if (NULL == (grid->grid_image = calloc(image->lAxes[0]*image->lAxes[1], sizeof(long)))) {
        err_msg("ERROR: problem allocating memory");
        return(1);        
    }

    
    // Loop through image to figure out which grid bin each
    // pixel belongs to.  Allocate buffers to match number of 
    // pixels per bin.
    
    int ii, jj;    
    for (jj=image->lAxes[1];jj--;) {
        double yterm = jj * grid->sin_angle;
        long out_y = jj * image->lAxes[0];

        for (ii=image->lAxes[0];ii--;) {
            long out_pix = ii+out_y;

            double pixval;
            pixval = get_image_value(image->data, image->dt, ii, jj,
                                     image->lAxes, image->mask);
            if (ds_dNAN(pixval)) {
                grid->grid_image[out_pix] = NULL_GRID_VALUE;
                continue;   // Skip pixel
            }

            double rotx; 
            rotx =  ii * grid->cos_angle + yterm;

            long grid_bin;
            grid_bin = floor((rotx - grid->x_min)/grid->binsize);  
            
            if (grid_bin < 0) {
                err_msg("I messed up somewhere\n");
                return(1);
            }

            if (grid_bin >= grid->num_bins) {
                err_msg("I messed up somewhere else\n");
                return(1);                
            }

            // calloc -> init to 0
            grid->buflen[grid_bin] += 1;

            grid->grid_image[out_pix] = grid_bin;
            
        } // end for ii
    } // end for jj


    int kk;
    for (kk=0;kk<grid->num_bins;kk++) {
        grid->buffers[kk] = calloc(grid->buflen[kk],sizeof(double));
    }
    
    return(0);
}


int get_coordinates(Image *image, Grid *grid)
{
    // Get coordinates of rotated x-axis.
    if (NULL == (grid->xmid = calloc(grid->num_bins,sizeof(double)))) {
        err_msg("ERROR: problem allocating memory");
        return(1);
    } 
    if (NULL == (grid->ymid = calloc(grid->num_bins,sizeof(double)))) {
        err_msg("ERROR: problem allocating memory");
        return(1);
    }
    
    int kk;
    for (kk=0;kk<grid->num_bins;kk++) {
        // The idea is to rotate and then project onto the X-axis
        // So in the rotated system we want coords for y=0.
        double xx = grid->x_min + kk*grid->binsize;
        double yy = 0;
        
        double rotx, roty;
        
        // cos(-a) = cos(a); sin(-a) = -sin(a)
        rotx = xx*grid->cos_angle - yy*grid->sin_angle;
        roty = xx*grid->sin_angle + yy*grid->cos_angle;
                
        double px, py;
        convert_coords(image, rotx, roty, &px, &py);
        grid->xmid[kk] = px;
        grid->ymid[kk] = py;
    } // end for kk

    return(0);    
}


int fill_buffers(Image *image, Grid *grid) 
{
    // Fill each grid's buffer with pixel values from image

    long *buffer_at;    
    if (NULL == (buffer_at = calloc(grid->num_bins,sizeof(long)))) {
        err_msg("ERROR: problem allocating memory");
        return(1);        
    }

    int ii, jj;    
    for (jj=image->lAxes[1];jj--;) {
        double yterm = jj * grid->sin_angle;

        for (ii=image->lAxes[0];ii--;) {
            double pixval;
            pixval = get_image_value(image->data, image->dt, ii, jj,
                                     image->lAxes, image->mask);
            if (ds_dNAN(pixval)) {
                continue;   // Skip pixel
            }

            double rotx; 
            rotx =  ii * grid->cos_angle + yterm;

            long grid_bin;
            grid_bin = floor((rotx - grid->x_min)/grid->binsize);  

            long kk = buffer_at[grid_bin];
            grid->buffers[grid_bin][kk] = pixval;
            buffer_at[grid_bin] += 1;
            
        } // end for ii
    } // end for jj

    free(buffer_at);

    return(0);    
}


Stat *compute_stats(Grid *grid)
{
    // Compute array of statistics for each grid
    
    char *list_of_stats[] = { "min", "max", "mean", "count", "sum", "median", "mode", "mid", 
                              "sigma", "extreme", "range", "q25", "q33", 
                              "q67", "q75", "olympic", NULL };

    // Count the number of statistcs
    int num_stats = 0; 
    while(list_of_stats[++num_stats]) {};
    
    // +1 so has NULL to terminate list
    Stat *stats = NULL;
    if (NULL == (stats = calloc(num_stats+1,sizeof(Stat)))) {
        err_msg("ERROR: Problem allocating memory");
        return(NULL);
    }

    
    int ii;
    for (ii=0;ii<num_stats;ii++) {
        // dmfilter.h
        double (*func)( double *vals, long nvals ) = NULL;
        if (NULL == (func = get_method(list_of_stats[ii]))) {
            err_msg("ERROR: Unknown filter %s", list_of_stats[ii]);
            return(NULL);
        }

        strcpy( stats[ii].name, list_of_stats[ii]);
        if (NULL == (stats[ii].values = calloc(grid->num_bins,sizeof(double)))) {
            err_msg("ERROR: problem allocating memory");
            return(NULL);
        }

        int gg;
        for (gg=0;gg<grid->num_bins;gg++) {
            double val;            
            val = func(grid->buffers[gg], grid->buflen[gg]);

            // Special case -- if 0 counts, set to 0 not NaN
            if ((strcmp(stats[ii].name, "count") == 0) && (ds_dNAN(val))) {
                val = 0;
            }

            stats[ii].values[gg] = val;            
        } // end for gg
        
    } // end for ii
    
    
    return(stats);
}

int write_output(char *outfile, Image *image, Grid *grid, Stat *stats) 
{
    // Write output results to table

    dmBlock *outBlock;
    if (NULL == (outBlock = dmTableCreate( outfile ) )) {
        err_msg("ERROR: Cannot create output file '%s'\n", outfile );
        return(1);
    }

    dmDataset *outDs = dmBlockGetDataset(outBlock);

    char units[100];
    memset( units, 0, sizeof(char)*100);
    dmGetUnit( dmImageGetDataDescriptor( image->block ), units, 99 );

    /* All the header stuff */
    Header_Type *hdr;
    hdr = getHdr( image->block, hdrDM_FILE );
    putHdr( outBlock, hdrDM_FILE, hdr, BASIC_STS, "dmimgproject2" );
    put_param_hist_info(outBlock, "dmimgproject2", NULL, 0);

    if (dmBlockGetNo(outBlock) != 1) {
        dmBlock *primary = dmDatasetMoveToBlock(outDs, 1);
        putHdr(primary, hdrDM_FILE, hdr, PRIMARY_STS, "dmimgproject2");
    }

    // Write data
    dmDescriptor *col;
    col = dmColumnCreate(outBlock, "X", dmDOUBLE, 0, "pix", "Rotated X coordinate");
    dmSetScalars_d(col, grid->xmid, 1, grid->num_bins);

    col = dmColumnCreate(outBlock, "Y", dmDOUBLE, 0, "pix", "Rotated Y coordinate");
    dmSetScalars_d(col, grid->ymid, 1, grid->num_bins);

    int ii = 0;
    while (strlen(stats[ii].name)) {        
        dmDataType dt;
        dt = strcmp(stats[ii].name, "count") ? dmDOUBLE : dmLONG;
        col = dmColumnCreate(outBlock, stats[ii].name, dt, 0, units, "Statistic Computed" );
        dmSetScalars_d(col, stats[ii].values, 1, grid->num_bins);
        ii++;
    }
    dmBlockClose(outBlock);


    // Save grid image (useful for debugging)
    outBlock = dmDatasetCreateImage(outDs, "GRID_IMG", dmLONG,
                    image->lAxes, 2);
    dmDescriptor *outDesc = dmImageGetDataDescriptor(outBlock);
    dmDescriptorSetNull_l(outDesc, NULL_GRID_VALUE);
    putHdr( outBlock, hdrDM_FILE, hdr, BASIC_STS, "dmimgproject2" );
    dmSetArray_l(outDesc, grid->grid_image, image->lAxes[0]*image->lAxes[1]);                    
    dmBlockClose(outBlock);
    dmDatasetClose(outDs);

    return(0);
}

Parameters *get_parameters()
{
    Parameters *pars;
    
    if (NULL == (pars = calloc(1,sizeof(Parameters)))) {
        err_msg("ERROR: problem allocating memory");
        return(NULL);
    }

    clgetstr("infile", pars->infile, DS_SZ_PATHNAME);
    clgetstr("outfile", pars->outfile, DS_SZ_PATHNAME);
    pars->rotang = clgetd("angle");
    pars->binsize = clgetd("binsize");
    pars->verbose = clgeti("verbose");
    pars->clobber = clgetb("clobber");

    if (pars->verbose >= 1){
        printf("dmimgproject2\n");
        printf("%15s = %-s\n", "infile", pars->infile);        
        printf("%15s = %-s\n", "outfile", pars->outfile);        
        printf("%15s = %g\n", "angle", pars->rotang);        
        printf("%15s = %g\n", "binsize", pars->binsize);        
        printf("%15s = %d\n", "verbose", pars->verbose);        
        printf("%15s = %-s\n", "clobber", pars->clobber ? "yes" : "no");                
    }
    
    return(pars);
}


int dmimgproject2()
{
    Parameters *pars;
    if (NULL == (pars = get_parameters())) {
        return(-9);
    }
    

    Image *image;    
    if (NULL == (image = load_infile(pars->infile))) {
        return(-1);
    }

    if ( ds_clobber(pars->outfile, pars->clobber, NULL ) != 0 ) {
        return(-2);
    }

    Grid *grid;
    if (NULL == (grid = get_output_grid(image, pars->rotang, pars->binsize))) {
        return(-3);
    }

    if (0 != setup_buffers(image, grid)) {
        return(-5);
    }

    if (0 != fill_buffers(image,grid)) {
        return(-7);
    }

    Stat *stats;
    if (NULL == (stats = compute_stats(grid))) {
        return(-4);
    }

    if (0 != get_coordinates(image, grid)) {
        return(-8);
    }

    if (0 != write_output(pars->outfile, image, grid, stats)) {
        return(-6);
    }

    return(0);

}








