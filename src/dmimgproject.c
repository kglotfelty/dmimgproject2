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


#include <dslib.h>
#include <dsnan.h>
#include <histlib.h>
#include <math.h>
#include <float.h>
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
    char name[30];
    double *values;    
} Stat;


typedef struct {
    double x_min;
    double x_max;
    double binsize;

    double cos_angle;
    double sin_angle;

    long num_bins;
    long *buflen;
    double **buffers;

    double *xmid;
    double *ymid;
    
} Grid;





Image* load_infile(char *infile);
int dmimgproject();
int convert_coords( Image *image, double x_in, double y_in, 
        double *x_out, double *y_out);
Grid *get_output_grid(Image *image, double rotang, double binsize);
int setup_buffers(Image *image, Grid *grid);
int fill_buffers(Image *image, Grid *grid);
Stat *compute_stats(Grid *grid);
int write_output(char *outfile, Image *image, Grid *grid, Stat *stats);
int get_coordinates(Image *image, Grid *grid);




Image* load_infile(char *infile){

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
    Grid *grid;
    if (NULL == (grid = calloc(1, sizeof(Grid)))) {
        err_msg("ERROR: Problem allocating memory");
        return(NULL);
    }

    grid->x_min = DBL_MAX;
    grid->x_max = -DBL_MAX;

    grid->cos_angle = cos(rotang*3.141592/180.0);
    grid->sin_angle = sin(rotang*3.141592/180.0);

    printf("cos %g\n",grid->cos_angle);
    printf("sin %g\n",grid->sin_angle);

    
    int ii;
    double x_corners[4] = {0,0,image->lAxes[0]-1,image->lAxes[0]-1};
    double y_corners[4] = {0,image->lAxes[1]-1,image->lAxes[1]-1,0};
    
    for (ii=0;ii<4;ii++) {
        double xx,yy;

        xx = x_corners[ii];
        yy = y_corners[ii];
        
        double rotx; //, roty;        
        rotx =  xx * grid->cos_angle + yy * grid->sin_angle;
        // roty = -xx * grid->sin_angle + yy * grid->cos_angle + y0;

        if (rotx < grid->x_min) {
            grid->x_min = rotx;
        }
        if (rotx > grid->x_max) {
            grid->x_max = rotx;
        }

        printf("lo = %g\n", grid->x_min);
        printf("hi = %g\n", grid->x_max);
    }   

    long grid_len = floor((grid->x_max - grid->x_min)/binsize)+1;
    printf("%ld\n",grid_len);

    grid->binsize = binsize;
    grid->num_bins = grid_len;
        
    return(grid);
}


int setup_buffers(Image *image, Grid *grid)
{

    grid->buflen = calloc(grid->num_bins, sizeof(long));
    grid->buffers = calloc(grid->num_bins, sizeof(double*));
    
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

            double rotx; //, roty;        
            rotx =  ii * grid->cos_angle + yterm;

            long grid_bin;
            grid_bin = floor((rotx - grid->x_min)/grid->binsize);  
            
            if (grid_bin < 0) {
                printf("I messed up somewehre\n");
                return(1);
            }

            if (grid_bin >= grid->num_bins) {
                printf("I messed up somewehre\n");
                return(1);                
            }

            grid->buflen[grid_bin] += 1;
            
        } // end for ii
    } // end for jj

    int kk;
    for (kk=0;kk<grid->num_bins;kk++) {
        grid->buffers[kk] = calloc(grid->buflen[kk],sizeof(double));
    }

    fill_buffers(image,grid);
    
    return(0);
}


int get_coordinates(Image *image, Grid *grid)
{
    grid->xmid = calloc(grid->num_bins,sizeof(double));
    grid->ymid = calloc(grid->num_bins,sizeof(double));
    

    int kk;
    for (kk=0;kk<grid->num_bins;kk++) {
        double xx = grid->x_min + (kk+0.5)*grid->binsize;
        double yy = 0;
        
        double rotx, roty;
        
        // cos(-a) = cos(a); sin(-a) = -sin(a)
        rotx = xx*grid->cos_angle - yy*grid->sin_angle;
        roty = xx*grid->sin_angle + yy*grid->cos_angle;
                
        double px, py;
        convert_coords(image, rotx, roty, &px, &py);
        grid->xmid[kk] = px;
        grid->ymid[kk] = py;
    }

    return(0);    
}


int fill_buffers(Image *image, Grid *grid) 
{

    long *buffer_at = calloc(grid->num_bins,sizeof(long));

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

            double rotx; //, roty;        
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
    Stat *stats = NULL;
    

    char *list_of_stats[] = { "min", "max", "mean", "count", "sum", "median", "mode", "mid", 
                              "sigma", "extreme", "range", "q25", "q33", 
                              "q67", "q75", "olympic", NULL };
    char *stat_name;
    int ii, num_stats;
    
    num_stats=0;
    while(list_of_stats[++num_stats]) {};
    
    stats = calloc(num_stats+1,sizeof(Stat));  // +1 so has NULL to terminate list
    
    for (ii=0;ii<num_stats;ii++) {
        double (*func)( double *vals, long nvals ) = NULL;
        double val;


        if (NULL == (func = get_method(list_of_stats[ii]))) {
            err_msg("ERROR: Unknown filter %s", list_of_stats[ii]);
            return(NULL);
        }

        printf("Working on stat %s\n", list_of_stats[ii]);
        strcpy( stats[ii].name, list_of_stats[ii]);
        stats[ii].values = calloc(grid->num_bins,sizeof(double));

        int gg;
        for (gg=0;gg<grid->num_bins;gg++) {
            double val;
            
            val = func(grid->buffers[gg], grid->buflen[gg]);
            stats[ii].values[gg] = val;
            
        } // end for gg
        
    } // end for ii
    
    
    return(stats);
    
}

int write_output(char *outfile, Image *image, Grid *grid, Stat *stats) 
{

    dmBlock *outBlock;
    if (NULL == (outBlock = dmTableCreate( outfile ) )) {
        err_msg("ERROR: Cannot create output file '%s'\n", outfile );
        return(1);
    }

    char units[100];
    memset( units, 0, sizeof(char)*100);
    dmGetUnit( dmImageGetDataDescriptor( image->block ), units, 99 );

    /* All the header stuff */
    Header_Type *hdr;
    hdr = getHdr( image->block, hdrDM_FILE );
    putHdr( outBlock, hdrDM_FILE, hdr, BASIC_STS, "dmimgproject" );
    put_param_hist_info(outBlock, "dmimgproject", NULL, 0);

    dmDescriptor *col;
    col = dmColumnCreate(outBlock, "X", dmDOUBLE, 0, "pix", "X coordinate ish");
    dmSetScalars_d(col, grid->xmid, 1, grid->num_bins);

    col = dmColumnCreate(outBlock, "Y", dmDOUBLE, 0, "pix", "Y coordinate ish");
    dmSetScalars_d(col, grid->ymid, 1, grid->num_bins);

    int ii = 0;
    while (strlen(stats[ii].name)) {
        
        col = dmColumnCreate(outBlock, stats[ii].name, dmDOUBLE, 0, units, stats[ii].name );
        dmSetScalars_d(col, stats[ii].values, 1, grid->num_bins);
        ii++;
    }

    dmTableClose(outBlock);

    
    return(0);
}


int dmimgproject()
{
    char infile[DS_SZ_PATHNAME];
    char outfile[DS_SZ_PATHNAME];
    double rotang;
    double binsize;
    short verbose;
    short clobber;
    
    clgetstr("infile", infile, DS_SZ_PATHNAME);
    clgetstr("outfile", outfile, DS_SZ_PATHNAME);
    rotang = clgetd("angle");
    binsize = clgetd("binsize");
    verbose = clgeti("verbose");
    clobber = clgetb("clobber");

    if (verbose >= 1){
        printf("dmimgproject\n");
        printf("%15s = %-s\n", "infile", infile);        
        printf("%15s = %-s\n", "outfile", outfile);        
        printf("%15s = %g\n", "angle", rotang);        
        printf("%15s = %g\n", "binsize", binsize);        
        printf("%15s = %d\n", "verbose", verbose);        
        printf("%15s = %-s\n", "clobber", clobber ? "yes" : "no");                
    }


    Image *image;    
    if (NULL == (image = load_infile(infile))) {
        return(-1);
    }

    if ( ds_clobber(outfile, clobber, NULL ) != 0 ) {
        return(-2);
    }

    Grid *grid;
    if (NULL == (grid = get_output_grid(image, rotang, binsize))) {
        return(-3);
    }

    if (0 != setup_buffers(image, grid)) {
        return(-5);
    }

    Stat *stats;
    if (NULL == (stats = compute_stats(grid))) {
        return(-4);
    }

    get_coordinates(image, grid);

    if (0 != write_output(outfile, image, grid, stats)) {
        return(-6);
    }

    return(0);

}








