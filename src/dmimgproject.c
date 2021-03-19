/*                                                                
**  Copyright (C) 2004-2008  Smithsonian Astrophysical Observatory 
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

typedef enum { PROJ_X=0, PROJ_Y=1 } Projection_Type;


double get_image_value( void *data, dmDataType dt, long *lAxes, 
			long xx, long yy, regRegion *dss, 
			dmDescriptor *xAxis, dmDescriptor *yAxis,
			long nullval, short has_null )
{

  long npix = xx + (yy * lAxes[0] );
  double retval;

  /* Okay, first get all the data from the different data types.  
     Cast everything to doubles */

  switch ( dt ) {
    
  case dmBYTE: {
    unsigned char *img = (unsigned char*)data;
    retval = img[npix];
    break;
  }
    
  case dmSHORT: {
    short *img = (short*)data;
    retval = img[npix];
    break;
  }
    
  case dmUSHORT: {
    unsigned short *img = (unsigned short*)data;
    retval = img[npix];
    break;
  }
    
  case dmLONG: {
    long *img = (long*)data;
    retval = img[npix];
    break;
  }
    
  case dmULONG: {
    unsigned long *img = (unsigned long*)data;
    retval = img[npix];
    break;
  }
    
  case dmFLOAT: {
    float *img = (float*)data;
    retval = img[npix];
    break;
  }
  case dmDOUBLE: {
    double *img = (double*)data;
    retval = img[npix];
    break;
  }
  default:
    ds_MAKE_DNAN( retval );

  }

  
  /* Now ... if it is an integer data type, it could possibly have a
     null value. Check for that */

  if ( has_null && ( retval == nullval ) ) 
    ds_MAKE_DNAN( retval );

  /* If the image has a data sub space (aka a region filter applied)
     then need to convert coords to physical and check */
  if ( dss && xAxis ) {
    double pos[2];
    double loc[2];
    pos[0]=xx+1;
    pos[1]=yy+1;

    if (yAxis) {  /* If no y axis, then xAxis has 2 components */
      dmCoordCalc_d( xAxis, pos, loc );
      dmCoordCalc_d( yAxis, pos+1, loc+1 );
    } else {
      dmCoordCalc_d( xAxis, pos, loc );
    }
    if ( !regInsideRegion( dss, loc[0], loc[1] ) )
      ds_MAKE_DNAN( retval );
  }

  return(retval);

}



/* Load the data into memory,  check for DSS, null values */
dmDataType get_image_data( dmBlock *inBlock, void **data, long **lAxes,
			   regRegion **dss, long *nullval, short *nullset )
{

  dmDescriptor *imgDesc;
  dmDataType dt;
  dmDescriptor *grp;
  dmDescriptor *imgdss;

  long naxes;
  long npix;
  char ems[1000];

  *nullval = INDEFL;
  *dss = NULL;
  *nullset = 0;
  
  imgDesc = dmImageGetDataDescriptor( inBlock );

  /* Sanity check, only 2D images */
  naxes = dmGetArrayDimensions( imgDesc, lAxes );
  if ( naxes != 2 ) {
    return( dmUNKNOWNTYPE );
  }
  npix = (*lAxes)[0] * (*lAxes)[1];
  dt = dmGetDataType( imgDesc );


  /* Okay, first lets get the image descriptor */
  grp = dmArrayGetAxisGroup( imgDesc, 1 );
  dmGetName( grp, ems, 1000);
  imgdss = dmSubspaceColOpen( inBlock, ems );
  if ( imgdss )
    *dss = dmSubspaceColGetRegion( imgdss);
  
  
  switch ( dt ) 
    {
    case dmBYTE:
      *data = ( void *)calloc( npix, sizeof(char ));
      dmGetArray_ub( imgDesc, (unsigned char*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
	*nullset=0;
      } else
	*nullset=1;
      break;
      
    case dmSHORT:
      *data = ( void *)calloc( npix, sizeof(short ));
      dmGetArray_s( imgDesc, (short*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
	*nullset=0;
      } else
	*nullset=1;
      break;
      
    case dmUSHORT:
      *data = ( void *)calloc( npix, sizeof(short ));
      dmGetArray_us( imgDesc, (unsigned short*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
	*nullset=0;
      } else
	*nullset=1;
      break;
      
    case dmLONG:
      *data = ( void *)calloc( npix, sizeof(long ));
      dmGetArray_l( imgDesc, (long*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
	*nullset=0;
      } else
	*nullset=1;
      break;
      
    case dmULONG:
      *data = ( void *)calloc( npix, sizeof(long ));
      dmGetArray_ul( imgDesc, (unsigned long*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
	*nullset=0;
      } else
	*nullset=1;
      break;
      
    case dmFLOAT:
      *data = ( void *)calloc( npix, sizeof(float ));
      dmGetArray_f( imgDesc, (float*) *data, npix );
      *nullset = 0;
      break;
      
    case dmDOUBLE:
      *data = ( void *)calloc( npix, sizeof(double ));
      dmGetArray_d( imgDesc, (double*) *data, npix );
      *nullset = 0;
      break;
      
    default:
      return( dmUNKNOWNTYPE );
    }

  return(dt);

}



/* Get the WCS descriptor */
short  get_image_wcs( dmBlock *imgBlock, dmDescriptor **xAxis, 
		      dmDescriptor **yAxis )
{
  

  dmDescriptor *imgData;
  long n_axis_groups;

  imgData = dmImageGetDataDescriptor( imgBlock );
  n_axis_groups = dmArrayGetNoAxisGroups( imgData );
  

  /* This is the usual trick ... can have 1 axis group w/ 
     dimensionality 2 (eg a vector column) or can have
     2 axis groups w/ dimensionaity 1 (eg 2 disjoint columns)*/

  if ( n_axis_groups == 1 ) {
    dmDescriptor *pos = dmArrayGetAxisGroup( imgData, 1 );
    dmDescriptor *xcol;
    long n_components;
    
    n_components = dmGetElementDim( pos );
    if ( n_components != 2 ) {
      err_msg("ERROR: could not find 2D image\n");
      return(-1);
    }
    
    xcol = dmGetCpt( pos, 1 );
    
    *xAxis = pos;
    *yAxis = NULL;
    
  } else if ( n_axis_groups == 2 ) {
    dmDescriptor *xcol;
    dmDescriptor *ycol;
  
    xcol = dmArrayGetAxisGroup( imgData, 1 );
    ycol = dmArrayGetAxisGroup( imgData, 2 );

    *xAxis = xcol;
    *yAxis = ycol;
    
  } else {
    err_msg("Invalid number of axis groups\n");
    *xAxis = NULL;
    *yAxis = NULL;
    return(-1);
  }

  return(0);

}


/* Routine to write out the data */
int write_outfile( dmBlock *inBlock,
		   dmBlock *outBlock,
		   long *bin,
		   double *grid,
		   double *counts,
		   double *area,
		   double *norm,
		   double *mins,
		   double *maxs,
		   double *modes,
		   double *ranges,
		   double *medians,
		   double *q25,
		   double *q33,
		   double *q67,
		   double *q75,
		   double *sigs,
		   double *nmodes,
		   double *loc_min,
		   double *loc_max,
		   double *loc_med,
		   long npoints,
		   Projection_Type pt
		   )
{
  Header_Type *hdr;
  dmDescriptor *bcol, *gcol, *ccol, *acol, *ncol;

  dmDescriptor *mins_col;
  dmDescriptor *maxs_col;
  dmDescriptor *modes_col;
  dmDescriptor *ranges_col;
  dmDescriptor *medians_col;
  dmDescriptor *q25_col;
  dmDescriptor *q33_col;
  dmDescriptor *q67_col;
  dmDescriptor *q75_col;
  dmDescriptor *sigs_col;
  dmDescriptor *nmodes_col;
  dmDescriptor *lmin_col;
  dmDescriptor *lmax_col;
  dmDescriptor *lmed_col;
  

  char *name;
  
  char units[100];

  
  if ( pt == PROJ_X ) 
    name = "X" ;
  else
    name = "Y";


  memset( units, 0, sizeof(char)*100);
  dmGetUnit( dmImageGetDataDescriptor( inBlock ), units, 99 );


  /* All the header stuff */
  hdr = getHdr( inBlock, hdrDM_FILE );
  putHdr( outBlock, hdrDM_FILE, hdr, BASIC_STS, "dmimgproject" );
  put_param_hist_info(outBlock, "dmimgproject", NULL, 0);

  if ( NULL == ( bcol = dmColumnCreate( outBlock, "BIN", dmLONG,
				0, NULL, "Bin number" ) ) ) {
    /* ERROR */
    return(-1);
  }

  if ( NULL == ( gcol = dmColumnCreate( outBlock, name, dmDOUBLE,
				0, NULL, "location along axis" ) ) ) {
    /* ERROR */
    return(-1);
  }

  if ( NULL == ( ccol = dmColumnCreate( outBlock, "SUM", dmDOUBLE,
				0, units, "Raw projection" ) ) ) {
    /* ERROR */
    return(-1);
  }
  
  if ( NULL == ( acol = dmColumnCreate( outBlock, "NUMBER", dmDOUBLE,
				0, "pixels", "Number of valid pixels" ) ) ) {
    /* ERROR */
    return(-1);
  }

  if ( NULL == ( ncol = dmColumnCreate( outBlock, "MEAN", dmDOUBLE,
					0, NULL, "/pixel" ) ) ) {
    /* ERROR */
    return(-1);
  }


  
  mins_col =dmColumnCreate(outBlock,"MIN",dmDOUBLE,0,units,"Min value");
  lmin_col = dmColumnCreate( outBlock, "MIN_LOC", dmDOUBLE, 0, NULL, 
			     "Location of Min value");

  maxs_col =dmColumnCreate(outBlock,"MAX",dmDOUBLE,0,units,"Max value");
  lmax_col = dmColumnCreate( outBlock, "MAX_LOC", dmDOUBLE, 0, NULL, 
			     "Location of Max value");

  medians_col =dmColumnCreate(outBlock,"MEDIAN",dmDOUBLE,0,units,
			      "Middle value (lower)");
  lmed_col = dmColumnCreate( outBlock, "MED_LOC", dmDOUBLE, 0, NULL, 
			     "Location of Median value");

  modes_col =dmColumnCreate(outBlock,"MODE",dmDOUBLE,0,units,"Mode=3*median-2*mean");
  ranges_col =dmColumnCreate(outBlock,"RANGE",dmDOUBLE,0,units,"Max - Min");
  q25_col =dmColumnCreate(outBlock,"QUANT_25",dmDOUBLE,0,units,"25% Quantile");
  q33_col =dmColumnCreate(outBlock,"QUANT_33",dmDOUBLE,0,units,"33% Quantile");
  q67_col =dmColumnCreate(outBlock,"QUANT_67",dmDOUBLE,0,units,"67% Quantile");
  q75_col =dmColumnCreate(outBlock,"QUANT_75",dmDOUBLE,0,units,"75% Quantile");
  sigs_col =dmColumnCreate(outBlock,"STDEV",dmDOUBLE,0,units,"Sigma/Standard Deviation");
  nmodes_col =dmColumnCreate(outBlock,"NMODE",dmDOUBLE,0,units,"Normalized mode=mode/median");
  



  dmSetScalars_l( bcol, bin, 1, npoints );
  dmSetScalars_d( gcol, grid, 1, npoints );
  dmSetScalars_d( ccol, counts, 1, npoints );
  dmSetScalars_d( acol, area, 1, npoints );
  dmSetScalars_d( ncol, norm, 1, npoints );


  dmSetScalars_d(mins_col, mins, 1, npoints);
  dmSetScalars_d(maxs_col, maxs, 1, npoints);
  dmSetScalars_d(modes_col, modes, 1, npoints);
  dmSetScalars_d(ranges_col, ranges, 1, npoints);
  dmSetScalars_d(medians_col, medians, 1, npoints);
  dmSetScalars_d(q25_col, q25, 1, npoints);
  dmSetScalars_d(q33_col, q33, 1, npoints);
  dmSetScalars_d(q67_col, q67, 1, npoints);
  dmSetScalars_d(q75_col, q75, 1, npoints);
  dmSetScalars_d(sigs_col, sigs, 1, npoints);
  dmSetScalars_d(nmodes_col, nmodes, 1, npoints);

  dmSetScalars_d(lmin_col, loc_min, 1, npoints );
  dmSetScalars_d(lmax_col, loc_max, 1, npoints );
  dmSetScalars_d(lmed_col, loc_med, 1, npoints );
  
  return(0);
}

long find_value( double *vals, long nvals, double which ) 
{
  long retval;
  long ii;

  for(ii=0;ii<nvals;ii++) {
    if ( vals[ii] == which ) break;
  }

  if ( ii == nvals ) 
    retval=-1;
  else
    retval=(nvals-ii); /*nvals because vals are ordered from top to bottom */

  return(retval);
}






int dmimgproject()
{

  char infile[DS_SZ_PATHNAME];
  char outfile[DS_SZ_PATHNAME];
  char paxis[10];
  short clobber;
  short verbose;

  void *data;
  long *lAxes;
  regRegion *dss;
  long null;
  short has_null;

  long ii, jj;
  double *counts;
  double *area;
  double *norm;
  double *grid;
  double *mins;
  double *maxs;
  double *modes;
  double *ranges;
  double *medians;
  double *q25;
  double *q33;
  double *q67;
  double *q75;
  double *sigs;
  /* double *extremes */
  double *nmodes;

  double *loc_min;
  double *loc_max;
  double *loc_med;



  long *bin;

  double *vals;

  long *xx, *yy;

  dmDataType dt;
  dmBlock *inBlock;
  dmDescriptor *xdesc, *ydesc;
  
  dmBlock *outBlock;
  
  char proj[10];

  Projection_Type indi, sum;


  /* Get the parameters */
  clgetstr( "infile", infile, DS_SZ_FNAME );
  clgetstr( "outfile", outfile, DS_SZ_FNAME );
  clgetstr( "axis", paxis, 10 );
  clobber = clgetb( "clobber" );
  verbose = clgeti( "verbose" );


  
  if ( NULL == ( inBlock = dmImageOpen( infile) ) ) {
    err_msg("ERROR: Cannot open image '%s'\n", infile );
    return(-1);
  }

  if ( dmUNKNOWNTYPE == ( dt = get_image_data( inBlock, &data,  &lAxes, 
					       &dss, &null, &has_null ) ) ) {
    err_msg("ERROR: Cannot get image data or unknown image data-type for "
	    "file '%s'\n", infile);
    return(-1);
  }

  if ( 0 != get_image_wcs( inBlock, &xdesc, &ydesc ) ) {
    err_msg("ERROR: Cannot load WCS for file '%s'\n", infile );
    return(-1);
  }

  if ( ds_clobber( outfile, clobber, NULL ) != 0 ) {
    return(-1);
  }

  if ( *paxis  == 'x' ) {
    indi = PROJ_X;  /* indi is the independent axis, get the x-axis */
    sum = PROJ_Y;   /* sum is the axis that is being summed */
    xx = &ii; /* makes things go a little faster */
    yy = &jj; /* rather than having another if in the for loops */

  } else {
    indi = PROJ_Y;
    sum = PROJ_X;
    xx = &jj;
    yy = &ii;
  }
  

  /* Allocate all the memeory */
  counts = (double*)calloc( lAxes[indi], sizeof(double));
  area = (double*)calloc( lAxes[indi], sizeof(double));
  norm = (double*)calloc( lAxes[indi], sizeof(double));
  grid = (double*)calloc( lAxes[indi], sizeof(double));
  bin = (long*)calloc(lAxes[indi],sizeof(long));


  mins= (double*)calloc( lAxes[indi], sizeof(double));
  maxs = (double*)calloc( lAxes[indi], sizeof(double));
  modes = (double*)calloc( lAxes[indi], sizeof(double));
  ranges = (double*)calloc( lAxes[indi], sizeof(double));
  medians = (double*)calloc( lAxes[indi], sizeof(double));
  q25 = (double*)calloc( lAxes[indi], sizeof(double));
  q33 = (double*)calloc( lAxes[indi], sizeof(double));
  q67 = (double*)calloc( lAxes[indi], sizeof(double));
  q75 = (double*)calloc( lAxes[indi], sizeof(double));
  sigs = (double*)calloc( lAxes[indi], sizeof(double));
  nmodes = (double*)calloc( lAxes[indi], sizeof(double));

  loc_min = (double*)calloc( lAxes[indi], sizeof(double));
  loc_max = (double*)calloc( lAxes[indi], sizeof(double));
  loc_med = (double*)calloc( lAxes[indi], sizeof(double));


  vals = (double*)calloc(lAxes[sum],sizeof(double));

  /* Loop through all the data */

  for (ii=lAxes[indi];ii--;) {
    long nvals;
    nvals = 0;
    for (jj=lAxes[sum];jj--;) {
      double dat;
      
      dat = get_image_value( data, dt, lAxes, *xx, *yy,
			     dss, xdesc, ydesc, null, has_null );
      
      if ( ds_dNAN( dat ) ) {
	/* do nothing */
      } else {
	vals[nvals]=dat;
	nvals++;

	/* counts[ii] += dat; */
	/* area[ii] +=1;*/
      }

      
    } /* end for jj */

    if ( nvals > 0 ) {
      area[ii]   = _filtCOUNT( vals, nvals );
      counts[ii] = _filtSUM( vals, nvals );
      norm[ii]   = _filtMEAN( vals, nvals );    
      mins[ii]   = _filtMIN( vals, nvals );
      maxs[ii]   = _filtMAX( vals, nvals );
      modes[ii]  = _filtMODE( vals, nvals );
      ranges[ii] = _filtRANGE(vals, nvals );
      medians[ii] = _filtMEDIAN(vals, nvals );
      q25[ii] = _filtQUANTILE_25( vals, nvals );
      q33[ii] = _filtQUANTILE_33( vals, nvals );
      q67[ii] = _filtQUANTILE_67( vals, nvals );
      q75[ii] = _filtQUANTILE_75( vals, nvals );
      sigs[ii] = _filtSIG( vals, nvals );
      nmodes[ii] = _filtNORM_MODE(vals, nvals );

      loc_min[ii]=find_value( vals, nvals, mins[ii]);
      loc_max[ii]=find_value( vals, nvals, maxs[ii]);
      loc_med[ii]=find_value( vals, nvals, medians[ii]);

    } else {
      area[ii] = 0;
      ds_MAKE_DNAN(( counts[ii] ));
      ds_MAKE_DNAN(( norm[ii] ));
      ds_MAKE_DNAN(( modes[ii] ));
      ds_MAKE_DNAN(( ranges[ii] ));
      ds_MAKE_DNAN(( medians[ii] ));
      ds_MAKE_DNAN(( q25[ii] ));
      ds_MAKE_DNAN(( q33[ii] ));
      ds_MAKE_DNAN(( q67[ii] ));
      ds_MAKE_DNAN(( q75[ii] ));
      ds_MAKE_DNAN(( sigs[ii] ));
      ds_MAKE_DNAN(( nmodes[ii] ));
      ds_MAKE_DNAN(( mins[ii] ));
      ds_MAKE_DNAN(( maxs[ii] ));

    }
    bin[ii] = ii+1;


    /* Fill in the grid */
    if (xdesc) {
      double loc[2];
      double pos[2];

      pos[indi] = bin[ii];
      pos[sum] = 1;

      if (ydesc) {
	dmCoordCalc_d( xdesc, pos, loc );
	dmCoordCalc_d( ydesc, pos+1, loc+1 );
      } else {
	dmCoordCalc_d( xdesc, pos, loc );
      }
      grid[ii] = loc[indi];

      /* repeat for loc's */
      pos[indi] = 1;
      pos[sum] = loc_min[ii];
      if (ydesc) {
	dmCoordCalc_d( xdesc, pos, loc );
	dmCoordCalc_d( ydesc, pos+1, loc+1 );
      } else {
	dmCoordCalc_d( xdesc, pos, loc );
      }
      loc_min[ii] = loc[sum];

      /* repeat for loc's */
      pos[indi] = 1;
      pos[sum] = loc_max[ii];
      if (ydesc) {
	dmCoordCalc_d( xdesc, pos, loc );
	dmCoordCalc_d( ydesc, pos+1, loc+1 );
      } else {
	dmCoordCalc_d( xdesc, pos, loc );
      }
      loc_max[ii] = loc[sum];

      /* repeat for loc's */
      pos[indi] = 1;
      pos[sum] = loc_med[ii];
      if (ydesc) {
	dmCoordCalc_d( xdesc, pos, loc );
	dmCoordCalc_d( ydesc, pos+1, loc+1 );
      } else {
	dmCoordCalc_d( xdesc, pos, loc );
      }
      loc_med[ii] = loc[sum];


    } else {
      grid[ii] = bin[ii];
    }

  } /* end for ii */


  if (NULL == (outBlock = dmTableCreate( outfile ) )) {

    err_msg("ERROR: Cannot create output file '%s'\n", outfile );
    return(-1);
  }
  if ( 0 != write_outfile( inBlock, outBlock, bin, grid, counts,
			   area, norm, mins, maxs, modes,
			   ranges, medians, q25, q33, q67, q75,
			   sigs, nmodes, loc_min, loc_max, loc_med,
			   lAxes[indi], indi ) ) {
    err_msg("ERROR: Problems writing to output file '%s'\n", outfile );
    return(-1);
  }


  dmTableClose(outBlock);

  /* Note: left all allocs unfreed and inBlock open, just going
     to free when program exits anyway, save a little time */

  return(0);

}








