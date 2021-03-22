# `dmimgproject2`

This is a re-write of the `dmimgproject` tool which projects 
the pixels in an image onto an arbitrary axis.  The original tool only
projected onto either the X or Y axis.

```bash
dmimgproject2 in_image.fits out_table.fits angle=0
```

## The output

The output file contains two blocks of data.

### The projection results

The first block is a TABLE with the projection results.

The X and Y columns are the coordinates of the X-axis (so Y=1) rotated by the
specified angle.

The remaining columns are different statistics computed from the values
along the projection. This includes

- sum
- mean (average)
- count (number of pixels)
- min (min pixel value)
- max
- median
- sigma (standard deviation)
- and several more.

The filter functions are described
[in the dmimgfilt ahelp file](https://cxc.harvard.edu/ciao/ahelp/dmimgfilt.html#Filter_functions)
(not all filters described are included in the dmimgreproject2 output).

### The output grid

The 2nd block contains an image showing which pixels are included in
each projection grid.  The pixel values match the BIN value in the 
table.

Pixels which are NaN or NULL in the input image are skipped.  They are
set to the NULL value (-999).

## Installation

This tool requires [CIAO 4.13](https://cxc.cfa.harvard.edu/ciao).

```bash

for d in dmimgio dmfilters dmimgproject2
do
  git clone https://github.com/kglotfelty/${d}
  (cd $d && ./configure --prefix=$ASCDS_INSTALL && make && make install)
done

```

You may optionally run the `dmimgproject2` test suite before installing:

```bash
cd dmimgproject
make check
```


