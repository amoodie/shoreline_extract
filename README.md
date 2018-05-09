# Shoreline_extract

Matlab scripts for extracting subaerial delta shoreline from a collection of Landsat images.

Simplified demonstration version of the code used in Moodie et al. to make the shoreline traces and determine the rates of delta development of the Yellow River delta.


## Dependencies

This code was written with Matlab 2017b on Ubuntu 16.04. 
It has intentionally been written to be OS agnostic, but does depend on Matlab 2017 or newer.

## Limitations

This code was written to process Landsat images of the Yellow River delta, China.
While most of the code is agnostic of the image/location it is processing, there are some operations that are specific to the location. 

These limitations mainly occur in the context of image cropping and shoreline extraction:

* the limits of cropping are explicitly defined in the code
* a padding operation assumes land is in the lower left of the image
* some size parameters in shoreline extraction may need to be tweeked

The code will probably not work out of the box when applied to another system because of these limitations, __BUT__ with a minimal amount of work (the authors are happy to help!) it could be made to work.


## Using the program

### Obtain the source code

The code should be downloaded by cloning this repository. 
__WARNING:__ the default branch of this repository is almost 1 GB in size, because it contains 4 Landsat Level 1 products (though these have also been stripped to only necessary parts to reduce their size).

This default branch can be cloned with:

```
git clone https://github.com/amoodie/shoreline_extract.git
```

A branch of the repository without the data is also available, for a reduced file size or if you wish to use your own data with the codebase.
This version can be cloned with:

```
git clone
```

### File/folder structure

A tree of the default branch, simplified for brevity, is produced below
```
.
├── data/
│   ├── empty.txt
│   └── qinshuigou_channelline.csv
├── output/
├── private/
├── source/
│   ├── build_shorelineset.m
│   └── explore_shorelineset.m
├── LICENSE.txt
└── README.md
```

The main functions live in the `source` folder.
This folder contains two files:

* `build_shorelineset.m`
* `explore_shorelineset.m`

These files make up the application, and should be run sequentially. 
The `build_shorelineset.m` file processes the raw data into a collection of data files that contain xy coordinates of the shoreline.
Then, `explore_shorelineset.m` manipulates the extracted shorelines to determine usable data from the shoreline trajectories.

The data used in developing the script are provided as part of this repository, but you can also obtain them (or others) for yourself from the Level 1 Landsat data at the [USGS Earth Explorer website](https://earthexplorer.usgs.gov/).
There were four scenes used to make the images and plots in this readme:

* LT51210341985329HAJ00
* LT41210341989028XXX02
* LT51210341992317HAJ00
* LT51210341995357CLT00

You can download these scenes explicitly, or any other scenes that fall into the `WRS_PATH = 121` and `WRS_ROW = 034` flight path from any time.

These files should be placed into the `data` folder as uncompressed folders (see the file tree above for where the folders should be located).
Actually, the folders could be placed anywhere, if the `meta.directory` parameter in `build_shorelineset.m` is changed to the appropriate path.




## License



## Disclaimer and acknowledgments 