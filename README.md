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


## Using the code

```
.
├── data/
│   ├── empty.txt
│   ├── LT41210341989028XXX02/
│   ├── LT51210341985329HAJ00/
│   ├── LT51210341992317HAJ00/
│   ├── LT51210341995357CLT00/
│   └── qinshuigou_channelline.csv
├── LICENSE.txt
├── output/
│   ├── meta_1985-11-25.mat
│   ├── meta_1989-01-28.mat
│   ├── meta_1992-11-12.mat
│   ├── meta_1995-12-23.mat
│   ├── shoreline_1985-11-25.csv
│   ├── shoreline_1985-11-25.mat
│   ├── shoreline_1989-01-28.csv
│   ├── shoreline_1989-01-28.mat
│   ├── shoreline_1992-11-12.csv
│   ├── shoreline_1992-11-12.mat
│   ├── shoreline_1995-12-23.csv
│   └── shoreline_1995-12-23.mat
├── private/
├── README.md
└── source/
    ├── build_shorelineset.m
    └── explore_shorelineset.m
```

The main functions live in the `source` folder. 
This folder contains two files:

* `build_shorelineset.m`
* `explore_shorelineset.m`

These files are the bulk of the application, and should be run sequentially in the order described. 
The `build_shorelineset.m` file processes the raw data into a collection of data files that contain xy coordinates of the shoreline.

__The Github repository does not contain the data used to run the script!!__
You will need to download the Level 1 Landsat data from the [USGS Earth Explorer website](https://earthexplorer.usgs.gov/).
There were four scenes used to make the images and plots in this readme:

* LT51210341985329HAJ00
* LT41210341989028XXX02
* LT51210341992317HAJ00
* LT51210341995357CLT00

You can download these scenes explicitly, or any other scenes that fall into the `WRS_PATH = 121` and `WRS_ROW = 034` flight path from any time.

These files should be placed into the `data` folder as uncompressed folders (see the file tree above for where the folders should be located).
Actually, the folders could be placed anywhere, if the `build_shorelineset.m` code parameter `meta.directory` is changed.



## License



## Disclaimer and acknowledgments 