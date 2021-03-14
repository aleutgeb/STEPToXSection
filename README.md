# STEPToXSection
The program STEPToXSection extracts the contour of a planar cross section of solids contained in STEP files.

# Description
The program STEPToXSection is a command line utility to export the contour of a planar cross section of solids contained in STEP files. The contour is a list of line segments. The supported output file formats are ply and xyz. In the case of xyz, two consecutive vertices belong to the same edge.  STEPToXSection is based on OpenCASCADE (https://www.opencascade.com). The program uses cxxops (https://github.com/jarro2783/cxxopts) for parsing the command line.

# Requirements
 * CMake installation (https://cmake.org)
 * Visual Studio C++ installation (https://visualstudio.microsoft.com)
 * OpenCASCADE installation (https://old.opencascade.com/content/latest-release, download needs registration)

# Usage
Listing the contents (solids) of a STEP file:
`STEPToXSection -c -i <step file>`

Computing the planar cross section contour of the overall file content (solids):

`STEPToXSection -i <step file> -o <output file> -d <deflection> -p <plane>`

The parameter `<deflection>` controls the resolution of the approximation with line segments.

Computing the planar cross section contour of selected solids of the file:

`STEPToXSection -i <step file> -o <output file> -d <deflection> -p <plane> -s <solid1>,<solid2>,<...>`

In order to change the default output format xyz the command line argument `-f ply`has to be specified.

Following the help text from the command line:
```
STEPToXSection.exe -h
Extracts the contour of a planar cross section of solids contained in the STEP file
Usage:
  STEPToXSection [OPTION...]

  -i, --in arg          Input file
  -o, --out arg         Output file
  -f, --format arg      Output file format (xyz or ply) (default: xyz)
  -c, --content         List content (solids)
  -s, --select arg      Select solids by name (comma seperated list)
  -d, --deflection arg  deflection
  -p, --plane arg       Plane (a,b,c,d), in which a*x + b*y + c*z + d = 0
  -h, --help            Print usage
```

# Examples
 * See `examples` directory
 
# Remarks
This code has been tested with an OpenCASCADE 7.4.0 prebuilt binary (`opencascade-7.4.0-vc14-64.exe`). With changes in the configuration section in the `CMakeLists.txt` file the build should also work with other OpenCASCADE versions. For building on Linux the `CMakeLists.txt` file needs further adaptions.
