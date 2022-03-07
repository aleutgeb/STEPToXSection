# STEPToXSection
The program STEPToXSection extracts the contour of a planar cross section of solids contained in STEP files. Additionally it supports offsetting of the cross section contour.

# Description
The program STEPToXSection is a command line utility to export the contour of a planar cross section of solids contained in STEP files. Additionally it supports offsetting of the cross section contour. The contour is a list of line segments. The supported output file formats are ply and xyz. In the case of xyz, two consecutive vertices belong to the same edge. A popular viewer for the different file types is MeshLab (https://www.meshlab.net). STEPToXSection is based on OpenCASCADE (https://www.opencascade.com). The program uses cxxops (https://github.com/jarro2783/cxxopts) for parsing the command line.

The following examples were generated using the examples spheres.stp and bone_bocket.stp with the same values for the plane (-p) and deflection (-d). For example spheres.stp the offset values range from -5.0 to +10.0 and for example bone_bocket.stp the offset values range from -4.0 to +5.0. Example calls:
```
STEPToXSection.exe -i spheres.stp -o out.ply -f ply -p 1.0,0.0,0.0,0.0 -d 0.01 -t 1.0
STEPToXSection.exe -i bone_bocket.stp -o out.ply -f ply -p 1.0,0.0,0.0,0.0 -d 0.01 -t 1.0
```

| Solids (example spheres) |
| :---: |
| ![Image Solids-Spheres](examples/spheres/solids.png) |

| Planar cross-section (example spheres) |
| :---: |
| ![Image Cross-Section-Spheres](examples/spheres/cross_section.png) |

| Positive offset curves (example spheres) |
| :---: |
| ![Image Positive-Offset-Curves-Spheres](examples/spheres/positive_offset_curves.png) |

| Negative offset curves (example spheres) |
| :---: |
| ![Image Negative-Offset-Curves-Spheres](examples/spheres/negative_offset_curves.png) |

| Solid (example bone_pocket) |
| :---: |
| ![Image Solid-Bone-Pocket](examples/bone_pocket/solid.png) |

| Planar cross-section (example bone_pocket) |
| :---: |
| ![Image Cross-Section-Bone-Pocket](examples/bone_pocket/cross_section.png) |

| Positive offset curves (example bone_pocket) |
| :---: |
| ![Image Positive-Offset-Curves-Bone-Pocket](examples/bone_pocket/positive_offset_curves.png) |

| Negative offset curves (example bone_pocket) |
| :---: |
| ![Image Negative-Offset-Curves-Bone-Pocket](examples/bone_pocket/negative_offset_curves.png) |

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
Extracts the contour of a planar cross section of solids contained in the STEP file. Additionally it supports offsetting of the cross section contour.
Usage:
  STEPToXSection [OPTION...]

  -i, --in arg          Input file
  -o, --out arg         Output file
  -f, --format arg      Output file format (xyz or ply) (default: xyz)
  -c, --content         List content (solids)
  -s, --select arg      Select solids by name or index (comma seperated list, index starts with 1)
  -d, --deflection arg  deflection
  -t, --offset arg      offset (default: 0.0)
  -p, --plane arg       Plane (a,b,c,d), in which a*x + b*y + c*z + d = 0
  -h, --help            Print usage
```

# Examples
 * See `examples` directory
 
# Remarks
This code has been tested with an OpenCASCADE 7.5.0 prebuilt binary (`opencascade-7.5.0-vc14-64.exe`) on Windows, as well as OpenCASCADE system packages on openSUSE Linux. With changes in the configuration section in the `CMakeLists.txt` file the build should also work with other OpenCASCADE versions.