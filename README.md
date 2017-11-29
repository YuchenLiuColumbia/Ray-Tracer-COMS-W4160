# Introduction

A ray tracing program, using the open source library OpenEXR (developed by IL&M), allow using specific input or existing OBJ files to create rendered image

# Environment

This program uses OpenEXR library. Before running codes in your computer, you have to setup OpenEXR and IlmBase, as well as necessary .h files first.
Mac OS and Linux are recommended. Environment set-up under windows is extremely difficult, but could work (with no further version above Visual Studio 2013).
For more information, visit http://www.openexr.com/

# Realization

This specific ray-tracing renderer realizes functions include:
// Triangle and Sphere shading
// Single point lights and area lights
// Diffuse contribution and specular contribution of lights
// Shadows
// True reflections
// Fast graphic processing - Bounding Box and BVH tree
// Monte-carlo method and soft shadows
// Readin from standard triangle-defined OBJ files
In the final version of my code, fast graphic processing (implement of BVH tree) is by default. It would make the image rendering much faster; however, you can choose not to use it by changing specific areas in 'raytra.cc' and 'camera.cc'.

# Structure

sphere.h
sphere.cc
// Defined storage and behaviors of spheres, as well as its intersection check
triangle.cc
triangle.cc
// Defined storage and behaviors of triangles, as well as its intersection check
bbox.h
bbox.cc
// Defined behaviors of bounding boxes, includes their initialization, combination, and calculation of whether intersecting or not
material.h
material.cc
// Defined specific materials, includes lumen of lights, coefficients of diffuse reflection, specular reflection, and true reflection, for all spheres & triangles
ray.h
// Defined usage of rays, include initial ray, reflect ray, and shadow ray
light.h
// Defined 3 kinds of lights - point light, area light, and background (ambient) light
camera.cc
camera.h
// The observation of whole image, include its processing and rendering
// The most important part of image rendering
// Tip: only one camera should by added
basicmath.h
// Defined 'vector' and 'point' as two basic element in every images
parser.h
parser.cc
// Image-definition file processing
raytra.cc
// main

# Sample input and output

input: chess.scn
output: chess.exr
