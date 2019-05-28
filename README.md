# Introduction

A ray tracing program, using the open source library OpenEXR (developed by IL&M), allow using specific input or existing OBJ files to create rendered image.<br>

# Purpose

This is the work-through class project of COMS W4160 Computer Graphics, taught by prof. Michael Reed. This is also my first project in computer science area, so as you can see from the code, it's not very organized, and did not follow the google code style. Also, I didn't use the default inherit logic. <br>

Please don't copy the code for any reason.

# Environment

This program uses OpenEXR library. Before running codes in your computer, you have to setup OpenEXR and IlmBase, as well as necessary .h files first.<br>
Mac OS and Linux are recommended. Environment set-up under windows is extremely difficult, but could work (with no further version above Visual Studio 2013).<br>
For more information, visit http://www.openexr.com/ <br>

# Realization

This specific ray-tracing renderer realizes functions include:<br>
>Triangle and Sphere shading<br>
>Single point lights and area lights<br>
>Diffuse contribution and specular contribution of lights<br>
>Shadows<br>
>True reflections<br>
>Fast graphic processing - Bounding Box and BVH tree<br>
>Monte-carlo method and soft shadows<br>
>Readin from standard triangle-defined OBJ files<br>

In the final version of my code, fast graphic processing (implement of BVH tree) is by default. It would make the image rendering much faster; however, you can choose not to use it by changing specific areas in 'raytra.cc' and 'camera.cc'.

# Structure

sphere.h<br>
sphere.cc<br>
> Defined storage and behaviors of spheres, as well as its intersection check<br>

triangle.cc<br>
triangle.cc<br>
> Defined storage and behaviors of triangles, as well as its intersection check<br>

bbox.h<br>
bbox.cc<br>
> Defined behaviors of bounding boxes, includes their initialization, combination, and calculation of whether intersecting or not<br>

material.h<br>
material.cc<br>
> Defined specific materials, includes lumen of lights, coefficients of diffuse reflection, specular reflection, and true reflection, for all spheres & triangles<br>

ray.h<br>
> Defined usage of rays, include initial ray, reflect ray, and shadow ray<br>

light.h<br>
> Defined 3 kinds of lights - point light, area light, and background (ambient) light<br>

camera.cc<br>
camera.h<br>
> The observation of whole image, include its processing and rendering<br>
> The most important part of image rendering<br>
>> Tip: only one camera should by added<br>

basicmath.h<br>
> Defined 'vector' and 'point' as two basic element in every images<br>

parser.h<br>
parser.cc<br>
> Image-definition file processing<br>

raytra.cc<br>
> main

# Sample input and output

input: ./raytra chess.scn chess.exr 2 2 <br>
output: chess.exr
