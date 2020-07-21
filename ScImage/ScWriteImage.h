// Adapted from the files of Dr. George Bebis
// Utilized for writing PPM files from a known image format
// Major change: added a third dimension to the writing for RGB values.

#pragma once
#include "ScImage.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>

#ifndef SCWRITEIMAGE_H
#define SCWRITEIMAGE_H

using namespace std;

namespace ScWriteImage
{
	int writePPMImage(const char fname[], ScImage& image);
};

#endif

