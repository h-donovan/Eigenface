// Adapted from the files of Dr. George Bebis
// Utilized for writing PPM files from a known image format
// Major change: added a third dimension to the writing for RGB values.

#pragma once
#include "ScImage.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>

#ifndef SCREADIMAGE_H
#define SCREADIMAGE_H

using namespace std;

namespace ScReadImage
{
	int readPPMImage(const char fname[], ScImage& image);
};

#endif