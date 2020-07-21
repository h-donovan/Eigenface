// Adapted from the files of Dr. George Bebis
// Utilized for writing PPM files from a known image format
// Major change: added a third dimension to the writing for RGB values.

#pragma once
#ifndef SCIMAGE_H
#define SCIMAGE_H

class ScImage
{
	int resolution_x;
	int resolution_y;
	int pixel_depth;
	int bytes_per_pixel;

	int*** pixels;

public:

	ScImage();

	ScImage(int x, int y, int depth, int bpp);

	ScImage(ScImage& otherimage);

	~ScImage();

	void setImageData(int x, int y, int bpp, int depth);

	void getImageData(int& x, int& y, int& bpp, int& depth);

	void setPixel(int x, int y, int depth, int val);

	int getPixel(int x, int y, int depth);

	int getPixel(int x, int y, char depth);
};

#endif
