#include "ScImage.h"

ScImage::ScImage()
{

	resolution_x = 0;
	resolution_y = 0;
	pixel_depth = 0;
	bytes_per_pixel = 0;

	pixels = nullptr;
}

ScImage::ScImage(int x, int y, int depth, int bpp)
{
	resolution_x = x;
	resolution_y = y;
	pixel_depth = depth;
	bytes_per_pixel = bpp;

	pixels = new int** [resolution_x];
	for (int i = 0; i < resolution_x; i++)
	{
		pixels[i] = new int* [resolution_y];
		for (int j = 0; j < resolution_y; j++)
		{
			pixels[i][j] = new int[bpp];
			for (int k = 0; k < bpp; k++)
			{
				pixels[i][j][k] = 0;
			}
		}
	}
}

ScImage::ScImage(ScImage& otherimage)
{
	resolution_x = otherimage.resolution_x;
	resolution_y = otherimage.resolution_y;
	pixel_depth = otherimage.pixel_depth;
	bytes_per_pixel = otherimage.bytes_per_pixel;

	pixels = new int** [resolution_x];
	for (int i = 0; i < resolution_x; i++)
	{
		pixels[i] = new int* [resolution_y];
		for (int j = 0; j < resolution_y; j++)
		{
			pixels[i][j] = new int[bytes_per_pixel];
			for (int k = 0; k < bytes_per_pixel; k++)
			{
				pixels[i][j][k] = otherimage.pixels[i][j][k];
			}
		}
	}
}

ScImage::~ScImage()
{
	if (!pixels)
		return;
	for (int i = 0; i < resolution_x; i++)
	{
		if (!pixels[i])
			break;
		for (int j = 0; j < resolution_y; j++)
		{
			if(pixels[i][j])
				delete[] pixels[i][j];
		}
		delete[] pixels[i];
	}
	delete[] pixels;
	
	pixels = nullptr;
}

void ScImage::setImageData(int x, int y, int bpp, int depth)
{
	int t = bytes_per_pixel;
	resolution_x = x;
	resolution_y = y;
	pixel_depth = depth;
	bytes_per_pixel = bpp;

	if ( t != bpp)
	{
		if (pixels != nullptr)
		{
			for (int i = 0; i < resolution_x; i++)
			{
				for (int j = 0; j < resolution_y; j++)
				{
					delete[] pixels[i][j];
				}
				delete[] pixels[i];
			}
			delete[] pixels;
		}

		pixels = new int** [resolution_x];
		for (int i = 0; i < resolution_x; i++)
		{
			pixels[i] = new int* [resolution_y];
			for (int j = 0; j < resolution_y; j++)
			{
				pixels[i][j] = new int[bytes_per_pixel];
				for (int k = 0; k < bytes_per_pixel; k++)
				{
					pixels[i][j][k] = 0;
				}
			}
		}
	}
}

void ScImage::getImageData(int& x, int& y, int& bpp, int& depth)
{
	x = resolution_x;
	y = resolution_y;
	depth = pixel_depth;
	bpp = bytes_per_pixel;
}

void ScImage::setPixel(int x, int y, int depth, int val)
{
	pixels[x][y][depth] = val;
}

int ScImage::getPixel(int x, int y, int depth)
{
	return pixels[x][y][depth];
}

int ScImage::getPixel(int x, int y, char depth)
{
	switch (depth)
	{
		case 'r':
		case 'R':
			return pixels[x][y][0];
		case 'g':
		case 'G':
			return depth > 1 ? pixels[x][y][1] : 0;
		case 'b':
		case 'B':
			return depth > 2 ? pixels[x][y][2] : 0;
		default:
			return 0;
	}
}