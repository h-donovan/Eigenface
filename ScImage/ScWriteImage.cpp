#include "ScWriteImage.h"

int ScWriteImage::writePPMImage(const char fname[], ScImage& image)
{
	int x, y, bpp, depth;
	unsigned char* charImage;
	ofstream ofp;

	image.getImageData(x, y, bpp, depth);

	charImage = (unsigned char*) new unsigned char[x * y * bpp];

	int val;

	for (int j = 0; j < y; j++)
	{
		for (int i = 0; i < x; i++)
		{
			for (int k = 0; k < bpp; k++)
			{
				val = image.getPixel(i, j, k);
				charImage[(j * x * bpp) + (i * bpp) + k] = (unsigned char)val;
			}
		}
	}

	ofp.open(fname, ios::out | ios::binary);

	if (!ofp) {
		cout << "Can't open file: " << fname << endl;
		exit(1);
	}

	ofp << "P6" << endl;
	ofp << x << " " << y << endl;
	ofp << depth << endl;

	ofp.write(reinterpret_cast<char*>(charImage), (x * y * bpp) * sizeof(unsigned char));

	if (ofp.fail())
	{
		cout << "Can't write image " << fname << endl;
		exit(0);
	}

	ofp.close();

	delete[] charImage;

	return 1;
}
