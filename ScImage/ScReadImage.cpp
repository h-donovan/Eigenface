#include "ScReadImage.h"

int ScReadImage::readPPMImage(const char fname[], ScImage& image)
{
	int x, y, depth, bpp;
	unsigned char* charImage;
	char header[100], * ptr;
	ifstream ifp;

	ifp.open(fname, ios::in | ios::binary);

	if (!ifp)
	{
		cout << "Can't read image: " << fname << endl;
		exit(1);
	}

	ifp.getline(header, 100, '\n');
	if ((header[0] != 80) ||	// 'P'
		header[1] != 54) {		// '6'
		cout << "Image " << fname << " cannot be read as PPM" << endl;
		exit(1);
	}

	ifp.getline(header, 100, '\n');
	while (header[0] == '#')
		ifp.getline(header, 100, '\n');
	
	x = strtol(header, &ptr, 0);
	y = atoi(ptr);

	ifp.getline(header, 100, '\n');
	depth = strtol(header, &ptr, 0);

	bpp = 3;

	image.setImageData(x, y, bpp, depth);
	cout << "Image size: " << x << ":" << y << endl;

	charImage = (unsigned char*) new unsigned char[x * y * bpp];

	ifp.read(reinterpret_cast<char*>(charImage), (x * y * bpp) * sizeof(unsigned char));

	if (ifp.fail())
	{
		cout << "Image " << fname << " has wrong size." << endl;
		exit(1);
	}

	ifp.close();

	int val;

	for (int j = 0; j < y; j++)
	{
		for (int i = 0; i < x; i++)
		{
			for (int k = 0; k < bpp; k++)
			{
				val = (int)charImage[ (j * x * bpp) + (i * bpp) + k];
				image.setPixel(i, j, k, val);
			}
		}
	}

	delete[] charImage;

	return 1;
}
