#include <assert.h>
#include <math.h>
#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include "features.h"
#include "ImageLib/FileIO.h"
#include "test.h"

void testGetFloatFromMatrix()
{
	CImageOf<double> img = GetImageFromMatrix((double *)sobelX, 3, 3);
	bool areSame = true;	
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			float imgPixel = (float)img.Pixel(j, i, 0);
			float sobelPixel = (float)sobelX[i * 3 + j];

			if (img.Pixel(j, i, 0) != sobelX[i * 3 + j])
			{
				areSame = false;
			}
		}
	}

	if (areSame)
		printf("Same\n");
	else
		printf("Not the same\n");
}

void testOrientation()
{
	float Am[] = {
		2.
	};
	float Bm[] = {
		1.
	};
	float Cm[] = {
		2.
	};
	CFloatImage A = GetImageFromMatrix(Am, 1, 1);
	CFloatImage B = GetImageFromMatrix(Bm, 1, 1);
	CFloatImage C = GetImageFromMatrix(Cm, 1, 1);

	CFloatImage partialX = GetImageFromMatrix(Bm, 1, 1);
	CFloatImage partialY = GetImageFromMatrix(Bm, 1, 1);

	double angle = GetCanonicalOrientation(0, 0, A, B, C, partialX, partialY);
	printf("Angle: %f\n", angle);
}

void testGet41x41()
{
}