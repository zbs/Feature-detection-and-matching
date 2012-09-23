#include <assert.h>
#include <math.h>
#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include "features.h"
#include "ImageLib/FileIO.h"

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
		printf("Same");
	else
		printf("Not the same");
}