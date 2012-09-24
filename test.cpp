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

void testRotation()
{
	
		CFloatImage matrixImage = GetImageFromMatrix((float *)featureMatrix, 10, 10);
		CTransform3x3 translationNegative;
		CTransform3x3 translationPositive;
		CTransform3x3 rotation;
		CFloatImage postHomography;

		Feature f;
		f.x = 6;
		f.y = 5;
		f.angleRadians = PI;

		translationNegative = translationNegative.Translation(f.x,f.y);
		translationPositive = translationPositive.Translation(-f.x,-f.y);

		rotation = rotation.Rotation(-f.angleRadians * 180/ PI);


		WarpGlobal(matrixImage, postHomography, translationNegative*rotation*translationPositive, eWarpInterpLinear, eWarpInterpNearest);
		for (int i = 0; i < postHomography.Shape().height; i++)
		{
			for (int j = 0; j < postHomography.Shape().width; j++)
			{
				printf("%.0f\t", postHomography.Pixel(j, i, 0));
			}
			printf("\n");
		}
}

void testGetWindow()
{
	CFloatImage matrixImage = GetImageFromMatrix((float *)featureMatrix, 10, 10);

	CFloatImage sample5x5 = GetXWindowAroundPixel(matrixImage, 5, 4, 5);

	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			printf("%f\t", sample5x5.Pixel(j, i, 0));
		}
		printf("\n");
	}

}