/* features.cpp */
//computeFeatures ./graf/img1.ppm ./graf/features.f 2 2
//computeFeatures ./Yosemite/Yosemite1.jpg ./Yosemite_1.f 2 2
//computeFeatures ./Yosemite/Yosemite2.jpg ./Yosemite_2.f 2 2
// Get intensity differences on a per-40X40-patch level vs on a entire-image level?
	//per-patch level
//roc ./graf_1.f ./graf_2.f ./graf/H1to2p 1 ./graf/roc1.txt ./graf/auc1.txt
//roc ./Yosemite_1.f ./Yosemite_2.f ./Yosemite/H1to2p 1 ./Yosemite/roc1.txt ./Yosemite/auc1.txt
// 
// change standard dev to something normal
//search for PROVISIONAL MEASURES in this code
#include <assert.h>
#include <math.h>
#include <hash_map>
#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include "features.h"
#include "ImageLib/FileIO.h"

#define PI 3.14159265358979323846

// Compute features of an image.
bool computeFeatures(CFloatImage &image, FeatureSet &features, int featureType, int descriptorType)
{
    // TODO: Instead of calling dummyComputeFeatures, implement
    // Harris feature detector.  This step fills in "features"
    // with information needed for descriptor computation.
    switch (featureType) {
    case 1:
        dummyComputeFeatures(image, features);
        break;
    case 2:
        ComputeHarrisFeatures(image, features);
        break;
    default:
        return false;
    }

    // TODO: You will implement two descriptors for this project
    // (see webpage).  This step fills in "features" with
    // descriptors.  The third "custom" descriptor is extra credit.
    switch (descriptorType) {
    case 1:
        ComputeSimpleDescriptors(image, features);
        break;
    case 2:
        ComputeMOPSDescriptors(image, features);
        break;
    case 3:
        ComputeCustomDescriptors(image, features);
        break;
    default:
        return false;
    }

    // This is just to make sure the IDs are assigned in order, because
    // the ID gets used to index into the feature array.
    for (unsigned int i=0; i<features.size(); i++) {
        features[i].id = i+1;
    }

    return true;
}

// Perform a query on the database.  This simply runs matchFeatures on
// each image in the database, and returns the feature set of the best
// matching image.
bool performQuery(const FeatureSet &f, const ImageDatabase &db, int &bestIndex, vector<FeatureMatch> &bestMatches, double &bestScore, int matchType) {
    // Here's a nice low number.
    bestScore = -1e100;

    vector<FeatureMatch> tempMatches;
    double tempScore;

    for (unsigned int i=0; i<db.size(); i++) {
        if (!matchFeatures(f, db[i].features, tempMatches, tempScore, matchType)) {
            return false;
        }

        if (tempScore > bestScore) {
            bestIndex = i;
            bestScore = tempScore;
            bestMatches = tempMatches;
        }
    }

    return true;
}

// Match one feature set with another.
bool matchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore, int matchType) {
    // TODO: We have given you the ssd matching function, you must write your own
    // feature matching function for the ratio test.
        
    printf("\nMatching features.......\n");

    switch (matchType) {
    case 1:
        ssdMatchFeatures(f1, f2, matches, totalScore);
        return true;
    case 2:
        ratioMatchFeatures(f1, f2, matches, totalScore);
        return true;
    default:
        return false;
    }
}

// Evaluate a match using a ground truth homography.  This computes the
// average SSD distance between the matched feature points and
// the actual transformed positions.
double evaluateMatch(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9]) {
    double d = 0;
    int n = 0;

    double xNew;
    double yNew;

    unsigned int num_matches = matches.size();
    for (unsigned int i=0; i<num_matches; i++) {
        int id1 = matches[i].id1;
        int id2 = matches[i].id2;
        applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);
        d += sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
        n++;
    }       

    return d / n;
}

void addRocData(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9],vector<bool> &isMatch,double threshold,double &maxD) {
    double d = 0;

    double xNew;
    double yNew;

    unsigned int num_matches = matches.size();
    for (unsigned int i=0; i<num_matches; i++) {
        int id1 = matches[i].id1;
        int id2 = matches[i].id2;
        applyHomography(f1[id1-1].x, f1[id1-1].y, xNew, yNew, h);

        // Ignore unmatched points.  There might be a better way to
        // handle this.
        d = sqrt(pow(xNew-f2[id2-1].x,2)+pow(yNew-f2[id2-1].y,2));
        if (d<=threshold)
            {
		isMatch.push_back(1);
            }
        else
            {
		isMatch.push_back(0);
            }

        if (matches[i].score>maxD)
            maxD=matches[i].score;
    }       
}

vector<ROCPoint> computeRocCurve(vector<FeatureMatch> &matches,
                                 vector<bool> &isMatch,
                                 vector<double> &thresholds)
{
    vector<ROCPoint> dataPoints;

    for (int i=0; i < (int)thresholds.size();i++) {
        //printf("Checking threshold: %lf.\r\n",thresholds[i]);
        int tp=0;
        int actualCorrect=0;
        int fp=0;
        int actualError=0;
        int total=0;

        int num_matches = (int) matches.size();
        for (int j=0;j < num_matches;j++)
            {
                if (isMatch[j])
                    {
                        actualCorrect++;
                        if (matches[j].score<thresholds[i])
                            {
                                tp++;
                            }
                    }
                else
                    {
                        actualError++;
                        if (matches[j].score<thresholds[i])
                            {
                                fp++;
                            }
                    }                           
                        
                total++;
            }

        ROCPoint newPoint;
        //printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);
        newPoint.trueRate=(double(tp)/actualCorrect);
        newPoint.falseRate=(double(fp)/actualError);
        //printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);

        dataPoints.push_back(newPoint);
    }

    return dataPoints;
}


// Compute silly example features.  This doesn't do anything
// meaningful.
void dummyComputeFeatures(CFloatImage &image, FeatureSet &features) {
    CShape sh = image.Shape();
    Feature f;

    for (int y=0; y<sh.height; y++) {
        for (int x=0; x<sh.width; x++) {
            double r = image.Pixel(x,y,0);
            double g = image.Pixel(x,y,1);
            double b = image.Pixel(x,y,2);

            if ((int)(255*(r+g+b)+0.5) % 100  == 1) {
		// If the pixel satisfies this meaningless criterion,
		// make it a feature.
                                
		f.type = 1;
		f.id += 1;
		f.x = x;
		f.y = y;
		f.angleRadians = 0; // default value
		features.push_back(f);
            }
        }
    }
}
void GetHarrisComponents(CFloatImage &srcImage, CFloatImage &A, CFloatImage &B, CFloatImage &C, CFloatImage *partialX, CFloatImage *partialY)
{
	int w = srcImage.Shape().width;
    int h = srcImage.Shape().height;

	CFloatImage *partialXPtr;
	CFloatImage *partialYPtr;

	if (partialX != nullptr && partialY != nullptr)
	{
		partialXPtr = partialX;
		partialYPtr = partialY;
	}
	else
	{
		partialXPtr = new CFloatImage(srcImage.Shape());
		partialYPtr = new CFloatImage(srcImage.Shape());
	}

	CFloatImage partialXX(srcImage.Shape());
	CFloatImage partialYY(srcImage.Shape());
	CFloatImage partialXY(srcImage.Shape());

	CFloatImage gaussianImage = GetImageFromMatrix((float *)gaussian5x5Float, 5, 5);

	Convolve(srcImage, *partialXPtr, ConvolveKernel_SobelX);
	Convolve(srcImage, *partialYPtr, ConvolveKernel_SobelY);
	
	for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
			float *xxPixel = &partialXX.Pixel(x, y, 0);
			float *yyPixel = &partialYY.Pixel(x, y, 0);
			float *xyPixel = &partialXY.Pixel(x, y, 0);
			
			// The 1/8 factor is to do the scaling inherent in sobel filtering
			*xxPixel = pow((double)(1./8. *8. * partialXPtr->Pixel(x, y, 0)), 2.);
			*yyPixel = pow((double)(1./8. *8. * partialYPtr->Pixel(x, y, 0)), 2.);
			*xyPixel = pow(1./8. *8., 2.) * partialXPtr->Pixel(x, y, 0) * partialYPtr->Pixel(x, y, 0);
		}
	}

	Convolve(partialXX, A, gaussianImage);
	Convolve(partialXY, B, gaussianImage);
	Convolve(partialYY, C, gaussianImage);
}

double GetCanonicalOrientation(int x, int y, CFloatImage A, CFloatImage B, CFloatImage C, CFloatImage partialX, CFloatImage partialY)
{
	float aPixel = A.Pixel(x, y, 0);	
	float bPixel = B.Pixel(x, y, 0);	
	float cPixel = C.Pixel(x, y, 0);	

	/*double a = 1;
	double b = -(aPixel+cPixel);
	double c = (aPixel * cPixel - pow((double)bPixel, 2.));*/

	//double lambda = (- b + sqrt(pow((double)b, 2.) - 4.*a*c)) / (2.*a);

	double lambda = 1./2. *((aPixel+cPixel) + sqrt(4.*pow((double)bPixel,2.) + pow(((double)aPixel - cPixel), 2.)));
	double yComponent = aPixel - lambda - bPixel;
	double xComponent = cPixel - lambda - bPixel;
	/*y = -b
	x = a - lambd
	*/
	if (xComponent == 0.)
	{
		return (partialY.Pixel(x, y, 0) > 0)? PI/2. : -PI/2.;
	}

	double first = aPixel*xComponent + bPixel*yComponent;
	double second = bPixel*xComponent + cPixel*yComponent;
	double first_precision = pow((first -lambda * xComponent), 2.);
	double second_precision = pow((second - lambda * yComponent), 2.);

	//Sanity check for eigenvalues: this confirms our given eigenvalue is
	//computed correctly
	if (first_precision > .0000001 || second_precision > .0000001 )
	{
		int z = 3;
	}

	return (partialX.Pixel(x, y, 0) > 0)? atan(-bPixel/(aPixel-lambda)) : atan(-bPixel/(aPixel-lambda)) + PI;

	//return (partialX.Pixel(x, y, 0) > 0)? atan(yComponent/xComponent) : atan(yComponent/xComponent) + PI;
}
void ComputeHarrisFeatures(CFloatImage &image, FeatureSet &features)
{
    //Create grayscale image used for Harris detection
    CFloatImage grayImage=ConvertToGray(image);

    //Create image to store Harris values
    CFloatImage harrisImage(image.Shape().width,image.Shape().height,1);

    //Create image to store local maximum harris values as 1, other pixels 0
    CByteImage harrisMaxImage(image.Shape().width,image.Shape().height,1);

    //compute Harris values puts harris values at each pixel position in harrisImage. 
    //You'll need to implement this function.
    computeHarrisValues(grayImage, harrisImage);
        
    // Threshold the harris image and compute local maxima.  You'll need to implement this function.
    computeLocalMaxima(harrisImage,harrisMaxImage);

    // Prints out the harris image for debugging purposes
    CByteImage tmp(harrisImage.Shape());
    convertToByteImage(harrisImage, tmp);
    WriteFile(tmp, "harris.tga");

    // TO DO--------------------------------------------------------------------
    //Loop through feature points in harrisMaxImage and fill in information needed for 
    //descriptor computation for each point above a threshold. We fill in id, type, 
    //x, y, and angle.
	CFloatImage A(grayImage.Shape());
	CFloatImage B(grayImage.Shape());
	CFloatImage C(grayImage.Shape());

	CFloatImage partialX(grayImage.Shape());
	CFloatImage partialY(grayImage.Shape());

	GetHarrisComponents(grayImage, A, B, C, &partialX, &partialY);
	int featureCount = 0;
    for (int y=0;y<harrisMaxImage.Shape().height;y++) {
        for (int x=0;x<harrisMaxImage.Shape().width;x++) {
                
            // Skip over non-maxima
            if (harrisMaxImage.Pixel(x, y, 0) == 0)
				continue;

            //TO DO---------------------------------------------------------------------
            // Fill in feature with descriptor data here. 
            Feature f;
			f.type = 2;
			f.id = featureCount++;
			f.x = x;
			f.y = y;
			f.angleRadians = GetCanonicalOrientation(x, y, A, B, C, partialX, partialY);
				//atan(partialY.Pixel(x, y, 0)/partialX.Pixel(x, y, 0));
			// Add the feature to the list of features
            features.push_back(f);
        }
    }
}

template <class T>
CImageOf<T> GetImageFromMatrix(T *matrix, int width, int height)
{
	// Allocate the new image
    CShape dShape(width, height, 1);
	CImageOf<T> dst(dShape);

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			// Allocate the new image
			T *pixel = &dst.Pixel(x,y,0);
			*pixel = matrix[y * width + x];
		}
	}
	return dst;
}

void test()
{
	GetImageFromMatrix((double *)sobelX, 3, 3);
}

//TO DO---------------------------------------------------------------------
//Loop through the image to compute the harris corner values as described in class
// srcImage:  grayscale of original image
// harrisImage:  populate the harris values per pixel in this image
void computeHarrisValues(CFloatImage &srcImage, CFloatImage &harrisImage)
{
	int h = srcImage.Shape().height;
	int w = srcImage.Shape().width;

	CFloatImage A(srcImage.Shape());
	CFloatImage B(srcImage.Shape());
	CFloatImage C(srcImage.Shape());

	GetHarrisComponents(srcImage, A, B, C);

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
			double determinant = A.Pixel(x, y, 0) * C.Pixel(x, y, 0) - B.Pixel(x, y, 0)* B.Pixel(x, y, 0);
			double trace = A.Pixel(x, y, 0) + C.Pixel(x, y, 0);
			
			float *pixel = &harrisImage.Pixel(x, y, 0);

			if (trace == 0)
			{
				*pixel = 0;
			}
			else
			{
				*pixel = determinant / trace;
			}
        }
    }
}

bool isLocalMax(CFloatImage srcImage, int x, int y)
{
	int width = srcImage.Shape().width;
	int height = srcImage.Shape().height;
	float centerPixel = srcImage.Pixel(x, y, 0);

	for (int row = 0; row < 5; row++)
	{
		for (int col = 0; col < 5; col++)
		{
			int xOffset = x - 2 + col;
			int yOffset = y - 2 + row;

			if (xOffset == x && yOffset == y)
			{
				continue;
			}

			float pixelAtOffset;
			if (xOffset < 0 || yOffset < 0 || xOffset >= width || yOffset >= height)
			{
				pixelAtOffset = 0.;
			}
			else
			{
				pixelAtOffset = srcImage.Pixel(xOffset, yOffset, 0);
			}
			if (pixelAtOffset >= centerPixel)
			{
				return false;
			}
		}
	}

	return true;
}

// TO DO---------------------------------------------------------------------
// Loop through the harrisImage to threshold and compute the local maxima in a neighborhood
// srcImage:  image with Harris values
// destImage: Assign 1 to a pixel if it is above a threshold and is the local maximum in 3x3 window, 0 otherwise.
//    You'll need to find a good threshold to use.
void computeLocalMaxima(CFloatImage &srcImage,CByteImage &destImage)
{
	int width = srcImage.Shape().width;
	int height = srcImage.Shape().height;

	double mean, stdDev;
	double sum = 0;
	double squareSum = 0;

	for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
			float pixel = srcImage.Pixel(x, y, 0);
			if (!(pixel >= 0 || pixel < 0))
			{
				auto error = "TRUE";
			}
			sum += srcImage.Pixel(x, y, 0);
		}
    }

	mean = sum / (float)(width * height);

	for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            squareSum += pow((srcImage.Pixel(x, y, 0) - mean), 2.);
        }
    }

	stdDev = sqrt(squareSum / (float)(width * height - 1));
	int count = 0;
	for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
			unsigned char *pixel = &destImage.Pixel(x, y, 0);
			if (srcImage.Pixel(x, y, 0) >= 3.*stdDev + mean && isLocalMax(srcImage, x, y))
			{
				count++;
				*pixel = 1;
			}
			else
			{
				*pixel = 0;
			}
        }
    }
}


// Compute Simple descriptors.
void ComputeSimpleDescriptors(CFloatImage &image, FeatureSet &features)
{
    //Create grayscale image used for Harris detection
    CFloatImage grayImage=ConvertToGray(image);

    vector<Feature>::iterator i = features.begin();
    while (i != features.end()) {
        Feature &f = *i;
		//these fields should already be set in the computeFeatures function
		int x = f.x;
		int y = f.y;

		// now get the 5x5 window surrounding the feature and store them in the features
		for(int row=(y-2); row<=(y+2); row++)
		{
			for(int col=(x-2); col<=(x+2); col++)
			{
				//if the pixel is out of bounds, assume it is black
				if(row<0 || row>=grayImage.Shape().height || col<0 || col>=grayImage.Shape().width)
				{
					f.data.push_back(0.0);
				}
				else
				{
					f.data.push_back(grayImage.Pixel(col,row,0));
				}
			}
		}
		printf("feature num %d\n", i->id);
        i++;
    }
}

CFloatImage GetXWindowAroundPixel(CFloatImage srcImage, int x, int y, int size)
{
	float *matrix = new float[size * size];

	for(int row=(y-(size-1)/2); row<=(y+(size-1)/2); row++)
		{
			for(int col=(x-(size-1)/2);col<=(x+(size-1)/2);col++)
			{
				if(row<0 || row>=srcImage.Shape().height || col<0 || col>=srcImage.Shape().width)
				{
					matrix[(row-(y-(size-1)/2))*size + (col-(x-(size-1)/2))] = 0.;
				}
				else
				{
					matrix[(row-(y-(size-1)/2))*size + (col-(x-(size-1)/2))] = srcImage.Pixel(col, row, 0);
				}
			}
		}
	return GetImageFromMatrix(matrix, size, size);
}

// Compute MOPs descriptors.
void ComputeMOPSDescriptors(CFloatImage &image, FeatureSet &features)
{
	CFloatImage grayImage=ConvertToGray(image);
	CFloatImage blurredImage;
	
	Convolve(grayImage, blurredImage, ConvolveKernel_7x7);

	CFloatImage postHomography = CFloatImage();
	CFloatImage gaussianImage = GetImageFromMatrix((float *)gaussian5x5Float, 5, 5);


	//first make the image invariant to changes in illumination by subtracting off the mean
	int grayHeight = grayImage.Shape().height;
	int grayWidth = grayImage.Shape().width;

	// now make this rotation invariant
    vector<Feature>::iterator featureIterator = features.begin();
    while (featureIterator != features.end()) {
		Feature &f = *featureIterator;

		CTransform3x3 scaleTransform = CTransform3x3();
		CTransform3x3 translationNegative;
		CTransform3x3 translationPositive;
		CTransform3x3 rotation;

		double scaleFactor = 41/8;
		scaleTransform[0][0] = scaleFactor;
		scaleTransform[1][1] = scaleFactor;

		translationNegative = translationNegative.Translation(f.x,f.y);
		translationPositive = translationPositive.Translation(-4, -4);
		rotation = rotation.Rotation(f.angleRadians * 180/ PI);
		
		CTransform3x3 finalTransformation = translationNegative * rotation * scaleTransform * translationPositive;
		//CFloatImage sample61x61Window = 
		//CFloatImage pixelWindow = GetXWindowAroundPixel(grayImage, f.x, f.y, 61);

		WarpGlobal(blurredImage, postHomography, finalTransformation, eWarpInterpLinear, 1.0f);

		//now we get the 41x41 box around the feature
		for(int row=0; row< 8; row++)
		{
			for(int col=0;col< 8;col++)
			{
				f.data.push_back(postHomography.Pixel(col, row, 0));
			}
		}
		/*
		// now we do the subsampling first round to reduce to a 20x20
		int imgSize = 41;
		subsample(&f, imgSize, gaussianImage);

		//second round of subsampling to get it to a 10x10
		imgSize = 20;
		subsample(&f, imgSize, gaussianImage);	

		imgSize = 10;
		CFloatImage img = featureToImage(f, imgSize, imgSize);
		CFloatImage blurredImg(img.Shape());
		Convolve(img, blurredImg, gaussianImage);
		featuresFromImage(&f,blurredImg,imgSize,imgSize);
		
		int count = 0;
		for(int y=0; y<imgSize; y++)
		{
			for(int x=0; x<imgSize; x++)
			{
				if(x == 3 || x == 7 || y == 3 || y == 7)
				{
					f.data.erase(f.data.begin() + count);
				}
				else
				{
					count++;
				}
			}
		}
		*/
		normalizeIntensities(&f, 8, 8);
		featureIterator++;

	}
}

void normalizeIntensities(Feature* f, int width, int height)
{
	double mean = 0.;
	vector<double, std::allocator<double>>::iterator it;
	it = f->data.begin();
	// calculate the mean
	for(int y=0; y<height; y++)
	{
		for(int x=0; x<width; x++)
		{
			mean+=*it;
			it++;
		}
	}
	mean = (mean/((double)width*height));
	
	// calculate the standard deviation
	it=f->data.begin();
	double stddev = 0.;
	for(int y=0; y<height; y++)
	{
		for(int x=0; x<width; x++)
		{
			stddev+=pow((*it-mean),2);
			it++;
		}
	}
	stddev = stddev/((double)width*height);
	stddev = sqrt(stddev);

	// subtract the mean and divide by the standard deviation
	Feature returnFeature;
	returnFeature.angleRadians = f->angleRadians;
	returnFeature.id = f->id;
	returnFeature.x = f->x;
	returnFeature.y = f->y;

	it = f->data.begin();
	for(int y=0; y<height; y++)
	{
		for(int x=0; x<width; x++)
		{
			double newVal = (*it - mean)/stddev;
			returnFeature.data.push_back(newVal);
			it++;
		}
	}
	*f = returnFeature;
}

void subsample(Feature* f, int imgSize, CFloatImage gaussianImage)
{
	vector<double, std::allocator<double>>::iterator it;
	CFloatImage img = featureToImage(*f, imgSize, imgSize);
	CFloatImage blurredImg(img.Shape());
	Convolve(img, blurredImg, gaussianImage);
	featuresFromImage(f,blurredImg,imgSize,imgSize);

	int count = 0;
	for(int y=0; y<imgSize; y++)
	{
		for(int x=0; x<imgSize; x++)
		{
			if(x%2 == 0 || y%2 == 0)
			{
				f->data.erase(f->data.begin() + count);
			}
			else
			{
				count++;
			}
		}
	}
}

CFloatImage featureToImage(Feature f, int width, int height)
{
	vector<double, std::allocator<double>>::iterator it;
	float *matrix = new float[width*height];
	int matIndex = 0;
	it = f.data.begin();
	while (it != f.data.end()) {
		double freq = *it;
		matrix[matIndex] = freq;
		matIndex++;
		it++;
	}
	CFloatImage img = GetImageFromMatrix(matrix, width, height);
	return img;
}

void featuresFromImage(Feature* f, CFloatImage img, int width, int height)
{
	vector<double, std::allocator<double>>::iterator it;
	f->data.clear();
	for(int y=0; y<height; y++)
	{
		for(int x=0; x<width; x++)
		{
			f->data.push_back(img.Pixel(x,y,0));
		}
	}
}

// Compute Custom descriptors (extra credit)
void ComputeCustomDescriptors(CFloatImage &image, FeatureSet &features)
{

}

// Perform simple feature matching.  This just uses the SSD
// distance between two feature vectors, and matches a feature in the
// first image with the closest feature in the second image.  It can
// match multiple features in the first image to the same feature in
// the second image.
void ssdMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) {
    int m = f1.size();
    int n = f2.size();

    matches.resize(m);
    totalScore = 0;

    double d;
    double dBest;
    int idBest;

    for (int i=0; i<m; i++) {
        dBest = 1e100;
        idBest = 0;

        for (int j=0; j<n; j++) {
            d = distanceSSD(f1[i].data, f2[j].data);

            if (d < dBest) {
		dBest = d;
		idBest = f2[j].id;
            }
        }

        matches[i].id1 = f1[i].id;
        matches[i].id2 = idBest;
        matches[i].score = dBest;
        totalScore += matches[i].score;
    }
}

// TODO: Write this function to perform ratio feature matching.  
// This just uses the ratio of the SSD distance of the two best matches as the score
// and matches a feature in the first image with the closest feature in the second image.
// It can match multiple features in the first image to the same feature in
// the second image.  (See class notes for more information, and the sshMatchFeatures function above as a reference)
void ratioMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) 
{
	int m = f1.size();
	int n = f2.size();

	matches.resize(m);
	totalScore = 0;

	double d;
	double dBest;
	double dSecondBest;
	int idBest;
	int idSecondBest;

	for (int i=0; i<m; i++) {
		dBest = 1e100;
	    idBest = 0;
		dSecondBest = 1e100;
		idSecondBest = 0;
		for (int j=0; j<n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);

			if (d < dBest) {
				dSecondBest = dBest;
				idSecondBest = idBest;
				dBest = d;
				idBest = f2[j].id;
			}
			else if(d >= dBest && d<dSecondBest)
			{
				dSecondBest = d;
				idSecondBest = f2[j].id;
			}
		}

		matches[i].id1 = f1[i].id;
		matches[i].id2 = idBest;
		matches[i].score = dBest/dSecondBest;
		totalScore += matches[i].score;
	}
}


// Convert Fl_Image to CFloatImage.
bool convertImage(const Fl_Image *image, CFloatImage &convertedImage) {
    if (image == NULL) {
        return false;
    }

    // Let's not handle indexed color images.
    if (image->count() != 1) {
        return false;
    }

    int w = image->w();
    int h = image->h();
    int d = image->d();

    // Get the image data.
    const char *const *data = image->data();

    int index = 0;

    for (int y=0; y<h; y++) {
        for (int x=0; x<w; x++) {
            if (d < 3) {
		// If there are fewer than 3 channels, just use the
		// first one for all colors.
		convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
		convertedImage.Pixel(x,y,1) = ((uchar) data[0][index]) / 255.0f;
		convertedImage.Pixel(x,y,2) = ((uchar) data[0][index]) / 255.0f;
            }
            else {
		// Otherwise, use the first 3.
		convertedImage.Pixel(x,y,0) = ((uchar) data[0][index]) / 255.0f;
		convertedImage.Pixel(x,y,1) = ((uchar) data[0][index+1]) / 255.0f;
		convertedImage.Pixel(x,y,2) = ((uchar) data[0][index+2]) / 255.0f;
            }

            index += d;
        }
    }
        
    return true;
}

// Convert CFloatImage to CByteImage.
void convertToByteImage(CFloatImage &floatImage, CByteImage &byteImage) {
    CShape sh = floatImage.Shape();

    assert(floatImage.Shape().nBands == byteImage.Shape().nBands);
    for (int y=0; y<sh.height; y++) {
        for (int x=0; x<sh.width; x++) {
            for (int c=0; c<sh.nBands; c++) {
		float value = floor(255*floatImage.Pixel(x,y,c) + 0.5f);

		if (value < byteImage.MinVal()) {
                    value = byteImage.MinVal();
		}
		else if (value > byteImage.MaxVal()) {
                    value = byteImage.MaxVal();
		}

		// We have to flip the image and reverse the color
		// channels to get it to come out right.  How silly!
		byteImage.Pixel(x,sh.height-y-1,sh.nBands-c-1) = (uchar) value;
            }
        }
    }
}

// Compute SSD distance between two vectors.
double distanceSSD(const vector<double> &v1, const vector<double> &v2) {
    int m = v1.size();
    int n = v2.size();

    if (m != n) {
        // Here's a big number.
        return 1e100;
    }

    double dist = 0;

        
    for (int i=0; i<m; i++) {
        dist += pow(v1[i]-v2[i], 2);
    }
        
        
    return sqrt(dist);
}

// Transform point by homography.
void applyHomography(double x, double y, double &xNew, double &yNew, double h[9]) {
    double d = h[6]*x + h[7]*y + h[8];

    xNew = (h[0]*x + h[1]*y + h[2]) / d;
    yNew = (h[3]*x + h[4]*y + h[5]) / d;
}

// Compute AUC given a ROC curve
double computeAUC(vector<ROCPoint> &results)
{
    double auc=0;
    double xdiff,ydiff;
    for (int i = 1; i < (int) results.size(); i++) {
        //fprintf(stream,"%lf\t%lf\t%lf\n",thresholdList[i],results[i].falseRate,results[i].trueRate);
        xdiff=(results[i].falseRate-results[i-1].falseRate);
        ydiff=(results[i].trueRate-results[i-1].trueRate);
        auc=auc+xdiff*results[i-1].trueRate+xdiff*ydiff/2;
    }
    return auc;
}
