#pragma once

#ifndef FEATURES_H
#define FEATURES_H

#include "ImageLib/ImageLib.h"
#include "ImageDatabase.h"

class Fl_Image;

//5x5 Gaussian
const double gaussian5x5[25] = {0.003663, 0.014652,  0.025641,  0.014652,  0.003663, 
0.014652, 0.0586081, 0.0952381, 0.0586081, 0.014652, 
0.025641,   0.0952381, 0.150183,  0.0952381, 0.025641, 
0.014652, 0.0586081, 0.0952381, 0.0586081, 0.014652, 
0.003663, 0.014652,  0.025641,  0.014652,  0.003663 };


//5x5 Gaussian
const float gaussian5x5Float[25] = {0.003663, 0.014652,  0.025641,  0.014652,  0.003663, 
0.014652, 0.0586081, 0.0952381, 0.0586081, 0.014652, 
0.025641,   0.0952381, 0.150183,  0.0952381, 0.025641, 
0.014652, 0.0586081, 0.0952381, 0.0586081, 0.014652, 
0.003663, 0.014652,  0.025641,  0.014652,  0.003663 };

static float derivativeX3x3[9] = {0.0, 0.0, 0.0,
1.0, 0.0, -1.0,
0.0, 0.0, 0.0};

static float derivativeY3x3[9] = {0.0, 1.0, 0.0,
0.0, 0.0, 0.0,
0.0, -1.0, 0.0};

//7x7 Gaussian
const double gaussian7x7[49] = {0.000896861, 0.003587444, 0.006278027, 0.00896861,  0.006278027, 0.003587444, 0.000896861,
0.003587444, 0.010762332, 0.023318386, 0.029596413, 0.023318386, 0.010762332, 0.003587444, 
0.006278027, 0.023318386, 0.049327354, 0.06367713,  0.049327354, 0.023318386, 0.006278027,
0.00896861,  0.029596413, 0.06367713,  0.08161435,  0.06367713,  0.029596413, 0.00896861,  
0.006278027, 0.023318386, 0.049327354, 0.06367713,  0.049327354, 0.023318386, 0.006278027,
0.003587444, 0.010762332, 0.023318386, 0.029596413, 0.023318386, 0.010762332, 0.003587444,
0.000896861, 0.003587444, 0.006278027, 0.00896861,  0.006278027, 0.003587444, 0.000896861};

const double sobelX[9] = 
{
	-1./8., 0., 1./8.,
	-2./8., 0., 2./8.,
	-1./8., 0., 1./8.
};

const double sobelY[9] = 
{
	-1./8., -2./8., -1./8.,
	0., 0., 0.,
	1./8., 2./8., 1./8.
};

struct ROCPoint
{
	double trueRate;
	double falseRate;
};

void test();

template <class T>
CImageOf<T> GetImageFromMatrix(T *matrix, int width, int height);

bool isLocalMax(CFloatImage srcImage, int x, int y);

CFloatImage GetXWindowAroundPixel(CFloatImage srcImage, int x, int y, int size);

void GetHarrisComponents(CFloatImage &srcImage, CFloatImage &A, CFloatImage &B, CFloatImage &C, CFloatImage *partialX=nullptr, CFloatImage *partialY=nullptr);

double GetCanonicalOrientation(int x, int y, CFloatImage A, CFloatImage B, CFloatImage C, CFloatImage partialX, CFloatImage partialY);

// Compute harris values of an image.
void computeHarrisValues(CFloatImage &srcImage,CFloatImage &destImage);

//  Compute local maximum of Harris values in an image.
void computeLocalMaxima(CFloatImage &srcImage,CByteImage &destImage);

// Compute features of an image.
bool computeFeatures(CFloatImage &image, FeatureSet &features, int featureType, int descriptorType);

// Perform a query on the database.
bool performQuery(const FeatureSet &f1, const ImageDatabase &db, int &bestIndex, vector<FeatureMatch> &bestMatches, double &bestScore, int matchType);

// Match one feature set with another.
bool matchFeatures(const FeatureSet &f, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore, int matchType);

// Add ROC curve data to the data vector
void addRocData(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9],vector<bool> &isMatch,double threshold,double &maxD);

// Evaluate a match using a ground truth homography.
double evaluateMatch(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9]);

// Silly example feature detector
void dummyComputeFeatures(CFloatImage &image, FeatureSet &features);

// Harris feature detector
void ComputeHarrisFeatures(CFloatImage &image, FeatureSet &features);

// Compute Simple descriptors
void ComputeSimpleDescriptors(CFloatImage &image, FeatureSet &features);

// Compute MOPS descriptors
void ComputeMOPSDescriptors(CFloatImage &image, FeatureSet &features);

void subsample(Feature* f, int imgSize, CFloatImage gaussianImage);

CFloatImage featureToImage(Feature f, int width, int height);

void featuresFromImage(Feature* f, CFloatImage img, int width, int height);

void normalizeIntensities(Feature* f, int width, int height);

// Compute Custom descriptors
void ComputeCustomDescriptors(CFloatImage &image, FeatureSet &features);

// Perform ssd feature matching.
void ssdMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore);

// Perform ratio feature matching.  You must implement this.
void ratioMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore);

// Convert Fl_Image to CFloatImage.
bool convertImage(const Fl_Image *image, CFloatImage &convertedImage);

// Convert CFloatImage to CByteImage.
void convertToByteImage(CFloatImage &floatImage, CByteImage &byteImage);

// Compute SSD distance between two vectors.
double distanceSSD(const vector<double> &v1, const vector<double> &v2);

// Transform point by homography.
void applyHomography(double x, double y, double &xNew, double &yNew, double h[9]);

// Computes points on the Roc curve
vector<ROCPoint> computeRocCurve(vector<FeatureMatch> &matches,vector<bool> &isMatch,vector<double> &thresholds);

// Compute AUC given a ROC curve
double computeAUC(vector<ROCPoint> &results);

#endif