#ifndef _COLORSEGMENTATION_H_
#define _COLORSEGMENTATION_H_

//ITK Includes
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

//STL and Other Includes
#include <iostream>
#include <algorithm>
#include <limits.h>
#include <float.h>
#include <math.h>

#include "dhEvalState.h"
#include "dhSlice.h"
#include "dhHistogram.h"
#include "RGB_Atype.h"
#include "dhSeedGrid.h"
#include "dhClassifiers.h"

typedef unsigned char UcharPixelType;
typedef unsigned short UshrtPixelType;
typedef float FloatPixelType;
typedef itk::RGBPixel< UcharPixelType > RGBPixelType;
typedef itk::Image< RGBPixelType, 3 > RGBImageType;
typedef itk::Image< UcharPixelType, 3 > UcharImageType;
typedef itk::Image< UshrtPixelType, 3 > UshrtImageType;
typedef itk::Image< FloatPixelType, 3 > FloatImageType;

enum Pixel_Class { UNKNOWN, RED_CELL, BLUE_CELL, BKGD_FIELD };

class ColorSegmentation
{
public:
	//Constructor
	ColorSegmentation(RGBImageType::Pointer input);
	~ColorSegmentation(){};//Destructor

	void SetTesting(bool t = true){ TESTING = t; };		//default is false
	void SetIgnoreBackground(bool i = true){ IGNORE_BACKGROUND = i; }; //default is false
	void SetLightBackground(bool d = true){ LIGHT_BACKGROUND = d; }; //default is false

	//Methods:
	void TransformToRLI();			//First step
	void FindArchetypalColors();	//Compute Archetypal Colors
	void SetArchetypalColors(dh::RLI r, dh::RLI b, dh::RLI w);
	void SetArchetypalColors(dh::_RGB r, dh::_RGB b, dh::_RGB w);
	void ComputeClassWeights();		//Get Grayscales Based On Distances From Atypes

	//Get Results:
	UcharImageType::Pointer ComputeBinary(int num_bins, int num_in_fg, bool fgrnd_dark = false);

	//A few parameters:
	bool IGNORE_BACKGROUND;
	bool LIGHT_BACKGROUND;
	bool TESTING;

	//void RunInitialBinarization();	//MOVED TO BOTTOM OF CPP

protected:
	//Image Pointers
	RGBImageType::Pointer rgb_input;

	UcharImageType::Pointer red_image;
	UcharImageType::Pointer lime_image;
	UcharImageType::Pointer intensity_image;

	UcharImageType::Pointer red_weights;
	UcharImageType::Pointer blue_weights;

	//Intermediate values
	// These actually contain RLI values (/2):
	dh::RLI archTypRED, archTypBLUE, archTypBACK; //1->Red-ish 2->Blue-ish

private:
	void go_best_dir( dh::Histogram *hist, dh::_RGB& ma, bool& moved, const dh::_RGB& sa, const dh::_RGB& bkgnd, const int r = 1, int res = 1 );
};

#endif
