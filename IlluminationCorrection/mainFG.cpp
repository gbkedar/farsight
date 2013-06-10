/* 
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include <vector>
#include <string>
#include <cstdio>
#include <sstream>
#include <iomanip>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkIntTypes.h"
#include "itkNumericTraits.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"

typedef unsigned short USPixelType;
typedef unsigned char  UCPixelType;
const unsigned int     Dimension3 = 3;
const unsigned int     Dimension2 = 2;
typedef itk::Image< USPixelType, Dimension3 > US3ImageType;
typedef itk::Image< USPixelType, Dimension2 > US2ImageType;

#define WSz 128

void usage( const char *funcName )
{
  std::cout << "USAGE:"
	    << " " << funcName << " InputImage OutputImage\n";
}

void GetTile( US2ImageType::Pointer &currentTile, US3ImageType::Pointer &readImage, unsigned i )
{
  typedef itk::ExtractImageFilter< US3ImageType, US2ImageType > DataExtractType;
  DataExtractType::Pointer deFilter = DataExtractType::New();
  US3ImageType::RegionType dRegion  = readImage->GetLargestPossibleRegion();
  dRegion.SetSize (2,0);
  dRegion.SetIndex(2,i);
  deFilter->SetExtractionRegion(dRegion);
  deFilter->SetDirectionCollapseToIdentity();
  deFilter->SetInput( readImage );
  try
  {
    deFilter->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught !" << excep << std::endl;
    exit (EXIT_FAILURE);
  }
  currentTile = deFilter->GetOutput();
  currentTile->Register();
}

US2ImageType::PixelType returnthresh( itk::SmartPointer<US2ImageType> input_image,
					int num_bin_levs, int num_in_fg )
{
  //Instantiate the different image and filter types that will be used
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIteratorType;
  typedef itk::Statistics::Histogram< float > HistogramType;
  typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;

  //Create a temporary histogram container:
  const int numBins = itk::NumericTraits<US2ImageType::PixelType>::max();
  double *tempHist;
  tempHist = (double*) malloc( sizeof(double) * numBins );
  for(US2ImageType::PixelType i=0; i<numBins; ++i)
     tempHist[i] = 0;

  US2ImageType::PixelType maxval = itk::NumericTraits<US2ImageType::PixelType>::ZeroValue();
  US2ImageType::PixelType minval = itk::NumericTraits<US2ImageType::PixelType>::max();
  //Populate the histogram (assume pixel type is actually is some integer type):
  ConstIteratorType it( input_image, input_image->GetRequestedRegion() );
  for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
  {
    US2ImageType::PixelType pix = it.Get();
    ++tempHist[pix];
    if( pix > maxval ) maxval = pix;
    if( pix < minval ) minval = pix;
  }
  //return max of type if there is no variation in the staining
  if( (maxval-minval)<3 ) return itk::NumericTraits<US2ImageType::PixelType>::max(); 
  const US2ImageType::PixelType numBinsPresent = maxval+1;
  
  //Find max value in the histogram
  double floatIntegerMax = itk::NumericTraits<US2ImageType::PixelType>::max();
  double max = 0.0;
  for(US2ImageType::PixelType i=0; i<numBinsPresent; ++i)
    if( tempHist[i] > max )
    	max = tempHist[i];

  double scaleFactor = 1;
  if(max >= floatIntegerMax)
    scaleFactor = floatIntegerMax / max;

  HistogramType::Pointer histogram = HistogramType::New() ;
  // initialize histogram
  HistogramType::SizeType size;
  HistogramType::MeasurementVectorType lowerBound;
  HistogramType::MeasurementVectorType upperBound;

  lowerBound.SetSize(1);
  upperBound.SetSize(1);
  size.SetSize(1);

  lowerBound.Fill(0.0);
  upperBound.Fill((double)maxval);
  size.Fill(numBinsPresent);

  histogram->SetMeasurementVectorSize(1);
  histogram->Initialize(size, lowerBound, upperBound ) ;

  US2ImageType::PixelType i=0;
  for (HistogramType::Iterator iter = histogram->Begin(); iter != histogram->End(); ++iter )
  {
    float norm_freq = (float)(tempHist[i] * scaleFactor);
    iter.SetFrequency(norm_freq);
    ++i;
  }
  free( tempHist );

  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetNumberOfThresholds( num_bin_levs );
  calculator->SetInputHistogram( histogram );
  calculator->Update();
  const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
  CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();

  float thresh;

  for(US2ImageType::PixelType i=0; i < num_in_fg; ++itNum, ++i)
    thresh = (static_cast<float>(*itNum));

  std::cout<<"Threshold computed: "<<thresh<<std::endl;

  return (US2ImageType::PixelType)(thresh+0.5);
}

int main(int argc, char *argv[])
{ 
  if( argc < 3 )
  {
    usage(argv[0]);
    std::cerr << "PRESS ENTER TO EXIT\n";
    getchar();
    return EXIT_FAILURE;
  }

  std::string inputImageName  = argv[0]; //Just In case..
  std::string outputImageName = argv[1]; //Name of the input image

  typedef itk::ImageFileReader< US3ImageType >    ReaderType;
  typedef itk::ImageFileWriter< US3ImageType >    WriterType;
  typedef itk::MedianImageFilter< US2ImageType, US2ImageType > MedianFilterType;
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  typedef itk::ImageRegionIteratorWithIndex< US2ImageType > IterType;
  typedef itk::ImageRegionIteratorWithIndex< US3ImageType > IterType3d;

  ReaderType::Pointer reader = ReaderType::New();
  reader = ReaderType::New();
  reader->SetFileName( inputImageName.c_str() );
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }
  US3ImageType::Pointer inputImage = reader->GetOutput();

  itk::IndexValueType numSlices = inputImage->GetLargestPossibleRegion().GetSize()[2];
  itk::IndexValueType numCol = inputImage->GetLargestPossibleRegion().GetSize()[1];
  itk::IndexValueType numRow = inputImage->GetLargestPossibleRegion().GetSize()[0];
  US2ImageType::PixelType lowT  = itk::NumericTraits< US2ImageType::PixelType >::max();
  US2ImageType::PixelType highT = 0;
  double meanT = 0;

  std::vector< itk::SmartPointer<US2ImageType> > medianImages;
  std::vector< itk::SmartPointer<US2ImageType> > contourImages;
  std::vector< US2ImageType::PixelType > thresholds;
  thresholds.resize( numSlices ); medianImages.resize( numSlices );
  contourImages.resize( numSlices );
  #pragma omp parallel for
  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    US2ImageType::Pointer currentSlice;
    GetTile( currentSlice, inputImage, (unsigned)i );

    //Median filter for each slice to remove thermal noise
    MedianFilterType::Pointer medFilter = MedianFilterType::New();
    medFilter->SetInput( currentSlice );
    medFilter->SetRadius( 3 );

    //Get Contours using a large median filter 
    MedianFilterType::Pointer medFilter1 = MedianFilterType::New();
    medFilter1->SetInput( currentSlice );
    medFilter1->SetRadius( 9 );
    try
    {
      medFilter->Update();
      medFilter1->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception caught !" << excep << std::endl;
      exit (EXIT_FAILURE);
    }
    US2ImageType::Pointer medFiltIm = medFilter->GetOutput();
    US2ImageType::Pointer contourIm = medFilter1->GetOutput();
    
    ConstIterType it( medFiltIm, medFiltIm->GetRequestedRegion() );
    IterType it1( contourIm, contourIm->GetRequestedRegion() );
    it1.GoToBegin();
    for( it.GoToBegin(); !it.IsAtEnd(); ++it, ++it1 )
    {
      if( it.Get()>it1.Get() )
        it1.Set( it.Get()-it1.Get() );
      else
	it1.Set( 0 );
    }
    US2ImageType::PixelType curThresh = returnthresh( contourIm, 1, 1 );
    for( it1.GoToBegin(); !it1.IsAtEnd(); ++it1 )
    {
      it1.Set( 0 );
      if( it1.Get()>curThresh )
	it1.Set( 1 );
    }
    medFiltIm->Register(); medianImages.at(i) = medFiltIm;
    contourIm->Register(); contourImages.at(i) = contourIm;
    thresholds.at(i) = returnthresh( medFiltIm, 1, 1 );
    if( thresholds.at(i)>highT )
    {
      #pragma omp critical
      {
        if( thresholds.at(i)>highT )
	  highT = thresholds.at(i);
      }
    }
    if( thresholds.at(i)<lowT )
    {
      #pragma omp critical
      {
        if( thresholds.at(i)<lowT )
	  lowT = thresholds.at(i);
      }
    }
    #pragma omp critical
    meanT += thresholds.at(i);
  }
  meanT /= numSlices;

  //Compute thresh for all slices from threshs from individual slices
  double sigmaT = 0, threshSlices = 0;
  {
    typedef itk::Statistics::Histogram< float > HistogramType;
    typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;
    std::vector< double > histVals( (highT-lowT+1), 0 );
    for( itk::SizeValueType i=0; i<thresholds.size(); ++i )
    {
      histVals[(thresholds.at(i)-1)] += 1;
      sigmaT += (((double)(meanT-thresholds.at(i)))*((double)(meanT-thresholds.at(i))));
    }
    sigmaT /= numSlices;
    sigmaT = std::sqrt(sigmaT);
    HistogramType::Pointer histogram = HistogramType::New() ;
    // initialize histogram
    HistogramType::SizeType size;
    HistogramType::MeasurementVectorType lowerBound;
    HistogramType::MeasurementVectorType upperBound;

    lowerBound.SetSize(1);
    upperBound.SetSize(1);
    size.SetSize(1);

    lowerBound.Fill((double)lowT);
    upperBound.Fill((double)highT);
    size.Fill((highT-lowT+1));

    histogram->SetMeasurementVectorSize(1);
    histogram->Initialize(size, lowerBound, upperBound ) ;
    US2ImageType::PixelType i = 0;
    for( HistogramType::Iterator iter = histogram->Begin(); iter != histogram->End(); ++iter )
    {
      iter.SetFrequency(histVals.at(i));
      ++i;
    }

    CalculatorType::Pointer calculator = CalculatorType::New();
    calculator->SetNumberOfThresholds( 1 );
    calculator->SetInputHistogram( histogram );
    calculator->Update();
    const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
    CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();
    for(US2ImageType::PixelType i=0; i < 1; ++itNum, ++i)
      threshSlices = (static_cast<float>(*itNum));
  }
  std::cout<<"Main threshold: "<<threshSlices<<"\t"
  	   <<"Sigma: "<<sigmaT<<"\t"
	   <<"Mod thr: "<<(threshSlices-(2*sigmaT))<<std::endl;
  threshSlices -= (2*sigmaT);
  //Windowed thresholding and get fg points for all the slices
  #pragma omp parallel for
  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    itk::IndexValueType numjParts = std::ceil( ((double)numRow)/((double)WSz));
    itk::IndexValueType numkParts = std::ceil( ((double)numCol)/((double)WSz));
    IterType  itFull( medianImages.at(i),   medianImages.at(i)->GetRequestedRegion() );
    IterType itCFull( contourImages.at(i), contourImages.at(i)->GetRequestedRegion() );
    for( itk::IndexValueType j=0; j<numjParts; ++j )
      for( itk::IndexValueType k=0; k<numkParts; ++k )
      {
	US2ImageType::Pointer windowIm = US2ImageType::New();
	//Allocate windowed image
	US2ImageType::PointType origin;
	origin[0] = 0; origin[1] = 0;
	windowIm->SetOrigin( origin );
	US2ImageType::IndexType start;
	start[0] = 0; start[1] = 0;
	US2ImageType::SizeType size;
	size[0] = (((j+1)*WSz-1)>=numRow) ? numRow-(j*WSz) : WSz;
	size[1] = (((k+1)*WSz-1)>=numCol) ? numCol-(k*WSz) : WSz;
	US2ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	windowIm->SetRegions( region );
	windowIm->Allocate();
	windowIm->FillBuffer(0);
	windowIm->Update();
	//Copy windowed image
	IterType itCrop( windowIm, windowIm->GetRequestedRegion() );
	for( itk::SizeValueType cVal=0; cVal<size[1]; ++cVal )
	  for( itk::SizeValueType rVal=0; rVal<size[0]; ++rVal )
	  {
	    US2ImageType::IndexType wIndex, fIndex;
	    wIndex[0] = rVal; wIndex[1] = cVal;
	    fIndex[0] = (j*WSz)+rVal; fIndex[1] = (j*WSz)+cVal;
	    itCrop.SetIndex( wIndex );
	    itFull.SetIndex( fIndex );
	    itCrop.Set( itFull.Get() );
	  }
	US2ImageType::PixelType curThr = returnthresh( windowIm, 1, 1 );
	if( curThr < threshSlices )
	  curThr = threshSlices;
	for( itk::SizeValueType cVal=0; cVal<size[1]; ++cVal )
	  for( itk::SizeValueType rVal=0; rVal<size[0]; ++rVal )
	  {
	    US2ImageType::IndexType fIndex,cIndex;
	    fIndex[0] = (j*WSz)+rVal; fIndex[1] = (j*WSz)+cVal;
	    cIndex[0] = (j*WSz)+rVal; cIndex[1] = (j*WSz)+cVal;
	    itFull.SetIndex( fIndex ); itCFull.SetIndex( cIndex );
	    if( itFull.Get() < curThr )
	    {
	      if( itCFull.Get() != 1 )
	      {
	        itCFull.Set( 0 );
	      }
	    }
	    else
	    {
	      itCFull.Set( 2 );
	    }
	  }
      }
  }

  //Copy into 3d image
  US3ImageType::Pointer outputImage = US3ImageType::New();
  US3ImageType::PointType origin;
  origin[0] = 0; origin[1] = 0; origin[2] = 0;
  outputImage->SetOrigin( origin );
  US3ImageType::IndexType start;
  start[0] = 0; start[1] = 0; start[2] = 0;
  US3ImageType::SizeType  size;
  size[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
  size[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
  size[2] = inputImage->GetLargestPossibleRegion().GetSize()[2];
  US3ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  outputImage->SetRegions( region );
  outputImage->Allocate();
  outputImage->FillBuffer(0);
  outputImage->Update();

  #pragma omp parallel for
  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    IterType3d itOut( outputImage, outputImage->GetRequestedRegion() );
    IterType itCFull( contourImages.at(i), contourImages.at(i)->GetRequestedRegion() );
    for( itk::IndexValueType j=0; j<numCol; ++j )
      for( itk::IndexValueType k=0; k<numRow; ++k )
      {
        US2ImageType::IndexType cIndex; cIndex[0] = k; cIndex[1] = j;
	US3ImageType::IndexType oIndex; oIndex[0] = k; oIndex[1] = j; oIndex[2] = k;
	itOut.SetIndex( oIndex ); itCFull.SetIndex( cIndex );
	itOut.Set( itCFull.Get() );
      }
  }

  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    medianImages.at(i)->UnRegister();
    contourImages.at(i)->UnRegister();
  }

  WriterType::Pointer writer = WriterType::New();
  writer = WriterType::New();
  writer->SetFileName( outputImageName.c_str() );
  writer->SetInput( outputImage );
  try
  {
    writer->Update();
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }
  exit( EXIT_SUCCESS );
}
