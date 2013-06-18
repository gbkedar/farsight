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
#include "itkMeanImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"

#define WinSz 50

typedef unsigned short	USPixelType;
typedef unsigned char	UCPixelType;
typedef double		CostPixelType;
const unsigned int	Dimension3 = 3;
const unsigned int	Dimension2 = 2;
typedef itk::Image< USPixelType, Dimension3 > US3ImageType;
typedef itk::Image< UCPixelType, Dimension3 > UC3ImageType;
typedef itk::Image< USPixelType, Dimension2 > US2ImageType;
typedef itk::Image< CostPixelType, Dimension2 > CostImageType;

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

itk::SizeValueType ComputeHistogram(
	std::vector< itk::SmartPointer<US2ImageType>  > &medFiltImages,
	std::vector< double > &histogram,
	US2ImageType::IndexType &start )
{
  itk::SizeValueType max = 0;
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  for( itk::SizeValueType i=0; i<medFiltImages.size(); ++i )
  {
    US2ImageType::SizeType size; size[0] = WinSz; size[1] = WinSz;
    US2ImageType::RegionType region;
    region.SetSize( size ); region.SetIndex( start );
    ConstIterType constIter ( medFiltImages.at(i), region );
    for( constIter.GoToBegin(); !constIter.IsAtEnd(); ++constIter )
      ++histogram[constIter.Get()];
  }
  for( itk::SizeValueType i=0; i<histogram.size(); ++i )
  {
    if( histogram[i] )
      max = i;
    histogram[i] /= (double)WinSz*(double)WinSz;
  }
  return max;
}

void ComputeCosts( std::vector< itk::SmartPointer<US2ImageType>  > &medFiltImages,
		   std::vector< itk::SmartPointer<CostImageType> > &autoFlourCosts,
		   std::vector< itk::SmartPointer<CostImageType> > &flourCosts )
{
  itk::IndexValueType numCol = medFiltImages.at(0)->GetLargestPossibleRegion().GetSize()[1];
  itk::IndexValueType numRow = medFiltImages.at(0)->GetLargestPossibleRegion().GetSize()[0];
  itk::IndexValueType WinSz2 = (itk::IndexValueType)std::round(((double)WinSz)/2-0.5);
#ifdef _OPENMP
  #pragma omp parallel for
#endif
#endif
  for( itk::IndexValueType i=0; i<numRow; ++i )
  {
    for( itk::IndexValueType j=0; j<numCol; ++j )
    {
      //Compute histogram at point i,j with window size define WinSz
      itk::SizeValueType histSize = itk::NumericTraits<US2ImageType::PixelType>::max()+1;
      std::vector< double > histogram( histSize, 0 );
      US2ImageType::IndexType start; start[0] = numRow-WinSz-1; start[1] = numCol-WinSz-1;
      if( (i-WinSz2)<0 ) start[0] = 0; if( (j-WinSz2)<0 ) start[1] = 0;
      if( start[0] && ((i+WinSz2+1)<numRow) ) start[0] = i-WinSz2;
      if( start[1] && ((j+WinSz2+1)<numCol) ) start[1] = j-WinSz2;
      itk::SizeValueType max = ComputeHistogram( medFiltImages, histogram, start );
      US2ImageType::IndexType curPoint; curPoint[0] = i; curPoint[1] = j;
      ComputeCostsFromHist( curPoint, &medFiltImages, autoFlourCosts, flourCosts );
    }
  }
  return;
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

  std::string inputImageName  = argv[1]; //Just In case..
  std::string outputImageName = argv[2]; //Name of the input image

  typedef itk::ImageFileReader< US3ImageType >    ReaderType;
  typedef itk::ImageFileWriter< UC3ImageType >    WriterType;
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

  std::cout<<"Number of slices:"<<numSlices<<std::endl;

  std::vector< itk::SmartPointer<US2ImageType> > medFiltImages;
  medFiltImages.resize( numSlices );
#ifdef _OPENMP
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
#if _OPENMP >= 200805L
  #pragma omp parallel for schedule(dynamic,1)
#else
  #pragma omp parallel for
#endif
#endif
  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    US2ImageType::Pointer currentSlice;
    GetTile( currentSlice, inputImage, (unsigned)i );
    //Median filter for each slice to remove thermal noise
    MedianFilterType::Pointer medFilter = MedianFilterType::New();
    medFilter->SetInput( currentSlice );
    medFilter->SetRadius( 3 );
    try
    {
      medFilter ->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception caught !" << excep << std::endl;
      exit (EXIT_FAILURE);
    }
    US2ImageType::Pointer medFiltIm = medFilter->GetOutput();
    medFiltIm->Register();
    medFiltImages.at(i) = medFiltIm;
    currentSlice->UnRegister();
  }

  std::cout<<"Done! Computing Costs\n"<<std::flush;

  //Allocate two images of double to store the auto-flour and the flour costs
  std::vector< itk::SmartPointer<CostImageType> > autoFlourCosts, flourCosts;
  autoFlourCosts.resize( numSlices );
  flourCosts.resize( numSlices );
#ifdef _OPENMP
#if _OPENMP >= 200805L
  #pragma omp parallel for schedule(dynamic,1)
#else
  #pragma omp parallel for
#endif
#endif
  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    CostImageType::Pointer costs1 = CostImageType::New();
    CostImageType::Pointer costs2 = CostImageType::New();
    CostImageType::PointType origin;	origin[0] = 0;	  origin[1] = 0;
    CostImageType::IndexType start;	start[0] = 0;	  start[1] = 0;
    CostImageType::SizeType size;	size[0] = numRow; size[1] = numCol;
    CostImageType::RegionType region;
    region.SetSize( size ); region.SetIndex( start );
    costs1->SetOrigin( origin );	costs2->SetOrigin( origin );
    costs1->SetRegions( region );	costs2->SetRegions( region );
    costs1->Allocate();			costs2->Allocate();
    costs1->FillBuffer(0);		costs2->FillBuffer(0);
    costs1->Update();			costs2->Update();
    costs1->Register();			costs2->Register();
    autoFlourCosts.at(i) = costs1; 	flourCosts.at(i) = costs2;
  }

  ComputeCosts( medFiltImages, autoFlourCosts, flourCosts );



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

#ifdef _OPENMP
#if _OPENMP >= 200805L
  #pragma omp parallel for  schedule(dynamic,1)
#else
  #pragma omp parallel for
#endif
#endif
  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    IterType3d itOut( outputImage, outputImage->GetRequestedRegion() );
    IterType itCFull( threshImages.at(i), threshImages.at(i)->GetRequestedRegion() );
    for( itk::IndexValueType j=0; j<numCol; ++j )
      for( itk::IndexValueType k=0; k<numRow; ++k )
      {
        US2ImageType::IndexType cIndex; cIndex[0] = k; cIndex[1] = j;
	US3ImageType::IndexType oIndex; oIndex[0] = k; oIndex[1] = j; oIndex[2] = i;
	itOut.SetIndex( oIndex ); itCFull.SetIndex( cIndex );
	itOut.Set( itCFull.Get() ); //Writing out binary type change should be ok
      }
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

  try
  {
    for( itk::IndexValueType i=0; i<numSlices; ++i )
    {
      threshImages.at(i)->UnRegister();
    }
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }

  exit( EXIT_SUCCESS );
}
