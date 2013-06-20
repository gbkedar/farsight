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
#define DBGG

#include <vector>
#include <string>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <float.h>

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
#define NumBins 1000

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
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  itk::SizeValueType histSize = itk::NumericTraits<US2ImageType::PixelType>::max()+1;
  std::vector< double > histogramInternal( histSize, 0 );
  itk::SizeValueType max = 0;
  for( itk::SizeValueType i=0; i<medFiltImages.size(); ++i )
  {
    US2ImageType::SizeType size; size[0] = WinSz; size[1] = WinSz;
    US2ImageType::RegionType region;
    region.SetSize( size ); region.SetIndex( start );
    ConstIterType constIter ( medFiltImages.at(i), region );
    for( constIter.GoToBegin(); !constIter.IsAtEnd(); ++constIter )
      ++histogramInternal[constIter.Get()];
    if( constIter.Get()>max )
      max = constIter.Get();
  }
  itk::SizeValueType valsPerBin = std::floor( ((double)max+1)/((double)NumBins) + 0.5 );
  itk::SizeValueType overFlow = 0;
  if( valsPerBin*NumBins >= (max+1) )
    overFlow = max+1-valsPerBin*NumBins;
#ifdef DBGG
  double sum = 0.0;
#endif //DBGG
  for( itk::SizeValueType i=0; i<NumBins; ++i )
  {
    for( itk::SizeValueType j=0; j<valsPerBin; ++j )
      histogram.at(i) += histogramInternal.at(i*valsPerBin+j);
    if( i==(NumBins-1) )
      for( itk::SizeValueType j=0; j<overFlow; ++j )
	histogram.at(i) += histogramInternal.at(i+1*valsPerBin+j);
    histogram.at(i) /= ((double)WinSz)*((double)WinSz)*((double)medFiltImages.size());
    sum += histogram.at(i);
  }
#ifdef DBGG
  std::cout<<"Sum:"<<sum<<std::endl;
#endif //DBGG
  if( valsPerBin )
    max = valsPerBin*NumBins - 1;
  return max;
}

void computePoissonParams( std::vector< double > &histogram,
			   std::vector< double > &parameters )
{
  itk::SizeValueType max = histogram.size()-1;
  std::cout<<"computing pos params with max " << max << "\n" << std::flush;
  //The three-level min error thresholding algorithm
  double min_J = DBL_MAX;
  double P0, U0, P1, U1, P2, U2, U, J;
  // Try this: we need to define a penalty term that depends on the number of parameters
  //The penalty term is given as 0.5*k*ln(n)
  //where k is the number of parameters of the model and n is the number of samples
  //In this case, k=6 and n=256
  double PenTerm3 = sqrt(6.0)*log(((double)max));
  for( itk::SizeValueType i=0; i<(max-1); ++i )//to set the first threshold
  {
    //compute the current parameters of the first component
    P0 = U0 = 0.0;
    for( itk::SizeValueType l=0; l<=i; l++ )
    {
      P0 += histogram.at(l);
      U0 += (l+1)*histogram.at(l);
    }
    U0 /= P0;

    for( itk::SizeValueType j=i+1; j<max; ++j )//to set the second threshold
    {
      //compute the current parameters of the second component
      P1 = U1 = 0.0;
      for( itk::SizeValueType l=i+1; l<=j; ++l )
      {
        P1 += histogram.at(l);
        U1 += (l+1)*histogram.at(l);
      }
      U1 /= P1;

      //compute the current parameters of the third component
      P2 = U2 = 0.0;
      for( itk::SizeValueType l=j+1; l<=max; ++l)
      {
        P2 += histogram.at(l);
        U2 += (l+1)*histogram.at(l);
      }
      U2 /= P2;

      //compute the overall mean
      U = P0*U0 + P1*U1 + P2*U2;

      //Compute the current value of the error criterion function
      J = U - (P0*(log(P0)+U0*log(U0))+ P1*(log(P1)+U1*log(U1)) + P2*(log(P2)+U2*log(U2)));
      //Add the penalty term
      J += PenTerm3;

      if( J<min_J )
      {
        min_J = J;
        parameters.at(0) = U0;
        parameters.at(1) = U1;
        parameters.at(2) = U2; //Just a negative number to let the program knows that two levels will be used
        parameters.at(3) = P0;
        parameters.at(4) = P1;
      }
    }
  }
#ifdef DBGG
  std::cout<<"Parameters1: ";
  for( itk::SizeValueType j=0; j<parameters.size(); ++j )
    std::cout<<parameters.at(j)<<"\t";
  std::cout<<"\n"<<std::flush;
#endif //DBGG

  //try this: see if using two components is better
  //The penalty term is given as sqrt(k)*ln(n)
  //In this case, k=4 and n=256
  double PenTerm2 = 2.0*log(((double)max));
  for( itk::SizeValueType i=0; i<(max-1); ++i )//to set the first threshold
  {
    //compute the current parameters of the first component
    P0 = U0 = 0.0;
    for( itk::SizeValueType l=0; l<=i; ++l )
    {
      P0 += histogram.at(l);
      U0 += (l+1)*histogram.at(l);
    }
    U0 /= P0;

    for( itk::SizeValueType j=i+1; j<max; ++j )//to set the second threshold
    {
      //compute the current parameters of the second component
      P1 = U1 = 0.0;
      for( itk::SizeValueType l=j; l<=max; ++l )
      {
        P1 += histogram.at(l);
        U1 += (l+1)*histogram.at(l);
      }
      U1 /= P1;

      //compute the overall mean
      U = P0*U0 + P1*U1;

      //Compute the current value of the error criterion function
      J = U - (P0*(log(P0)+U0*log(U0))+ P1*(log(P1)+U1*log(U1)));
      //Add the penalty term
      J += PenTerm2;
      if( J<min_J )
      {
        std::cout<<"This image does not need 3 level separation!\n";
        throw;
      }
    }
  }
#ifdef DBGG
  std::cout<<"Parameters2: ";
  for( itk::SizeValueType j=0; j<parameters.size(); ++j )
    std::cout<<parameters.at(j)<<"\t";
  std::cout<<"\n"<<std::flush;
#endif //DBGG
  return;
}

void ComputeCostsFromHist( US2ImageType::IndexType &curPoint,
	std::vector< itk::SmartPointer<US2ImageType>  > &medFiltImages,
	std::vector< itk::SmartPointer<CostImageType> > &autoFlourCosts,
	std::vector< itk::SmartPointer<CostImageType> > &flourCosts,
	std::vector< itk::SmartPointer<CostImageType> > &autoFlourCostsBG,
	std::vector< itk::SmartPointer<CostImageType> > &flourCostsBG,
	std::vector< double > &histogram, itk::SizeValueType &max )
{
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  typedef itk::ImageRegionIterator< CostImageType > CostIterType;
#ifdef DBGG
  std::cout<<"computing costs from hist\n"<<std::flush;
#endif //DBGG
  std::vector< double > parameters( 5, 0 );
  computePoissonParams( histogram, parameters );
  for( itk::SizeValueType i=0; i<medFiltImages.size(); ++i )
  {
    US2ImageType::SizeType size; size[0] = 1; size[1] = 1;
    US2ImageType::RegionType region;
    region.SetSize( size ); region.SetIndex( curPoint );
    ConstIterType constIter ( medFiltImages.at(i), region );
    CostIterType costIterFlour		( flourCosts.at(i),	region );
    CostIterType costIterFlourBG	( flourCostsBG.at(i),	region );
    CostIterType costIterAutoFlour	( autoFlourCosts.at(i),	region );
    CostIterType costIterAutoFlourBG	( autoFlourCostsBG.at(i),region );
    constIter.GoToBegin(); costIterFlour.GoToBegin(); costIterFlourBG.GoToBegin();
    costIterAutoFlour.GoToBegin(); costIterAutoFlourBG.GoToBegin();
    US2ImageType::PixelType curPix = constIter.Get();
    
  }

  return;
}

void ComputeCosts( std::vector< itk::SmartPointer<US2ImageType>  > &medFiltImages,
		   std::vector< itk::SmartPointer<CostImageType> > &autoFlourCosts,
		   std::vector< itk::SmartPointer<CostImageType> > &flourCosts
		   std::vector< itk::SmartPointer<CostImageType> > &autoFlourCostsBG,
		   std::vector< itk::SmartPointer<CostImageType> > &flourCostsFG )
{
  itk::IndexValueType numCol = medFiltImages.at(0)->GetLargestPossibleRegion().GetSize()[1];
  itk::IndexValueType numRow = medFiltImages.at(0)->GetLargestPossibleRegion().GetSize()[0];
  itk::IndexValueType WinSz2 = (itk::IndexValueType)floor(((double)WinSz)/2+0.5);
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for( itk::IndexValueType i=0; i<numRow; ++i )
  {
    for( itk::IndexValueType j=0; j<numCol; ++j )
    {
      //Compute histogram at point i,j with window size define WinSz
      std::vector< double > histogram( NumBins, 0 );
      US2ImageType::IndexType start; start[0] = numRow-WinSz-1; start[1] = numCol-WinSz-1;
      if( (i-WinSz2)<0 ) start[0] = 0; if( (j-WinSz2)<0 ) start[1] = 0;
      if( start[0] && ((i+WinSz2+1)<numRow) ) start[0] = i-WinSz2;
      if( start[1] && ((j+WinSz2+1)<numCol) ) start[1] = j-WinSz2;
      itk::SizeValueType max = ComputeHistogram( medFiltImages, histogram, start );
      US2ImageType::IndexType curPoint; curPoint[0] = i; curPoint[1] = j;
#ifdef DBGG
      std::cout<<"histogram computed for " << curPoint << "\n" << std::flush;
#endif //DBGG
      ComputeCostsFromHist( curPoint, medFiltImages, autoFlourCosts,
      			    flourCosts, histogram, max );
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

  std::cout<<"Number of slices:"<<numSlices<<std::endl;

  std::vector< itk::SmartPointer<US2ImageType> > medFiltImages;
  std::vector< itk::SmartPointer<CostImageType> > autoFlourCosts, flourCosts;
  std::vector< itk::SmartPointer<CostImageType> > autoFlourCostsBG, flourCostsBG;
  medFiltImages.resize( numSlices );
  flourCosts.resize( numSlices );   autoFlourCosts.resize( numSlices );
  flourCostsBG.resize( numSlices ); autoFlourCostsBG.resize( numSlices );
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
    //Allocate space for costs
      CostImageType::Pointer costs1 = CostImageType::New();
      CostImageType::Pointer costs2 = CostImageType::New();
      CostImageType::Pointer costs3 = CostImageType::New();
      CostImageType::Pointer costs4 = CostImageType::New();
      CostImageType::PointType origin;	origin[0] = 0;	  origin[1] = 0;
      CostImageType::IndexType start;	start[0] = 0;	  start[1] = 0;
      CostImageType::SizeType size;	size[0] = numRow; size[1] = numCol;
      CostImageType::RegionType region;
      region.SetSize( size ); region.SetIndex( start );
      costs1->SetOrigin( origin );	costs2->SetOrigin( origin );
      costs3->SetOrigin( origin );	costs4->SetOrigin( origin );
      costs1->SetRegions( region );	costs2->SetRegions( region );
      costs3->SetRegions( region );	costs4->SetRegions( region );
      costs1->Allocate();		costs2->Allocate();
      costs3->Allocate();		costs4->Allocate();
      costs1->FillBuffer(0);		costs2->FillBuffer(0);
      costs3->FillBuffer(0);		costs4->FillBuffer(0);
      costs1->Update();			costs2->Update();
      costs3->Update();			costs4->Update();
      costs1->Register();		costs2->Register();
      costs3->Register();		costs4->Register();
      autoFlourCosts.at(i)   = costs1; 	flourCosts.at(i)   = costs2;
      autoFlourCostsBG.at(i) = costs3; 	flourCostsBG.at(i) = costs4;
  }

  std::cout<<"Done! Starting to compute costs\n"<<std::flush;

  ComputeCosts( medFiltImages, autoFlourCosts, flourCosts,
  				autoFlourCostsBG, flourCostsBG );

  std::cout<<"Done! Computing Cuts\n"<<std::flush;

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
/*
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
*/
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
      medFiltImages.at(i)->UnRegister();
    }
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }

  exit( EXIT_SUCCESS );
}
