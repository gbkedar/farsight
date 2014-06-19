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

//#define DEBUG_RESCALING_N_COST_EST
//#define NOISE_THR_DEBUG
//#define DEBUG_THREE_LEVEL_LABELING
//#define DEBUG_MEAN_PROJECTIONS
#define DEBUG_ELASTIC_NETS
//#define DEBUG_CORRECTION_SURFACES
#define NO_GRAPH_CUTS

#include <vector>
#include <algorithm>
#include <string>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <float.h>
#include <math.h>
#include <limits.h>
#include "new_graph.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkIntTypes.h"
#include "itkNumericTraits.h"
#include "itkExtractImageFilter.h"
#include "itkMinimumProjectionImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkSumProjectionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"
#include "itkImageDuplicator.h"
#include "itkSCIFIOImageIO.h"

#include <mlpack/core.hpp>
#include <mlpack/methods/lars/lars.hpp>
#include <mlpack/methods/linear_regression/linear_regression.hpp>

#define WinSz 128	//Histogram computed on this window
#define CWin  16	//This is half the inner window and must be a divisor of WinSz
#define NumBins 1024	//Downsampled to these number of bins
#define NN 10.0		//The bottom NN percent are used to estimate the BG
#define NNBF 75.0	//The bottom NNBF percent are used to estimate the BG
#define MinMax 100	//Number of pixels used to compute the min/max
#define HistPCMin 20.0	//Min percentage of histogram to search while 
			//updating poisson params

#define ORDER 4		//Order of the polynomial 2-4
#define numCoeffs 14	//((ORDER+1)*(ORDER+2)/2)-1 Compute and enter!
#define useCholesky true//Use Cholesky decomposition in ARMA
#define MaxIter 3	//Maximum number of iterations
#define IterThresh 0.01 //The smallest max-min that is needed to run an iteration
#define LowerNoiseThr 14000 //Noise threshold should be at least this value
#define MY_EPS 1.8E-10  //Defining a a very loose threshold change in values for diff L1 constraints
#define L2LowThr 1.0E-4 //Very little changes in the solution below this threshold in terms of sparsity
#define L1LowThr 1.0E-4 //Very little changes in the solution below this threshold in terms of coeff
#define L2IterThr 3	//Max number of times to refine L2 soln

typedef unsigned short	USPixelType;
typedef unsigned char	UCPixelType;
typedef double		CostPixelType;
const unsigned int	Dimension3 = 3;
const unsigned int	Dimension2 = 2;
typedef itk::Image< USPixelType, Dimension3 > US3ImageType;
typedef itk::Image< UCPixelType, Dimension3 > UC3ImageType;
typedef itk::Image< USPixelType, Dimension2 > US2ImageType;
typedef itk::Image< CostPixelType, Dimension2 > CostImageType;
typedef itk::Image< CostPixelType, Dimension3 > CostImageType3d;

//Output name includes directory wiz full path with the extension stripped
std::string nameTemplate;
//Stores the current iter as a string for debugging outputs
std::string iterTemplate; 

double lambda1start = 1;//Seed for lambda 1 in the elastic nets 
double L2THR = 1.0;	//Minimum L2 norm to be achieved 
double L2MINTHR = 0.1;	//Min accuracy with which L2 condition in Regression should be satisfied
itk::SizeValueType poissonCompsRepeated = 0;
itk::SizeValueType poissonCompsFailed   = 0;

void usage( const char *funcName )
{
  std::cout << "USAGE:"
	    << " " << funcName << " InputImage NumberOfColsInMeanderScan(0 to skip) "
	    << "NumberOfThreads(Optional-default=24) UseSignleLevel(Optional-default=0) "
	    << "UseNoiseThr(For max of input pixel type use 0)\n";
}

template<typename InputImageType> void WriteITKImage
  ( typename InputImageType::Pointer inputImagePointer,
    std::string outputName )
{
  typedef typename itk::ImageFileWriter< InputImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputName.c_str() );
  writer->SetInput( inputImagePointer );
  try
  {
    writer->Update();
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }
  return;
}

template<typename InputImageType> typename InputImageType::Pointer
  ReadITKImageScifio( std::string inputName )
{
  typedef typename itk::ImageFileReader< InputImageType > ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  itk::SCIFIOImageIO::Pointer io = itk::SCIFIOImageIO::New();
  reader->SetImageIO( io );
  reader->SetFileName( inputName.c_str() );
  try
  {
    reader->Update();
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }
  typename InputImageType::Pointer inputImagePointer = reader->GetOutput();
  return inputImagePointer;
}

template<typename InputImageType> typename InputImageType::Pointer
  ReadITKImage( std::string inputName )
{
  typedef typename itk::ImageFileReader< InputImageType > ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputName.c_str() );
  try
  {
    reader->Update();
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }
  typename InputImageType::Pointer inputImagePointer = reader->GetOutput();
  return inputImagePointer;
}

template<typename InputImageType, typename OutputImageType> 
  typename OutputImageType::Pointer  SumProject3dImageTo2d
  ( typename InputImageType::Pointer inputImage )
{
  typedef typename itk::SumProjectionImageFilter< InputImageType, OutputImageType >
	SumProjFilterType;
  typename SumProjFilterType::Pointer sumProjFiltFlIm = SumProjFilterType::New();
  sumProjFiltFlIm->SetInput( inputImage );
  sumProjFiltFlIm->SetProjectionDimension( 2 );
#ifdef _OPENMP
itk::MultiThreader::SetGlobalDefaultNumberOfThreads
		(itk::MultiThreader::GetGlobalMaximumNumberOfThreads());
#endif
  try
  {
    sumProjFiltFlIm->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught sum projection filter!" << excep << std::endl;
    exit (EXIT_FAILURE);
  }
#ifdef _OPENMP
itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
#endif
  typename OutputImageType::Pointer outputImage = sumProjFiltFlIm->GetOutput();
  return outputImage;
}

template<typename InputImageType, typename OutputImageType> 
  typename OutputImageType::Pointer CastImage
    ( typename InputImageType::Pointer inputImage )
{
  typedef typename itk::CastImageFilter<InputImageType,
  					OutputImageType> CastFilterType;
  if( InputImageType::ImageDimension!=OutputImageType::ImageDimension )
  {
    std::cout<<"Error! This cast function needs equal input and output dimensions";
  }
  typename CastFilterType::Pointer cast = CastFilterType::New();
  cast->SetInput( inputImage );
  try
  {
    cast->Update();
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }
  return cast->GetOutput();
}

template<typename InputImageType, typename OutputImageType> void CastNWriteImage
  ( typename InputImageType::Pointer inputImage,
    std::string &outFileName )
{
  typedef typename itk::CastImageFilter<InputImageType,
  					OutputImageType> CastFilterType;
  if( InputImageType::ImageDimension!=OutputImageType::ImageDimension )
  {
    std::cout<<"This function needs equal input and output dimensions";
    return;
  }
  typename OutputImageType::Pointer cast = 
		CastImage< InputImageType, OutputImageType >( inputImage );
  WriteITKImage< OutputImageType >( cast, outFileName );
  return;
}

template<typename InputImageType, typename OutputImageType> void
RescaleCastNWriteImage
( typename InputImageType::Pointer inputImage,
    std::string &outFileName )
{
  typedef itk::RescaleIntensityImageFilter< InputImageType, OutputImageType > RescaleFilterType;
  if( InputImageType::ImageDimension!=OutputImageType::ImageDimension )
  {
    std::cout<<"This function needs equal input and output dimensions";
    return;
  }
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(inputImage);
  rescaleFilter->SetOutputMinimum( itk::NumericTraits< typename OutputImageType::PixelType >::min() );
  rescaleFilter->SetOutputMaximum( itk::NumericTraits< typename OutputImageType::PixelType >::max() );
  WriteITKImage< OutputImageType >( rescaleFilter->GetOutput(), outFileName );
  return;
}

template<typename InputImageType> typename InputImageType::Pointer
  CreateDefaultCoordsNAllocateSpace( typename InputImageType::SizeType size )
{
  typename InputImageType::Pointer inputImagePointer = InputImageType::New();
  typename InputImageType::PointType origin;
  typename InputImageType::IndexType start;
  const int imDims = InputImageType::ImageDimension;
  for( itk::IndexValueType i=0;
       i<imDims; ++i )
  {
    origin[i] = 0; start[i] = 0;
  }
  typename InputImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  inputImagePointer->SetOrigin( origin );
  inputImagePointer->SetRegions( region );
  inputImagePointer->Allocate();
  inputImagePointer->FillBuffer(0);
  try
  {
    inputImagePointer->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught in allocating space for image!" << excep << std::endl;
    exit (EXIT_FAILURE);
  }
  return inputImagePointer;
}

template<typename InputImageType, typename OutputImageType>
  void CastNWriteImage2DStackOfVecsTo3D
   ( std::vector< typename InputImageType::Pointer > &inputImage,
     std::string &outFileName )
{
  typedef typename itk::ImageRegionConstIterator< InputImageType > ConstIterType2d;
  typedef typename itk::ImageRegionIteratorWithIndex< OutputImageType > IterType3d;
  typedef typename itk::ImageFileWriter< OutputImageType > WriterType;

  //Check if input vector is 2D and output is 3D n print error if not
  if( InputImageType::ImageDimension!=2 || OutputImageType::ImageDimension!=3 )
  {
    std::cout<<"ERROR: Utility expects a 2D image stack and outputs a 3D image.\n";
    return;
  }

  //Allocate space
  typename OutputImageType::SizeType size;
  typename OutputImageType::IndexType start;
  typename OutputImageType::RegionType region;
  start[0] = 0; start[1] = 0;
  size[0]  = inputImage.at(0)->GetLargestPossibleRegion().GetSize()[0];
  size[1]  = inputImage.at(0)->GetLargestPossibleRegion().GetSize()[1];
  size[2]  = inputImage.size();
  typename OutputImageType::Pointer outputImage =
    CreateDefaultCoordsNAllocateSpace<OutputImageType>( size );

  //Cast n Write Values
  itk::SizeValueType typeMax =
		itk::NumericTraits<typename OutputImageType::PixelType>::max();
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for( itk::IndexValueType i=0; i<inputImage.size(); ++i )
  {
    //Start index
    start[2] = i; size[2] = 1; //Reset for writing out the slices
    typename InputImageType::IndexType start2d; start2d[0] = 0;      start2d[1] = 0;
    typename InputImageType::SizeType  size2d;   size2d[0] = size[0]; size2d[1] = size[1];
    region.SetSize( size ); region.SetIndex( start );
    typename InputImageType::RegionType region2d;
    region2d.SetSize( size2d ); region2d.SetIndex( start2d );
    ConstIterType2d iter2d ( inputImage.at(i), region2d );
    IterType3d iter3d( outputImage, region );
    for( ; !iter2d.IsAtEnd(); ++iter2d, ++iter3d )
      iter3d.Set( (typename OutputImageType::PixelType)std::floor( iter2d.Get()+0.5 ) );
  }
  WriteITKImage< OutputImageType >( outputImage, outFileName );
  return;
}

US2ImageType::Pointer GetTile( US3ImageType::Pointer &readImage,
		itk::SizeValueType i )
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
    std::cerr << "Exception caught in slice extraction filter!" << excep << std::endl;
    exit (EXIT_FAILURE);
  }
  US2ImageType::Pointer currentTile = deFilter->GetOutput();
  return currentTile;
}

void ComputeGlobalHistogram(
	US3ImageType::Pointer inpImage,
	std::vector< double > &histogram,
	US3ImageType::PixelType valsPerBin )
{
  typedef itk::ImageRegionConstIterator< US3ImageType > ConstIterType;
  US3ImageType::SizeType size;
  size[0] = inpImage->GetLargestPossibleRegion().GetSize()[0];
  size[1] = inpImage->GetLargestPossibleRegion().GetSize()[1];
  size[2] = inpImage->GetLargestPossibleRegion().GetSize()[2];
  US3ImageType::IndexType start; start[0] = 0; start[1] = 0; start[2] = 0;
  US3ImageType::RegionType region;
  region.SetSize( size ); region.SetIndex( start );
  ConstIterType constIter ( inpImage, region );
  for( constIter.GoToBegin(); !constIter.IsAtEnd(); ++constIter )
    ++histogram[(itk::SizeValueType)std::floor((double)constIter.Get()/(double)valsPerBin)];
  double normalizeFactor = size[0]*size[1]*size[2];
  for( itk::SizeValueType j=0; j<histogram.size(); ++j )
  {
    histogram.at(j) /= normalizeFactor;
  }
  return;
}

void ComputeHistogram(
	std::vector< US2ImageType::Pointer  > &medFiltImages,
	std::vector< double > &histogram,
	US2ImageType::IndexType &start, US3ImageType::PixelType valsPerBin )
{
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  for( itk::SizeValueType i=0; i<medFiltImages.size(); ++i )
  {
    US2ImageType::SizeType size; size[0] = WinSz; size[1] = WinSz;
    US2ImageType::RegionType region;
    region.SetSize( size ); region.SetIndex( start );
    ConstIterType constIter ( medFiltImages.at(i), region );
    for( constIter.GoToBegin(); !constIter.IsAtEnd(); ++constIter )
      ++histogram[(itk::SizeValueType)std::floor((double)constIter.Get()/(double)valsPerBin)];
  }
  double normalizeFactor = ((double)WinSz)*((double)WinSz)*((double)medFiltImages.size());
  for( itk::SizeValueType j=0; j<histogram.size(); ++j )
  {
    histogram.at(j) /= normalizeFactor;
  }
  return;
}

void UpdateHistogram(
	std::vector< US2ImageType::Pointer  > &medFiltImages,
	std::vector< double > &histogram,
	US2ImageType::IndexType &start, US2ImageType::IndexType &prevStart,
	US3ImageType::PixelType valsPerBin )
{
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  double pixelContrib = 1.00/(((double)WinSz)*((double)WinSz)*((double)medFiltImages.size()));
  //Remove the previous col from the histogram
  US2ImageType::IndexType startRem, startAdd;
  startRem[0] = start[0]; startRem[1] = (start[1]-CWin)<0 ? 0 : (start[1]-CWin);
  US2ImageType::SizeType size;
  size[0] = WinSz; size[1] = start[1]-prevStart[1];//(start[1]-CWin)<0 ? CWin-start[1]+CWin : CWin;
  US2ImageType::RegionType RemRegion, AddRegion;
  RemRegion.SetSize( size ); RemRegion.SetIndex( startRem );
  for( itk::SizeValueType i=0; i<medFiltImages.size(); ++i )
  {
    ConstIterType constIter ( medFiltImages.at(i), RemRegion );
    for( constIter.GoToBegin(); !constIter.IsAtEnd(); ++constIter )
      histogram[(itk::SizeValueType)std::floor((double)constIter.Get()/(double)valsPerBin)]
	-=pixelContrib;
  }
  startAdd[0] = start[0]; startAdd[1] = start[1]+WinSz-size[1];
  AddRegion.SetSize( size ); AddRegion.SetIndex( startAdd );
  for( itk::SizeValueType i=0; i<medFiltImages.size(); ++i )
  {
    ConstIterType constIter ( medFiltImages.at(i), AddRegion );
    for( constIter.GoToBegin(); !constIter.IsAtEnd(); ++constIter )
      histogram[(itk::SizeValueType)std::floor((double)constIter.Get()/(double)valsPerBin)]
	+=pixelContrib;
  }
  for( itk::SizeValueType i=0; i<histogram.size(); ++i ) 
    if( histogram.at(i) < pixelContrib ) histogram.at(i) = 0;
  return;
}

void computePoissonParams( std::vector< double > &histogram,
			   std::vector< double > &parameters, bool firstPass )
{
  itk::SizeValueType max = histogram.size()-1;
  std::vector< double > parametersCopy( parameters );
  //The three-level min error thresholding algorithm
  double min_J = DBL_MAX;
  double P0, U0, P1, U1, P2, U2, U, J;
  // Try this: we need to define a penalty term that depends on the number of parameters
  //The penalty term is given as 0.5*k*ln(n)
  //where k is the number of parameters of the model and n is the number of samples
  //In this case, k=6 and n=256
  double PenTerm3 = sqrt(6.0)*log(((double)max));
  double pc5=0, pc10=0, pc15=0, pc85=0, pc90=0, histPc=0;
  for( itk::SizeValueType i=0; i<histogram.size(); ++i )
  {
    histPc += histogram.at(i);
    if( histPc > 0.05 && !pc5 ) pc5 = i ? i-1 : i;
    if( histPc > 0.10 && !pc10 ) pc10 = i ? i-1 : i;
    if( histPc > 0.15 && !pc15 ) pc15 = i ? i-1 : i;
    if( histPc > 0.85 && !pc85 ) pc85 = i ? i-1 : i;
    if( histPc > 0.90 && !pc90 ) pc90 = i ? i-1 : i;
  }
  itk::SizeValueType min_i, max_i;
  if( firstPass )
  {
    min_i = pc10;
    max_i = pc85;
  }
  else
  {
    double histPcChange = 2.0/WinSz*100;
    if( histPcChange<HistPCMin ) histPcChange=HistPCMin;
    histPcChange = std::ceil(histPcChange/100.0*((double)histogram.size()));
    min_i = (parameters.at(0)-histPcChange)<pc5 ? pc5 : (parameters.at(0)-histPcChange);
    max_i = (parameters.at(0)+histPcChange)>pc85 ?
		pc85 : (parameters.at(0)+histPcChange);
  }
  for( itk::SizeValueType i=min_i; i<max_i; ++i )//to set the first threshold
  {
    //compute the current parameters of the first component
    P0 = U0 = 0.0;
    for( itk::SizeValueType l=0; l<=i; l++ )
    {
      P0 += histogram.at(l);
      U0 += (l+1)*histogram.at(l);
    }
    U0 /= P0;

    for( itk::SizeValueType j=i+1>pc10?i+1:pc10; j<pc90; ++j )//to set the second threshold
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
        parameters.at(0) = U0; //Lowest mean
        parameters.at(1) = U1; //Intermediate mean
        parameters.at(2) = U2; //Highest mean
        parameters.at(3) = P0; //Prior for the lowest
        parameters.at(4) = P1; //Prior for the intermediate, highest will be 1-(P0+P1)
      }
    }
  }
#ifdef DEBUG_POS_EST
  std::cout<<"Parameters1: ";
  for( itk::SizeValueType j=0; j<parameters.size(); ++j )
    std::cout<<parameters.at(j)<<"\t";
  std::cout<<"\n"<<std::flush;

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
//        throw;
      }
    }
  }
  std::cout<<"Parameters2: ";
  for( itk::SizeValueType j=0; j<parameters.size(); ++j )
    std::cout<<parameters.at(j)<<"\t";
  std::cout<<"\n"<<std::flush;
#endif //DEBUG_POS_EST
  double delU0Pc = firstPass?0.0:(parameters.at(0)-parametersCopy.at(0))/ 
  			((double)histogram.size());
  double delU1Pc = firstPass?0.0:(parameters.at(1)-parametersCopy.at(1))/
  			((double)histogram.size());
  double delU2Pc = firstPass?0.0:(parameters.at(2)-parametersCopy.at(2))/
  			((double)histogram.size());

  if( parameters.at(0) >= histogram.size() ||
      parameters.at(1) >= histogram.size() ||
      parameters.at(2) >= histogram.size() || 
      delU0Pc > HistPCMin ||
      delU1Pc > HistPCMin ||
      delU2Pc > HistPCMin )
  {
#ifdef _OPENMP
#pragma omp critical
#endif //_OPENMP
    ++poissonCompsRepeated;

    if( !firstPass )
      computePoissonParams( histogram, parameters, true );
    else
    {
      if( parameters.at(0)>=histogram.size() )
      {
	parameters.at(0)=histogram.size()-3;
	parameters.at(1)=histogram.size()-2;
	parameters.at(2)=histogram.size()-1;
      }
      else if( parameters.at(1)>=histogram.size() )
      {
	parameters.at(1)=histogram.size()-2;
	parameters.at(2)=histogram.size()-1;
      }
      else parameters.at(2)=histogram.size()-1;

#ifdef _OPENMP
#pragma omp critical
#endif //_OPENMP
      ++poissonCompsFailed;

    }
  }
  return;
}

//Intialize pdf vector with max_intensity+1, 1
void ComputePoissonProbability( double &alpha, std::vector<double> &pdf )
{
  double A;
  A = exp(-alpha);
  pdf.at(0) = 1.0;
  double epsThresh = std::numeric_limits<double>::epsilon()*2;
  for( itk::SizeValueType i=1; i<pdf.size(); ++i )
  {
    pdf.at(i) = pdf.at(i-1);
    pdf.at(i) = pdf.at(i) * (alpha/((double)i));
    if( pdf.at(i) < epsThresh )
    {
      //Sine the full interval should not be more than 0-10^3
      pdf.at(i) = std::numeric_limits<double>::epsilon();
      for( itk::SizeValueType j=i+1; j<pdf.size(); ++j )
	pdf.at(j) = std::numeric_limits<double>::epsilon();
      break;
    }
  }
  for( itk::SizeValueType i=0; i<pdf.size(); ++i )
  {
    pdf.at(i) *= A;
    if( pdf.at(i) < std::numeric_limits<double>::epsilon() )
      pdf.at(i) = std::numeric_limits<double>::epsilon();
  }
  return;
}

void returnthresh
	( US2ImageType::Pointer input_image, int num_bin_levs,
	  std::vector< US2ImageType::PixelType > &returnVec )
{
  //Instantiate the different image and filter types that will be used
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIteratorType;
  typedef itk::Statistics::Histogram< float > HistogramType;
  typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;

  std::cout<<"Starting threshold computation\n";

  //Create a temporary histogram container:
  const itk::SizeValueType numBins = itk::NumericTraits<US3ImageType::PixelType>::max()+1;
  std::vector< double > tempHist( numBins, 0 );

  US3ImageType::PixelType maxval = itk::NumericTraits<US3ImageType::PixelType>::ZeroValue();
  US3ImageType::PixelType minval = itk::NumericTraits<US3ImageType::PixelType>::max();
  //Populate the histogram (assume pixel type is actually is some integer type):
  ConstIteratorType it( input_image, input_image->GetLargestPossibleRegion() );
  for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
  {
    US3ImageType::PixelType pix = it.Get();
    ++tempHist.at(pix);
    if( pix > maxval ) maxval = pix;
    if( pix < minval ) minval = pix;
  }
  //return max of type if there is no variation in the staining
  if( (maxval-minval)<3 )
  {
    for( unsigned i=0; i<returnVec.size(); ++i )
      returnVec.at(i) = itk::NumericTraits<US3ImageType::PixelType>::max();
    return;
  }
  const itk::SizeValueType numBinsPresent = maxval-minval+1;
  
  //Find max value in the histogram
  double floatIntegerMax = (double)itk::NumericTraits<US3ImageType::PixelType>::max()/2.0;
  double max = 0.0;
  for( itk::SizeValueType i=minval; i<=maxval; ++i )
    if( tempHist.at(i) > max )
      max = tempHist.at(i);

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
  upperBound.Fill((double)numBinsPresent);
  size.Fill(numBinsPresent);

  histogram->SetMeasurementVectorSize(1);
  histogram->Initialize( size, lowerBound, upperBound ) ;

  itk::SizeValueType i=minval;
  for( HistogramType::Iterator iter = histogram->Begin(); iter != histogram->End(); ++iter, ++i )
  {
    float norm_freq = (float)(tempHist.at(i) * scaleFactor);
    iter.SetFrequency(norm_freq);
  }

  std::cout<<"Histogram computed\n";

  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetNumberOfThresholds( num_bin_levs );
  calculator->SetInputHistogram( histogram );
  calculator->Update();
  std::cout<<"Threshold computed: ";
  const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
  CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();

  for(US3ImageType::PixelType i=0; i < num_bin_levs; ++itNum, ++i)
  {
    returnVec.at(i) = (static_cast<float>(*itNum))+minval;
    std::cout<<returnVec.at(i)<<std::endl;
  }
  return;
}

US3ImageType::PixelType SetSaturatedFGPixelsToMin( US3ImageType::Pointer InputImage,
				int numThreads,	US3ImageType::PixelType lowNoiseThr )
{
  typedef itk::MaximumProjectionImageFilter< US3ImageType, US2ImageType > MaxProjFilterType;
  typedef itk::MinimumProjectionImageFilter< US3ImageType, US2ImageType > MinProjFilterType;
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  typedef itk::ImageRegionIterator< US3ImageType > IterTypeUS3d;

  MaxProjFilterType::Pointer maxIntProjFilt = MaxProjFilterType::New();
  maxIntProjFilt->SetInput( InputImage );
  maxIntProjFilt->SetProjectionDimension( 2 );

  MinProjFilterType::Pointer minIntProjFilt = MinProjFilterType::New();
  minIntProjFilt->SetInput( InputImage );
  minIntProjFilt->SetProjectionDimension( 2 );

  try
  {
    minIntProjFilt->Update();
    maxIntProjFilt->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught in min/max intensity projection!" << excep << std::endl;
    exit (EXIT_FAILURE);
  }

#ifdef NOISE_THR_DEBUG
  std::string minIntPFileName = nameTemplate + "minIntProj.tif";
  WriteITKImage< US2ImageType >( minIntProjFilt->GetOutput(), minIntPFileName );
  std::string maxIntPFileName = nameTemplate + "maxIntProj.tif";
  WriteITKImage< US2ImageType >( maxIntProjFilt->GetOutput(), maxIntPFileName );
  itk::SizeValueType countCorrected = 0;
#endif //NOISE_THR_DEBUG

  std::cout<<"Size: "<< minIntProjFilt->GetOutput()->GetLargestPossibleRegion().GetSize()[0]
  	<< " " << minIntProjFilt->GetOutput()->GetLargestPossibleRegion().GetSize()[1]
	<< std::endl;

  std::vector< US2ImageType::PixelType > thresholdVec(1,0);
  returnthresh( maxIntProjFilt->GetOutput(), 1, thresholdVec );
  
  double size = InputImage->GetLargestPossibleRegion().GetSize()[0] *
  		InputImage->GetLargestPossibleRegion().GetSize()[1];
  double meanMin = 0;
  
  ConstIterType minIter( minIntProjFilt->GetOutput(),
  			 minIntProjFilt->GetOutput()->GetLargestPossibleRegion() );
  minIter.GoToBegin();
  for( ; !minIter.IsAtEnd(); ++minIter )
    meanMin += (double)minIter.Get()/size;
  meanMin = std::ceil( meanMin );
  US3ImageType::PixelType noiseThr = (US3ImageType::PixelType) thresholdVec.at(0);
  bool thresholdForSmallSpecs = false;
  bool thresholdForNoise      = false;
  if( noiseThr < lowNoiseThr ) //Don't do any noise thresholding
  {
    ConstIterType maxIter( maxIntProjFilt->GetOutput(),
			   maxIntProjFilt->GetOutput()->GetLargestPossibleRegion() );
    maxIter.GoToBegin();
    US3ImageType::PixelType maxPix = 0;
    for( ; !maxIter.IsAtEnd(); ++maxIter )
      if( maxIter.Get()>maxPix )
        maxPix = maxIter.Get();
    if( maxPix < (2*noiseThr) )
      noiseThr = maxPix;
    else
    {
      noiseThr = 2*noiseThr;
      thresholdForSmallSpecs = true;
    }
  }
  else
    thresholdForNoise = true;

  if( thresholdForNoise || thresholdForSmallSpecs )
   //Do thresholding
{
  std::cout<<"Noise threshold is: "<<noiseThr<<"\tAverage min is: "<<meanMin<<std::endl;

  itk::IndexValueType numSlices = InputImage->GetLargestPossibleRegion().GetSize()[2];
#ifdef _OPENMP
  #pragma omp parallel for num_threads(numThreads)
#endif
  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    US3ImageType::SizeType size;
    size[0] = InputImage->GetLargestPossibleRegion().GetSize()[0];
    size[1] = InputImage->GetLargestPossibleRegion().GetSize()[1];
    size[2] = 1;
    US3ImageType::IndexType start;
    start[0] = 0; start[1] = 0; start[2] = i;
    US3ImageType::RegionType region;
    region.SetSize( size ); region.SetIndex( start );
    IterTypeUS3d iter( InputImage, region );
    iter.GoToBegin();
    for( ; !iter.IsAtEnd(); ++iter )
      if( iter.Get()>noiseThr )
      {
        if(thresholdForNoise) iter.Set( meanMin );
	else iter.Set( noiseThr );
#ifdef NOISE_THR_DEBUG
	++countCorrected;
#endif
      }
  }
#ifdef NOISE_THR_DEBUG
  std::cout<<"countCorrected = "<<countCorrected<<std::endl;
  std::string noiseThrFileName = nameTemplate + "noisecorrected.tif";
  WriteITKImage< US3ImageType >( InputImage, noiseThrFileName );
#endif //NOISE_THR_DEBUG
}
  return (US3ImageType::PixelType) noiseThr;
}

void ComputeCosts( int numThreads,
		   std::vector< US2ImageType::Pointer  > &medFiltImages,
		   std::vector< double > &globalParameters,
		   std::vector< CostImageType::Pointer > &autoFlourCosts,
		   std::vector< CostImageType::Pointer > &flourCosts,
		   std::vector< CostImageType::Pointer > &autoFlourCostsBG,
		   std::vector< CostImageType::Pointer > &flourCostsBG,
		   US3ImageType::PixelType valsPerBin
#ifdef DEBUG_RESCALING_N_COST_EST
, std::vector< US2ImageType::Pointer > &resacledImages
#endif //DEBUG_RESCALING_N_COST_EST
		   )
{
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  typedef itk::ImageRegionIterator< CostImageType > CostIterType;
  typedef itk::ImageRegionIterator< US3ImageType > IterTypeUS3d;
  itk::IndexValueType numCol =  medFiltImages.at(0)->
				GetLargestPossibleRegion().GetSize()[1];
  itk::IndexValueType numRow =  medFiltImages.at(0)->
				GetLargestPossibleRegion().GetSize()[0];
  itk::IndexValueType WinSz2 = (itk::IndexValueType)floor(((double)WinSz)/2+0.5)
			      -(itk::IndexValueType)floor(((double)CWin)/2+0.5);
  itk::IndexValueType WinSzHalf = (itk::IndexValueType)floor(((double)WinSz)/2+0.5);

#ifdef DEBUG_RESCALING_N_COST_EST
  unsigned count = 0;
  clock_t start_time = clock();
#endif //DEBUG_RESCALING_N_COST_EST

#ifdef _OPENMP
  #pragma omp parallel for num_threads(numThreads)
#endif
  for( itk::IndexValueType i=0; i<numRow; i+=CWin )
  {
    std::vector< double > parameters( globalParameters );//(5,0);
    US2ImageType::IndexType prevStart; prevStart[0] = 0; prevStart[1] = 0;
    for( itk::IndexValueType j=0; j<numCol; j+=CWin )
    {
      //Compute histogram at point i,j with window size define WinSz
      US2ImageType::IndexType curPoint; curPoint[0] = i; curPoint[1] = j;
      US2ImageType::IndexType start; start[0] = i-WinSz2; start[1] = j-WinSz2;
      if( start[0]<0 ) start[0] = 0; if( start[1]<0 ) start[1] = 0;
      if( (start[0]+WinSz)>=numRow ) start[0] =  numRow-WinSz-1;
      if( (start[1]+WinSz)>=numCol ) start[1] =  numCol-WinSz-1;
      std::vector< double > histogram( NumBins, 0 ); //Move outside for loop after
      //debugging UpdateHistogram
      /*if( j==0 )
      {
        ComputeHistogram( medFiltImages, histogram, start, valsPerBin );
	computePoissonParams( histogram, parameters, true );
      }
      else if( start[1]>0 && j<(numRow-WinSzHalf-1) )
      {*/
        ComputeHistogram( medFiltImages, histogram, start, valsPerBin );
	//UpdateHistogram( medFiltImages, histogram, start, prevStart, valsPerBin );
	computePoissonParams( histogram, parameters, false );
      /*}*/
      prevStart[0] = start[0]; prevStart[1] = start[1];

      std::vector< double > pdf0( histogram.size()+1, 1 );
      ComputePoissonProbability( parameters.at(0), pdf0 );
      std::vector< double > pdf1( histogram.size()+1, 1 );
      ComputePoissonProbability( parameters.at(1), pdf1 );
      std::vector< double > pdf2( histogram.size()+1, 1 );
      ComputePoissonProbability( parameters.at(2), pdf2 );

#ifdef DEBUG_RESCALING_N_COST_EST
      if( !omp_get_thread_num() )
      {
	std::cout << "histogram computed for " << curPoint << "\t" ;
	std::cout << "computing costs from hist\t";
	std::cout << "Time: " << (clock()-start_time)/((float)CLOCKS_PER_SEC) << "\n";
	start_time = clock();
      }
#endif //DEBUG_RESCALING_N_COST_EST
      for( itk::SizeValueType k=0; k<medFiltImages.size(); ++k )
      {
	//Declare iterators for the four images
	US2ImageType::SizeType size; size[0] = CWin; size[1] = CWin;
	if( (i+CWin)>=numRow ) size[0] = numRow-i-1;
	if( (j+CWin)>=numCol ) size[1] = numCol-j-1;
	US2ImageType::RegionType region;
	region.SetSize( size ); region.SetIndex( curPoint );
	ConstIterType constIter ( medFiltImages.at(k), region );
	CostIterType costIterFlour	( flourCosts.at(k),	region );
	CostIterType costIterFlourBG	( flourCostsBG.at(k),	region );
	CostIterType costIterAutoFlour	( autoFlourCosts.at(k),	region );
	CostIterType costIterAutoFlourBG( autoFlourCostsBG.at(k),region );
	constIter.GoToBegin(); costIterFlour.GoToBegin(); costIterFlourBG.GoToBegin();
	costIterAutoFlour.GoToBegin(); costIterAutoFlourBG.GoToBegin();
#ifdef DEBUG_RESCALING_N_COST_EST
//	US2ImageType::SizeType size; size[0] = CWin; size[1] = CWin;
//	if( (i+CWin)>=numRow ) size[0] = numRow-i-1;
//	if( (j+CWin)>=numCol ) size[1] = numCol-j-1;
//	US2ImageType::RegionType region;
//	ConstIterType constIter ( medFiltImages.at(k), region );
//	constIter.GoToBegin();
	typedef itk::ImageRegionIterator< US2ImageType > IterType;
	IterType rescaleIter ( resacledImages.at(k), region );
	rescaleIter.GoToBegin();
//	for( ; !rescaleIter.IsAtEnd(); ++rescaleIter )
#endif //DEBUG_RESCALING_N_COST_EST

	for( ; !constIter.IsAtEnd(); ++constIter, ++costIterFlour, ++costIterFlourBG,
#ifdef DEBUG_RESCALING_N_COST_EST
		++costIterAutoFlour, ++costIterAutoFlourBG, ++rescaleIter )
#else  //DEBUG_RESCALING_N_COST_EST
		++costIterAutoFlour, ++costIterAutoFlourBG )
#endif //DEBUG_RESCALING_N_COST_EST
	{
	  itk::SizeValueType currentPixel = (US2ImageType::PixelType)std::floor
						( ((double)constIter.Get())/((double)valsPerBin) );
#ifdef DEBUG_RESCALING_N_COST_EST
	  rescaleIter.Set( constIter.Get() );
#endif //DEBUG_RESCALING_N_COST_EST
	  if( currentPixel >= pdf0.size() )
	    currentPixel = pdf0.size()-1;

	  //Compute node costs for each type
	  double AF, AFBG, F, FBG;
	  if( currentPixel >= parameters.at(2) )
	    F  =  ( 1 - ( parameters.at(3) + parameters.at(4) ) ) * 
		  pdf2.at( std::floor( parameters.at(2)+0.5 ) );
	  else
	    F  =  ( 1-(parameters.at(3)+parameters.at(4))) *
		  pdf2.at( currentPixel );
	  if( currentPixel >= parameters.at(1) )
	    AF =  F + ( parameters.at(4) ) *	//Easier to estimate AF+F and sub F after cuts
		  pdf1.at( std::floor( parameters.at(1)+0.5 ) );
	  else
	    AF =  F + ( parameters.at(4) ) *
		  pdf1.at( currentPixel );
	  if( currentPixel <= parameters.at(0) )
	    AFBG = parameters.at(3) *
	  	  pdf0.at( std::floor( parameters.at(0)+0.5 ) );
	  else
	    AFBG = parameters.at(3) *
		  pdf0.at( currentPixel );
	  if( currentPixel <= parameters.at(1) )
	    FBG = AFBG + parameters.at(4) *
		  pdf1.at( std::floor( parameters.at(1)+0.5 ) );
	  else
	    FBG = AFBG + parameters.at(4) *
		  pdf1.at( currentPixel );
	  if( currentPixel < 1 )
	  {
	    FBG = AFBG = 10000.0;
	    F   =  AF  = 0;
	  }
	  else
	  {
	    F    = -log( F    ); if( F    > 10000.0 ) F    = 10000.0;
	    AF   = -log( AF   ); if( AF   > 10000.0 ) AF   = 10000.0;
	    FBG  = -log( FBG  ); if( FBG  > 10000.0 ) FBG  = 10000.0;
	    AFBG = -log( AFBG ); if( AFBG > 10000.0 ) AFBG = 10000.0;
	  }
	  costIterFlour.Set( F ); 	costIterAutoFlour.Set( AF );
	  costIterFlourBG.Set( FBG );	costIterAutoFlourBG.Set( AFBG );
        }
      }
    }
  }
  std::cout<<"Number of Poisson estimation steps repeated:"<<poissonCompsRepeated
	<<std::endl;
  std::cout<<"Number of Poisson estimation steps failed:"<<poissonCompsFailed
	<<std::endl;
  return;
}

void ComputeCut( itk::IndexValueType slice,
		 std::vector< US2ImageType::Pointer  > &medFiltImages,
		 std::vector< CostImageType::Pointer > &flourCosts,
		 std::vector< CostImageType::Pointer > &flourCostsBG,
		 UC3ImageType::Pointer outputImage,
		 UC3ImageType::PixelType foregroundValue
		)
{
  double sigma = 25.0; //What! A hard coded constant check Boykov's paper!! Also check 20 in weights
  typedef itk::ImageRegionIteratorWithIndex< CostImageType > CostIterType;
  typedef itk::ImageRegionIteratorWithIndex< US2ImageType > US2IterType;
  typedef itk::ImageRegionIteratorWithIndex< UC3ImageType > UC3IterType;
  //Compute the number of nodes and edges
  itk::SizeValueType numRow   = medFiltImages.at(0)->GetLargestPossibleRegion().GetSize()[0];
  itk::SizeValueType numCol   = medFiltImages.at(0)->GetLargestPossibleRegion().GetSize()[1];
  itk::SizeValueType numNodes = numCol*numRow;
  itk::SizeValueType numEdges = 3*numCol*numRow //Down, right and diagonal
  				+ 1 - 2*numCol	//No Down At Bottom
				- 2*numRow;	//No Right At Edge
  //Declare some common iterators
  US2IterType medianIter( medFiltImages.at(slice),
  			  medFiltImages.at(slice)->GetLargestPossibleRegion() );
  CostIterType AFCostIter( flourCosts.at(slice), 
  			   flourCosts.at(slice)->GetLargestPossibleRegion() );
  CostIterType AFBGCostIter( flourCostsBG.at(slice), 
  			     flourCostsBG.at(slice)->GetLargestPossibleRegion() );
  UC3ImageType::IndexType start;start[0] = 0;	  start[1] = 0;	    start[2] = slice;
  UC3ImageType::SizeType size;   size[0] = numRow; size[1] = numCol; size[2] = 1;
  UC3ImageType::RegionType region; region.SetSize( size ); region.SetIndex( start );
  UC3IterType outputIter( outputImage, region );
#ifndef NO_GRAPH_CUTS
  typedef Graph_B < double, double, double > GraphType;
  GraphType *graph = new GraphType( numNodes, numEdges );
  //Iterate and add terminal weights
  for( itk::SizeValueType i=0; i<numRow; ++i )
  {
    for( itk::SizeValueType j=0; j<numCol; ++j )
    {
      CostImageType::IndexType index; index[0] = i; index[1] = j;
      AFCostIter.SetIndex( index ); AFBGCostIter.SetIndex( index );
      itk::SizeValueType indexCurrentNode = i*numCol+j; 
      graph->add_node();
      graph->add_tweights( indexCurrentNode, AFCostIter.Get(), AFBGCostIter.Get() );
    }
  }
  for( itk::SizeValueType i=0; i<numRow-1; ++i )
  {
    for( itk::SizeValueType j=0; j<numCol-1; ++j )
    {
      US2ImageType::IndexType index; index[0] = i; index[1] = j; medianIter.SetIndex( index );
      double currentVal = medianIter.Get();
      //Intensity discontinuity terms as edges as done in Yousef's paper
      itk::SizeValueType indexCurrentNode  = i*numCol+j;
      itk::SizeValueType indexRightNode    = indexCurrentNode+1;
      itk::SizeValueType indexBelowNode    = indexCurrentNode+numCol;
      itk::SizeValueType indexDiagonalNode = indexBelowNode+1;
      //Right
      index[0] = i; index[1] = j+1; medianIter.SetIndex( index );
      double rightCost = 20*exp(-pow(currentVal-medianIter.Get(),2)/(2*pow(sigma,2)));
      graph->add_edge( indexCurrentNode, indexRightNode, rightCost, rightCost );
      //Below
      index[0] = i+1; index[1] = j; medianIter.SetIndex( index );
      double downCost = 20*exp(-pow(currentVal-medianIter.Get(),2)/(2*pow(sigma,2)));
      graph->add_edge( indexCurrentNode, indexBelowNode, downCost, downCost );
      //Diagonal
      index[0] = i+1; index[1] = j+1; medianIter.SetIndex( index );
      double diagonalCost = 20*exp(-pow(currentVal-medianIter.Get(),2)/(2*pow(sigma,2)));
      graph->add_edge( indexCurrentNode, indexDiagonalNode, diagonalCost, diagonalCost );
    }
  }
  //Max flow:
  graph->maxflow();

  //Iterate and write out the pixels
  for( itk::SizeValueType i=0; i<numRow-1; ++i )
  {
    for( itk::SizeValueType j=0; j<numCol-1; ++j )
    {
      UC3ImageType::IndexType index; index[0] = i; index[1] = j; index[2] = slice;
      outputIter.SetIndex( index );
      itk::SizeValueType indexCurrentNode = i*numCol+j;
      if( graph->what_segment( indexCurrentNode ) != GraphType::SOURCE )
	outputIter.Set( foregroundValue );
    }
  }
  delete graph;
#else //NO_GRAPH_CUTS
  AFCostIter.GoToBegin(); AFBGCostIter.GoToBegin(); outputIter.GoToBegin();
  for( ; !AFCostIter.IsAtEnd(); ++AFCostIter, ++AFBGCostIter, ++outputIter )
    if( AFBGCostIter.Get()>AFCostIter.Get() )
      outputIter.Set( foregroundValue );
#endif //NO_GRAPH_CUTS
}

std::vector< CostImageType::Pointer >
  ComputeMeanImages ( UC3ImageType::Pointer labelImage,
	std::vector< US2ImageType::Pointer > &medFiltIms, int numThreads,
	int useSingleLev )
{
  typedef itk::ImageRegionIteratorWithIndex< CostImageType > CostIterType;
  typedef itk::ImageRegionIteratorWithIndex< CostImageType3d > CostIterType3d;
  typedef itk::ImageRegionIteratorWithIndex< US2ImageType > US2IterType;
  typedef itk::ImageRegionIteratorWithIndex< UC3ImageType > UC3IterType;

  US2ImageType::SizeType size;
  size[0] = medFiltIms.at(0)->GetLargestPossibleRegion().GetSize()[0];
  size[1] = medFiltIms.at(0)->GetLargestPossibleRegion().GetSize()[1];
  CostImageType::Pointer BGAvgIm =
		CreateDefaultCoordsNAllocateSpace<CostImageType>( size );

  US3ImageType::SizeType size3dd;
  size3dd[0] = size[0]; size3dd[1] = size[1]; size3dd[2] = numThreads;
  CostImageType3d::Pointer flAvgCounts =
		CreateDefaultCoordsNAllocateSpace<CostImageType3d>( size3dd );
  CostImageType3d::Pointer AFAvgCounts =
		CreateDefaultCoordsNAllocateSpace<CostImageType3d>( size3dd );
  CostImageType3d::Pointer flAvgIms =
		CreateDefaultCoordsNAllocateSpace<CostImageType3d>( size3dd );
  CostImageType3d::Pointer AFAvgIms =
		CreateDefaultCoordsNAllocateSpace<CostImageType3d>( size3dd );

  //The background needs a vector pixels on a 2D grid to be sorted
  std::vector< std::vector<US2ImageType::PixelType> > pixelVectForBG;
  pixelVectForBG.resize(size[0]*size[1]);
  for( itk::IndexValueType k=0; k<(size[0]*size[1]); ++k )
    pixelVectForBG.at(k).resize( medFiltIms.size() );

#ifdef _OPENMP
#if _OPENMP >= 200805L
  #pragma omp parallel for schedule(dynamic,1) num_threads(numThreads)
#else
  #pragma omp parallel for
#endif
#endif
  for( itk::IndexValueType k=0; k<medFiltIms.size(); ++k )
  {
    int tid = omp_get_thread_num();
    //Declare iterators for each image
    CostImageType3d::SizeType size3d; CostImageType3d::IndexType start3d;
    CostImageType3d::RegionType region;
    size3d[0]  = size[0];  size3d[1]  = size[1]; size3d[2]  = 1;
    start3d[0] = 0;        start3d[1] = 0;       start3d[2] = tid;
    region.SetSize( size3d ); region.SetIndex( start3d );
    //Average and count iterators
    CostIterType3d flAvgImIter( flAvgIms, region ), flAvgCountIter( flAvgCounts, region );
    CostIterType3d AFAvgImIter( AFAvgIms, region ), AFAvgCountIter( AFAvgCounts, region );
    flAvgImIter.GoToBegin(); flAvgCountIter.GoToBegin(); AFAvgImIter.GoToBegin();
    AFAvgCountIter.GoToBegin();

    //Median image iterator
    US2ImageType::SizeType size2d; US2ImageType::IndexType start2d;
    size2d[0]  = size[0]; size2d[1]  = size[1];
    start2d[0] = 0;       start2d[1] = 0;
    US2ImageType::RegionType region2d;
    region2d.SetSize( size2d ); region2d.SetIndex( start2d );
    US2IterType medFiltImIter( medFiltIms.at(k), region2d ); medFiltImIter.GoToBegin();

    //Label image iterator
    US3ImageType::IndexType startlab; startlab[0] = 0; startlab[1] = 0; startlab[2] = k;
    US3ImageType::RegionType regionlab;
    regionlab.SetSize( size3d ); regionlab.SetIndex( startlab );
    UC3IterType labelImageIter( labelImage, regionlab ); labelImageIter.GoToBegin();

    for( ; !flAvgImIter.IsAtEnd(); ++flAvgImIter, ++flAvgCountIter, ++AFAvgImIter,
				++AFAvgCountIter, ++medFiltImIter, ++labelImageIter )
    {
      itk::IndexValueType indexValue2d = labelImageIter.GetIndex()[0]*size[1] + 
						labelImageIter.GetIndex()[1];
      if( !labelImageIter.Get() )
	pixelVectForBG.at( indexValue2d ).at(k) = medFiltImIter.Get();
      else
	pixelVectForBG.at( indexValue2d ).at(k) = 
				std::numeric_limits<unsigned short>::max();

      if( labelImageIter.Get()==1 )
      {
	double tempValue = AFAvgImIter.Get()+medFiltImIter.Get();
	AFAvgImIter.Set( tempValue );
	double tempCount = AFAvgCountIter.Get()+1;
	AFAvgCountIter.Set( tempCount );
      }
      else if( labelImageIter.Get()==2 )
      {
	double tempValue = flAvgImIter.Get()+medFiltImIter.Get();
	flAvgImIter.Set( tempValue );
	double tempCount = flAvgCountIter.Get()+1;
	flAvgCountIter.Set( tempCount );
      }
    }
  }

  CostImageType::Pointer flAvgCount =
	SumProject3dImageTo2d<CostImageType3d,CostImageType>( flAvgCounts );
  CostImageType::Pointer AFAvgCount =
	SumProject3dImageTo2d<CostImageType3d,CostImageType>( AFAvgCounts );
  CostImageType::Pointer flAvgIm =
	SumProject3dImageTo2d<CostImageType3d,CostImageType>( flAvgIms );
  CostImageType::Pointer AFAvgIm =
	SumProject3dImageTo2d<CostImageType3d,CostImageType>( AFAvgIms );

  //Average and count iterators
  CostIterType  flAvgImIter( flAvgIm, flAvgIm->GetLargestPossibleRegion() ),
		AFAvgImIter( AFAvgIm, AFAvgIm->GetLargestPossibleRegion() ),
		flAvgCountIter( flAvgCount, flAvgCount->GetLargestPossibleRegion() ),
		AFAvgCountIter( AFAvgCount, AFAvgCount->GetLargestPossibleRegion() );
  flAvgImIter.GoToBegin(); flAvgCountIter.GoToBegin();
  AFAvgImIter.GoToBegin(); AFAvgCountIter.GoToBegin();
  //Divide by counts to get average
  for( ; !flAvgImIter.IsAtEnd(); ++flAvgImIter, ++flAvgCountIter, ++AFAvgImIter,
  				 ++AFAvgCountIter )
  {
    if( flAvgCountIter.Get() )
      flAvgImIter.Set( flAvgImIter.Get()/flAvgCountIter.Get() );
    else
      flAvgImIter.Set( std::numeric_limits<unsigned short>::max() );
    if( AFAvgCountIter.Get() )
      AFAvgImIter.Set( AFAvgImIter.Get()/AFAvgCountIter.Get() );
    else
      AFAvgImIter.Set( std::numeric_limits<unsigned short>::max() );
  }

  //Sort vectors to get the min of stacks at each pixel
#ifdef _OPENMP
#if _OPENMP >= 200805L
  #pragma omp parallel for schedule(dynamic,1) num_threads(numThreads)
#else
  #pragma omp parallel for
#endif
#endif
  for( itk::IndexValueType i=0; i<pixelVectForBG.size(); ++i )
    std::sort( pixelVectForBG.at(i).begin(), pixelVectForBG.at(i).end() );

  //Take the average of the bottom NN percent of the non-flour pixels
  itk::SizeValueType NNPcIndex = std::floor(((double)pixelVectForBG.size())*NN/100.0+0.5);
  if( useSingleLev ) //Use different percentage in case of brightfield
    NNPcIndex = std::floor(((double)pixelVectForBG.size())/NNBF+0.5);
  CostIterType BGAvgImIter( BGAvgIm, BGAvgIm->GetLargestPossibleRegion() ); 
  for( BGAvgImIter.GoToBegin(); !BGAvgImIter.IsAtEnd(); ++BGAvgImIter )
  {
    itk::IndexValueType indexValue2d = BGAvgImIter.GetIndex()[0]*size[1] + 
						BGAvgImIter.GetIndex()[1];
    std::vector<US2ImageType::PixelType>::iterator low;
    low = std::lower_bound( pixelVectForBG.at(indexValue2d).begin(),
			    pixelVectForBG.at(indexValue2d).end(),
			    std::numeric_limits<unsigned short>::max() );
    itk::IndexValueType validLength = low-pixelVectForBG.at(indexValue2d).begin();
    if( !validLength ) BGAvgImIter.Set( std::numeric_limits<unsigned short>::max() );
    else
    {
      double average = 0;
      for( itk::IndexValueType i=0; i<validLength; ++i )
	average += pixelVectForBG.at(indexValue2d).at(i);
      average /= validLength;
      BGAvgImIter.Set( average );
    }
  }
  std::vector< CostImageType::Pointer > returnVec; returnVec.resize(3);
  flAvgIm->Register(); returnVec.at(0) = flAvgIm;
  AFAvgIm->Register(); returnVec.at(1) = AFAvgIm;
  BGAvgIm->Register(); returnVec.at(2) = BGAvgIm;

  return returnVec;
}

std::vector< CostImageType::Pointer >
  ComputeMeanImages ( UC3ImageType::Pointer labelImage,
	US3ImageType::Pointer medFiltIms, int numThreads,
	int useSingleLev, US3ImageType::PixelType upperThreshold )
{
  typedef itk::ImageRegionIteratorWithIndex< CostImageType > CostIterType;
  typedef itk::ImageRegionIteratorWithIndex< CostImageType3d > CostIterType3d;
  typedef itk::ImageRegionIteratorWithIndex< US3ImageType > US3IterType;
  typedef itk::ImageRegionIteratorWithIndex< UC3ImageType > UC3IterType;

  US2ImageType::SizeType size;
  size[0] = medFiltIms->GetLargestPossibleRegion().GetSize()[0];
  size[1] = medFiltIms->GetLargestPossibleRegion().GetSize()[1];
  CostImageType::Pointer BGAvgIm =
		CreateDefaultCoordsNAllocateSpace<CostImageType>( size );

  US3ImageType::SizeType size3dd;
  size3dd[0] = size[0]; size3dd[1] = size[1]; size3dd[2] = numThreads;
  CostImageType3d::Pointer flAvgCounts =
		CreateDefaultCoordsNAllocateSpace<CostImageType3d>( size3dd );
  CostImageType3d::Pointer AFAvgCounts =
		CreateDefaultCoordsNAllocateSpace<CostImageType3d>( size3dd );
  CostImageType3d::Pointer flAvgIms =
		CreateDefaultCoordsNAllocateSpace<CostImageType3d>( size3dd );
  CostImageType3d::Pointer AFAvgIms =
		CreateDefaultCoordsNAllocateSpace<CostImageType3d>( size3dd );

  //The background needs a vector pixels on a 2D grid to be sorted
  std::vector< std::vector<US2ImageType::PixelType> > pixelVectForBG;
  pixelVectForBG.resize(size[0]*size[1]);
  for( itk::IndexValueType k=0; k<(size[0]*size[1]); ++k )
    pixelVectForBG.at(k).resize( medFiltIms->GetLargestPossibleRegion().GetSize()[2] );

#ifdef _OPENMP
#if _OPENMP >= 200805L
  #pragma omp parallel for schedule(dynamic,1) num_threads(numThreads)
#else
  #pragma omp parallel for
#endif
#endif
  for( itk::IndexValueType k=0; k<medFiltIms->GetLargestPossibleRegion().GetSize()[2]; ++k )
  {
    int tid = omp_get_thread_num();
    //Declare iterators for each image
    CostImageType3d::SizeType size3d; CostImageType3d::IndexType start3d;
    CostImageType3d::RegionType region;
    size3d[0]  = size[0];  size3d[1]  = size[1]; size3d[2]  = 1;
    start3d[0] = 0;        start3d[1] = 0;       start3d[2] = tid;
    region.SetSize( size3d ); region.SetIndex( start3d );
    //Average and count iterators
    CostIterType3d flAvgImIter( flAvgIms, region ), flAvgCountIter( flAvgCounts, region );
    CostIterType3d AFAvgImIter( AFAvgIms, region ), AFAvgCountIter( AFAvgCounts, region );
    flAvgImIter.GoToBegin(); flAvgCountIter.GoToBegin(); AFAvgImIter.GoToBegin();
    AFAvgCountIter.GoToBegin();

    //Label image iterator
    US3ImageType::IndexType startlab; startlab[0] = 0; startlab[1] = 0; startlab[2] = k;
    US3ImageType::RegionType regionlab;
    regionlab.SetSize( size3d ); regionlab.SetIndex( startlab );
    UC3IterType labelImageIter( labelImage, regionlab ); labelImageIter.GoToBegin();
    US3IterType medFiltImIter( medFiltIms, regionlab ); medFiltImIter.GoToBegin();

    for( ; !flAvgImIter.IsAtEnd(); ++flAvgImIter, ++flAvgCountIter, ++AFAvgImIter,
				++AFAvgCountIter, ++medFiltImIter, ++labelImageIter )
    {
      itk::IndexValueType indexValue2d = labelImageIter.GetIndex()[0]*size[1] + 
						labelImageIter.GetIndex()[1];
      if( !labelImageIter.Get() )
	pixelVectForBG.at( indexValue2d ).at(k) = medFiltImIter.Get();
      else
	pixelVectForBG.at( indexValue2d ).at(k) = 
				std::numeric_limits<unsigned short>::max();

      if( medFiltImIter.Get() < upperThreshold )
      {
      if( labelImageIter.Get()==1 )
      {
	double tempValue = AFAvgImIter.Get()+medFiltImIter.Get();
	AFAvgImIter.Set( tempValue );
	double tempCount = AFAvgCountIter.Get()+1;
	AFAvgCountIter.Set( tempCount );
      }
      else if( labelImageIter.Get()==2 )
      {
	double tempValue = flAvgImIter.Get()+medFiltImIter.Get();
	flAvgImIter.Set( tempValue );
	double tempCount = flAvgCountIter.Get()+1;
	flAvgCountIter.Set( tempCount );
      }
      }
    }
  }

  CostImageType::Pointer flAvgCount =
	SumProject3dImageTo2d<CostImageType3d,CostImageType>( flAvgCounts );
  CostImageType::Pointer AFAvgCount =
	SumProject3dImageTo2d<CostImageType3d,CostImageType>( AFAvgCounts );
  CostImageType::Pointer flAvgIm =
	SumProject3dImageTo2d<CostImageType3d,CostImageType>( flAvgIms );
  CostImageType::Pointer AFAvgIm =
	SumProject3dImageTo2d<CostImageType3d,CostImageType>( AFAvgIms );

  //Average and count iterators
  CostIterType  flAvgImIter( flAvgIm, flAvgIm->GetLargestPossibleRegion() ),
		AFAvgImIter( AFAvgIm, AFAvgIm->GetLargestPossibleRegion() ),
		flAvgCountIter( flAvgCount, flAvgCount->GetLargestPossibleRegion() ),
		AFAvgCountIter( AFAvgCount, AFAvgCount->GetLargestPossibleRegion() );
  flAvgImIter.GoToBegin(); flAvgCountIter.GoToBegin();
  AFAvgImIter.GoToBegin(); AFAvgCountIter.GoToBegin();
  //Divide by counts to get average
  for( ; !flAvgImIter.IsAtEnd(); ++flAvgImIter, ++flAvgCountIter, ++AFAvgImIter,
  				 ++AFAvgCountIter )
  {
    if( flAvgCountIter.Get() )
      flAvgImIter.Set( flAvgImIter.Get()/flAvgCountIter.Get() );
    else
      flAvgImIter.Set( std::numeric_limits<unsigned short>::max() );
    if( AFAvgCountIter.Get() )
      AFAvgImIter.Set( AFAvgImIter.Get()/AFAvgCountIter.Get() );
    else
      AFAvgImIter.Set( std::numeric_limits<unsigned short>::max() );
  }

  //Sort vectors to get the min of stacks at each pixel
#ifdef _OPENMP
#if _OPENMP >= 200805L
  #pragma omp parallel for schedule(dynamic,1) num_threads(numThreads)
#else
  #pragma omp parallel for
#endif
#endif
  for( itk::IndexValueType i=0; i<pixelVectForBG.size(); ++i )
    std::sort( pixelVectForBG.at(i).begin(), pixelVectForBG.at(i).end() );

  //Take the average of the bottom NN percent of the non-flour pixels
  itk::SizeValueType NNPcIndex = std::floor(((double)pixelVectForBG.size())/NN+0.5);
  if( useSingleLev ) //Use different percentage in case of brightfield
    NNPcIndex = std::floor(((double)pixelVectForBG.size())/NNBF+0.5);
  CostIterType BGAvgImIter( BGAvgIm, BGAvgIm->GetLargestPossibleRegion() ); 
  for( BGAvgImIter.GoToBegin(); !BGAvgImIter.IsAtEnd(); ++BGAvgImIter )
  {
    itk::IndexValueType indexValue2d = BGAvgImIter.GetIndex()[0]*size[1] + 
						BGAvgImIter.GetIndex()[1];
    std::vector<US2ImageType::PixelType>::iterator low;
    low = std::lower_bound( pixelVectForBG.at(indexValue2d).begin(),
			    pixelVectForBG.at(indexValue2d).end(),
			    std::numeric_limits<unsigned short>::max() );
    itk::IndexValueType validLength = low-pixelVectForBG.at(indexValue2d).begin();
    if( !validLength ) BGAvgImIter.Set( std::numeric_limits<unsigned short>::max() );
    else
    {
      double average = 0;
      for( itk::IndexValueType i=0; i<validLength; ++i )
	average += pixelVectForBG.at(indexValue2d).at(i);
      average /= validLength;
      BGAvgImIter.Set( average );
    }
  }
  std::vector< CostImageType::Pointer > returnVec; returnVec.resize(3);
  flAvgIm->Register(); returnVec.at(0) = flAvgIm;
  AFAvgIm->Register(); returnVec.at(1) = AFAvgIm;
  BGAvgIm->Register(); returnVec.at(2) = BGAvgIm;

  return returnVec;
}

void Regresss( arma::mat matX, arma::mat matY, std::vector<double> &outCoeffs, int numThreads,
		double lambda1, double lambda2, int task, double spacing = 0, unsigned iter = 0 )
{
  //Run many regression problems in parallel start by allocating space for output coeffs
  std::vector< std::vector< double > > coeffsParallel;
  int numCoeffsPar = numThreads<numCoeffs ? numCoeffs : numThreads;
  std::vector< double > normsPar( numCoeffsPar, 0 );
  for( int i=0; i<numCoeffsPar; ++i )
  {
    std::vector< double > tempArr( numCoeffs, 0 );
    coeffsParallel.push_back( tempArr );
  }
#ifdef _OPENMP
  #pragma omp parallel for num_threads(numThreads)
#endif //_OPENMP
  for( int i=0; i<numCoeffsPar; ++i )
  {
    double lambda1Cur = lambda1;
    i1( task==1 ) lambda1Cur /= pow( 2, ((double)i) );
    double lambda2Cur = lambda2;
    if( task==2 ) lambda2Cur /= pow( 2, ((double)i) );
    if( task==3 ) lambda2Cur = lambda2+spacing*i;
    mlpack::regression::LARS lars( useCholesky, lambda1Cur, lambda2Cur );
    arma::vec beta;
    lars.Regress( matX, matY, beta, true );
    normsPar.at(i) = arma::norm( beta, 2);
    for( itk::SizeValueType j=0; j<numCoeffs; ++j )
      coeffsParallel.at(i).at(j)=beta(j,0);
#ifdef DEBUG_ELASTIC_NETS
#pragma omp critical
  if( task>0 )
{
  std::cout<< i << "\t";
  std::cout<<"L1="<<lambda1Cur<<"\t"<<"L2="<<lambda2Cur<<"\t"
	<<"L2Norm="<<normsPar.at(i)<<"\n";
  beta.t().print();
}
#endif //DEBUG_ELASTIC_NETS
  }
  int firstGoodVal = -1;
  if( task==1 )
  { //Loosen L1 constraints till all the variables are included and only L2 constraint is active
    for( itk::SizeValueType i=0;i<coeffsParallel.size(); ++i )
    {
      unsigned zeroCount = 0;
      for( itk::SizeValueType j=0; j<coeffsParallel.at(i).size(); ++j )
      {
	if( std::abs(coeffsParallel.at(i).at(j))<MY_EPS )
	  ++zeroCount;
      }
      if( !zeroCount || ( ( lambda1/pow( 2, ((double)i) ) ) < L1LowThr ) )
      {
	bool vecCompareSame = true;
	if( i<(coeffsParallel.size()-1) )
	{//If only L2 constraint is active the coefficients should remain same with loosening L1
	  for( itk::SizeValueType j=0; j<coeffsParallel.at(i).size(); ++j )
	    if( std::abs( ( coeffsParallel.at(i).at(j)-coeffsParallel.at(i+1).at(j) )/
	        coeffsParallel.at(i).at(j) ) >= 0.1 )
	    {
	      vecCompareSame = false;
	      break;
	    }
	}
	else vecCompareSame = false; //If it is the last in the current parallel set continue
	if( vecCompareSame || ( ( lambda1/pow( 2, ((double)i) ) ) < L1LowThr ) )
	{ //Task 1 is done
	  firstGoodVal = i;
	  lambda1 /= pow( 2, ((double)i) );
	  task=2;
#ifdef DEBUG_ELASTIC_NETS
	  std::cout<<i<<" L1 set to "<<lambda1<<std::endl;
#endif //DEBUG_ELASTIC_NETS
	  break;
	}
      }
    }
    if( firstGoodVal<0 )
    {
      lambda1 /= pow( 2, ((double)(numCoeffsPar-2)) );
#ifdef DEBUG_ELASTIC_NETS
      std::cout<<"Loosening L1 constraints to Lambda 1 ="<<lambda1<<"\n";
#endif //DEBUG_ELASTIC_NETS
      Regresss( matX, matY, outCoeffs, numThreads, lambda1, lambda2, 1 );
      return;
    }
  }
  if( firstGoodVal>-1 )
  { //Adjust L2 constraints
    if( normsPar.at(firstGoodVal)>L2THR )
    {
      lambda2 *= pow( 2, ((double)(numCoeffsPar-2)) );
#ifdef DEBUG_ELASTIC_NETS
      std::cout<<"Starting with L2 constraints at Lambda 2 ="<<lambda2<<"\n";
#endif //DEBUG_ELASTIC_NETS
      Regresss( matX, matY, outCoeffs, numThreads, lambda1, lambda2, 2 );
    }
    else
    {
#ifdef DEBUG_ELASTIC_NETS
      std::cout<<"Starting with L2 constraints at Lambda 2 ="<<lambda2<<"\n";
#endif //DEBUG_ELASTIC_NETS
      Regresss( matX, matY, outCoeffs, numThreads, lambda1, lambda2, 2 );
    }
    return;
  }
  if( task==2 )
  { //Adjust L2 Parmeters till  L2norm of coeffs are around 1.0
    for( itk::SizeValueType i=0;i<coeffsParallel.size(); ++i )
    {
      unsigned zeroCount = 0;
      for( itk::SizeValueType j=0; j<coeffsParallel.at(i).size(); ++j )
      {
	if( std::abs(coeffsParallel.at(i).at(j))<MY_EPS )
	  ++zeroCount;
      }
      if( ( zeroCount && lambda1>L1LowThr ) &&
	  ( ( lambda2/pow( 2,((double)(i))) ) > L2LowThr ) )
      {
	Regresss( matX, matY, outCoeffs, numThreads, lambda1,
			lambda2/pow(2,((double)(i))), 1 );
	return;
      }
    }
    if( (normsPar.at(0)<L2THR && normsPar.at(numCoeffsPar-1)>L2THR) )
    { 
      for( itk::SizeValueType i=1; i<normsPar.size(); ++i )
      {
	if( normsPar.at(i-1)<L2THR && normsPar.at(i)>L2THR )
	{
	  lambda2 /= pow( 2, ((double)i) );
#ifdef DEBUG_ELASTIC_NETS
	  std::cout<<"L2 is near Lambda 2 ="<<lambda2<<"\n";
#endif //DEBUG_ELASTIC_NETS
	}
      }
    }
    else
    {
      if( (lambda2/pow(2,((double)(numCoeffsPar-1)))) > L2LowThr )
      {
      if( normsPar.at(0)>L2THR )
      {
	lambda2 *= pow( 2, ((double)(numCoeffsPar-2)) );
#ifdef DEBUG_ELASTIC_NETS
	std::cout<<"Tightening L2 constraints to Lambda 2 ="<<lambda2<<"\n";
#endif //DEBUG_ELASTIC_NETS
	Regresss( matX, matY, outCoeffs, numThreads, lambda1, lambda2, 2 );
      }
      else
      {
	lambda2 /= pow( 2, ((double)(numCoeffsPar-2)) );
#ifdef DEBUG_ELASTIC_NETS
	std::cout<<"Loosening L2 constraints to Lambda 2 ="<<lambda2<<"\n";
#endif //DEBUG_ELASTIC_NETS
	Regresss( matX, matY, outCoeffs, numThreads, lambda1, lambda2, 2 );
      }
      return;
      }
      else
	lambda2 /= pow( 2, ((double)(numCoeffsPar-1)) );
    }
#ifdef DEBUG_ELASTIC_NETS
    std::cout<<"Searching with Lambda 1:"<<lambda1<<"\tLambda 2 lower bound is:"<<lambda2<<std::endl;
#endif //DEBUG_ELASTIC_NETS
    spacing = (lambda2*2-lambda2)/(numCoeffsPar-1);
    if( lambda2>L2LowThr )
    {
      Regresss( matX, matY, outCoeffs, numThreads, lambda1, lambda2, 3, spacing, 1 );
      return;
    }
    else task = 3;
  }
  if( task==3 ) //The L2 constrained solution is still not close enough to the threshold tighten bounds
    for( itk::SizeValueType i=1; i<coeffsParallel.size(); ++i )
      if( (normsPar.at(i)<L2THR) || (lambda2<L2LowThr) )
      {
	if( ( (normsPar.at(i-1)-L2THR) > L2MINTHR ) && ( spacing>MY_EPS ) && ( lambda2>L2LowThr ) )
	{
	  lambda2 = lambda2+spacing*(i-1);
	  spacing /= (numCoeffsPar-1);
#ifdef DEBUG_ELASTIC_NETS
    std::cout<<"Searching with Lambda 1:"<<lambda1<<"\tLambda 2 lower bound is:"<<lambda2
	<<"\t Current Spacing is:"<<spacing<<std::endl;
#endif //DEBUG_ELASTIC_NETS
	  Regresss( matX, matY, outCoeffs, numThreads, lambda1, lambda2, 3, spacing );
	  return;
	}
	else
	{ //Accepatable solution found
	  for( int j=0;j<numCoeffs; ++j )
	    outCoeffs.at(j) = coeffsParallel.at(i).at(j);
#ifdef DEBUG_ELASTIC_NETS
	  std::cout<<"Lambda 2 is:"<<lambda2<<" L2norm is:"<<normsPar.at(i-1)<<"\n";
#endif //DEBUG_ELASTIC_NETS
	  break;
	}
      }

  //Compute plain regression error estimate
//  mlpack::regression::LinearRegression lr( matX, matY );
//  arma::vec parameters = lr.Parameters();
//  arma::vec predictions;
//  lr.Predict(matX, predictions);
//  arma::mat errorVec = predictions-matY;
//  double regressionErr = arma::norm( errorVec, 2);
//  std::cout<<"Least squares regression error:"<<regressionErr<<std::endl;
//  //parameters.t().print();

  //Compute current elastic net estimate
  arma::mat coeffs( outCoeffs.size(), 1 );
  for( int j=0;j<numCoeffs; ++j )
    coeffs(j,0) = outCoeffs.at(j);
  bool repeat = true;
  while( repeat )
  {
    //Tighten L1 and pick a reduction is RSS?--Check if theorems from regression apply
#ifdef DEBUG_ELASTIC_NETS
    std::cout<< matX.n_rows << " " << matX.n_cols << " " << coeffs.n_rows << " "<< coeffs.n_cols
	<< " are the matrix and col sizes for RSS computation\n"
	<< matY.n_rows << " " << matY.n_cols << "\n" << std::flush;
#endif //DEBUG_ELASTIC_NETS
    arma::mat predictions = matX.t()*coeffs;
    predictions -= matY;
    double currentRSS = arma::norm( predictions, 2);
    std::cout<<"Current RSS:"<<currentRSS<<"\n";
    //Add more code here
    repeat = false;
  }
#ifdef DEBUG_ELASTIC_NETS
  std::cout << "The output coefficients are: "; coeffs.t().print();
#endif //DEBUG_ELASTIC_NETS
  lambda1start = pow( 2, ( std::ceil( std::log( lambda1 ) / std::log( 2 ) ) ) );
  return;
}


typedef std::vector<double> DblVec;
double RunRegression( DblVec X, DblVec Y, DblVec X2, DblVec Y2, DblVec XY, DblVec &ImVals,
	DblVec &normConstants, unsigned normIndex //Idex to store the normalization consts
#if ORDER>2
		, DblVec X3, DblVec X2Y, DblVec XY2, DblVec Y3
#endif
#if ORDER>3
		, DblVec X4, DblVec X3Y, DblVec X2Y2, DblVec XY3, DblVec Y4
#endif
		, DblVec &outCoeffs, int numThreads )
{
//Take log of ImVals and at the same time ignore vals and indices where they are 
//equal to std::numeric_limits<unsigned short>::max()
  itk::SizeValueType countFail=0;
  //Compute valid indices and means & 2-norms
  double XMean=0, YMean=0, X2Mean=0, Y2Mean=0, XYMean=0, ImValsMean=0,
	XNorm=0, YNorm=0, X2Norm=0, Y2Norm=0, XYNorm=0, ImValsNorm=0
#if ORDER>2
	, X3Mean=0, X2YMean=0, XY2Mean=0, Y3Mean=0
	, X3Norm=0, X2YNorm=0, XY2Norm=0, Y3Norm=0
#endif
#if ORDER>3
	, X4Mean=0, X3YMean=0, X2Y2Mean=0, XY3Mean=0, Y4Mean=0
	, X4Norm=0, X3YNorm=0, X2Y2Norm=0, XY3Norm=0, Y4Norm=0
#endif
	;

  for( itk::IndexValueType i=0; i<ImVals.size(); ++i )
  {
    if( ImVals.at(i)<std::numeric_limits<unsigned short>::max() )
    {
	ImVals.at(i) = ImVals.at(i) ? log( ImVals.at(i) ) : 0; //Changing to log later
	XMean+=X.at(i); YMean+=Y.at(i); X2Mean+=X2.at(i); Y2Mean+=Y2.at(i); XYMean+=XY.at(i);
	ImValsMean+=ImVals.at(i);
#if ORDER>2
	X3Mean+=X3.at(i); X2YMean+=X2Y.at(i); XY2Mean+=XY2.at(i); Y3Mean+=Y3.at(i);
#endif
#if ORDER>3
	X4Mean+=X4.at(i); X3YMean+=X3Y.at(i); X2Y2Mean+=X2Y2.at(i); XY3Mean+=XY3.at(i);	Y4Mean+=Y4.at(i);
#endif
    }
    else
      ++countFail;
  }

  double countIndex = ImVals.size()-countFail;
  XMean/=countIndex; YMean/=countIndex;
  X2Mean/=countIndex; Y2Mean/=countIndex; XYMean/=countIndex; ImValsMean/=countIndex;
  X3Mean/=countIndex; X2YMean/=countIndex; XY2Mean/=countIndex; Y3Mean/=countIndex;
  X4Mean/=countIndex; X3YMean/=countIndex; X2Y2Mean/=countIndex; XY3Mean/=countIndex; Y4Mean/=countIndex;

  for( itk::IndexValueType i=0; i<ImVals.size(); ++i )
  {
    if( ImVals.at(i)<std::numeric_limits<unsigned short>::max() )
    {
	X.at(i)-=XMean; Y.at(i)-=YMean; X2.at(i)-=X2Mean; Y2.at(i)-=Y2Mean; XY.at(i)-=XYMean;
	ImVals.at(i) -= ImValsMean;
	XNorm+=(X.at(i)*X.at(i)); YNorm+=(Y.at(i)*Y.at(i)); X2Norm+=(X2.at(i)*X2.at(i));
	Y2Norm+=(Y2.at(i)*Y2.at(i)); XYNorm+=(XY.at(i)*XY.at(i));
	ImValsNorm+=(ImVals.at(i)*ImVals.at(i));
#if ORDER>2
	X3.at(i)-=X3Mean; X2Y.at(i)-=X2YMean; XY2.at(i)-=XY2Mean; Y3.at(i)-=Y3Mean;
	X3Norm+=(X3.at(i)*X3.at(i)); X2YNorm+=(X2Y.at(i)*X2Y.at(i)); XY2Norm+=(XY2.at(i)*XY2.at(i));
	Y3Norm+=(Y3.at(i)*Y3.at(i));
#endif
#if ORDER>3
	X4.at(i)-=X4Mean; X3Y.at(i)-=X3YMean; X2Y2.at(i)-=X2Y2Mean; XY3.at(i)-=XY3Mean;	Y4.at(i)-=Y4Mean;
	X4Norm+=(X4.at(i)*X4.at(i)); X3YNorm+=(X3Y.at(i)*X3Y.at(i)); X2Y2Norm+=(X2Y2.at(i)*X2Y2.at(i));
	XY3Norm+=(XY3.at(i)*XY3.at(i)); Y4Norm+=(Y4.at(i)*Y4.at(i));
#endif
    }
  }

  XNorm=sqrt(XNorm); YNorm=sqrt(YNorm); X2Norm=sqrt(X2Norm); Y2Norm=sqrt(Y2Norm); XYNorm=sqrt(XYNorm);
  ImValsNorm=sqrt(ImValsNorm);
#if ORDER>2
  X3Norm=sqrt(X3Norm); X2YNorm=sqrt(X2YNorm); XY2Norm=sqrt(XY2Norm); Y3Norm=sqrt(Y3Norm);
#endif
#if ORDER>3
  X4Norm=sqrt(X4Norm); X3YNorm=sqrt(X3YNorm); X2Y2Norm=sqrt(X2Y2Norm); XY3Norm=sqrt(XY3Norm);
  Y4Norm=sqrt(Y4Norm);
#endif

  //Store the normalization constants
  normConstants.at(normIndex)=XMean; normConstants.at(normIndex+1)=YMean;
  normConstants.at(normIndex+2)=X2Mean; normConstants.at(normIndex+3)=Y2Mean;
  normConstants.at(normIndex+4)=XYMean;
  normConstants.at(normIndex+numCoeffs)=XNorm; normConstants.at(normIndex+numCoeffs+1)=YNorm;
  normConstants.at(normIndex+numCoeffs+2)=X2Norm; normConstants.at(normIndex+numCoeffs+3)=Y2Norm;
  normConstants.at(normIndex+numCoeffs+4)=XYNorm;
#if ORDER>2
  normConstants.at(normIndex+5)=X3Mean; normConstants.at(normIndex+6)=X2YMean;
  normConstants.at(normIndex+7)=XY2Mean; normConstants.at(normIndex+8)=Y3Mean;
  normConstants.at(normIndex+numCoeffs+5)=X3Norm; normConstants.at(normIndex+numCoeffs+6)=X2YNorm;
  normConstants.at(normIndex+numCoeffs+7)=XY2Norm; normConstants.at(normIndex+numCoeffs+8)=Y3Norm;
#endif
#if ORDER>3
  normConstants.at(normIndex+9)=X4Mean; normConstants.at(normIndex+10)=X3YMean;
  normConstants.at(normIndex+11)=X2Y2Mean; normConstants.at(normIndex+12)=XY3Mean;
  normConstants.at(normIndex+13)=Y4Mean;
  normConstants.at(normIndex+numCoeffs+9)=X4Norm; normConstants.at(normIndex+numCoeffs+10)=X3YNorm;
  normConstants.at(normIndex+numCoeffs+11)=X2Y2Norm; normConstants.at(normIndex+numCoeffs+12)=XY3Norm;
  normConstants.at(normIndex+numCoeffs+13)=Y4Norm;
#endif

  arma::mat matX( numCoeffs, ImVals.size()-countFail );
  arma::mat matY( ImVals.size()-countFail, 1 );
  std::cout<<"Norm of predictor before scaling:"<<ImValsNorm<<std::endl;

  for( itk::SizeValueType i=0,j=0; i<ImVals.size(); ++i )
  {
    if( ImVals.at(i)<std::numeric_limits<unsigned short>::max() )
    {
	X.at(i)/=XNorm; Y.at(i)/=YNorm; X2.at(i)/=X2Norm; Y2.at(i)/=Y2Norm; XY.at(i)/=XYNorm;
	matX(0,j)=X.at(i); matX(1,j)=Y.at(i); matX(2,j)=X2.at(i); matX(3,j)=Y2.at(i); matX(4,j)=XY.at(i);
	ImVals.at(i)/=ImValsNorm;
	matY(j,0)=ImVals.at(i);
#if ORDER>2
	X3.at(i)/=X3Norm; X2Y.at(i)/=X2YNorm; XY2.at(i)/=XY2Norm; Y3.at(i)/=Y3Norm;
	matX(5,j)=X3.at(i); matX(6,j)=X2Y.at(i); matX(7,j)=XY2.at(i); matX(8,j)=Y3.at(i);
#endif
#if ORDER>3
	X4.at(i)/=X4Norm; X3Y.at(i)/=X3YNorm; X2Y2.at(i)/=X2Y2Norm; XY3.at(i)/=XY3Norm;	Y4.at(i)/=Y4Norm;
	matX(9,j)=X4.at(i); matX(10,j)=X3Y.at(i); matX(11,j)=X2Y2.at(i); matX(12,j)=XY3.at(i);
	matX(13,j)=Y4.at(i);
#endif
	++j;
    }
  }
  Regresss( matX, matY, outCoeffs, numThreads, lambda1start, 1.0, 1 );
  return ImValsNorm;
}

void ComputePolynomials( CostImageType::Pointer flAvgIm, CostImageType::Pointer AFAvgIm,
  CostImageType::Pointer BGAvgIm, std::vector<double> &flPolyCoeffs,
  std::vector<double> &AFPolyCoeffs, std::vector<double> &BGPolyCoeffs,
  std::vector<double> &normConstants, std::vector<double> &imNorms,
  int numThreads, int useSingleLev )
{
  typedef itk::ImageRegionIteratorWithIndex< CostImageType > CostIterType;
  std::vector<double>
	X, Y, X2, Y2, XY, FlVals, AFVals, BGVals
#if ORDER>2
	, X3, X2Y, XY2, Y3
#endif
#if ORDER>3
	, X4, X3Y, X2Y2, XY3, Y4
#endif
	;
  //Iterate through the images and pick out the indices and values
  CostIterType  flAvgImIter( flAvgIm, flAvgIm->GetLargestPossibleRegion() ),
		BGAvgImIter( BGAvgIm, BGAvgIm->GetLargestPossibleRegion() ),
		AFAvgImIter( AFAvgIm, AFAvgIm->GetLargestPossibleRegion() );
  flAvgImIter.GoToBegin(); BGAvgImIter.GoToBegin(); AFAvgImIter.GoToBegin();
  for( ; !flAvgImIter.IsAtEnd(); ++flAvgImIter, ++BGAvgImIter, ++AFAvgImIter )
  {
    X.push_back( flAvgImIter.GetIndex()[1] ); Y.push_back( flAvgImIter.GetIndex()[0] );
    X2.push_back( X.back()*X.back() ); Y2.push_back( Y.back()*Y.back() );
    XY.push_back( X.back()*Y.back() );
    FlVals.push_back( flAvgImIter.Get() ); AFVals.push_back( AFAvgImIter.Get() );
    BGVals.push_back( BGAvgImIter.Get() );
#if ORDER>2
    X3.push_back( X2.back()*X.back() ); X2Y.push_back( X2.back()*Y.back() );
    XY2.push_back( X.back()*Y2.back() ); Y3.push_back( Y.back()*Y2.back() );
#endif
#if ORDER>3
    X4.push_back( X2.back()*X2.back() ); X3Y.push_back( X3.back()*Y.back() );
    X2Y2.push_back( X2.back()*Y2.back() ); XY3.push_back( X.back()*Y3.back() );
    Y4.push_back( Y2.back()*Y2.back() );
#endif
  }
  if( !useSingleLev )
{
  imNorms.at(0) = RunRegression( X, Y, X2, Y2, XY, FlVals, normConstants,
  				0 //Idex to store the normalization consts
#if ORDER>2
				, X3, X2Y, XY2, Y3
#endif
#if ORDER>3
				, X4, X3Y, X2Y2, XY3, Y4
#endif
		, flPolyCoeffs, numThreads );
  imNorms.at(1) = RunRegression( X, Y, X2, Y2, XY, AFVals, normConstants,
  				2*numCoeffs //Index to store the normalization consts
#if ORDER>2
				, X3, X2Y, XY2, Y3
#endif
#if ORDER>3
				, X4, X3Y, X2Y2, XY3, Y4
#endif
				, AFPolyCoeffs, numThreads );
}
  imNorms.at(2) = RunRegression( X, Y, X2, Y2, XY, BGVals, normConstants,
  				4*numCoeffs //Idex to store the normalization consts
#if ORDER>2
				, X3, X2Y, XY2, Y3
#endif
#if ORDER>3
				, X4, X3Y, X2Y2, XY3, Y4
#endif
				, BGPolyCoeffs, numThreads );
  return;
}

struct IndexStructType
{ double X, Y, FlVals, AFVals, BGVals;
};

void GetSurfaceForIndices( US3ImageType::Pointer inputImage,
	std::vector< IndexStructType > &IndexVector, std::vector<double> &normConstants,
	std::vector<double> &flPolyCoeffs, std::vector<double> &AFPolyCoeffs,
	std::vector<double> &BGPolyCoeffs, std::vector<double> &imNorms )
{
  typedef itk::ImageRegionIteratorWithIndex< US3ImageType > IterType3d;
  US3ImageType::SizeType size;
  size[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
  size[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];

  US3ImageType::IndexType start3d;
  US3ImageType::SizeType  size3d;
  US3ImageType::RegionType region3d;
  start3d[0] = 0;      start3d[1] = 0;       start3d[2] = 0;
   size3d[0] = size[0]; size3d[1] = size[1];  size3d[2] = 1;
  region3d.SetSize( size3d ); region3d.SetIndex( start3d );
  IterType3d iter3d( inputImage, region3d );
  for( iter3d.GoToBegin(); !iter3d.IsAtEnd(); ++iter3d )
  {
    US3ImageType::IndexType index = iter3d.GetIndex();
    IndexStructType currentIndex; //Check order in the estimation code
    currentIndex.X = index[1];
    currentIndex.Y = index[0];
    double X2 = currentIndex.X*currentIndex.X;
    double Y2 = currentIndex.Y*currentIndex.Y;
    double XY = currentIndex.X*currentIndex.Y;
#if ORDER>2
    double X3  = X2*currentIndex.X;
    double X2Y = X2*currentIndex.Y;
    double XY2 = currentIndex.X*Y2;
    double Y3  = currentIndex.Y*Y2;
#endif
#if ORDER>3
    double X4   = X2*X2;
    double X3Y  = X3*currentIndex.Y;
    double X2Y2 = X2*Y2;
    double XY3  = currentIndex.X*Y3;
    double Y4   = Y2*Y2;
#endif
    currentIndex.FlVals = (
    			flPolyCoeffs.at(0)*((currentIndex.X-normConstants.at(0))
			/normConstants.at(numCoeffs))
+flPolyCoeffs.at(1)*((currentIndex.Y -normConstants.at(1))/normConstants.at(numCoeffs+1))
+flPolyCoeffs.at(2)*((X2-normConstants.at(2))/normConstants.at(numCoeffs+2))
+flPolyCoeffs.at(3)*((Y2-normConstants.at(3))/normConstants.at(numCoeffs+3))
+flPolyCoeffs.at(4)*((XY-normConstants.at(4))/normConstants.at(numCoeffs+4))
#if ORDER>2
+flPolyCoeffs.at(5)*((X3 -normConstants.at(5))/normConstants.at(numCoeffs+5))
+flPolyCoeffs.at(6)*((X2Y-normConstants.at(6))/normConstants.at(numCoeffs+6))
+flPolyCoeffs.at(7)*((XY2-normConstants.at(7))/normConstants.at(numCoeffs+7))
+flPolyCoeffs.at(8)*((Y3 -normConstants.at(8))/normConstants.at(numCoeffs+8))
#endif
#if ORDER>3
+flPolyCoeffs.at(9) *((X4  -normConstants.at(9))/normConstants.at(numCoeffs+9))
+flPolyCoeffs.at(10)*((X3Y -normConstants.at(10))/normConstants.at(numCoeffs+10))
+flPolyCoeffs.at(11)*((X2Y2-normConstants.at(11))/normConstants.at(numCoeffs+11))
+flPolyCoeffs.at(12)*((XY3 -normConstants.at(12))/normConstants.at(numCoeffs+12))
+flPolyCoeffs.at(13)*((Y4  -normConstants.at(13))/normConstants.at(numCoeffs+13))
) *
imNorms.at(0)
#endif
	;
    currentIndex.AFVals = (
			AFPolyCoeffs.at(0)*((currentIndex.X-normConstants.at(2*numCoeffs))
			/normConstants.at(3*numCoeffs))
+AFPolyCoeffs.at(1)*((currentIndex.Y -normConstants.at(2*numCoeffs+1))/normConstants.at(3*numCoeffs+1))
+AFPolyCoeffs.at(2)*((X2-normConstants.at(2*numCoeffs+2))/normConstants.at(3*numCoeffs+2))
+AFPolyCoeffs.at(3)*((Y2-normConstants.at(2*numCoeffs+3))/normConstants.at(3*numCoeffs+3))
+AFPolyCoeffs.at(4)*((XY-normConstants.at(2*numCoeffs+4))/normConstants.at(3*numCoeffs+4))
#if ORDER>2
+AFPolyCoeffs.at(5)*((X3 -normConstants.at(2*numCoeffs+5))/normConstants.at(3*numCoeffs+5))
+AFPolyCoeffs.at(6)*((X2Y-normConstants.at(2*numCoeffs+6))/normConstants.at(3*numCoeffs+6))
+AFPolyCoeffs.at(7)*((XY2-normConstants.at(2*numCoeffs+7))/normConstants.at(3*numCoeffs+7))
+AFPolyCoeffs.at(8)*((Y3 -normConstants.at(2*numCoeffs+8))/normConstants.at(3*numCoeffs+8))
#endif
#if ORDER>3
+AFPolyCoeffs.at(9) *((X4  -normConstants.at(2*numCoeffs+9))/normConstants.at(3*numCoeffs+9))
+AFPolyCoeffs.at(10)*((X3Y -normConstants.at(2*numCoeffs+10))/normConstants.at(3*numCoeffs+10))
+AFPolyCoeffs.at(11)*((X2Y2-normConstants.at(2*numCoeffs+11))/normConstants.at(3*numCoeffs+11))
+AFPolyCoeffs.at(12)*((XY3 -normConstants.at(2*numCoeffs+12))/normConstants.at(3*numCoeffs+12))
+AFPolyCoeffs.at(13)*((Y4  -normConstants.at(2*numCoeffs+13))/normConstants.at(3*numCoeffs+13))
) *
imNorms.at(1)
#endif
	;
    currentIndex.BGVals = (
			BGPolyCoeffs.at(0)*((currentIndex.X-normConstants.at(4*numCoeffs))
			/normConstants.at(5*numCoeffs))
+BGPolyCoeffs.at(1)*((currentIndex.Y -normConstants.at(4*numCoeffs+1))/normConstants.at(5*numCoeffs+1))
+BGPolyCoeffs.at(2)*((X2-normConstants.at(4*numCoeffs+2))/normConstants.at(5*numCoeffs+2))
+BGPolyCoeffs.at(3)*((Y2-normConstants.at(4*numCoeffs+3))/normConstants.at(5*numCoeffs+3))
+BGPolyCoeffs.at(4)*((XY-normConstants.at(4*numCoeffs+4))/normConstants.at(5*numCoeffs+4))
#if ORDER>2
+BGPolyCoeffs.at(5)*((X3 -normConstants.at(4*numCoeffs+5))/normConstants.at(5*numCoeffs+5))
+BGPolyCoeffs.at(6)*((X2Y-normConstants.at(4*numCoeffs+6))/normConstants.at(5*numCoeffs+6))
+BGPolyCoeffs.at(7)*((XY2-normConstants.at(4*numCoeffs+7))/normConstants.at(5*numCoeffs+7))
+BGPolyCoeffs.at(8)*((Y3 -normConstants.at(4*numCoeffs+8))/normConstants.at(5*numCoeffs+8))
#endif
#if ORDER>3
+BGPolyCoeffs.at(9) *((X4   -normConstants.at(4*numCoeffs+9))/normConstants.at(5*numCoeffs+9))
+BGPolyCoeffs.at(10)*((X3Y -normConstants.at(4*numCoeffs+10))/normConstants.at(5*numCoeffs+10))
+BGPolyCoeffs.at(11)*((X2Y2-normConstants.at(4*numCoeffs+11))/normConstants.at(5*numCoeffs+11))
+BGPolyCoeffs.at(12)*((XY3 -normConstants.at(4*numCoeffs+12))/normConstants.at(5*numCoeffs+12))
+BGPolyCoeffs.at(13)*((Y4  -normConstants.at(4*numCoeffs+13))/normConstants.at(5*numCoeffs+13))
) *
imNorms.at(2)
#endif
    ;
    IndexVector.push_back( currentIndex );
  }
  return;
}

bool IndMin( const IndexStructType &a, const IndexStructType &b)
{ if( a.Y==b.Y ) return (a.X < b.X); else return (a.Y < b.Y); }

double CorrectImages( std::vector<double> &flPolyCoeffs,
	std::vector<double> &AFPolyCoeffs, std::vector<double> &BGPolyCoeffs,
	std::vector<double> &normConstants, US3ImageType::Pointer inputImage,
	UC3ImageType::Pointer labelImage, CostImageType::Pointer flAvgIm,
	CostImageType::Pointer AFAvgIm,	CostImageType::Pointer BGAvgIm,
	std::vector<double> &imNorms, int numThreads, int useSingleLev,
	US3ImageType::PixelType noiseThr )
{
  typedef itk::ImageRegionIteratorWithIndex< CostImageType > CostIterType;
  typedef itk::ImageRegionIteratorWithIndex< US3ImageType  > InputIterType;
  typedef itk::ImageRegionIteratorWithIndex< UC3ImageType  > LabelIterType;
  US3ImageType::SizeType size;
  size[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
  size[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
  size[2] = inputImage->GetLargestPossibleRegion().GetSize()[2];
  itk::SizeValueType numSlices = size[2];
  std::vector< IndexStructType > IndexVector;
  GetSurfaceForIndices( inputImage, IndexVector, normConstants, flPolyCoeffs, AFPolyCoeffs,
  			BGPolyCoeffs, imNorms );

  std::sort( IndexVector.begin(), IndexVector.end(), IndMin );

  //Rescale to between 0 and max-min for all three surfaces
  double flMinSurface, flMaxSurface, AFMaxSurface, AFMinSurface, BGMaxSurface, BGMinSurface;
  flMinSurface = AFMinSurface = BGMinSurface = std::numeric_limits< double >::max();
  flMaxSurface = AFMaxSurface = BGMaxSurface = -std::numeric_limits< double >::max();
  for( itk::IndexValueType i=0; i<IndexVector.size(); ++i )
  {
    if( flMinSurface>IndexVector.at(i).FlVals ) flMinSurface=IndexVector.at(i).FlVals;
    if( flMaxSurface<IndexVector.at(i).FlVals ) flMaxSurface=IndexVector.at(i).FlVals;
    if( AFMinSurface>IndexVector.at(i).AFVals ) AFMinSurface=IndexVector.at(i).AFVals;
    if( AFMaxSurface<IndexVector.at(i).AFVals ) AFMaxSurface=IndexVector.at(i).AFVals;
    if( BGMinSurface>IndexVector.at(i).BGVals ) BGMinSurface=IndexVector.at(i).BGVals;
    if( BGMaxSurface<IndexVector.at(i).BGVals ) BGMaxSurface=IndexVector.at(i).BGVals;
  }
  double delta = flMaxSurface-flMinSurface;
  if( (AFMaxSurface-AFMinSurface)>delta ) delta = AFMaxSurface-AFMinSurface;
  if( (BGMaxSurface-BGMinSurface)>delta ) delta = BGMaxSurface-BGMinSurface;
  if( delta<IterThresh )
    return delta;
  CostImageType::SizeType size2d; size2d[0] = size[0]; size2d[1] = size[1];
  CostImageType::Pointer flSurf = CreateDefaultCoordsNAllocateSpace<CostImageType>( size2d );
  CostImageType::Pointer AFSurf = CreateDefaultCoordsNAllocateSpace<CostImageType>( size2d );
  CostImageType::Pointer BGSurf = CreateDefaultCoordsNAllocateSpace<CostImageType>( size2d );
  CostIterType flIter( flSurf, flSurf->GetLargestPossibleRegion() );
  CostIterType AFIter( AFSurf, AFSurf->GetLargestPossibleRegion() );
  CostIterType BGIter( BGSurf, BGSurf->GetLargestPossibleRegion() );
  for( itk::IndexValueType i=0; i<IndexVector.size(); ++i )
  {
    CostImageType::IndexType index; index[0] = IndexVector.at(i).Y; index[1] = IndexVector.at(i).X;
    flIter.SetIndex( index ); AFIter.SetIndex( index ); BGIter.SetIndex( index );
    flIter.Set( (IndexVector.at(i).FlVals-flMinSurface) );
    AFIter.Set( (IndexVector.at(i).AFVals-AFMinSurface) );
    BGIter.Set( (IndexVector.at(i).BGVals-BGMinSurface) );
  }

#ifdef DEBUG_CORRECTION_SURFACES
  std::cout<<"The min and maxes for the surfaces are:"<<iterTemplate<<"\n"<<
  "flMin="<<flMinSurface<<"\tflMax"<<flMaxSurface<<"\n"<<
  "AFMin="<<AFMinSurface<<"\tAFMax"<<AFMaxSurface<<"\n"<<
  "BGMin="<<BGMinSurface<<"\tBGMax"<<BGMaxSurface<<"\n";
  std::string flSurfName = nameTemplate+iterTemplate+"flSurface.tif";
  std::string AFSurfName = nameTemplate+iterTemplate+"AFSurface.tif";
  std::string BGSurfName = nameTemplate+iterTemplate+"BGSurface.tif";
  std::cout<<"Writing correction surfaces\n";
if( !useSingleLev )
{
  RescaleCastNWriteImage<CostImageType,US2ImageType>(flSurf,flSurfName);
  RescaleCastNWriteImage<CostImageType,US2ImageType>(AFSurf,AFSurfName);
}
  RescaleCastNWriteImage<CostImageType,US2ImageType>(BGSurf,BGSurfName);
#endif //DEBUG_CORRECTION_SURFACES

  size[2] = 1; //Process one slice at a time
#ifdef _OPENMP
#if _OPENMP >= 200805L
  #pragma omp parallel for schedule(dynamic,1) num_threads(numThreads)
#else
  #pragma omp parallel for
#endif
#endif
  for( itk::IndexValueType k=0; k<numSlices; ++k )
  {
    US3ImageType::IndexType index; index[0]=0; index[1]=0; index[2]=k;
    US3ImageType::RegionType region; region.SetSize( size ); region.SetIndex( index );
    InputIterType inputIter( inputImage, region );
    LabelIterType labelIter( labelImage, region );
    CostIterType flIterPerThr( flSurf, flSurf->GetLargestPossibleRegion() );
    CostIterType AFIterPerThr( AFSurf, AFSurf->GetLargestPossibleRegion() );
    CostIterType BGIterPerThr( BGSurf, BGSurf->GetLargestPossibleRegion() );
    inputIter.GoToBegin(); labelIter.GoToBegin();
    flIterPerThr.GoToBegin(); AFIterPerThr.GoToBegin(); BGIterPerThr.GoToBegin();
    for( ; !flIterPerThr.IsAtEnd(); ++inputIter, ++labelIter,
    				++flIterPerThr, ++AFIterPerThr, ++BGIterPerThr )
    {
      double curPix = log(inputIter.Get());
      double minusVal = 0;
      if( !labelIter.Get() ) minusVal = BGIterPerThr.Get();
      if( labelIter.Get()==1 ) minusVal = AFIterPerThr.Get();
      if( labelIter.Get()==2 || (noiseThr<inputIter.Get()) ) minusVal = flIterPerThr.Get();
      curPix -= minusVal;
      curPix = floor(exp(curPix)+0.5);
      inputIter.Set( (US3ImageType::PixelType)curPix );
    }
  }
  return delta;
}

void DivideImageByScanDirection( US3ImageType::Pointer Input3dImage, int numColsForDiv,
		std::vector< US3ImageType::Pointer  > &dividedOutputs )
{
  typedef itk::ImageRegionIteratorWithIndex< US3ImageType  > OutputIterType;
  typedef itk::ImageRegionConstIteratorWithIndex< US3ImageType  > InputIterConstType;
  itk::SizeValueType numSlicesFull = Input3dImage->GetLargestPossibleRegion().GetSize()[2];
  itk::SizeValueType numCol = Input3dImage->GetLargestPossibleRegion().GetSize()[1];
  itk::SizeValueType numRow = Input3dImage->GetLargestPossibleRegion().GetSize()[0];
  US3ImageType::SizeType size; size[0] = numRow; size[1] = numCol;
  if( (numSlicesFull/numColsForDiv)%2 )
    size[2] = std::floor((((double)numSlicesFull)/(2.0*(double)numColsForDiv))+0.5)
		*numColsForDiv;
  else size[2] = numSlicesFull/2;
  US3ImageType::Pointer dir1Im = CreateDefaultCoordsNAllocateSpace<US3ImageType>( size );
  if( (numSlicesFull/numColsForDiv)%2 )
    size[2] = std::floor(((double)numSlicesFull)/(2.0*(double)numColsForDiv))*numColsForDiv;
  else size[2] = numSlicesFull/2;
  US3ImageType::Pointer dir2Im = CreateDefaultCoordsNAllocateSpace<US3ImageType>( size );
  dividedOutputs.at(0) = dir1Im; dir1Im->Register();
  dividedOutputs.at(1) = dir2Im; dir2Im->Register();
  bool first = true;
  itk::IndexValueType firstStackIndex = 0;
  itk::IndexValueType secondStackIndex = 0;
  //Copy slices in alternate rows
  for( unsigned i=0; i<numSlicesFull; i+=numColsForDiv )
  {
    for( unsigned j=i; j<(i+numColsForDiv); ++j )
    {
      US3ImageType::IndexType inpIndex; inpIndex[0]=0; inpIndex[1]=0; inpIndex[2]=j;
      US3ImageType::IndexType outIndex; outIndex[0]=0; outIndex[1]=0;
      if( first ) outIndex[2] = firstStackIndex++;
      else outIndex[2] = secondStackIndex++;
      US3ImageType::SizeType sliceSize;
      sliceSize[0] = numRow; sliceSize[1] = numCol; sliceSize[2] = 1;
      US3ImageType::RegionType regionInput, regionOutput;
      regionInput.SetSize( sliceSize ); regionInput.SetIndex( inpIndex );
      regionOutput.SetSize( sliceSize ); regionOutput.SetIndex( outIndex );
      InputIterConstType inputIter( Input3dImage, regionInput );
      if( first )
      {
	OutputIterType outputIter( dir1Im, regionOutput );
	outputIter.GoToBegin(); inputIter.GoToBegin();
	for( ; !inputIter.IsAtEnd() && !outputIter.IsAtEnd(); ++inputIter, ++outputIter )
	  outputIter.Set( inputIter.Get() );
      }
      else
      {
	OutputIterType outputIter( dir2Im, regionOutput );
	outputIter.GoToBegin(); inputIter.GoToBegin();
	for( ; !inputIter.IsAtEnd() && !outputIter.IsAtEnd(); ++inputIter, ++outputIter )
	  outputIter.Set( inputIter.Get() );
      }
    }
    first = !first;
  }
  return;
}

template<typename InputImageType> void MergeAlternateScanRows(
	typename InputImageType::Pointer Output3dImage, int numColsForDiv,
	std::vector< typename InputImageType::Pointer  > &dividedOutputs )
{
  typedef typename itk::ImageRegionIteratorWithIndex< InputImageType  > OutputIterType;
  typedef typename itk::ImageRegionConstIteratorWithIndex< InputImageType  > InputIterConstType;
  itk::SizeValueType numSlicesFull = Output3dImage->GetLargestPossibleRegion().GetSize()[2];
  itk::SizeValueType numCol = Output3dImage->GetLargestPossibleRegion().GetSize()[1];
  itk::SizeValueType numRow = Output3dImage->GetLargestPossibleRegion().GetSize()[0];
  bool first = true;
  itk::IndexValueType firstStackIndex = 0;
  itk::IndexValueType secondStackIndex = 0;
  //Copy slices from alternate rows
  for( unsigned i=0; i<numSlicesFull; i+=numColsForDiv )
  {
    for( unsigned j=i; j<(i+numColsForDiv); ++j )
    {
      typename InputImageType::IndexType inpIndex; inpIndex[0]=0; inpIndex[1]=0;
      if( first ) inpIndex[2] = firstStackIndex++;
      else inpIndex[2] = secondStackIndex++;
      typename InputImageType::IndexType outIndex; outIndex[0]=0; outIndex[1]=0; outIndex[2]=j;
      typename InputImageType::SizeType sliceSize;
      sliceSize[0] = numRow; sliceSize[1] = numCol; sliceSize[2] = 1;
      typename InputImageType::RegionType regionInput, regionOutput;
      regionInput.SetSize( sliceSize ); regionInput.SetIndex( inpIndex );
      regionOutput.SetSize( sliceSize ); regionOutput.SetIndex( outIndex );
      OutputIterType outputIter( Output3dImage, regionOutput );
      if( first )
      {
	InputIterConstType inputIter( dividedOutputs.at(0), regionInput );
	outputIter.GoToBegin(); inputIter.GoToBegin();
	for( ; !inputIter.IsAtEnd() && !outputIter.IsAtEnd(); ++inputIter, ++outputIter )
	  outputIter.Set( inputIter.Get() );
      }
      else
      {
	InputIterConstType inputIter( dividedOutputs.at(1), regionInput );
	outputIter.GoToBegin(); inputIter.GoToBegin();
	for( ; !inputIter.IsAtEnd() && !outputIter.IsAtEnd(); ++inputIter, ++outputIter )
	  outputIter.Set( inputIter.Get() );
      }
    }
    first = !first;
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

  std::string inputImageName = argv[1]; //Name of the input image
  unsigned found = inputImageName.find_last_of(".");
  nameTemplate = inputImageName.substr(0,found) + "_";
  int numColsForDiv = atoi(argv[2]);
  if( numColsForDiv )
    std::cout<<"Number of tiles in each row for meander scan "<<numColsForDiv<<std::endl;
  else
    std::cout<<"Meander scan code not used\n";
  int numThreads = 24;
  int useSingleLev = 0;
  US3ImageType::PixelType lowNoiseThr = 0;
  if( argc > 3 )
    numThreads = atoi( argv[3] );
  if( argc > 4 )
  {
    useSingleLev = atoi( argv[4] );
    std::cout<<"Single level flag set to "<<useSingleLev<<std::endl;
  }
  if( argc > 5 ) lowNoiseThr = (US3ImageType::PixelType) atoll( argv[5] );
  else lowNoiseThr = (US3ImageType::PixelType) LowerNoiseThr;
  if( useSingleLev )
    lowNoiseThr = itk::NumericTraits< US3ImageType::PixelType >::max();
  std::cout<<"Lower noise threshold set to "<<lowNoiseThr<<std::endl;

  double reducedThreadsDbl = std::floor((double)numThreads*0.95);
  int reducedThreads9 = 1>reducedThreadsDbl? 1 : (int)reducedThreadsDbl;
  std::cout<<"Using "<<numThreads<<" and "<<reducedThreads9<<" threads\n";
  US3ImageType::Pointer inputImageFull = ReadITKImageScifio<US3ImageType>( inputImageName );
  itk::SizeValueType numSlicesFull = inputImageFull->GetLargestPossibleRegion().GetSize()[2];
  itk::SizeValueType numCol = inputImageFull->GetLargestPossibleRegion().GetSize()[1];
  itk::SizeValueType numRow = inputImageFull->GetLargestPossibleRegion().GetSize()[0];

  //Dupulicate for noise thresholding
  typedef itk::ImageDuplicator<US3ImageType> DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( inputImageFull );
  try
  {
    duplicator->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught duplicator filter!" << excep << std::endl;
    exit (EXIT_FAILURE);
  }
  US3ImageType::Pointer clonedImageFull = duplicator->GetModifiableOutput();

  //Noise threshold images if needed
  US3ImageType::PixelType upperThreshold = itk::NumericTraits< US3ImageType::PixelType >::max();
  if( lowNoiseThr && !useSingleLev ) 
    upperThreshold = SetSaturatedFGPixelsToMin( inputImageFull, numThreads, lowNoiseThr );

  //Compute global poisson parameters
  US3ImageType::PixelType valsPerBin = 1;
  while( ((double)(upperThreshold+1)/(double)valsPerBin) > NumBins )
  {
    ++valsPerBin;
  }
  std::vector< double > globalHistogram( NumBins, 0 );
  std::vector< double > globalParameters( 5, 0 );
  if( !useSingleLev )
  {
    ComputeGlobalHistogram( inputImageFull, globalHistogram, valsPerBin );
    computePoissonParams( globalHistogram, globalParameters, true );
  }
  std::cout<< globalParameters.at(0) << " Lowest mean\t"
	<< globalParameters.at(1) << " Intermediate mean\t"
	<< globalParameters.at(2) << " Highest mean\t"
	<< globalParameters.at(3) << " Prior for the lowest\t"
	<< globalParameters.at(4) << " Prior for the intermediate\t"
	<< valsPerBin << " vals per bin\n";

  std::cout<<"Number of slices:"<<numSlicesFull<<" Rows[0]:"<<numRow<<" Cols: "<<numCol<<"\n";
  std::vector< US3ImageType::Pointer  > dividedInputs;
  std::vector< US3ImageType::Pointer > dividedOutputs;
  
  if( numColsForDiv )
  {
    dividedInputs.resize(2); dividedOutputs.resize(2);
    DivideImageByScanDirection( inputImageFull, numColsForDiv, dividedInputs );
    DivideImageByScanDirection( clonedImageFull, numColsForDiv, dividedOutputs );
  }
  else
  {
    dividedInputs.push_back( inputImageFull );
    dividedOutputs.push_back( clonedImageFull );
  }

#ifdef DEBUG_THREE_LEVEL_LABELING
  std::vector< UC3ImageType::Pointer  > dividedLabels;
  dividedLabels.resize( dividedInputs.size() );
#endif //DEBUG_THREE_LEVEL_LABELING

/*WriteITKImage<US3ImageType>( dividedInputs.at(0), im11nm );
WriteITKImage<US3ImageType>( dividedOutputs.at(0), im12nm );
if( numColsForDiv )
{
  std::string im21nm = nameTemplate + "image21.nrrd";
  std::string im22nm = nameTemplate + "image22.nrrd";
  WriteITKImage<US3ImageType>( dividedInputs.at(1), im21nm );
  WriteITKImage<US3ImageType>( dividedOutputs.at(1), im22nm );
}*/

for( unsigned numMeanders=0; numMeanders<dividedInputs.size(); ++numMeanders )
{//Start scoping for for loop on scan dir division
  US3ImageType::Pointer inputImage = dividedInputs.at(numMeanders);
  US3ImageType::Pointer clonedImage = dividedOutputs.at(numMeanders);
  itk::SizeValueType numSlices =
	dividedInputs.at(numMeanders)->GetLargestPossibleRegion().GetSize()[2];
  std::vector< US2ImageType::Pointer  > medFiltImages;
  std::vector< CostImageType::Pointer > autoFlourCosts, flourCosts;
  std::vector< CostImageType::Pointer > autoFlourCostsBG, flourCostsBG;
  medFiltImages.resize( numSlices );
  flourCosts.resize( numSlices );   autoFlourCosts.resize( numSlices );
  flourCostsBG.resize( numSlices ); autoFlourCostsBG.resize( numSlices );

#ifdef DEBUG_RESCALING_N_COST_EST
  std::vector< US2ImageType::Pointer > resacledImages;
  resacledImages.resize( numSlices );
#endif //DEBUG_RESCALING_N_COST_EST

//{ //Scoping for Noise thresholded image
  typedef itk::MedianImageFilter< US2ImageType, US2ImageType > MedianFilterType;

#ifdef _OPENMP
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
#if _OPENMP >= 200805L
  #pragma omp parallel for schedule(dynamic,1) num_threads(numThreads) \
	shared(inputImage, medFiltImages, autoFlourCosts, flourCosts, autoFlourCostsBG, flourCostsBG)
#else
  #pragma omp parallel for
#endif
#endif
  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    US2ImageType::Pointer currentSlice = GetTile( inputImage, (itk::SizeValueType)i );
    //Median filter for each slice to remove thermal noise
    //MedianFilterType::Pointer medFilter = MedianFilterType::New();
    //medFilter->SetInput( currentSlice );
    //medFilter->SetRadius( 5 );
    //try
    //{
    //  medFilter ->Update(); 
    //}
    //catch( itk::ExceptionObject & excep )
    //{
    // std::cerr << "Exception caught median filter!" << excep << std::endl;
    //  exit (EXIT_FAILURE);
    //}
    //US2ImageType::Pointer medFiltIm = medFilter->GetOutput();
    US2ImageType::Pointer medFiltIm = currentSlice; //No median filter

    medFiltIm->Register();
    medFiltImages.at(i) = medFiltIm;
//  } //Remove brace after debug
    //Allocate space for costs
    CostImageType::SizeType size;
    size[0] = numRow; size[1] = numCol;
    CostImageType::Pointer costs1 = CreateDefaultCoordsNAllocateSpace<CostImageType>( size );
    CostImageType::Pointer costs2 = CreateDefaultCoordsNAllocateSpace<CostImageType>( size );
    CostImageType::Pointer costs3 = CreateDefaultCoordsNAllocateSpace<CostImageType>( size );
    CostImageType::Pointer costs4 = CreateDefaultCoordsNAllocateSpace<CostImageType>( size );
    costs1->Register(); costs2->Register(); costs3->Register(); costs4->Register();
    autoFlourCosts.at(i)   = costs1; flourCosts.at(i)   = costs2;
    autoFlourCostsBG.at(i) = costs3; flourCostsBG.at(i) = costs4;
#ifdef DEBUG_RESCALING_N_COST_EST
    US2ImageType::Pointer rescIm = CreateDefaultCoordsNAllocateSpace<US2ImageType>( size );
    rescIm->Register();
    resacledImages.at(i) = rescIm;
#endif //DEBUG_RESCALING_N_COST_EST
  }
#ifdef DEBUG_RESCALING_N_COST_EST
  std::string OutFiles = nameTemplate + "costImageMedFilt.nrrd";
  CastNWriteImage2DStackOfVecsTo3D<US2ImageType,US3ImageType>( medFiltImages, OutFiles );
#endif //DEBUG_RESCALING_N_COST_EST

  std::cout<<"Done! Starting to compute costs\n"<<std::flush;


  if( !useSingleLev )
  {
    ComputeCosts( numThreads, medFiltImages, globalParameters, autoFlourCosts, flourCosts,
#ifdef DEBUG_RESCALING_N_COST_EST
			autoFlourCostsBG, flourCostsBG, valsPerBin, resacledImages );
#else
			autoFlourCostsBG, flourCostsBG, valsPerBin );
#endif //DEBUG_RESCALING_N_COST_EST
  }
#ifdef DEBUG_RESCALING_N_COST_EST
  std::string OutFiles1 = nameTemplate + "costImageF.nrrd";
  std::string OutFiles2 = nameTemplate + "costImageFBG.nrrd";
  std::string OutFiles3 = nameTemplate + "costImageAF.nrrd";
  std::string OutFiles4 = nameTemplate + "costImageAFBG.nrrd";
  std::string OutFiles5 = nameTemplate + "costImageInputResc.nrrd";
  CastNWriteImage2DStackOfVecsTo3D<CostImageType,US3ImageType>( flourCosts,	  OutFiles1 );
  CastNWriteImage2DStackOfVecsTo3D<CostImageType,US3ImageType>( flourCostsBG,	  OutFiles2 );
  CastNWriteImage2DStackOfVecsTo3D<CostImageType,US3ImageType>( autoFlourCosts,   OutFiles3 );
  CastNWriteImage2DStackOfVecsTo3D<CostImageType,US3ImageType>( autoFlourCostsBG, OutFiles4 );
  CastNWriteImage2DStackOfVecsTo3D<US2ImageType ,US3ImageType>( resacledImages,   OutFiles5 );
#endif //DEBUG_RESCALING_N_COST_EST

  //Copy into 3d image
  UC3ImageType::SizeType  size;
  size[0] = numRow;
  size[1] = numCol;
  size[2] = numSlices;
  UC3ImageType::Pointer labelImage = CreateDefaultCoordsNAllocateSpace<UC3ImageType>( size );

  std::cout<<"Done! Starting Cuts\n"<<std::flush;
  unsigned count = 0;
  if( !useSingleLev )
{
#ifdef _OPENMP
#if _OPENMP >= 200805L
  #pragma omp parallel for schedule(dynamic,1) num_threads(reducedThreads9) \
	shared( medFiltImages, autoFlourCosts, autoFlourCostsBG, labelImage )
#else
  #pragma omp parallel for num_threads(reducedThreads9)
#endif
#endif
  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    ComputeCut( i, medFiltImages, autoFlourCosts, autoFlourCostsBG, labelImage, 1 );
    ComputeCut( i, medFiltImages, flourCosts, flourCostsBG, labelImage, 2 );
    autoFlourCosts.at(i)->UnRegister();
    autoFlourCostsBG.at(i)->UnRegister();
    flourCosts.at(i)->UnRegister();
    flourCostsBG.at(i)->UnRegister();
    #pragma omp critical
    {
      ++count;
      if( !( (unsigned)std::floor((double)count*100.0/(double)numSlices)%(unsigned)5 ) )
	std::cout<<(unsigned)((double)count*100.0/(double)numSlices)<<"\% Done\r"<<std::flush;
    }
  }
}
  autoFlourCosts.clear();
  autoFlourCostsBG.clear();
  flourCosts.clear();
  flourCostsBG.clear();

  std::cout<<std::endl<<std::flush;
// ************************ //Removed for debugging
/*std::string labelImageName = "/data/kedar/sim/ht/sim_label.tif";
  UC3ImageType::Pointer labelImage = ReadITKImage<UC3ImageType>( labelImageName );
  std::cout<<"Label Image Size: "<<labelImage->GetLargestPossibleRegion().GetSize()[0]
    <<" "<<labelImage->GetLargestPossibleRegion().GetSize()[1]
    <<" "<<labelImage->GetLargestPossibleRegion().GetSize()[2]<<"\n";

  labelImage->Register();
  std::cout<<"Label Image read\n";*/

#ifdef DEBUG_THREE_LEVEL_LABELING
  dividedLabels.at( numMeanders ) = labelImage;
  labelImage->Register();
#endif //DEBUG_THREE_LEVEL_LABELING

//Debugging code to replace function ComputeMeanImages
/*std::vector< CostImageType::Pointer > avgImsVec;
  US2ImageType::Pointer BG_US = ReadITKImage<US2ImageType>( "/data/kedar/sim/ht/sim_BGAvg.tif" );
  US2ImageType::Pointer AF_US = ReadITKImage<US2ImageType>( "/data/kedar/sim/ht/sim_AFAvg.tif" );
  US2ImageType::Pointer FG_US = ReadITKImage<US2ImageType>( "/data/kedar/sim/ht/sim_flAvg.tif" );
  avgImsVec.push_back( CastImage<US2ImageType,CostImageType>( FG_US ) );
  avgImsVec.push_back( CastImage<US2ImageType,CostImageType>( AF_US ) );
  avgImsVec.push_back( CastImage<US2ImageType,CostImageType>( BG_US ) );
*/
  std::cout<<"Computing mean Images\n"<<std::flush;
  std::vector< CostImageType::Pointer > avgImsVec 
	= ComputeMeanImages( labelImage, medFiltImages, numThreads, useSingleLev );

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
  medFiltImages.clear();

//******************
//  std::string labelImageName = nameTemplate + "label.tif";
//  std::string flAvgName = nameTemplate + "flAvg.tif";
//  std::string AFAvgName = nameTemplate + "AFAvg.tif";
//  std::string BGAvgName = nameTemplate + "BGAvg.tif";
//  UC3ImageType::Pointer labelImage = //ReadITKImage<UC3ImageType>( labelImageName );
//  CostImageType::Pointer flAvgIm = //ReadITKImage<CostImageType>( flAvgName );
//  CostImageType::Pointer AFAvgIm = //ReadITKImage<CostImageType>( AFAvgName );
//  CostImageType::Pointer BGAvgIm = //ReadITKImage<CostImageType>( BGAvgName );
//******************
  std::vector<double> flPolyCoeffs(numCoeffs,0), AFPolyCoeffs(numCoeffs,0),
			BGPolyCoeffs(numCoeffs,0), normConstants(2*3*numCoeffs,0);
  std::vector<double> imNorms(3,0);

  //Make pointers for the flour, autoflour n bg avg images
  CostImageType::Pointer flAvgIm, AFAvgIm, BGAvgIm;
  flAvgIm = avgImsVec.at(0); AFAvgIm = avgImsVec.at(1);
  BGAvgIm = avgImsVec.at(2);

#ifdef DEBUG_MEAN_PROJECTIONS
  std::string flAvgName = nameTemplate + "flAvg.tif";
  std::string AFAvgName = nameTemplate + "AFAvg.tif";
  std::string BGAvgName = nameTemplate + "BGAvg.tif";
  if( !useSingleLev)
{
  CastNWriteImage<CostImageType,US2ImageType>(flAvgIm,flAvgName);
  CastNWriteImage<CostImageType,US2ImageType>(AFAvgIm,AFAvgName);
}
  CastNWriteImage<CostImageType,US2ImageType>(BGAvgIm,BGAvgName);
#endif //DEBUG_MEAN_PROJECTIONS

  std::cout<<"Mean Images computed! Estimating polynomials\n"<<std::flush;
  ComputePolynomials( flAvgIm, AFAvgIm, BGAvgIm, flPolyCoeffs, AFPolyCoeffs, BGPolyCoeffs,
			normConstants, imNorms, numThreads, useSingleLev );
  std::cout<<"Polynomials estimated\n"<<std::flush;
  for( itk::SizeValueType i=0; i<numCoeffs; ++i )
    std::cout<<flPolyCoeffs.at(i)<<"\t"<<AFPolyCoeffs.at(i)<<"\t"<<BGPolyCoeffs.at(i)<<"\n"<<std::flush;

  double delta = CorrectImages( flPolyCoeffs, AFPolyCoeffs, BGPolyCoeffs, normConstants,
	clonedImage, labelImage, flAvgIm, AFAvgIm, BGAvgIm, imNorms, numThreads, useSingleLev,
	upperThreshold );
//  std::vector< double > deltaVec;
//  deltaVec.push_back( delta );
  flAvgIm->UnRegister(); AFAvgIm->UnRegister(); BGAvgIm->UnRegister();
  avgImsVec.clear();

  unsigned iterCount = 0;

  std::cout<<"Iteration:"<<++iterCount<<"\tCorrection delta:"<<delta<<std::endl;

  while( delta>IterThresh && iterCount<MaxIter )
  {
    L2THR *= 2;
    L2MINTHR *= 2;
    std::vector< CostImageType::Pointer > avgImsVecIter
	= ComputeMeanImages( labelImage, clonedImage, numThreads, useSingleLev, upperThreshold );
    std::cout<<"Mean Images computed! Estimating polynomials\n"<<std::flush;

    //Make pointers for the flour, autoflour n bg avg images
    CostImageType::Pointer flAvgImIt, AFAvgImIt, BGAvgImIt;
    flAvgImIt = avgImsVecIter.at(0); AFAvgImIt = avgImsVecIter.at(1);
    BGAvgImIt = avgImsVecIter.at(2);

    ComputePolynomials( flAvgImIt, AFAvgImIt, BGAvgImIt, flPolyCoeffs, AFPolyCoeffs, BGPolyCoeffs,
			normConstants, imNorms, numThreads, useSingleLev );
    std::cout<<"Polynomials estimated\n"<<std::flush;
    for( itk::SizeValueType i=0; i<numCoeffs; ++i )
      std::cout<<flPolyCoeffs.at(i)<<"\t"<<AFPolyCoeffs.at(i)<<"\t"<<BGPolyCoeffs.at(i)<<"\n"<<std::flush;

    delta = CorrectImages( flPolyCoeffs, AFPolyCoeffs, BGPolyCoeffs, normConstants,
	clonedImage, labelImage, flAvgImIt, AFAvgImIt, BGAvgImIt, imNorms, numThreads, useSingleLev,
	upperThreshold );

    std::stringstream ss;
    ss << iterCount;
    iterTemplate = ss.str();

#ifdef DEBUG_MEAN_PROJECTIONS
  std::string flAvgName = nameTemplate + iterTemplate + "flAvg.tif";
  std::string AFAvgName = nameTemplate + iterTemplate + "AFAvg.tif";
  std::string BGAvgName = nameTemplate + iterTemplate + "BGAvg.tif";
  if( !useSingleLev)
{
  CastNWriteImage<CostImageType,US2ImageType>(flAvgIm,flAvgName);
  CastNWriteImage<CostImageType,US2ImageType>(AFAvgIm,AFAvgName);
}
  CastNWriteImage<CostImageType,US2ImageType>(BGAvgIm,BGAvgName);
#endif //DEBUG_MEAN_PROJECTIONS

    flAvgImIt->UnRegister(); AFAvgImIt->UnRegister(); BGAvgImIt->UnRegister();
    avgImsVecIter.clear();

    std::cout<<"Iteration:"<<++iterCount<<"\tCorrection delta:"<<delta<<std::endl;

/*  deltaVec.push_back( delta );
    if( delta>IterThresh )
    {
      if( deltaVec.size()>4 )
      {
	double total=0;
	for( unsigned i=deltaVec.size()-5; i<deltaVec.size(); ++i )
	  total+=deltaVec.at(i);
	if( total<FiveIterThr )
	  break;
      }
    }
*/
  }
}//End scoping for for loop on scan dir division
  if( numColsForDiv )
  {
    std::cout<<"Combining images!\n"<<std::flush;
    MergeAlternateScanRows<US3ImageType>( clonedImageFull, numColsForDiv, dividedOutputs );
    dividedOutputs.at(1)->UnRegister();
    dividedInputs.at(1)->UnRegister();
  }

#ifdef DEBUG_THREE_LEVEL_LABELING
  UC3ImageType::Pointer labelImage;
  std::cout<<"Working on debug label images..\n"<<std::flush;
  if( numColsForDiv )
  {
    UC3ImageType::SizeType  sizeFull;
    sizeFull[0] = numRow;
    sizeFull[1] = numCol;
    sizeFull[2] = numSlicesFull;
    labelImage = CreateDefaultCoordsNAllocateSpace<UC3ImageType>( sizeFull );
  std::cout<<"Working on debug label images..\n"<<std::flush;
    MergeAlternateScanRows<UC3ImageType>( labelImage, numColsForDiv, dividedLabels );
    dividedLabels.at(1)->UnRegister();
  }
  else labelImage = dividedLabels.at(0);
  std::cout<<"Writing three level separation image\n"<<std::flush;
  std::string labelImageName = nameTemplate + "label.tif";
  WriteITKImage<UC3ImageType>( labelImage, labelImageName );
  dividedLabels.at(0)->UnRegister(); dividedLabels.clear();
#endif //DEBUG_THREE_LEVEL_LABELING

  std::string correctedImageName = nameTemplate + "IlluminationCorrected.nrrd";

  std::cout<<"Writing corrected image! "<<correctedImageName<<"\n"<<std::flush;

  WriteITKImage<US3ImageType>( clonedImageFull, correctedImageName );

  dividedOutputs.at(0)->UnRegister();
  dividedInputs.at(0)->UnRegister();
  dividedOutputs.clear(); dividedInputs.clear();  

  exit( EXIT_SUCCESS );
}
