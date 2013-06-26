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
#define DBGGG

#include <vector>
#include <string>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <float.h>
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
#include "itkRescaleIntensityImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOtsuMultipleThresholdsCalculator.h"

#define WinSz 100	//Histogram computed on this window
#define CWin  26	//This is the inner window
#define NumBins 1024	//Downsampled to these number of bins

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
	    << " " << funcName << " InputImage OutputImage NumberOfThreads(Optional-default=24)\n";
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
//    if( i==(NumBins-1) )
//      for( itk::SizeValueType j=0; j<overFlow; ++j )
//	histogram.at(i) += histogramInternal.at(i+1*valsPerBin+j);
    histogram.at(i) /= ((double)WinSz)*((double)WinSz)*((double)medFiltImages.size());
#ifdef DBGG
    sum += histogram.at(i);
  }
  std::cout<<"Sum:"<<sum<<std::endl;
#else
  }
#endif //DBGG
  if( valsPerBin )
    max = valsPerBin*NumBins - 1;
  return max;
}

void computePoissonParams( std::vector< double > &histogram,
			   std::vector< double > &parameters )
{
  itk::SizeValueType max = histogram.size()-1;
#ifdef DBGG
  std::cout<<"computing pos params with max " << max << "\n" << std::flush;
#endif
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
        parameters.at(0) = U0; //Lowest mean
        parameters.at(1) = U1; //Intermediate mean
        parameters.at(2) = U2; //Highest mean
        parameters.at(3) = P0; //Prior for the lowest
        parameters.at(4) = P1; //Prior for the intermediate, highest will be 1-(P0+P1)
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

#ifdef DBGG
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
#endif
  }
#ifdef DBGG
  std::cout<<"Parameters2: ";
  for( itk::SizeValueType j=0; j<parameters.size(); ++j )
    std::cout<<parameters.at(j)<<"\t";
  std::cout<<"\n"<<std::flush;
#endif //DBGG
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

void ComputeCosts( int numThreads,
		   std::vector< itk::SmartPointer<US2ImageType>  > &medFiltImages,
		   std::vector< itk::SmartPointer<CostImageType> > &autoFlourCosts,
		   std::vector< itk::SmartPointer<CostImageType> > &flourCosts,
		   std::vector< itk::SmartPointer<CostImageType> > &autoFlourCostsBG,
		   std::vector< itk::SmartPointer<CostImageType> > &flourCostsBG
#ifdef DBGGG
, std::vector< itk::SmartPointer<US2ImageType> > &resacledImages
#endif //DBGGG
		   )
{
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  typedef itk::ImageRegionIterator< CostImageType > CostIterType;
  itk::IndexValueType numCol = medFiltImages.at(0)->GetLargestPossibleRegion().GetSize()[1];
  itk::IndexValueType numRow = medFiltImages.at(0)->GetLargestPossibleRegion().GetSize()[0];
  itk::IndexValueType WinSz2 = (itk::IndexValueType)floor(((double)WinSz)/2+0.5)
			      -(itk::IndexValueType)floor(((double)CWin)/2+0.5);

#ifdef DBGGG
  unsigned count = 0;
  clock_t start_time = clock();
  double ratioMax=0, ratioMin=DBL_MAX;
#endif //DBGGG

#ifdef _OPENMP
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
  #pragma omp parallel for num_threads(numThreads)
#endif
  for( itk::IndexValueType i=0; i<numRow; i+=CWin )
  {
    for( itk::IndexValueType j=0; j<numCol; j+=CWin )
    {
      //Compute histogram at point i,j with window size define WinSz
      std::vector< double > histogram( NumBins, 0 );
      US2ImageType::IndexType start; start[0] = i-WinSz2; start[1] = j-WinSz2;
      if( start[0]<0 ) start[0] = 0; if( start[1]<0 ) start[1] = 0;
      if( (start[0]+WinSz)>=numRow ) start[0] =  numRow-WinSz-1;
      if( (start[1]+WinSz)>=numCol ) start[1] =  numCol-WinSz-1;
      itk::SizeValueType max = ComputeHistogram( medFiltImages, histogram, start );
      US2ImageType::IndexType curPoint; curPoint[0] = i; curPoint[1] = j;
      std::vector< double > parameters( 5, 0 );
      computePoissonParams( histogram, parameters );
      std::vector< double > pdf0( histogram.size()+1, 1 ); ComputePoissonProbability( parameters.at(0), pdf0 );
      std::vector< double > pdf1( histogram.size()+1, 1 ); ComputePoissonProbability( parameters.at(1), pdf1 );
      std::vector< double > pdf2( histogram.size()+1, 1 ); ComputePoissonProbability( parameters.at(2), pdf2 );
//      //Fix the params
//      double ratio = ((double)max)/((double)histogram.size());
//      parameters.at(0) *= ratio; parameters.at(2) *= ratio; parameters.at(2) *= ratio;
      double ratio = ((double)histogram.size())/((double)max);
#ifdef DBGGG
      if( ratio > ratioMax )
        ratioMax = ratio;
      if( ratio < ratioMin )
	ratioMin = ratio;
      if( !omp_get_thread_num() )
      {
	std::cout << "histogram computed for " << curPoint << "\t" ;
	std::cout << "computing costs from hist\t";
	std::cout << "Time: " << (clock()-start_time)/((float)CLOCKS_PER_SEC) << "\n";
	start_time = clock();
      }
#endif //DBGGG
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
#ifdef DBGGG
	typedef itk::ImageRegionIterator< US2ImageType > IterType;
	IterType rescaleIter ( resacledImages.at(k), region );
	rescaleIter.GoToBegin();
#endif //DBGGG

	for( ; !constIter.IsAtEnd(); ++constIter, ++costIterFlour, ++costIterFlourBG,
#ifdef DBGGG
		++costIterAutoFlour, ++costIterAutoFlourBG, ++rescaleIter )
#else  //DBGGG
		++costIterAutoFlour, ++costIterAutoFlourBG )
#endif //DBGGG
	{
	  itk::SizeValueType currentPixel = std::floor( ((double)constIter.Get())*ratio + 0.5 );
	  if( currentPixel >= pdf0.size() )
	    currentPixel = pdf0.size()-1;
#ifdef DBGGG
	  rescaleIter.Set( currentPixel );
#endif //DBGGG

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
	  else if( currentPixel > max )
	  {
	    FBG = AFBG = 0;
	    F   =  AF  = 10000.0;
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
#ifdef DBGG
    #pragma omp critical
    {
      std::cout<<++count<<"\t";
    }
#endif //DBGG
  }
#ifdef DBGG
  std::cout<<"\n";
#endif //DBGG
#ifdef DBGGG
  std::cout<<"Ratio max:"<<ratioMax<<"\t\tmin:"<<ratioMin<<std::endl<<std::flush;
#endif //DBGGG
  return;
}

#ifdef DBGGG
void CastNWriteImage( std::vector< itk::SmartPointer<CostImageType> > &inputImage, std::string &outFileName )
{
  typedef itk::ImageRegionConstIterator< CostImageType > ConstIterType2d;
  typedef itk::ImageRegionIteratorWithIndex< US3ImageType > IterType3d;
  typedef itk::ImageFileWriter< US3ImageType > WriterType;
  //Allocate space
  US3ImageType::Pointer outputImage = US3ImageType::New();
  US3ImageType::PointType origin;
  origin[0] = 0; origin[1] = 0; origin[2] = 0;
  outputImage->SetOrigin( origin );
  US3ImageType::IndexType start;
  start[0] = 0; start[1] = 0; start[2] = 0;
  US3ImageType::SizeType size;
  size[0] = inputImage.at(0)->GetLargestPossibleRegion().GetSize()[0];
  size[1] = inputImage.at(0)->GetLargestPossibleRegion().GetSize()[1];
  size[2] = inputImage.size();
  US3ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  outputImage->SetRegions( region );
  outputImage->Allocate();
  outputImage->FillBuffer(0);
  outputImage->Update();

  //Cast n Write Values
  bool warningWritten = false;
  itk::SizeValueType typeMax = itk::NumericTraits<US3ImageType::PixelType>::max();
  for( itk::SizeValueType i=0; i<inputImage.size(); ++i )
  {
    //Start index
    start[2] = i; size[2] = 1; //Reset for writing out the slices
    CostImageType::IndexType start2d; start2d[0] = 0;      start2d[1] = 0;
    CostImageType::SizeType  size2d;   size2d[0] = size[0]; size2d[1] = size[1];
    region.SetSize( size ); region.SetIndex( start );
    CostImageType::RegionType region2d; region2d.SetSize( size2d ); region2d.SetIndex( start2d );
    ConstIterType2d iter2d ( inputImage.at(i), region2d );
    IterType3d iter3d( outputImage, region );
    for( ; !iter2d.IsAtEnd(); ++iter2d, ++iter3d )
      iter3d.Set( (itk::SizeValueType)std::floor( iter2d.Get()+0.5 ) );
  }
  WriterType::Pointer writer = WriterType::New();
  writer = WriterType::New();
  writer->SetFileName( outFileName.c_str() );
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
  return;
}
void CastNWriteScaling( std::vector< itk::SmartPointer<US2ImageType> > &inputImage, std::string &outFileName )
{
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType2d;
  typedef itk::ImageRegionIteratorWithIndex< US3ImageType > IterType3d;
  typedef itk::ImageFileWriter< US3ImageType > WriterType;
  //Allocate space
  US3ImageType::Pointer outputImage = US3ImageType::New();
  US3ImageType::PointType origin;
  origin[0] = 0; origin[1] = 0; origin[2] = 0;
  outputImage->SetOrigin( origin );
  US3ImageType::IndexType start;
  start[0] = 0; start[1] = 0; start[2] = 0;
  US3ImageType::SizeType size;
  size[0] = inputImage.at(0)->GetLargestPossibleRegion().GetSize()[0];
  size[1] = inputImage.at(0)->GetLargestPossibleRegion().GetSize()[1];
  size[2] = inputImage.size();
  US3ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  outputImage->SetRegions( region );
  outputImage->Allocate();
  outputImage->FillBuffer(0);
  outputImage->Update();

  //Cast n Write Values
  bool warningWritten = false;
  itk::SizeValueType typeMax = itk::NumericTraits<US2ImageType::PixelType>::max();
  for( itk::SizeValueType i=0; i<inputImage.size(); ++i )
  {
    //Start index
    start[2] = i; size[2] = 1; //Reset for writing out the slices
    US2ImageType::IndexType start2d; start2d[0] = 0;      start2d[1] = 0;
    US2ImageType::SizeType  size2d;   size2d[0] = size[0]; size2d[1] = size[1];
    region.SetSize( size ); region.SetIndex( start );
    US2ImageType::RegionType region2d; region2d.SetSize( size2d ); region2d.SetIndex( start2d );
    ConstIterType2d iter2d ( inputImage.at(i), region2d );
    IterType3d iter3d( outputImage, region );
    for( ; !iter2d.IsAtEnd(); ++iter2d, ++iter3d )
      iter3d.Set( (itk::SizeValueType)std::floor( iter2d.Get()+0.5 ) );
  }
  WriterType::Pointer writer = WriterType::New();
  writer = WriterType::New();
  writer->SetFileName( outFileName.c_str() );
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
  return;
}
#endif //DBGGG

void ComputeCut( itk::IndexValueType slice,
		 std::vector< itk::SmartPointer<US2ImageType>  > &medFiltImages,
		 std::vector< itk::SmartPointer<CostImageType> > &flourCosts,
		 std::vector< itk::SmartPointer<CostImageType> > &flourCostsBG,
		 US3ImageType::Pointer outputImage,
		 US3ImageType::PixelType foregroundValue
		)
{
  double sigma = 25.0; //What! A hard coded constant check Boykov's paper!! Also check 20 in weights
  typedef itk::ImageRegionIteratorWithIndex< CostImageType > CostIterType;
  typedef itk::ImageRegionIteratorWithIndex< US2ImageType > US2IterType;
  typedef itk::ImageRegionIteratorWithIndex< US3ImageType > US3IterType;
  //Compute the number of nodes and edges
  itk::SizeValueType numRow   = medFiltImages.at(0)->GetLargestPossibleRegion().GetSize()[0];
  itk::SizeValueType numCol   = medFiltImages.at(0)->GetLargestPossibleRegion().GetSize()[1];
  itk::SizeValueType numNodes = numCol*numRow;
  itk::SizeValueType numEdges = 3*numCol*numRow /*Down, right and diagonal*/ + 1
  				- 2*numCol/*No Down At Bottom*/ - 2*numRow/*No Right At Edge*/;

  typedef Graph_B < unsigned, unsigned, unsigned > GraphType;
  GraphType *graph = new GraphType( numNodes, numEdges );
  US2IterType medianIter( medFiltImages.at(slice),
  			  medFiltImages.at(slice)->GetLargestPossibleRegion() );
  CostIterType AFCostIter( flourCosts.at(slice), 
  			   flourCosts.at(slice)->GetLargestPossibleRegion() );
  CostIterType AFBGCostIter( flourCostsBG.at(slice), 
  			     flourCostsBG.at(slice)->GetLargestPossibleRegion() );
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
  US3ImageType::IndexType start;start[0] = 0;	  start[1] = 0;	    start[2] = slice;
  US3ImageType::SizeType size;   size[0] = numRow; size[1] = numCol; size[2] = 1;
  US3ImageType::RegionType region; region.SetSize( size ); region.SetIndex( start );
  US3IterType outputIter( outputImage, region );
  for( itk::SizeValueType i=0; i<numRow-1; ++i )
  {
    for( itk::SizeValueType j=0; j<numCol-1; ++j )
    {
      US3ImageType::IndexType index; index[0] = i; index[1] = j; index[2] = slice;
      outputIter.SetIndex( index );
      itk::SizeValueType indexCurrentNode = i*numCol+j;
      if( graph->what_segment( indexCurrentNode ) != GraphType::SOURCE )
        outputIter.Set( foregroundValue );
    }
  }
  delete graph;
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

  std::string inputImageName  = argv[1]; //Name of the input image
  std::string outputImageName = argv[2];
  int numThreads = 24;
  if( argc == 4 )
    numThreads = atoi( argv[3] );

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

#ifdef DBGGG
  std::vector< itk::SmartPointer<US2ImageType> > resacledImages;
  resacledImages.resize( numSlices );
#endif //DBGGG

#ifdef _OPENMP
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
#if _OPENMP >= 200805L
  #pragma omp parallel for schedule(dynamic,1) num_threads(numThreads)
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
    costs1->Allocate();			costs2->Allocate();
    costs3->Allocate();			costs4->Allocate();
    costs1->FillBuffer(0);		costs2->FillBuffer(0);
    costs3->FillBuffer(0);		costs4->FillBuffer(0);
    try
    {
      costs1->Update();	costs2->Update();
      costs3->Update(); costs4->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception caught !" << excep << std::endl;
      exit (EXIT_FAILURE);
    }
    costs1->Register(); costs2->Register();
    costs3->Register(); costs4->Register();
    autoFlourCosts.at(i)   = costs1; flourCosts.at(i)   = costs2;
    autoFlourCostsBG.at(i) = costs3; flourCostsBG.at(i) = costs4;
#ifdef DBGGG
    US2ImageType::Pointer rescIm = US2ImageType::New();
    rescIm->SetOrigin( origin ); rescIm->SetRegions( region );
    rescIm->Allocate();		 rescIm->FillBuffer(0);
    try
    {
      rescIm->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception caught !" << excep << std::endl;
      exit (EXIT_FAILURE);
    }
    rescIm->Register();
    resacledImages.at(i) = rescIm;
#endif //DBGGG
  }
#ifdef DBGGG
#endif //DBGGG

  std::cout<<"Done! Starting to compute costs\n"<<std::flush;

  ComputeCosts( numThreads, medFiltImages, autoFlourCosts, flourCosts,
#ifdef DBGGG
  				autoFlourCostsBG, flourCostsBG, resacledImages );
#else
  				autoFlourCostsBG, flourCostsBG );
#endif //DBGGG

#ifdef DBGGG
  std::string OutFiles1 = "costImageF.nrrd";
  std::string OutFiles2 = "costImageFBG.nrrd";
  std::string OutFiles3 = "costImageAF.nrrd";
  std::string OutFiles4 = "costImageAFBG.nrrd";
  std::string OutFiles5 = "costImageInputResc.nrrd";
  CastNWriteImage( flourCosts,		OutFiles1 );
  CastNWriteImage( flourCostsBG,	OutFiles2 );
  CastNWriteImage( autoFlourCosts,	OutFiles3 );
  CastNWriteImage( autoFlourCostsBG,	OutFiles4 );
  CastNWriteScaling( resacledImages,	OutFiles5 );
#endif //DBGGG

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

  std::cout<<"Done! Computing Cuts\n"<<std::flush;

#ifdef _OPENMP
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
  #pragma omp parallel for num_threads(numThreads)
#endif
  for( itk::IndexValueType i=0; i<numSlices; ++i )
  {
    ComputeCut( i, medFiltImages, autoFlourCosts, autoFlourCostsBG, outputImage, 1 );
    ComputeCut( i, medFiltImages, flourCosts, flourCostsBG, outputImage, 2 );
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
