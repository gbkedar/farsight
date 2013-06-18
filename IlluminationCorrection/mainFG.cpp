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

#include "adaptive_binarization.h"

typedef unsigned short USPixelType;
typedef unsigned char  UCPixelType;
const unsigned int     Dimension3 = 3;
const unsigned int     Dimension2 = 2;
typedef itk::Image< USPixelType, Dimension3 > US3ImageType;
typedef itk::Image< UCPixelType, Dimension3 > UC3ImageType;
typedef itk::Image< USPixelType, Dimension2 > US2ImageType;

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

void AdaptivelyBinarizeTile( US2ImageType::Pointer InputImage,
			     US2ImageType::Pointer BinaryImage )
{
  typedef itk::RegionOfInterestImageFilter< US2ImageType, US2ImageType > ROIFilterType;
  typedef itk::ConnectedComponentImageFilter< US2ImageType, US2ImageType >
  								ConnectedComponentFilterType;
  typedef itk::ImageRegionConstIterator< US2ImageType > BinOutConstIteratorType;
  typedef itk::ImageRegionIteratorWithIndex< US2ImageType > BinMontageIteratorType;

  //Redefine some constants
  itk::SizeValueType TileSizeX = InputImage->GetLargestPossibleRegion().GetSize()[0];
  itk::SizeValueType TileSizeY = InputImage->GetLargestPossibleRegion().GetSize()[1];

  //Crop Input Image
  unsigned char *DataPtr;
  DataPtr = (unsigned char*)malloc( sizeof(unsigned char)*TileSizeX*TileSizeY );
{ //Scoping for temp cropimage
    US2ImageType::RegionType CroppedRegion;
    CroppedRegion.SetSize ( Size );
    CroppedRegion.SetIndex( Start );
    ROIFilterType::Pointer CropImageFilter = ROIFilterType::New();
    CropImageFilter->SetInput( InputImage );
    CropImageFilter->SetRegionOfInterest( CroppedRegion );
    try{ CropImageFilter->Update(); }
    catch( itk::ExceptionObject & excp )
    {
      std::cerr <<  "Extraction for binarization failed" << excp << std::endl;
    }
#if 0
  typedef itk::ImageFileWriter< US2ImageType > BinaryWriterType1;
  BinaryWriterType1::Pointer binwriter1 = BinaryWriterType1::New();
  binwriter1->SetInput( CropImageFilter->GetOutput() );
  std::stringstream filess1;
  filess1 << Start[0] << "_" <<Start[1] << "_" << Start[2];
  std::string OutFile1 = TempFolder + "/Temp_" + filess1.str() + "_crop.tif" ;
  binwriter1->SetFileName( OutFile1.c_str() );
  try{ binwriter1->Update(); }
  catch( itk::ExceptionObject & excp ){ std::cerr << excp << std::endl; }
#endif
    BinMontageIteratorType PixBuf(  CropImageFilter->GetOutput(),
    				CropImageFilter->GetOutput()->GetLargestPossibleRegion() );
    itk::SizeValueType Index = 0;
    for( PixBuf.GoToBegin(); !PixBuf.IsAtEnd(); ++PixBuf )
    {
      DataPtr[Index] = PixBuf.Get();
      ++Index;
    }
} //End scoping for temp cropimage

  //Create binarization object
  ftk::NuclearSegmentation * newNucSeg = new ftk::NuclearSegmentation();

  //Read Parameters
  for(int i=0; i<(int)definition->nuclearParameters.size(); ++i)
    newNucSeg->SetParameter(	definition->nuclearParameters.at(i).name,
				int(definition->nuclearParameters.at(i).value) );

  //Convert input image to FTKImage
  ftk::Image::Pointer FTKIImage = ftk::Image::New();
  std::vector<unsigned char> color;
  color.assign(3,255);
  FTKIImage->AppendChannelFromData3D( (void*)DataPtr,
					itk::ImageIOBase::UCHAR, sizeof(unsigned char),
					Size[0], Size[1], Size[2], "nuc", color, false );
  newNucSeg->SetInput( FTKIImage, "nuc", 0 );

  //Run Binarization
  newNucSeg->Binarize(true);

  //Get Output
  ftk::Image::Pointer FTKOImage = newNucSeg->GetLabelImage();
  US2ImageType::Pointer BinIm = FTKOImage->
					GetItkPtr< US2ImageType::PixelType >( 0, 0, ftk::Image::DEFAULT );

  //Check if the tile has more than 3 connected components
  ConnectedComponentFilterType::Pointer connectedComponentFilter = ConnectedComponentFilterType::New();
  connectedComponentFilter->SetInput( BinIm );
  connectedComponentFilter->FullyConnectedOn();
  try
  {
    connectedComponentFilter->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr <<  "CC filter for initial binarization failed" << excp << std::endl;
  }

  //Tiles with only background usually have 1-3 CCs
  if( connectedComponentFilter->GetObjectCount() > 3  )
  {
    BinOutConstIteratorType BinOutConstIter( BinIm, BinIm->GetLargestPossibleRegion() );
    BinMontageIteratorType BinMontagIter( BinaryImage, BinaryImage->GetLargestPossibleRegion() );

    //Copy output of tile into montage
    for( itk::SizeValueType k=Start[2]; k<numStacks; ++k )
      for( itk::SizeValueType l=Start[1]; l<(Start[1]+TileSize); ++l )
        for( itk::SizeValueType m=Start[0]; m<(Start[0]+TileSize); ++m )
        {
	  US2ImageType::IndexType FtkBinIndex;
	  FtkBinIndex[0] = m-Start[0];
	  FtkBinIndex[1] = l-Start[1];
	  FtkBinIndex[2] = k;
	  US2ImageType::IndexType BinaryImIndex;
	  BinaryImIndex[0] = m;
	  BinaryImIndex[1] = l;
	  BinaryImIndex[2] = k;
	  BinOutConstIter.SetIndex( FtkBinIndex );
	  BinMontagIter.SetIndex( BinaryImIndex );
	  if( BinOutConstIter.Get() )
	    BinMontagIter.Set( 255 );
        }
  }
#if 0
  typedef itk::ImageFileWriter< US2ImageType > BinaryWriterType;
  BinaryWriterType::Pointer binwriter = BinaryWriterType::New();
  binwriter->SetInput( BinIm );
  std::stringstream filess;
  filess << Start[0] << "_" <<Start[1] << "_" << Start[2];
  std::string OutFile = TempFolder + "/Temp_" + filess.str() + "_Bin.tif" ;
  binwriter->SetFileName( OutFile.c_str() );
  try{ binwriter->Update(); }
  catch( itk::ExceptionObject & excp ){ std::cerr << excp << std::endl; }
#endif
  delete newNucSeg;
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

  std::vector< itk::SmartPointer<US2ImageType> > threshImages;
  threshImages.resize( numSlices );
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
    medFilter->SetRadius( 5 );

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
    US2ImageType::Pointer thresholdIm;

    AdaptivelyBinarizeTile( medFiltIm, thresholdIm );
    thresholdIm->Register();
    
    threshImages.at(i) = thresholdIm;
  }

  std::cout<<"Done. Copying image."
	   <<std::endl<<std::flush;

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
