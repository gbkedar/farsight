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

#include "itkSCIFIOImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkMetaDataObject.h"
#include "itkStreamingImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

typedef itk::Image< unsigned short, 5 > InputImageType;
typedef itk::Image< unsigned short, 4 > InterimImageType;
typedef itk::Image< unsigned short, 3 > OutputImageType;
typedef itk::Image< unsigned short, 2 > TileImageType;

template<typename InputImageType> void WriteITKImage
  ( itk::SmartPointer<InputImageType> inputImagePointer,
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

template <typename InImageType> itk::SmartPointer<InImageType>
AllocateSpaceWithDefaultCoords( typename InImageType::SizeType &size )
{
  typename InImageType::PointType origin;
  typename InImageType::IndexType start;
  const int imDims = InImageType::ImageDimension;
  for( itk::IndexValueType i=0; i<imDims; ++i )
  {
    origin[i] = 0;
    start[i] = 0;
  }

  typename InImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  typename InImageType::Pointer inputImagePointer = InImageType::New();
  inputImagePointer->SetOrigin( origin );
  inputImagePointer->SetRegions( region );
  inputImagePointer->Allocate();
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

TileImageType::Pointer GetTile( OutputImageType::Pointer readImage, itk::IndexValueType i )
{
  typedef itk::ExtractImageFilter< OutputImageType, TileImageType > DataExtractType;
  DataExtractType::Pointer deFilter = DataExtractType::New();
  OutputImageType::RegionType dRegion  = readImage->GetLargestPossibleRegion();
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
  TileImageType::Pointer currentTile = deFilter->GetOutput();
  return currentTile;
}

//Reads a 5D 16-bit czi image and writes out 3D nrrd images with the redundant dimensions collapsed
void RunTempNrrdFileWriter ( std::string inputFileName )
{
  typedef itk::StreamingImageFilter< InputImageType, InputImageType > StreamingFilter;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  typedef itk::ExtractImageFilter< InputImageType, InterimImageType > DataExtractTypeInp2Interim;
  typedef itk::ExtractImageFilter< InterimImageType, OutputImageType > DataExtractTypeInterim2Output;
  typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutputIterType;
  typedef itk::ImageRegionIteratorWithIndex< TileImageType > TileIterType;

  unsigned found = inputFileName.find_last_of(".");
  std::string nameTemplate = inputFileName.substr(0,found) + "_";

  itk::SCIFIOImageIO::Pointer io = itk::SCIFIOImageIO::New();
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetImageIO( io );
  reader->SetFileName( inputFileName.c_str() );

  StreamingFilter::Pointer streamer = StreamingFilter::New();
  streamer->SetInput( reader->GetOutput() );
  streamer->SetNumberOfStreamDivisions( 1 );

  reader->UpdateOutputInformation();
  io->SetSeries(0);
  reader->Modified();

  itk::SizeValueType seriesEnd = io->GetSeriesCount();

  //Allocate space for all channels and read them
  InputImageType::SizeType size;
  size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  itk::IndexValueType numChannels = size[4];

  std::cout<<"Size of each series in the czi file: "
	<<size[0]<<"\t"<<size[1]<<"\t"<<size[2]<<"\t"<<size[3]<<"\t"<<size[4]<<"\n"
	<<"Number of series is: "<<seriesEnd<<std::endl
	<<"Number of channels is:"<<numChannels<<std::endl<<std::flush;

  OutputImageType::SizeType sizeOut;
  sizeOut[0] = size[0];
  sizeOut[1] = size[1];
  sizeOut[2] = seriesEnd;

  std::vector< itk::SmartPointer<OutputImageType> > outputPointers;
  for( itk::IndexValueType i=0; i<numChannels; ++i )
  {
    OutputImageType::Pointer tempPointer = AllocateSpaceWithDefaultCoords<OutputImageType>(sizeOut);
    outputPointers.push_back( tempPointer );
  }
  std::cout<<"Processing series: ";
  for( itk::SizeValueType i=0; i<seriesEnd; ++i )
  {
    //Start by collapsing the un-necessary dimensions
    if( i )
    {
      io->SetSeries(i);
      reader->Modified();
    }
    DataExtractTypeInp2Interim::Pointer toInterim = DataExtractTypeInp2Interim::New();
    InputImageType::RegionType dRegion1 = reader->GetOutput()->GetLargestPossibleRegion();
    dRegion1.SetSize (3,0);
    toInterim->SetExtractionRegion(dRegion1);
    toInterim->SetDirectionCollapseToIdentity();
    toInterim->SetInput( reader->GetOutput() );
    toInterim->UpdateOutputInformation();

    DataExtractTypeInterim2Output::Pointer toOutput = DataExtractTypeInterim2Output::New();
    InterimImageType::RegionType dRegion2 = toInterim->GetOutput()->GetLargestPossibleRegion();
    dRegion2.SetSize (2,0);
    toOutput->SetExtractionRegion(dRegion2);
    toOutput->SetDirectionCollapseToIdentity();
    toOutput->SetInput( toInterim->GetOutput() );

    try
    {
      toOutput->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception caught in slice extraction filter!" << excep << std::endl;
      exit (EXIT_FAILURE);
    }
    OutputImageType::Pointer currentSeries = toOutput->GetOutput();

    OutputImageType::SizeType size1;
    size1 = toOutput->GetOutput()->GetLargestPossibleRegion().GetSize();
    std::cout<<i<<"\t";

    //Write slices to images
#ifdef _OPENMP
  #pragma omp parallel for
#endif
    for( itk::IndexValueType j=0; j<numChannels; ++j )
    {
      TileImageType::Pointer currentTile;
#ifdef _OPENMP
      #pragma omp critical
#endif
      currentTile = GetTile( currentSeries, j );
      OutputImageType::RegionType currentSliceRegion = outputPointers.at(j)->GetLargestPossibleRegion();
      currentSliceRegion.SetSize(2,1);
      currentSliceRegion.SetIndex(2,i);
      OutputIterType outputIter( outputPointers.at(j), currentSliceRegion ); outputIter.GoToBegin();
      TileIterType tileIter( currentTile, currentTile->GetLargestPossibleRegion() ); tileIter.GoToBegin();
      for( ; !tileIter.IsAtEnd(); ++tileIter, ++outputIter )
	outputIter.Set( tileIter.Get() );
    }
  }
  std::cout<<std::endl;

  for( itk::IndexValueType i=0; i<numChannels; ++i )
  {
    std::stringstream ss; ss<<i<<".nrrd";
    std::string outputFilename = nameTemplate + "C" + ss.str();
    WriteITKImage<OutputImageType>( outputPointers.at(i), outputFilename );
  }

  return;
}

int main(int argc, char *argv[])
{
  if( argc < 2 )
  {
    std::cerr<<"Usage:"<<argv[0]<<" InputCziFile\n";
    std::cerr << "PRESS ENTER TO EXIT\n";
    getchar();
    return EXIT_FAILURE;
  }

  std::string inputFilename = argv[1];
  
  RunTempNrrdFileWriter( inputFilename );

  return EXIT_SUCCESS;
}
