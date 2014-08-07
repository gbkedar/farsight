#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCropImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"

//Extracts the N percent of the image along the diagonal
template<typename InputImageType > typename InputImageType::Pointer
  ExtractCenterNPercentOfImage( typename InputImageType::Pointer inputImage,
  double percentage )
{
  typedef typename itk::RegionOfInterestImageFilter
  			< InputImageType, InputImageType > FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  typename InputImageType::IndexType start;
  typename InputImageType::SizeType size, inSize;
  inSize = inputImage->GetLargestPossibleRegion().GetSize();
  double diagonal = std::sqrt( (double)inSize[0]*(double)inSize[0]
			      +(double)inSize[1]*(double)inSize[1] );
  double angle = std::asin( (double)inSize[1]/diagonal );
  diagonal *= (percentage/100.0);
  size[0] = std::floor( std::cos(angle)*diagonal );
  size[1] = std::floor( std::sin(angle)*diagonal );
  start[0] = std::floor(((double)(inSize[0]-size[0]))/2.0);
  start[1] = std::floor(((double)(inSize[1]-size[1]))/2.0);
  typename InputImageType::RegionType desiredRegion;
  desiredRegion.SetSize(size); desiredRegion.SetIndex(start);
  filter->SetRegionOfInterest(desiredRegion);
  filter->SetInput(inputImage);
  try
  {
    filter->Update();
  }
  catch(itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }
  typename InputImageType::Pointer outputImagePointer = filter->GetOutput();
  return outputImagePointer;
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

void usage( const char *funcName )
{
  std::cout << "USAGE:"
	    << " " << funcName << " InputImage MaxDim(Optional Set to 2000)\n";
}

typedef itk::Image< unsigned short, 2 > Ushort2DImageType;
typedef itk::Image< unsigned char, 2 > Uchar2DImageType;
typedef itk::RescaleIntensityImageFilter< Ushort2DImageType, Uchar2DImageType > RescaleIntensityType;
typedef itk::ImageFileWriter< Uchar2DImageType > ImageFileWriterType;

int main(int argc, char *argv[])
{
  if( argc < 2 )
  {
    usage(argv[0]);
    std::cerr << "PRESS ENTER TO EXIT\n";
    getchar();
    return EXIT_FAILURE;
  }

  int largestDim = 2000;
  if( argc > 2 )
    largestDim = atoi( argv[2] );

  std::string inputImageName = argv[1]; //Name of the input image
  unsigned found = inputImageName.find_last_of(".");
  std::string nameTemplate = inputImageName.substr(0,found) + "_";
  Ushort2DImageType::Pointer inputImage = ReadITKImage<Ushort2DImageType>( inputImageName );
  itk::SizeValueType imageDim = inputImage->GetLargestPossibleRegion().GetSize()[0];
  if( imageDim < inputImage->GetLargestPossibleRegion().GetSize()[1] )
    imageDim = inputImage->GetLargestPossibleRegion().GetSize()[1];
  double percentage = ((double)largestDim)/((double)imageDim)*100.0;
  Ushort2DImageType::Pointer croppedImage =
  	ExtractCenterNPercentOfImage<Ushort2DImageType>( inputImage,percentage );
  RescaleIntensityType::Pointer rescale = RescaleIntensityType::New();
  rescale->SetInput( croppedImage );
  rescale->SetOutputMaximum( itk::NumericTraits<Uchar2DImageType::PixelType>::max() );
  rescale->SetOutputMinimum( itk::NumericTraits<Uchar2DImageType::PixelType>::min() );

  std::string opstring = nameTemplate + "cropped.tif";
  ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
  writer->SetFileName( opstring.c_str() );
  writer->SetInput( rescale->GetOutput() );
  std::cout<<"Writing downsampled image:"<<opstring<<std::endl;
  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception in subsampling:\n"
      	<< excep << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
