#include <iostream>

#include "itkImage.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTIFFImageIO.h"

typedef itk::Image< unsigned short, 2 > Ushort2DImageType;
typedef itk::Image< unsigned char, 2 > Uchar2DImageType;

typedef itk::RescaleIntensityImageFilter< Ushort2DImageType, Uchar2DImageType > RescaleIntensityType;
typedef itk::ImageFileWriter< Uchar2DImageType > ImageFileWriterType;

void usage( const char *funcName )
{
  std::cout << "USAGE:"
	    << " " << funcName << " InputImage MaxDim(Optional Set to 2000)\n";
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

int main(int argc, char *argv[])
{
  if( argc < 2 )
  {
    usage(argv[0]);
    std::cerr << "PRESS ENTER TO EXIT\n";
    getchar();
    return EXIT_FAILURE;
  }

  typedef  itk::TIFFImageIO TIFFIOType;
  TIFFIOType::Pointer tiffIO1 = TIFFIOType::New();

  std::string inputImageName = argv[1]; //Name of the input image
  unsigned found = inputImageName.find_last_of(".");
  std::string nameTemplate = inputImageName.substr(0,found) + "_";
  Ushort2DImageType::Pointer inputImage = ReadITKImage<Ushort2DImageType>( inputImageName );

  RescaleIntensityType::Pointer rescale = RescaleIntensityType::New();
  rescale->SetInput( inputImage );
  rescale->SetOutputMaximum( itk::NumericTraits<Uchar2DImageType::PixelType>::max() );
  rescale->SetOutputMinimum( itk::NumericTraits<Uchar2DImageType::PixelType>::min() );

  std::string opstring = nameTemplate + "8bit.tif";
  ImageFileWriterType::Pointer writer = ImageFileWriterType::New();
  writer->SetFileName( opstring.c_str() );
  writer->SetInput( rescale->GetOutput() );
  writer->SetImageIO(tiffIO1);
  std::cout<<"Writing 8-bit image:"<<opstring<<std::endl;
  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception in rescaling:\n"
      	<< excep << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
