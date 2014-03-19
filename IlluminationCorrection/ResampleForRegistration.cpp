#include <iostream>

#include "itkImage.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

typedef itk::Image< unsigned short, 2 > Ushort2DImageType;
typedef itk::Image< unsigned char, 2 > Uchar2DImageType;
typedef itk::Image< float, 2 > Float3DImageType;
typedef itk::RescaleIntensityImageFilter< Ushort2DImageType, Uchar2DImageType > RescaleIntensityType;
typedef itk::CastImageFilter< Ushort2DImageType, Float3DImageType > CastFilterType;
typedef itk::RecursiveGaussianImageFilter< Float3DImageType, Float3DImageType > GaussianFilterType;
typedef itk::ResampleImageFilter< Float3DImageType, Ushort2DImageType > ResampleFilterType;
typedef itk::IdentityTransform< double, 2 >  TransformType;
typedef itk::LinearInterpolateImageFunction< Float3DImageType, double > InterpolatorType;
typedef itk::ImageFileWriter< Uchar2DImageType > ImageFileWriterType;

void usage( const char *funcName )
{
  std::cout << "USAGE:"
	    << " " << funcName << " InputImage MaxDim(Optional Set to 10000)\n";
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

  int largestDim = 2000;
  if( argc > 2 )
    largestDim = atoi( argv[2] );

  std::string inputImageName = argv[1]; //Name of the input image
  unsigned found = inputImageName.find_last_of(".");
  std::string nameTemplate = inputImageName.substr(0,found) + "_";
  Ushort2DImageType::Pointer inputImage = ReadITKImage<Ushort2DImageType>( inputImageName );

  double scaleFactor = inputImage->GetLargestPossibleRegion().GetSize()[1] >
			inputImage->GetLargestPossibleRegion().GetSize()[0] ?
			inputImage->GetLargestPossibleRegion().GetSize()[1] :
			inputImage->GetLargestPossibleRegion().GetSize()[0];
  scaleFactor = scaleFactor/(double)largestDim;
  const Ushort2DImageType::SpacingType& inputSpacing = inputImage->GetSpacing();

  CastFilterType::Pointer caster = CastFilterType::New();
  caster->SetInput( inputImage );

  GaussianFilterType::Pointer smoother = GaussianFilterType::New();
  smoother->SetInput( caster->GetOutput() );
  smoother->SetSigma( inputSpacing[0]*(scaleFactor/2) );
  smoother->SetNormalizeAcrossScale( true );

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  resampler->SetTransform( transform );

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  resampler->SetInterpolator( interpolator );

  resampler->SetDefaultPixelValue( 0 );

  Ushort2DImageType::SpacingType spacing;
  spacing[0] = inputSpacing[0] * scaleFactor;
  spacing[1] = inputSpacing[1] * scaleFactor;

  resampler->SetOutputSpacing( spacing );
  resampler->SetOutputOrigin( inputImage->GetOrigin() );
  resampler->SetOutputDirection( inputImage->GetDirection() );

  Ushort2DImageType::SizeType size;
  size[0] = ((double)inputImage->GetLargestPossibleRegion().GetSize()[0])/scaleFactor;
  size[1] = ((double)inputImage->GetLargestPossibleRegion().GetSize()[1])/scaleFactor;
  resampler->SetSize( size );

  resampler->SetInput( smoother->GetOutput() );

  RescaleIntensityType::Pointer rescale = RescaleIntensityType::New();
  rescale->SetInput( resampler->GetOutput() );
  rescale->SetOutputMaximum( itk::NumericTraits<Uchar2DImageType::PixelType>::max() );
  rescale->SetOutputMinimum( itk::NumericTraits<Uchar2DImageType::PixelType>::min() );

  std::string opstring = nameTemplate + "subsample.tif";
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
