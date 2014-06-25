#include "itkTranslationTransform.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
// #include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredTransformInitializer.h"
// #include "itkAffineTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCropImageFilter.h"
#include "itkCommand.h"
#include "itkTIFFImageIO.h"

const    unsigned int    Dimension = 2;
typedef  float           PixelType;

typedef itk::Image< PixelType, Dimension >  FloatImageType;
typedef itk::CropImageFilter< FloatImageType, FloatImageType > CropImageFilterType;

typedef  itk::TIFFImageIO TIFFIOType;

typedef itk::ImageRegistrationMethod< FloatImageType, FloatImageType > RegistrationType;
typedef RegistrationType::ParametersType ParametersType;
typedef itk::TranslationTransform< double, Dimension > TranslationTransformType;

class TranslateCommandIterationUpdate : public itk::Command
{
public:
  typedef  TranslateCommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

protected:
  TranslateCommandIterationUpdate() {};

public:

  typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
  typedef const OptimizerType                         *OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer affineOptimizer =
                         dynamic_cast< OptimizerPointer >( object );

    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }

    std::cout << affineOptimizer->GetCurrentIteration() << " = ";
    std::cout << affineOptimizer->GetValue() << " : ";
    std::cout << affineOptimizer->GetCurrentPosition() << std::endl;
  }
};

FloatImageType::Pointer ResampleImage( FloatImageType::Pointer inputImage, double scaleFactor )
{
  typedef itk::RecursiveGaussianImageFilter< FloatImageType, FloatImageType > GaussianFilterType;
  typedef itk::ResampleImageFilter< FloatImageType, FloatImageType > ResampleFilterType;
  typedef itk::IdentityTransform< double, 2 >  TransformType;
  typedef itk::LinearInterpolateImageFunction< FloatImageType, double > InterpolatorType;

  GaussianFilterType::Pointer smoother = GaussianFilterType::New();
  smoother->SetInput( inputImage );
  smoother->SetSigma( scaleFactor/2.0 );
  smoother->SetNormalizeAcrossScale( true );

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  resampler->SetTransform( transform );

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  resampler->SetInterpolator( interpolator );
  resampler->SetDefaultPixelValue( 0 );

  FloatImageType::SpacingType spacing;
  spacing[0] = inputImage->GetSpacing()[0];
  spacing[1] = inputImage->GetSpacing()[1];

  resampler->SetOutputSpacing( spacing );
  resampler->SetOutputOrigin( inputImage->GetOrigin() );
  resampler->SetOutputDirection( inputImage->GetDirection() );

  FloatImageType::SizeType size;
  size[0] = ((double)inputImage->GetLargestPossibleRegion().GetSize()[0])/scaleFactor;
  size[1] = ((double)inputImage->GetLargestPossibleRegion().GetSize()[1])/scaleFactor;
  resampler->SetSize( size );

  resampler->SetInput( smoother->GetOutput() );

  std::cout<<"Resampling at scale:"<<scaleFactor<<std::endl<<std::flush;
  try
  {
    resampler->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception in subsampling:\n"
      	<< excep << std::endl;
  }
  FloatImageType::Pointer ouputImage = resampler->GetOutput();
  return ouputImage;
}

RegistrationType::Pointer ResampleAndRegisterWithTranslations
			( FloatImageType::Pointer fixedIm, FloatImageType::Pointer movingIm )
{

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::MeanSquaresImageToImageMetric< FloatImageType, FloatImageType > MetricType;
  typedef itk::LinearInterpolateImageFunction< FloatImageType, double > InterpolatorType;

  double largestDim = fixedIm->GetLargestPossibleRegion().GetSize()[0];
  if( fixedIm->GetLargestPossibleRegion().GetSize()[1] > largestDim )
    largestDim = fixedIm->GetLargestPossibleRegion().GetSize()[1];
  
  double numIter = 0;
  while( (largestDim/std::pow(10.0,numIter)) > 100.00 )
     numIter++;

  TranslationTransformType::Pointer translationTransform = TranslationTransformType::New();
  RegistrationType::Pointer   OutputTranslationRegistration;
  ParametersType initialParameters( translationTransform->GetNumberOfParameters() );
  initialParameters[0] = 0.0;  // Initial offset in mm along X
  initialParameters[1] = 0.0;  // Initial offset in mm along Y

  for( int i=numIter; i>-1; i-- )
  {
    //Resample image
    double scaleFactorCurrentIter = std::pow( 10.0, ((double)i) );

    FloatImageType::Pointer resampledFixed;
    FloatImageType::Pointer resampledMoving;
    if( i )
    {
      resampledFixed  = ResampleImage( fixedIm,  scaleFactorCurrentIter );
      resampledMoving = ResampleImage( movingIm, scaleFactorCurrentIter );
    }
    else
    {
      resampledFixed  = fixedIm;
      resampledMoving = movingIm;
    }
    initialParameters[0] *= 10;
    initialParameters[1] *= 10;
    
    MetricType::Pointer         translationMetric		= MetricType::New();
    TranslationTransformType::Pointer translationTransform	= TranslationTransformType::New();
    OptimizerType::Pointer      translationOptimizer		= OptimizerType::New();
    InterpolatorType::Pointer   translationInterpolator		= InterpolatorType::New();
    RegistrationType::Pointer   translationRegistration		= RegistrationType::New();
    translationRegistration->SetMetric(        translationMetric        );
    translationRegistration->SetOptimizer(     translationOptimizer     );
    translationRegistration->SetTransform(     translationTransform     );
    translationRegistration->SetInterpolator(  translationInterpolator  );

    translationRegistration->SetFixedImage (    fixedIm    );
    translationRegistration->SetMovingImage(   movingIm    );

    translationRegistration->SetFixedImageRegion( fixedIm->GetBufferedRegion() );
    translationRegistration->SetInitialTransformParameters( initialParameters );
    translationOptimizer->SetMaximumStepLength( 4.00 );
    translationOptimizer->SetMinimumStepLength( 0.01 );
    translationOptimizer->SetNumberOfIterations( 200 );
    TranslateCommandIterationUpdate::Pointer translateObserver = TranslateCommandIterationUpdate::New();
    translationOptimizer->AddObserver( itk::IterationEvent(), translateObserver );
    try
    {
      translationRegistration->Update();
    }
    catch( itk::ExceptionObject & err )
    {
      std::cerr << "ExceptionObject caught in registration:\n";
      std::cerr << err << std::endl;
    }
    initialParameters = translationRegistration->GetLastTransformParameters();
    const double TranslationAlongX = initialParameters[0];
    const double TranslationAlongY = initialParameters[1];
    const unsigned int numberOfIterationsTranslation = translationOptimizer->GetCurrentIteration();
    const double bestValueTranslate = translationOptimizer->GetValue();
    std::cout << "Result: Translation X = " << TranslationAlongX  << "\tTranslation Y = " << TranslationAlongY
	<< "\tIterations = " << numberOfIterationsTranslation << "\tMetric value  = " << bestValueTranslate
	<< std::endl;
    if( !i )
       OutputTranslationRegistration = translationRegistration;
  }
  return  OutputTranslationRegistration;
}

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << "   fixedImageFile  movingImageFile1 [movingImageFile2] [movingImageFile3] ... [movingImageFileN]"<< std::endl;
    std::cerr << "This executable registers File1 and then uses the transform on all subsequent files\n";
    return EXIT_FAILURE;
    }
  typedef itk::ImageFileReader< FloatImageType  > FloatImageReaderType;
  typedef itk::ImageFileReader< FloatImageType > FloatImageReaderType;
  typedef itk::ResampleImageFilter< FloatImageType, FloatImageType > ResampleFilterType;

  FloatImageReaderType::Pointer  fixedImageReader  = FloatImageReaderType::New();
  TIFFIOType::Pointer tiffIO1 = TIFFIOType::New();
  // tiffIO1->SetPixelType(itk::ImageIOBase::USHORT);
  FloatImageReaderType::Pointer movingImageReader = FloatImageReaderType::New();
  TIFFIOType::Pointer tiffIO2 = TIFFIOType::New();
  // tiffIO2->SetPixelType(itk::ImageIOBase::USHORT);

  fixedImageReader->SetFileName(  argv[1] );
  fixedImageReader->SetImageIO(tiffIO1);
  movingImageReader->SetFileName( argv[2] );
  movingImageReader->SetImageIO(tiffIO2);

  try
  {
    fixedImageReader->Update();
    movingImageReader->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught in image reader:\n";
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  FloatImageType::Pointer fixedImage  = fixedImageReader->GetOutput();
  FloatImageType::Pointer movingImage = movingImageReader->GetOutput();

  RegistrationType::Pointer outputTranslationRegistration =
    ResampleAndRegisterWithTranslations( fixedImage, movingImage );
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( movingImageReader->GetOutput() );
  resampler->SetTransform( outputTranslationRegistration->GetOutput()->Get() );
  resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( 100 );
  try
  {
    resampler->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught in translation resampler:\n";
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  typedef  unsigned short  OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter< FloatImageType, OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  TIFFIOType::Pointer tiffIO3 = TIFFIOType::New();
  // tiffIO3->SetPixelType(itk::ImageIOBase::USHORT);
  CastFilterType::Pointer  caster =  CastFilterType::New();

  std::string name = argv[2];
  unsigned found = name.find_last_of(".");
  name = name.substr(0,found) + "_" + "registered.tif";
  writer->SetFileName( name.c_str() );

  caster->SetInput( resampler->GetOutput() ); //finalResampler->GetOutput()
  writer->SetInput( caster->GetOutput()   );
  writer->SetImageIO( tiffIO3 );
  writer->Update();

  for( int i=3; i<argc; ++i )
  {
    FloatImageReaderType::Pointer repeatMoveImageReader = FloatImageReaderType::New();
    repeatMoveImageReader->SetFileName( argv[i] );
	TIFFIOType::Pointer tiffIO4 = TIFFIOType::New();
    // tiffIO4->SetPixelType(itk::ImageIOBase::USHORT);
	repeatMoveImageReader->SetImageIO( tiffIO4 );
    try
    {
    repeatMoveImageReader->Update();
    }
    catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

    ResampleFilterType::Pointer repeatResampler = ResampleFilterType::New();
    repeatResampler->SetTransform( outputTranslationRegistration->GetOutput()->Get() );
    repeatResampler->SetInput( repeatMoveImageReader->GetOutput() );

    repeatResampler->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
    repeatResampler->SetOutputOrigin(  fixedImage->GetOrigin() );
    repeatResampler->SetOutputSpacing( fixedImage->GetSpacing() );
    repeatResampler->SetOutputDirection( fixedImage->GetDirection() );
    repeatResampler->SetDefaultPixelValue( 100 );

    WriterType::Pointer      repeatWriter =  WriterType::New();
    CastFilterType::Pointer  repeatCaster =  CastFilterType::New();

    std::string repeatName = argv[i];
    unsigned found = repeatName.find_last_of(".");
    repeatName = repeatName.substr(0,found) + "_" + "registered.tif";

    repeatWriter->SetFileName( repeatName.c_str() );
    repeatCaster->SetInput( repeatResampler->GetOutput() );
    repeatWriter->SetInput( repeatCaster->GetOutput()   );
    TIFFIOType::Pointer tiffIO5 = TIFFIOType::New();
    // tiffIO5->SetPixelType(itk::ImageIOBase::USHORT);
    repeatWriter->SetImageIO( tiffIO5 );

    try
    {
      repeatWriter->Update();
    }
    catch( itk::ExceptionObject & err )
    {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
