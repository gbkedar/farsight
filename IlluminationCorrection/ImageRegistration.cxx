#include "itkTranslationTransform.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
//#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkRegionOfInterestImageFilter.h"
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

#define DEBUG_SCALES 1

template<typename InputImageType> void WriteITKImage
  ( typename InputImageType::Pointer inputImagePointer,
    std::string outputName )
{
  typedef  itk::TIFFImageIO TIFFIOType;
  typedef typename itk::ImageFileWriter< InputImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  TIFFIOType::Pointer tiffIO = TIFFIOType::New();
  writer->SetFileName( outputName.c_str() );
  writer->SetInput( inputImagePointer );
  writer->SetImageIO( tiffIO );
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

template<typename InputImageType> typename InputImageType::Pointer
  ReadITKImage( std::string inputName )
{
  typedef  itk::TIFFImageIO TIFFIOType;
  typedef typename itk::ImageFileReader< InputImageType > ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  TIFFIOType::Pointer tiffIO = TIFFIOType::New();
  reader->SetFileName( inputName.c_str() );
  reader->SetImageIO( tiffIO );
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

const    unsigned int    Dimension = 2;
typedef  float           PixelType;
typedef  unsigned short  OutputPixelType;

typedef itk::Image< PixelType, Dimension >  FloatImageType;
typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;
typedef itk::CropImageFilter< FloatImageType, FloatImageType > CropImageFilterType;
typedef itk::ImageRegistrationMethod< FloatImageType, FloatImageType > RegistrationType;
typedef RegistrationType::ParametersType ParametersType;
typedef itk::TranslationTransform< double, Dimension > TranslationTransformType;
typedef itk::ResampleImageFilter< FloatImageType, FloatImageType > ResampleFilterType;

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

FloatImageType::Pointer ResampleByTranslating( FloatImageType::Pointer movingImage,
	FloatImageType::Pointer fixedImage,
	RegistrationType::Pointer outputTranslationRegistration )
{
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( movingImage );
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
  }
  FloatImageType::Pointer outputImage = resampler->GetOutput();
  return outputImage;
}


FloatImageType::Pointer ResampleImageByScaling
				( FloatImageType::Pointer inputImage, double scaleFactor )
{
  typedef itk::RecursiveGaussianImageFilter< FloatImageType, FloatImageType > GaussianFilterType;
  typedef itk::IdentityTransform< double, 2 >  TransformType;
  typedef itk::LinearInterpolateImageFunction< FloatImageType, double > InterpolatorType;

  GaussianFilterType::Pointer smoother = GaussianFilterType::New();
  smoother->SetInput( inputImage );
  smoother->SetSigma( ((double)inputImage->GetSpacing()[0])*scaleFactor/2.0 );
  smoother->SetNormalizeAcrossScale( true );

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  resampler->SetTransform( transform );

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  resampler->SetInterpolator( interpolator );
  resampler->SetDefaultPixelValue( 0 );

  FloatImageType::SpacingType spacing;
  spacing[0] = inputImage->GetSpacing()[0]*scaleFactor;
  spacing[1] = inputImage->GetSpacing()[1]*scaleFactor;

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
			( FloatImageType::Pointer fixedIm, FloatImageType::Pointer movingIm,
			  double xTranslate, double yTranslate )
{

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::MeanSquaresImageToImageMetric< FloatImageType, FloatImageType > MSQMetricType;
//  typedef itk::NormalizedCorrelationImageToImageMetric< FloatImageType, FloatImageType >
//  									       NCRMetricType;
  typedef itk::LinearInterpolateImageFunction< FloatImageType, double > InterpolatorType;

  double largestDim = fixedIm->GetLargestPossibleRegion().GetSize()[0];
  if( fixedIm->GetLargestPossibleRegion().GetSize()[1] > largestDim )
    largestDim = fixedIm->GetLargestPossibleRegion().GetSize()[1];
  
  double numIter = 0;
//  while( (largestDim/std::pow(10.0,numIter)) > 1000.00 )
//     numIter++;

  TranslationTransformType::Pointer translationTransform = TranslationTransformType::New();
  RegistrationType::Pointer   OutputTranslationRegistration;
  ParametersType initialParameters( translationTransform->GetNumberOfParameters() );
  initialParameters[0] = xTranslate;  // Initial offset in mm along X
  initialParameters[1] = yTranslate;  // Initial offset in mm along Y

  for( int i=numIter; i>-1; i-- )
  {
    //Resample image
    double scaleFactorCurrentIter = std::pow( 10.0, ((double)i) );

    FloatImageType::Pointer resampledFixed;
    FloatImageType::Pointer resampledMoving;
    if( i )
    {
      resampledFixed  = ResampleImageByScaling( fixedIm,  scaleFactorCurrentIter );
      resampledMoving = ResampleImageByScaling( movingIm, scaleFactorCurrentIter );
    }
    else
    {
      resampledFixed  = fixedIm;
      resampledMoving = movingIm;
    }

#ifdef DEBUG_SCALES
    if( i )
    {
      std::stringstream ss;
      ss << scaleFactorCurrentIter;
      std::string fxNm = "fixedImScale" + ss.str() + ".tif";
      CastNWriteImage<FloatImageType,OutputImageType>(resampledFixed,fxNm);
    }
#endif //DEBUG_SCALES
    if( i!=numIter )
    {
      initialParameters[0] *= 10;
      initialParameters[1] *= 10;
    }
    
    MSQMetricType::Pointer      translationMetricMSQ		= MSQMetricType::New();
//    NCRMetricType::Pointer      translationMetricNCR		= NCRMetricType::New();
    TranslationTransformType::Pointer translationTransform	= TranslationTransformType::New();
    OptimizerType::Pointer      translationOptimizer		= OptimizerType::New();
    InterpolatorType::Pointer   translationInterpolator		= InterpolatorType::New();
    RegistrationType::Pointer   translationRegistration		= RegistrationType::New();
//    if( 0 )
//    {
//      translationMetricNCR->SubtractMeanOn();
//      translationRegistration->SetMetric(        translationMetricNCR        );
//    }
//    else
//    {
      translationRegistration->SetMetric(        translationMetricMSQ        );
//    }
    translationRegistration->SetOptimizer(     translationOptimizer     );
    translationRegistration->SetTransform(     translationTransform     );
    translationRegistration->SetInterpolator(  translationInterpolator  );

    translationRegistration->SetFixedImage (    resampledFixed    );
    translationRegistration->SetMovingImage(   resampledMoving    );

    double maxStepLength = 1.0;//*std::pow(10.0,numIter);
    double minStepLength = 0.001;//0.01/std::pow(10.0,numIter);
    translationRegistration->SetFixedImageRegion( fixedIm->GetBufferedRegion() );
    translationRegistration->SetInitialTransformParameters( initialParameters );
    translationOptimizer->SetMaximumStepLength( maxStepLength );
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
#ifdef DEBUG_SCALES
    if( i )
    {
      std::stringstream ss;
      ss << scaleFactorCurrentIter;
      std::string mvNmTr = "movingImScaleTr" + ss.str() + ".tif";
      std::string mvNm = "movingImScale" + ss.str() + ".tif";
      FloatImageType::Pointer TransIm = ResampleByTranslating( resampledMoving, resampledFixed, translationRegistration );
      CastNWriteImage<FloatImageType,OutputImageType>(TransIm,mvNmTr);
      CastNWriteImage<FloatImageType,OutputImageType>(resampledMoving,mvNm);
    }
#endif //DEBUG_SCALES
  }
  return  OutputTranslationRegistration;
}

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << "   fixedImageFile  movingImageFile1 [TranslationX] [TranslationY] [movingImageFile2] [movingImageFile3]";
    std::cerr << "... [movingImageFileN]"<< std::endl;
    std::cerr << "This executable registers File1 and then uses the transform on all subsequent files\n";
    return EXIT_FAILURE;
    }

  std::string fixedImageName =  argv[1];
  std::string movingImageName = argv[2];

  double initialX = 0 , initialY = 0;

  if( argc > 4 )
  {
    initialX = atof(argv[3]); initialY = atof(argv[4]);
  }

  FloatImageType::Pointer fixedImage  = ReadITKImage<FloatImageType>(fixedImageName);
  FloatImageType::Pointer movingImage = ReadITKImage<FloatImageType>(movingImageName);

  FloatImageType::Pointer fixedImageCrop = ExtractCenterNPercentOfImage<FloatImageType>
  		( fixedImage, 60.0 );
  FloatImageType::Pointer movingImageCrop = ExtractCenterNPercentOfImage<FloatImageType>
  		( movingImage, 60.0 );
  RegistrationType::Pointer outputTranslationRegistration =
    ResampleAndRegisterWithTranslations( fixedImageCrop, movingImageCrop, initialX, initialY );
  FloatImageType::Pointer TranslatedImage = ResampleByTranslating( movingImage,
		fixedImage, outputTranslationRegistration );

  unsigned found = movingImageName.find_last_of(".");
  std::string outputFilename = movingImageName.substr(0,found) + "_" + "registered.tif";
  CastNWriteImage<FloatImageType,OutputImageType>(TranslatedImage,outputFilename);

  for( int i=5; i<argc; ++i )
  {
    movingImageName = argv[i];
    found = movingImageName.find_last_of(".");
    outputFilename = movingImageName.substr(0,found) + "_" + "registered.tif";
    FloatImageType::Pointer repeatMovingIm = ReadITKImage<FloatImageType>(movingImageName);
    FloatImageType::Pointer repeatTranslatedImage = ResampleByTranslating( repeatMovingIm,
		fixedImage, outputTranslationRegistration );
    CastNWriteImage<FloatImageType,OutputImageType>(repeatTranslatedImage,outputFilename);
  }
  return EXIT_SUCCESS;
}
