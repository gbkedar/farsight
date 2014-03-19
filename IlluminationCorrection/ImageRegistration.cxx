#include "itkTranslationTransform.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkAffineTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCommand.h"

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

class AffineCommandIterationUpdate : public itk::Command
{
public:
  typedef  AffineCommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  AffineCommandIterationUpdate() {};

public:
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef   const OptimizerType *                  OptimizerPointer;

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
      std::cout << affineOptimizer->GetCurrentIteration() << "   ";
      std::cout << affineOptimizer->GetValue() << "   ";
      std::cout << affineOptimizer->GetCurrentPosition();

      // Print the angle for the trace plot
      vnl_matrix<double> p(2, 2);
      p[0][0] = (double) affineOptimizer->GetCurrentPosition()[0];
      p[0][1] = (double) affineOptimizer->GetCurrentPosition()[1];
      p[1][0] = (double) affineOptimizer->GetCurrentPosition()[2];
      p[1][1] = (double) affineOptimizer->GetCurrentPosition()[3];
      vnl_svd<double> svd(p);
      vnl_matrix<double> r(2, 2);
      r = svd.U() * vnl_transpose(svd.V());
      double angle = vcl_asin(r[1][0]);
      std::cout << " AffineAngle: " << angle * 180.0 / vnl_math::pi << std::endl;
    }
};


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

  const    unsigned int    Dimension = 2;
  typedef  float           PixelType;

  typedef itk::Image< PixelType, Dimension >  FixedImageType;
  typedef itk::Image< PixelType, Dimension >  MovingImageType;

  typedef itk::AffineTransform< double, Dimension  > AffineTransformType;
  typedef itk::TranslationTransform< double, Dimension > TranslationTransformType;

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::NormalizedCorrelationImageToImageMetric< FixedImageType, MovingImageType > NCRMetricType;
  typedef itk::MeanSquaresImageToImageMetric< FixedImageType, MovingImageType > MetricType;
  typedef itk::LinearInterpolateImageFunction< MovingImageType, double > InterpolatorType;
  typedef itk::ImageRegistrationMethod< FixedImageType, MovingImageType > RegistrationType;

  MetricType::Pointer         translationMetric        = MetricType::New();
  TranslationTransformType::Pointer      translationTransform     = TranslationTransformType::New();
  OptimizerType::Pointer      translationOptimizer     = OptimizerType::New();
  InterpolatorType::Pointer   translationInterpolator  = InterpolatorType::New();
  RegistrationType::Pointer   translationRegistration  = RegistrationType::New();
  translationRegistration->SetMetric(        translationMetric        );
  translationRegistration->SetOptimizer(     translationOptimizer     );
  translationRegistration->SetTransform(     translationTransform     );
  translationRegistration->SetInterpolator(  translationInterpolator  );

  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );

  translationRegistration->SetFixedImage(    fixedImageReader->GetOutput()    );
  translationRegistration->SetMovingImage(   movingImageReader->GetOutput()   );
  try
    {
    fixedImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  translationRegistration->SetFixedImageRegion( fixedImageReader->GetOutput()->GetBufferedRegion() );
  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( translationTransform->GetNumberOfParameters() );

  initialParameters[0] = 0.0;  // Initial offset in mm along X
  initialParameters[1] = 0.0;  // Initial offset in mm along Y

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
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  ParametersType finalParametersTranslation = translationRegistration->GetLastTransformParameters();
  const double TranslationAlongX = finalParametersTranslation[0];
  const double TranslationAlongY = finalParametersTranslation[1];
  const unsigned int numberOfIterationsTranslation = translationOptimizer->GetCurrentIteration();
  const double bestValueTranslate = translationOptimizer->GetValue();
  std::cout << "Result = " << std::endl;
  std::cout << " Translation X = " << TranslationAlongX  << std::endl;
  std::cout << " Translation Y = " << TranslationAlongY  << std::endl;
  std::cout << " Iterations    = " << numberOfIterationsTranslation << std::endl;
  std::cout << " Metric value  = " << bestValueTranslate << std::endl;

  typedef itk::ResampleImageFilter< MovingImageType, FixedImageType > ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput( movingImageReader->GetOutput() );
  resampler->SetTransform( translationRegistration->GetOutput()->Get() );
  FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
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
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  //AFFINE STARTS HERE
  NCRMetricType::Pointer         affineMetric        = NCRMetricType::New();
  OptimizerType::Pointer      affineOptimizer     = OptimizerType::New();
  InterpolatorType::Pointer   affineInterpolator  = InterpolatorType::New();
  RegistrationType::Pointer   affineRegistration  = RegistrationType::New();

  affineRegistration->SetMetric(        affineMetric        );
  affineRegistration->SetOptimizer(     affineOptimizer     );
  affineRegistration->SetInterpolator(  affineInterpolator  );

  AffineTransformType::Pointer  affineTransform = AffineTransformType::New();
  affineRegistration->SetTransform( affineTransform );

  affineRegistration->SetFixedImage(    fixedImageReader->GetOutput()    );
  affineRegistration->SetMovingImage(   resampler->GetOutput()   );

  affineRegistration->SetFixedImageRegion(
     fixedImageReader->GetOutput()->GetBufferedRegion() );

  typedef itk::CenteredTransformInitializer< AffineTransformType, FixedImageType,
					MovingImageType >  TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform(   affineTransform );
  initializer->SetFixedImage(  fixedImageReader->GetOutput() );
  initializer->SetMovingImage( resampler->GetOutput() );
  //initializer->MomentsOn();
  initializer->InitializeTransform();

  affineRegistration->SetInitialTransformParameters( affineTransform->GetParameters() );
  double translationScale = 
    fixedImage->GetLargestPossibleRegion().GetSize()[0] > 
    fixedImage->GetLargestPossibleRegion().GetSize()[1] ?
    fixedImage->GetLargestPossibleRegion().GetSize()[0] :
    fixedImage->GetLargestPossibleRegion().GetSize()[1]; //1.0 / 10000.0;
  translationScale = 1.0/(translationScale*10.0);

  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType affineOptimizerScales( affineTransform->GetNumberOfParameters() );

  affineOptimizerScales[0] =  0.001;
  affineOptimizerScales[1] =  0.001;
  affineOptimizerScales[2] =  0.001;
  affineOptimizerScales[3] =  0.001;
  affineOptimizerScales[4] =  translationScale;
  affineOptimizerScales[5] =  translationScale;

  affineOptimizer->SetScales( affineOptimizerScales );

  double steplength = 0.01;
  unsigned int maxNumberOfIterations = 100;

  affineOptimizer->SetMaximumStepLength( steplength );
  affineOptimizer->SetMinimumStepLength( 0.0001 );
  affineOptimizer->SetNumberOfIterations( maxNumberOfIterations );
  affineOptimizer->MinimizeOn();

  // Create the Command observer and register it with the affineOptimizer.
  //
  AffineCommandIterationUpdate::Pointer affineObserver = AffineCommandIterationUpdate::New();
  affineOptimizer->AddObserver( itk::IterationEvent(), affineObserver );

  try
    {
    affineRegistration->Update();
    std::cout << "Optimizer stop condition: "
              << affineRegistration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  OptimizerType::ParametersType finalParametersAffine =
                    affineRegistration->GetLastTransformParameters();

  const double finalRotationCenterX = affineTransform->GetCenter()[0];
  const double finalRotationCenterY = affineTransform->GetCenter()[1];
  const double finalTranslationX    = finalParametersAffine[4];
  const double finalTranslationY    = finalParametersAffine[5];

  const unsigned int numberOfIterationsAffine = affineOptimizer->GetCurrentIteration();
  const double bestValueAffine = affineOptimizer->GetValue();

  std::cout << "Result = " << std::endl;
  std::cout << " Center X      = " << finalRotationCenterX  << std::endl;
  std::cout << " Center Y      = " << finalRotationCenterY  << std::endl;
  std::cout << " Translation X = " << finalTranslationX  << std::endl;
  std::cout << " Translation Y = " << finalTranslationY  << std::endl;
  std::cout << " Iterations    = " << numberOfIterationsAffine << std::endl;
  std::cout << " Metric value  = " << bestValueAffine    << std::endl;

  vnl_matrix<double> p(2, 2);
  p[0][0] = (double) finalParametersAffine[0];
  p[0][1] = (double) finalParametersAffine[1];
  p[1][0] = (double) finalParametersAffine[2];
  p[1][1] = (double) finalParametersAffine[3];
  vnl_svd<double> svd(p);
  vnl_matrix<double> r(2, 2);
  r = svd.U() * vnl_transpose(svd.V());
  double angle = vcl_asin(r[1][0]);

  const double angleInDegrees = angle * 180.0 / vnl_math::pi;

  std::cout << " Scale 1         = " << svd.W(0)        << std::endl;
  std::cout << " Scale 2         = " << svd.W(1)        << std::endl;
  std::cout << " Angle (degrees) = " << angleInDegrees  << std::endl;


  AffineTransformType::Pointer finalTransform = AffineTransformType::New();

  finalTransform->SetParameters( finalParametersAffine );
  finalTransform->SetFixedParameters( affineTransform->GetFixedParameters() );

  ResampleFilterType::Pointer finalResampler = ResampleFilterType::New();

  finalResampler->SetTransform( finalTransform );
  finalResampler->SetInput( resampler->GetOutput() );

  finalResampler->SetSize(    fixedImage->GetLargestPossibleRegion().GetSize() );
  finalResampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  finalResampler->SetOutputSpacing( fixedImage->GetSpacing() );
  finalResampler->SetOutputDirection( fixedImage->GetDirection() );
  finalResampler->SetDefaultPixelValue( 100 );

  typedef  unsigned short  OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  std::string name = argv[2];
  unsigned found = name.find_last_of(".");
  name = name.substr(0,found) + "_" + "registered.tif";
  writer->SetFileName( name.c_str() );

  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  writer->Update();

  //Add initial translation to the affine transform
  finalParametersAffine[4] += finalParametersTranslation[0];
  finalParametersAffine[5] += finalParametersTranslation[1];

  for( int i=3; i<argc; ++i )
  {
    MovingImageReaderType::Pointer repeatMoveImageReader = MovingImageReaderType::New();
    repeatMoveImageReader->SetFileName( argv[i] );
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
    AffineTransformType::Pointer combinedTransform = AffineTransformType::New();

    combinedTransform->SetParameters( finalParametersAffine );
    combinedTransform->SetFixedParameters( affineTransform->GetFixedParameters() );

    ResampleFilterType::Pointer repeatResampler = ResampleFilterType::New();

    repeatResampler->SetTransform( combinedTransform );
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
