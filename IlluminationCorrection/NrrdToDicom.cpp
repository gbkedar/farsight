#include <vector>
#include <string>
#include <cmath>

#include <tinyxml/tinyxml.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPasteImageFilter.h"
#include "itkExtractImageFilter.h"

typedef itk::IndexValueType IIVT;
typedef itk::SizeValueType  ISVT;
typedef unsigned short	USPixelType;
typedef double		CostPixelType;
const unsigned int	Dimension3 = 3;
const unsigned int	Dimension2 = 2;
typedef itk::Image< USPixelType, Dimension3 > US3ImageType;
typedef itk::Image< USPixelType, Dimension2 > US2ImageType;

std::string nameTemplate;

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

template<typename InputImageType> itk::SmartPointer<InputImageType>
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

template<typename InputImageType> itk::SmartPointer<InputImageType>
CreateDefaultCoordsNAllocateSpace
  ( typename InputImageType::SizeType size )
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

struct TilePositionStructType
{  IIVT X, Y;
   unsigned index;
   bool operator<(const TilePositionStructType& a) const
   {
     if( Y==a.Y )
       return X<a.X;
     else
       return Y<a.Y;
   }
};

void usage( const char *funcName )
{
  std::cout << "USAGE:"
	    << " " << funcName << " InputImage InputMetaData.xml "
	    << "NumberOfThreads(Optional-default=24)\n";
}

double FindPixelPitch( TiXmlDocument doc )
{
  double returnValue = -1.0;
  TiXmlElement* rootElement = doc.FirstChildElement();
  const char* docname = rootElement->Value();
  if ( strcmp( docname, "OME" ) != 0 )
  {
    std::cout<<"OME xml file required. Wrong xml input.\n";
    return returnValue;
  }

  TiXmlElement *parentElement = rootElement->FirstChildElement();
  while (parentElement)
  {
    const char *parent = parentElement->Value();
    if ( strcmp( parent, "Image" ) == 0 )
    {
      TiXmlElement *imageElement = parentElement->FirstChildElement();
      while (imageElement)
      {
        const char *image = imageElement->Value();
	if ( strcmp( image, "Pixels" ) == 0 )
	{
	  returnValue = atof( imageElement->Attribute("PhysicalSizeX") );
	  if( atof(imageElement->Attribute("PhysicalSizeY"))!=returnValue )
	    std::cerr<<"Warning. Pixel size X and Y do not match using pixel size X only.\n";
	  return returnValue;
	}
        imageElement = imageElement->NextSiblingElement();;
      }
    }
    parentElement = parentElement->NextSiblingElement();
  }
  return returnValue;
}

void GetTilePositions( TiXmlDocument doc, std::vector< TilePositionStructType > &positionVec,
			double div )
{
  TiXmlElement* rootElement = doc.FirstChildElement();
  TiXmlElement *parentElement = rootElement->FirstChildElement();
  while( parentElement )
  {
    const char *parent = parentElement->Value();
    if( strcmp( parent, "Image" ) == 0 )
    {
      TiXmlElement *imageElement = parentElement->FirstChildElement();
      bool planeNotFound = true;
      while( imageElement && planeNotFound )
      {
        const char *imageAttr = imageElement->Value();
	if ( strcmp( imageAttr, "Pixels" ) == 0 )
	{
	  TiXmlElement *pixelsElement = imageElement->FirstChildElement();
	  while( pixelsElement )
	  {
	    const char *pixelsAttr = pixelsElement->Value();
	    if ( strcmp( pixelsAttr, "Plane" )== 0 )
	    {
	      planeNotFound = false;
	      TilePositionStructType currentTile;
	      currentTile.X = std::floor( atof( pixelsElement->Attribute("PositionX") )/div + 0.5 );
	      currentTile.Y = std::floor( atof( pixelsElement->Attribute("PositionY") )/div + 0.5 );
	      if( positionVec.empty() )
		currentTile.index = 0;
	      else
		currentTile.index = positionVec.size();
	      positionVec.push_back( currentTile );
	      break;
	    }
	    pixelsElement = pixelsElement->NextSiblingElement();
	  }
	}
        imageElement = imageElement->NextSiblingElement();
      }
    }
    parentElement = parentElement->NextSiblingElement();
  }
  return;
}

void FlipCoordsForITK( IIVT xMin, IIVT xMax, IIVT yMin, IIVT yMax,
			std::vector< TilePositionStructType > &positionVec )
{
  //CHECK IF FUNCTION IS NEEDED
  return;
}

void GetDelta( ISVT &overlapX, ISVT &overlapY, std::vector< TilePositionStructType > &positionVec )
{
  ISVT i = 0;
  //Use different code L-R scan vs R-L scan
  double count=0, total=0;
  if( positionVec.at(0).X<positionVec.at(1).X )
  while( positionVec.at(i).X<positionVec.at(i+1).X )
  {
    total += ( positionVec.at(i+1).X - positionVec.at(i).X );
    ++i; ++count;
  }
  else
  while( positionVec.at(i).X>positionVec.at(i+1).X )
  {
    ++i;
  }
  overlapX = std::floor( total/count + 0.5 );
  if( overlapX%2 ) --overlapX;

  //Calculate Y overlap
  count=0; total=0; IIVT prevY = positionVec.at(0).Y;
  for( i=1; i<positionVec.size(); ++i )
  {
    if( prevY!= positionVec.at(i).Y )
    {
      total += std::abs( prevY - positionVec.at(i).Y );
      ++count;
      prevY= positionVec.at(i).Y;
    }
  }
  overlapY = std::floor( total/count + 0.5 );
  if( overlapY%2 ) --overlapY;
  return;
}

itk::SmartPointer<US2ImageType> GetTile( US3ImageType::Pointer readImage, itk::SizeValueType i )
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

US2ImageType::Pointer PasteNonOverLappingSections( US2ImageType::Pointer outputImage, US3ImageType::Pointer inputImage,
	ISVT overlapX, ISVT overlapY, std::vector< TilePositionStructType > &positionVec,
	IIVT xMin, IIVT yMin, IIVT xMax, IIVT yMax, int numThreads )
{
  typedef itk::PasteImageFilter <US2ImageType, US2ImageType> PasteImageFilterType;
/*#ifdef _OPENMP
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);
  #pragma omp parallel for
#endif*/
  for( ISVT i=0; i<positionVec.size(); ++i )
  {
    //Paste the one with the lower index first bled with higher later
    US2ImageType::Pointer currentTile = GetTile( inputImage, positionVec.at(i).index );
    //Compute region in montage
    US2ImageType::IndexType destIndex, srcIndex;
    US2ImageType::SizeType srcSize;
    if( !(positionVec.at(i).X - xMin) )
    {
      srcSize[0]   = currentTile->GetLargestPossibleRegion().GetSize()[0];
      destIndex[0] = 0; srcIndex[0] = 0;
    }
    else
    {
      srcSize[0]   = currentTile->GetLargestPossibleRegion().GetSize()[0]-overlapX;
      srcIndex[0]  = overlapX;
      destIndex[0] = positionVec.at(i).X - xMin + overlapX;
    }
    if( !(positionVec.at(i).Y - yMin) )
    {
      srcSize[1]   = currentTile->GetLargestPossibleRegion().GetSize()[1];
      destIndex[1] = 0; srcIndex[1] = 0;
    }
    else
    {
      srcSize[1]   = currentTile->GetLargestPossibleRegion().GetSize()[1]-overlapY;
      srcIndex[1]  = overlapY;
      destIndex[1] = positionVec.at(i).Y - yMin + overlapY;
    }
    US2ImageType::RegionType srcRegion; srcRegion.SetSize( srcSize ); srcRegion.SetIndex( srcIndex );
    PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New();
    pasteFilter->SetSourceImage( currentTile );
    pasteFilter->SetSourceRegion( srcRegion );
    pasteFilter->SetDestinationImage( outputImage );
    pasteFilter->SetDestinationIndex( destIndex );
    try
    {
//    #pragma omp critical
      pasteFilter->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception caught in slice extraction filter!" << excep << std::endl;
      exit (EXIT_FAILURE);
    }
    outputImage = pasteFilter->GetOutput();
  }
  return outputImage;
}

int main( int argc, char *argv[] )
{
  if( argc < 3 )
  {
    usage(argv[0]);
    std::cerr << "PRESS ENTER TO EXIT\n";
    getchar();
    return EXIT_FAILURE;
  }
  std::string inputImageName = argv[1]; //Name of the input image
  std::string inputMetaName  = argv[2]; //Name of the input metadata file
  unsigned found = inputImageName.find_last_of(".");
  nameTemplate = inputImageName.substr(0,found) + "_stitched.nrrd";

  int numThreads = 24;
  if( argc == 4 )
    numThreads = atoi( argv[2] );

  TiXmlDocument doc;
  if ( !doc.LoadFile( inputMetaName.c_str() ) )
  {
    return EXIT_FAILURE;
    std::cout<<"Couldn't find xml file\n";
  }


  double pixelPitch = FindPixelPitch( doc );
  std::cout<<"Pitch found "<< pixelPitch << std::endl << std::flush;
  if( pixelPitch<0 )
  {
    std::cout<<"Negative pixel pitch. Aborting!\n"<< std::flush;;
    return EXIT_FAILURE;
  }
  
  std::vector< TilePositionStructType > positionVec;
  GetTilePositions( doc, positionVec, pixelPitch );

  IIVT xMin=positionVec.at(0).X, yMin=positionVec.at(0).Y, xMax=positionVec.at(0).X, yMax=positionVec.at(0).Y;

  for( ISVT i=1; i<positionVec.size(); ++i )
  {
    if( positionVec.at(i).X < xMin ) xMin=positionVec.at(i).X;
    if( positionVec.at(i).Y < yMin ) yMin=positionVec.at(i).Y;
    if( positionVec.at(i).X > xMax ) xMax=positionVec.at(i).X;
    if( positionVec.at(i).Y > yMax ) yMax=positionVec.at(i).Y;
  }

  //FlipCoordsForITK( xMin, xMax, yMin, yMax, positionVec );
  //Read tiles and get size
  US3ImageType::Pointer inputImage = ReadITKImage<US3ImageType>( inputImageName );
  US3ImageType::SizeType nrrdSize;
  nrrdSize[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
  nrrdSize[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
  nrrdSize[2] = inputImage->GetLargestPossibleRegion().GetSize()[2];

  //Calculate size and create blank outuput image
  US2ImageType::SizeType size; size[0] = xMax-xMin+nrrdSize[0]; size[1] = yMax-yMin+nrrdSize[1];
  US2ImageType::Pointer outputImage = CreateDefaultCoordsNAllocateSpace<US2ImageType>(size);

  ISVT overlapX, overlapY;

  GetDelta( overlapX, overlapY, positionVec );
  overlapX = nrrdSize[0]-overlapX;
  overlapY = nrrdSize[1]-overlapY;

  std::cout<<"X="<<xMin<<"\t"<<xMax<<"\t"<<outputImage->GetLargestPossibleRegion().GetSize()[1]
  	<<"\t\tY="<<yMin<<"\t"<<yMax<<"\t"<<outputImage->GetLargestPossibleRegion().GetSize()[0]
	<<"\t\toverlap"<<overlapX<<"\t"<<overlapY<<std::endl;

  std::sort( positionVec.begin(), positionVec.end() );

  outputImage = PasteNonOverLappingSections( outputImage, inputImage, overlapX, overlapY, positionVec, xMin, yMin, xMax, yMax, numThreads );

  WriteITKImage<US2ImageType>( outputImage, nameTemplate ); 

  return EXIT_SUCCESS;
}
