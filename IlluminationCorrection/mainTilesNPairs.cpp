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

//#define DEBUG_GenerateRegistrationPairs

#include <tinyxml/tinyxml.h>

#include <vector>
#include <string>
#include <cstdio>
#include <sstream>
#include <iomanip>

#include "boost/filesystem/path.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/lexical_cast.hpp"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkShiftScaleImageFilter.h"

typedef unsigned short USPixelType;
typedef unsigned short UCPixelType;
const unsigned int     Dimension3 = 3;
const unsigned int     Dimension2 = 2;
typedef itk::Image< USPixelType, Dimension3 > US3ImageType;
typedef itk::Image< USPixelType, Dimension2 > US2ImageType;
typedef itk::Image< UCPixelType, Dimension2 > UC2ImageType;

void usage( const char *funcName )
{
  std::cout << "USAGE:"
	    << " " << funcName << " InputImage InputMetadataXML "
	    << "OutputTemplate <PairList>.txt NumChannelForRegistration(default 0)"
	    << " No_Ignore_Diagonal_Pairs(1 or 0-default 0) Scaling_comp_tile(default 0)\n";
  std::cout << "First three inputs are filenames\n";
}

unsigned GetNumberOfFilesFromXML( std::string inputXml )
{
  TiXmlDocument doc;
  if ( !doc.LoadFile( inputXml.c_str() ) )
  {
    std::cout<<"Failed to load the xml file\n";
    return 0;
  }
  TiXmlElement* rootElement = doc.FirstChildElement();
  const char* docname = rootElement->Value();
  if ( strcmp( docname, "OME" ) != 0 )
  {
    std::cout<<"Program needs an OME xml file\n";
    return 0;
  }
  //Count number of Image elements
  unsigned countImages = 0;
  TiXmlElement * parentElement = rootElement->FirstChildElement();
  while( parentElement )
  {
    const char * parent = parentElement->Value();
    if( strcmp( parent, "Image" ) == 0 )
      ++countImages;
    parentElement = parentElement->NextSiblingElement();
  }
  if( !countImages )
    std::cout<<"Did not find any Image elements in XML\n";
#ifdef DEBUG_GenerateRegistrationPairs
  else
    std::cout<<"Found "<<countImages<<" files in XML\n";
#endif //DEBUG_GenerateRegistrationPairs
  return countImages;
}

typedef struct{
  std::string ID;
  double physicalSizeX, physicalSizeY, positionX, positionY;
  unsigned sizeC, sizeX, sizeY;
} TileInfo;

void GetTileInfo( std::vector< TileInfo > &tilesInfo,
       std::string inputXml )
{
  //Populate tile info
  TiXmlDocument doc;
  doc.LoadFile( inputXml.c_str() );
  TiXmlElement* rootElement = doc.FirstChildElement(); //Element OME
  TiXmlElement* parentElement = rootElement->FirstChildElement();
  while( parentElement )
  {
    const char * parent = parentElement->Value();
    if( strcmp( parent, "Image" ) == 0 )
    {
      TileInfo currentTileInfo;
      currentTileInfo.ID = parentElement->Attribute("ID");

      //Get the child named Pixels
      TiXmlElement * pixelElement = parentElement->FirstChildElement();
      while( ( strcmp( pixelElement->Value(), "Pixels" ) != 0 ) &&
	     ( pixelElement ) )
	pixelElement = pixelElement->NextSiblingElement();
      if( !pixelElement )
      {
	std::cout<< "Could not find pixel info for the xml element " << currentTileInfo.ID
		 << std::endl;
	exit( EXIT_FAILURE );
      }
      //Populate data
      currentTileInfo.physicalSizeX = atof( pixelElement->Attribute( "PhysicalSizeX" ) );
      currentTileInfo.physicalSizeY = atof( pixelElement->Attribute( "PhysicalSizeY" ) );
      currentTileInfo.sizeC         = atoi( pixelElement->Attribute( "SizeC" ) );
      currentTileInfo.sizeX         = atoi( pixelElement->Attribute( "SizeX" ) );
      currentTileInfo.sizeY         = atoi( pixelElement->Attribute( "SizeY" ) );

      //Get the first child element named plane
      TiXmlElement * planeElement = pixelElement->FirstChildElement();
      while( ( strcmp( planeElement->Value(), "Plane" ) != 0 ) &&
	     ( planeElement ) )
	planeElement = planeElement->NextSiblingElement();
      if( !planeElement )
      {
	std::cout<< "Could not find plane info for the xml element " << currentTileInfo.ID
		 << std::endl;
	exit( EXIT_FAILURE );
      }
      //Populate data and assume rest of the channels use the same X and Y position
      currentTileInfo.positionX = atof( planeElement->Attribute( "PositionX" ) );
      currentTileInfo.positionY = atof( planeElement->Attribute( "PositionY" ) );
      tilesInfo.push_back( currentTileInfo );
    }
    parentElement = parentElement->NextSiblingElement();
  }
#ifdef DEBUG_GenerateRegistrationPairs
  std::cout<< "Got info for "<<tilesInfo.size()<<" images\n";
  std::cout<< "First two elements read:\n"
  	   << "ID\tPhysicalSizeX\tPhysicalSizeY\tSizeC\tSizeX\tSizeY\t"
  	   << "PositionX\tPositionY\n" << tilesInfo.at(0).ID << "\t" << tilesInfo.at(0).physicalSizeX
	   << "\t\t"<<tilesInfo.at(0).physicalSizeY << "\t\t" << tilesInfo.at(0).sizeC
	   << "\t"<<tilesInfo.at(0).sizeX<<"\t"     << tilesInfo.at(0).sizeY
	   << "\t"<<tilesInfo.at(0).positionX       << "\t\t"<< tilesInfo.at(0).positionY << "\n"
  	   << "ID\tPhysicalSizeX\tPhysicalSizeY\tSizeC\tSizeX\tSizeY\t"
  	   << "PositionX\tPositionY\n" << tilesInfo.at(1).ID << "\t" << tilesInfo.at(1).physicalSizeX
	   << "\t\t"<<tilesInfo.at(1).physicalSizeY << "\t\t" << tilesInfo.at(1).sizeC
	   << "\t"<<tilesInfo.at(1).sizeX<<"\t"     << tilesInfo.at(1).sizeY
	   << "\t"<<tilesInfo.at(1).positionX       << "\t"<< tilesInfo.at(1).positionY << "\n";
#endif //DEBUG_GenerateRegistrationPairs
  return;
}

bool CheckOverlap( std::vector< TileInfo > &tilesInfo, unsigned i, unsigned j, bool ignoreDiagonals )
{
  //A diagonal match is one whose center is greater than +/- 20% of the size of the tile in x&y
  if( !ignoreDiagonals )
  {
    double twentyPcX = tilesInfo.at(i).sizeX*tilesInfo.at(i).physicalSizeX*0.2;
    double twentyPcY = tilesInfo.at(i).sizeY*tilesInfo.at(i).physicalSizeY*0.2;
    if( std::abs( tilesInfo.at(i).positionX-tilesInfo.at(j).positionX ) > twentyPcX &&
	std::abs( tilesInfo.at(i).positionY-tilesInfo.at(j).positionY ) > twentyPcY )
    {
#ifdef DEBUG_GenerateRegistrationPairs
//      std::cout << "Skipping tiles X1:"<< tilesInfo.at(i).positionX << "\tX2:"
//		<< tilesInfo.at(j).positionX << "\tY1:" << tilesInfo.at(i).positionY << "\tY2:"
//		<< tilesInfo.at(j).positionY <<"\n";
#endif //DEBUG_GenerateRegistrationPairs
      return false;
    }
  }

  double iTileCornerX[4], jTileCornerX[4], iTileCornerY[4], jTileCornerY[4];
  for( unsigned k=0; k<4; ++k )
  {
    if( k<2 )
    {
      iTileCornerX[k] = tilesInfo.at(i).positionX - ((double)tilesInfo.at(i).sizeX)/2 * tilesInfo.at(i).physicalSizeX;
      jTileCornerX[k] = tilesInfo.at(j).positionX - ((double)tilesInfo.at(j).sizeX)/2 * tilesInfo.at(j).physicalSizeX;
    }
    else
    {
      iTileCornerX[k] = tilesInfo.at(i).positionX + ((double)tilesInfo.at(i).sizeX)/2 * tilesInfo.at(i).physicalSizeX;
      jTileCornerX[k] = tilesInfo.at(j).positionX + ((double)tilesInfo.at(j).sizeX)/2 * tilesInfo.at(j).physicalSizeX;
    }
    if( k%2==0 )
    {
      iTileCornerY[k] = tilesInfo.at(i).positionY - ((double)tilesInfo.at(i).sizeY)/2 * tilesInfo.at(i).physicalSizeY;
      jTileCornerY[k] = tilesInfo.at(j).positionY - ((double)tilesInfo.at(j).sizeY)/2 * tilesInfo.at(j).physicalSizeY;
    }
    else
    {
      iTileCornerY[k] = tilesInfo.at(i).positionY + ((double)tilesInfo.at(i).sizeY)/2 * tilesInfo.at(i).physicalSizeY;
      jTileCornerY[k] = tilesInfo.at(j).positionY + ((double)tilesInfo.at(j).sizeY)/2 * tilesInfo.at(j).physicalSizeY;
    }
  }

  //Compare each point in j to the bounds of i
  for( unsigned k=0; k<4; ++k )
  {
    if( ( jTileCornerX[k] >= iTileCornerX[0] ) && ( jTileCornerX[k] <= iTileCornerX[2] ) &&
	( jTileCornerY[k] >= iTileCornerY[0] ) && ( jTileCornerY[k] <= iTileCornerY[1] ) )
	return true;
  }
  return false;
}

void GenerateRegistrationPairs( std::vector< TileInfo > &tilesInfo,
				std::vector< std::pair< unsigned, unsigned > > &registrationPairs,
				bool ignoreDiagonals )
{
  //Assume PositionX n PositionY are in the middle of each tile and that they are rectangles
  for( unsigned i=0; i<(tilesInfo.size()-1); ++i )
  {
    for( unsigned j=i+1; j<tilesInfo.size(); ++j )
    {
      bool overLap = CheckOverlap( tilesInfo, i, j, ignoreDiagonals );
      if( !overLap )
	overLap = CheckOverlap( tilesInfo, j, i, ignoreDiagonals ); //One may be contained entirely/partially in the other
      if( overLap )
      {
	std::pair <int,int> registerPair = std::make_pair( i, j );
	registrationPairs.push_back( registerPair );
      }
    }
  }
#ifdef DEBUG_GenerateRegistrationPairs
  for( unsigned i=0; i<registrationPairs.size(); ++i )
  {
    std::cout << registrationPairs.at(i).first << " " << registrationPairs.at(i).second << "\n";
  }
#endif //DEBUG_GenerateRegistrationPairs
}

std::string CheckWritePermissionsNCreateTempFolder()
{
  //Get cwd
  boost::filesystem::path getcwd = boost::filesystem::current_path();
  static const char alphanum[] = "012345ABCDEabcde";
  std::string s = "1111";
  std::string temp_dir_str = getcwd.string() + "/Temp_nuc_seg_" + s;
  boost::filesystem::path temp_dir = temp_dir_str.c_str();
  //Generate a new temp dir
  while( boost::filesystem::is_directory( temp_dir ) )
  {
    for( unsigned i = 0; i<4; ++i )
      s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    temp_dir_str = getcwd.string() + "/Temp_nuc_seg_" + s;
    temp_dir = temp_dir_str.c_str();
  }
  try
  {
    boost::filesystem::create_directory( temp_dir );
  }
  catch( const boost::filesystem::filesystem_error& ex )
  {
    std::cout << ex.what() << std::endl;
    exit (EXIT_FAILURE);
  }
  temp_dir_str = temp_dir_str + "/";

#ifdef DEBUG_GenerateRegistrationPairs
  std::cout<<"The current working dir is: "<<getcwd<<std::endl;
  std::cout<<"The current temp dir string is: "<<temp_dir_str<<std::endl;
#endif //DEBUG_GenerateRegistrationPairs

  return temp_dir_str;
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

bool RescaleNCastTile( US2ImageType::Pointer &currentTile, UC2ImageType::Pointer &currentTileUC2,
	US2ImageType::PixelType scaleLower, US2ImageType::PixelType scaleHigher )
{
  typedef itk::ShiftScaleImageFilter< US2ImageType, UC2ImageType > ShiftRescaleUS2UCType;
  double scaling = itk::NumericTraits<UC2ImageType::PixelType>::max()/(scaleHigher-scaleLower);
  ShiftRescaleUS2UCType::Pointer shiftRescaleUS2UC = ShiftRescaleUS2UCType::New();
  shiftRescaleUS2UC->SetShift( scaleLower );
  shiftRescaleUS2UC->SetScale( scaling );
  shiftRescaleUS2UC->SetInput( currentTile );
  try
  {
    shiftRescaleUS2UC->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught !" << excep << std::endl;
    exit (EXIT_FAILURE);
  }
  //Give warning if over/underflow too high
  double overUnderFlowThresh = 0.1* currentTile->GetLargestPossibleRegion().GetSize()[0] *
  				    currentTile->GetLargestPossibleRegion().GetSize()[1];
  currentTileUC2 = shiftRescaleUS2UC->GetOutput();
  currentTileUC2->Register();
  if( shiftRescaleUS2UC->GetUnderflowCount()>overUnderFlowThresh ||
      shiftRescaleUS2UC->GetOverflowCount()>overUnderFlowThresh )
  {
    std::cout<<"Over/underflow greather than 10 percent in current tile\n";
    return true;
  }
  return false;
}

std::string GenerateFileNameString( std::string tempFolder, std::string templateName, unsigned i, unsigned j )
{
  std::string outString;
  std::string OutFile;
  std::stringstream filess;
  filess << i << "_" << j;
  outString = tempFolder + templateName + "_" +  filess.str();
  return outString;
}

template <typename ImageType>
void WriteChannel( std::string &fileName, typename ImageType::Pointer &writeFile )
{
  typedef typename itk::ImageFileWriter< ImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( fileName.c_str() );
  writer->SetInput( writeFile );
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }
}

void ComputeScalingConstsFromTileZero( US2ImageType::Pointer &tileZero,
	US2ImageType::PixelType &scaleLower, US2ImageType::PixelType &scaleHigher )
{
  typedef itk::MedianImageFilter< US2ImageType, US2ImageType > MedianFilterType;
  typedef itk::ImageRegionConstIterator< US2ImageType > ConstIterType;
  
  //Median filter to remove noise
  MedianFilterType::Pointer medFilter = MedianFilterType::New();
  medFilter->SetInput( tileZero );
  medFilter->SetRadius( 3 );
  try
  {
    medFilter->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught !" << excep << std::endl;
    exit (EXIT_FAILURE);
  }
  US2ImageType::Pointer filteredIm = medFilter->GetOutput();

  //Setup and initialize histogram
  const itk::SizeValueType numBins = itk::NumericTraits<US2ImageType::PixelType>::max();
  itk::SizeValueType *tempHist;
  tempHist = (itk::SizeValueType*) malloc( sizeof(itk::SizeValueType) * numBins );
  for( US2ImageType::PixelType i=0; i<numBins; ++i )
    tempHist[i] = 0;

  //Populate histogram
  ConstIterType it( filteredIm, filteredIm->GetRequestedRegion() );
  for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
  {
    US2ImageType::PixelType pix = it.Get();
    ++tempHist[pix];
  }

  double totalPixels, totalIter, threshold;
  threshold = 0.02; //Change if you want to saturate more pixels
  totalIter = 0;
  totalPixels = filteredIm->GetLargestPossibleRegion().GetSize()[0] *
  		filteredIm->GetLargestPossibleRegion().GetSize()[1];
  //Iterate till the bin that has 2% of pixels from each end and those
  //are the lower and higher values for rescaling/casting
  US2ImageType::PixelType i;
  for( i=0; i<numBins; ++i )
  {
    totalIter += tempHist[i];
    double ratio = totalIter/totalPixels;
    if( ratio >= threshold )
    {
      if( i>0 )
	--i;
      break;
    }
  }
  scaleLower = i;
  totalIter = 0;
  for( i=numBins-1; i>0; --i )
  {
    totalIter += tempHist[i];
    double ratio = totalIter/totalPixels;
    if( ratio >= threshold )
    {
      if( i<(numBins-1) )
	++i;
      break;
    }
  }
  scaleHigher = i;

  if( scaleHigher==scaleLower || scaleHigher<scaleLower )
  {
    std::cout<<"Failed to compute histogram. Check input tile.\n";
    exit( EXIT_FAILURE );
  }
  return;
}

void CreateTempFolderNWriteInputChannelTiles
	( std::string inputImage, std::vector< TileInfo > &tilesInfo,
	  std::vector< std::string > &registerPairFileNames, std::string &tempFolder,
	  std::string &templateName, unsigned registrationChannel, unsigned histogramTile
	)
{
  typedef itk::ImageFileReader< US3ImageType >    ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader = ReaderType::New();
  reader->SetFileName( inputImage );
  reader->SetUseStreaming( true );
  try
  {
    reader->UpdateOutputInformation();//Update();
  }
  catch (itk::ExceptionObject &e)
  {
    std::cerr << e << std::endl;
    exit( EXIT_FAILURE );
  }
  US3ImageType::Pointer    readImage = reader->GetOutput();
  US3ImageType::RegionType regionUS3 = readImage->GetLargestPossibleRegion();
  US3ImageType::SizeType     sizeUS3 = regionUS3.GetSize();

#ifdef DEBUG_GenerateRegistrationPairs
  std::cout<<sizeUS3<<std::endl;
#endif

  if( sizeUS3[2] != tilesInfo.at(0).sizeC*tilesInfo.size() )
  {
    std::cout<< "Metadata size "
	     << tilesInfo.at(0).sizeC*tilesInfo.size()
	     << " and image size " << sizeUS3[2]
	     << " do not match\n";
    exit( EXIT_FAILURE );
  }

  tempFolder = CheckWritePermissionsNCreateTempFolder();
  std::string templateNameUC = templateName + "_UC";

  //Get rescaling constants
  US2ImageType::PixelType scaleLower, scaleHigher;
  US2ImageType::Pointer tileZero;
  GetTile( tileZero, readImage, (histogramTile*tilesInfo.at(0).sizeC+registrationChannel) );
  ComputeScalingConstsFromTileZero( tileZero, scaleLower, scaleHigher );
#ifdef DEBUG_GenerateRegistrationPairs
	std::cout<<"Rescaling constants: L="<<scaleLower<<"\tH= "<<scaleHigher<<"\n";
#endif

  for( unsigned i=0; i<tilesInfo.size(); ++i )
    for( unsigned j=0; j<tilesInfo.at(i).sizeC; ++j )
    {
      //Write each tile for transofrmation and stitching
      US2ImageType::Pointer currentTile;
      GetTile( currentTile, readImage, (i*tilesInfo.at(0).sizeC+j) );
      std::string fileName = GenerateFileNameString( tempFolder, templateName, j, i ) + ".tif";
      WriteChannel< US2ImageType >( fileName, currentTile );
      //Write 8-bit images for register pair
      if( j==registrationChannel )
      {
#ifdef DEBUG_GenerateRegistrationPairs
	std::cout<<(i*tilesInfo.at(0).sizeC+j)<<"\t";
#endif
	UC2ImageType::Pointer currentTileUC2;
	if( RescaleNCastTile( currentTile, currentTileUC2, scaleLower, scaleHigher ) )
	  std::cout<<"Path to tile "<<currentTileUC2<<std::endl;
	fileName = GenerateFileNameString( tempFolder, templateNameUC, j, i ) + ".tif";
	registerPairFileNames.push_back( fileName );
	WriteChannel< UC2ImageType >( fileName, currentTileUC2 );
	currentTileUC2->UnRegister();
      }
      currentTile->UnRegister();
    }
#ifdef DEBUG_GenerateRegistrationPairs
    std::cout<<std::endl;
#endif
}

void WritePairsFile( std::vector< std::string > &registerPairFileNames,
	std::vector< std::pair< unsigned, unsigned > > &registrationPairs, std::string &registrationFile )
{
  std::ofstream ofs(registrationFile.c_str(), std::ofstream::out);
  for( unsigned i=0; i<registrationPairs.size(); ++i )
  {
    size_t findPath1, findPath2;
    findPath1 = registerPairFileNames.at(registrationPairs.at(i).first).find_last_of("/\\");
    findPath2 = registerPairFileNames.at(registrationPairs.at(i).second).find_last_of("/\\");
    ofs << registerPairFileNames.at(registrationPairs.at(i).first).substr(findPath1+1)
    	<< " " << registerPairFileNames.at(registrationPairs.at(i).second).substr(findPath2+1) << "\n";
  }
  ofs.close();
}

int main(int argc, char *argv[])
{ 
  if( argc < 5 || argc > 8 )
  {
    usage(argv[0]);
    std::cerr << "PRESS ENTER TO EXIT\n";
    getchar();
    return EXIT_FAILURE;
  }

  std::string myName	       = argv[0]; //Just In case..
  std::string inputImage       = argv[1]; //Name of the input image
  std::string inputXml         = argv[2]; //Name of the xml file with the metadata
  std::string outputTemplate   = argv[3]; //Template filename for the tiles
  std::string registrationFile = argv[4]; //Registration filename
  unsigned numChannel = 0;		  //Number of the channel to be used to run registration
  bool noIgnoreDiagonal = false;	  //Don't write diagonal pairs to registration
  unsigned histogramTile = 0;		  //Tile on which the histogram comp for rescaling is performed
  if( argc == 6 )
    numChannel = atoi(argv[5]);
  if( argc == 7 && atoi(argv[6]) == 1 )
    noIgnoreDiagonal = true;
  if( argc == 8 )
    histogramTile = atoi(argv[7]);

  //Read the xml file and get number of files
  unsigned numberOfFiles = GetNumberOfFilesFromXML( inputXml );
  if( !numberOfFiles )
    return EXIT_FAILURE;

  std::vector< TileInfo > tilesInfo;
  GetTileInfo( tilesInfo, inputXml );

  std::vector< std::pair< unsigned, unsigned > > registrationPairs;
  GenerateRegistrationPairs( tilesInfo, registrationPairs, noIgnoreDiagonal );
  if( registrationPairs.empty() )
  {
    std::cout<<"Found no overlapping tiles in the metadata\n";
    return EXIT_FAILURE;
  }

  std::vector< std::string > registerPairFileNames;
  std::string tempFolder;
  CreateTempFolderNWriteInputChannelTiles( inputImage, tilesInfo, registerPairFileNames, tempFolder,
	outputTemplate, numChannel, histogramTile );
  registrationFile = tempFolder + registrationFile;
  WritePairsFile( registerPairFileNames, registrationPairs, registrationFile ); 

  return EXIT_SUCCESS;
}
