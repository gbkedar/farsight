// ############################################################################################################################################################################

//#includes goes here

#include "ftkNuclearSegmentationNic.h"


#include <iostream>


// ############################################################################################################################################################################

int main( int argc, char * argv[] ){


	ftk::ftkNucSecNic::ftkNuclearSegmentationNic<int> nucSecNic;
	nucSecNic.ResetRealeaseAll();
	nucSecNic.ReleaseMemory();


	std::cout << std::endl;
	return 0;
};

// ############################################################################################################################################################################