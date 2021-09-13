/** \file
  * 
  * This code demonstrates the use of unsigned data types.
  * These can be used to set flags using bitwise operators
  *
  * An minimal example is given with generic operations functionals
  **/

#include "palabos2D.h"
#include "palabos2D.hh"


using namespace plb;
typedef double T;
typedef unsigned int unsignedType;

// For example purposes, simple functions to set and check for flags are provided
// These can be generalized for any location 
// with unsigned long, 64 flag positions are availibale
unsignedType setFlag1(unsignedType fieldValue){
    unsigned long loc = 0b1 << 1;
    return (fieldValue | loc);
}

unsignedType setFlag2(unsignedType fieldValue){
    unsigned long loc = 0b1 << 2;
    return (fieldValue | loc);
}

bool checkFlag1(unsignedType fieldValue){
    unsignedType loc = 0b1 << 1;
    return (fieldValue & loc) != 0u;
}


int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    makeDirectory("tmp",false);
    
    // Create a scalar field with unsigned data type
    MultiScalarField2D<unsignedType> flagField(100,100,0);

    Box2D square1(40,60,40,60);
    Box2D square2(20,80,20,80);

    // set Flags in the specified regions
    apply( setFlag2, flagField, square2);
    apply( setFlag1, flagField, square1);

    // evaluate returns a unqiue_ptr ScalarField with the same data type as the provided ScalarField
    // thus it needs to be converted to boolean.
    std::unique_ptr<MultiScalarField2D<bool>> boolFlag1 = 
        copyConvert<unsignedType, bool> (
            *evaluate( checkFlag1, flagField )
        );

    // Write images to show that only where the Flag1 was set, the evaluate returns true
    ImageWriter<unsignedType> imageWriter("leeloo.map");
    imageWriter.writeScaledGif("flagField", flagField, 1000, 1000);
    imageWriter.writeScaledGif("flag1", *copyConvert<bool, unsignedType>(*boolFlag1), 1000, 1000);
	
    // Write VTK to show that only where the Flag1 was set, the evaluate returns true
	VtkImageOutput2D<T> vtkOut("vtkOut", (T)1);
    vtkOut.writeData<unsignedType>(*copyConvert<unsignedType, T>(flagField), "flagField", (T)1);
    vtkOut.writeData<float>(*copyConvert<bool, T>(*boolFlag1), "flag1", (T)1);
    
	pcout << "Finished running unsignedTypes" << std::endl;
}
