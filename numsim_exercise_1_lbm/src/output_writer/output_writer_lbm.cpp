#include "output_writer/output_writer_lbm.h"

#include <iostream>

OutputWriterLBM::OutputWriterLBM(Lattice_boltzmann lattice_boltzmann)
 : lattice_boltzmann_(lattice_boltzmann), fileNo_(0)
{
  // create "out" subdirectory if it does not yet exist
  int returnValue = system("mkdir -p out");
  if (returnValue != 0)
    std::cout << "Could not create subdirectory \"out\"." << std::endl;
}
