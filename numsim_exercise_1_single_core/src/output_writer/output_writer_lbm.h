#pragma once

#include "lattice_boltzmann/lattice_boltzmann.h"

#include <memory>

/** Inteface class for writing simulation data output.
 */
class OutputWriterLBM
{
public:
  //! constructor
  //! @param lattice boltzmann object that will contain all the data to be written to the file
  OutputWriterLBM(Lattice_boltzmann lattice_boltzmann);

  //! write current velocities to file, filename is output_<count>.vti
  virtual void writeFile(double currentTime, const Lattice_boltzmann &lattice_boltzmann) = 0;

protected:

  Lattice_boltzmann lattice_boltzmann_;  //< a pointer to the lattice boltzmann method which contains all data that will be written to the file
  int fileNo_;   //< a counter that increments for every file, this number is part of the file name of output files
};
