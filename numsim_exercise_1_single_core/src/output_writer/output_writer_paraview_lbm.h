#pragma once

#include "output_writer/output_writer_lbm.h"

#include <array>

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include <memory>

/** Write *.vti files that can be viewed with ParaView.
 *  The mesh that can be visualized in ParaView corresponds to the mesh of the computational domain.
 *  All values are given for the nodes of the mesh, i.e., the corners of each cell.
 *  This means, values will be interpolated because the values are stored at positions given by the staggered grid.
 */
class OutputWriterParaviewLBM :
  public OutputWriterLBM
{
public:
  //! constructor
  //! @param discretization shared pointer to the discretization object that will contain all the data to be written to the file
	OutputWriterParaviewLBM(Lattice_boltzmann lattice_boltzmann, std::array<double, 2> meshwidth, std::array<int, 2> nCells);

  //! write current velocities to file, filename is output_<count>.vti
  void writeFile(double currentTime, const Lattice_boltzmann &lattice_boltzmann);

private:

  vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_;   //< vtk writer to write ImageData

  const double dx;
  const double dy;
  const double dz{1.0};
  std::array<int, 2> size;
};
