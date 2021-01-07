#pragma once

#include "array2d/array2d.h"
#include "output_writer/output_writer.h"

#include <array>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include <memory>

/** Write *.vti files that can be viewed with ParaView.
 *  The mesh that can be visualized in ParaView corresponds to the mesh of the computational domain.
 *  All values are given for the nodes of the mesh, i.e., the corners of each cell.
 *  This means, values will be interpolated because the values are stored at positions given by the staggered grid.
 */
class OutputWriterParaview : public OutputWriter
{
public:
  //! constructor
  //! @param discretization shared pointer to the discretization object that will contain all the data to be written to the file
  OutputWriterParaview(std::shared_ptr<Discretization> discretization, std::array<int, 2> nCells, std::array<int, 2> prcs, std::array<int, 2> local_nCells);

  //! write current velocities to file, filename is output_<count>.vti
  void writeFile(double currentTime, std::vector<Array2D> &u_gathered, std::vector<Array2D> &v_gathered, std::vector<Array2D> &p_gathered);

private:
  vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_; //< vtk writer to write ImageData

  const double dx;
  const double dy;
  const double dz{1.0};
  std::array<int, 2> size;
  std::array<int, 2> prcs;
  std::array<int, 2> local_nCells; //< lower left

  enum Directions
  {
    HORIZONTAL,
    VERTICAL,
    HORIZONTAL_VERTICAL,
  };
};
