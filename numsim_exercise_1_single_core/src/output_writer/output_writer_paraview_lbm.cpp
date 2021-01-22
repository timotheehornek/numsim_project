#include "output_writer/output_writer_paraview_lbm.h"

#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

OutputWriterParaviewLBM::OutputWriterParaviewLBM(Lattice_boltzmann lattice_boltzmann, std::array<double, 2> meshwidth, std::array<int, 2> nCells)
    : OutputWriterLBM(lattice_boltzmann), dx{ meshwidth[0] }, dy{ meshwidth[1] }, size{ nCells[0], nCells[1] }
{
  // Create a vtkWriter_
  vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
}

void OutputWriterParaviewLBM::writeFile(double currentTime, Lattice_boltzmann lattice_boltzmann)
{
  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << setfill('0') << fileNo_ << "." << vtkWriter_->GetDefaultFileExtension();

  // increment file no.
  fileNo_++;

  // assign the new file name to the output vtkWriter_
  vtkWriter_->SetFileName(fileName.str().c_str());

  // initialize data set that will be output to the file
  vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
  dataSet->SetOrigin(0, 0, 0);

  // set spacing of mesh
  dataSet->SetSpacing(dx, dy, dz);

  // set number of points in each dimension, 1 cell in z direction
  dataSet->SetDimensions(size[0], size[1], 1);

  // add pressure field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayPressure = vtkDoubleArray::New();

  // the pressure is a scalar which means the number of components is 1
  arrayPressure->SetNumberOfComponents(1);

  // Set the number of pressure values and allocate memory for it. We already know the number, it has to be the same as there are nodes in the mesh.
  arrayPressure->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  arrayPressure->SetName("pressure");


  // add velocity field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayVelocity = vtkDoubleArray::New();

  // here we have two components (u,v), but ParaView will only allow vector glyphs if we have an ℝ^3 vector,
  // therefore we use a 3-dimensional vector and set the 3rd component to zero
  arrayVelocity->SetNumberOfComponents(3);

  // set the number of values
  arrayVelocity->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  arrayVelocity->SetName("velocity");

  // loop over the nodes of the mesh and assign the interpolated p values in the vtk data structure
  // we only consider the cells that are the actual computational domain, not the helper values in the "halo"

  int index = 0;   // index for the vtk data structure, will be incremented in the inner loop
  //for (int j = 0; j < nCells[1]+1; j++)
  for (int j = 0; j < size[1]; j++)
  {
    //for (int i = 0; i < nCells[0]+1; i++, index++)
    for (int i = 0; i < size[0]; i++, index++)
    {
      arrayPressure->SetValue(index, lattice_boltzmann.p(i,j));

      std::array<double,3> velocityVector;
      velocityVector[0] = lattice_boltzmann.u(i,j);
      velocityVector[1] = lattice_boltzmann.v(i,j);

      velocityVector[2] = 0.0;    // z-direction is 0

      arrayVelocity->SetTuple(index, velocityVector.data());
    }
  }

  // now, we should have added as many values as there are points in the vtk data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arrayPressure);

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arrayVelocity);

  // add current time
  vtkSmartPointer<vtkDoubleArray> arrayTime = vtkDoubleArray::New();
  arrayTime->SetName("TIME");
  arrayTime->SetNumberOfTuples(1);
  arrayTime->SetTuple1(0, currentTime);
  dataSet->GetFieldData()->AddArray(arrayTime);

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
  vtkWriter_->SetInputData(dataSet);

  //vtkWriter_->SetDataModeToAscii();     // comment this in to get ascii text files: those can be checked in an editor
  vtkWriter_->SetDataModeToBinary();      // set file mode to binary files: smaller file sizes

  // finally write out the data
  vtkWriter_->Write();
}
