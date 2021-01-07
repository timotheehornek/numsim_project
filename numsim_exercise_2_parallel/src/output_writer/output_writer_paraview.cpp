#include "output_writer/output_writer_paraview.h"

#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

OutputWriterParaview::OutputWriterParaview(std::shared_ptr<Discretization> discretization, std::array<int, 2> nCells, std::array<int, 2> prcs, std::array<int, 2> local_nCells)
    : OutputWriter(discretization), dx{ discretization->dx() }, dy{ discretization->dy() }, size{ nCells[0] + 2,nCells[1] + 2 }, prcs{prcs}, local_nCells{local_nCells}
{
  // Create a vtkWriter_
  vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
}

void OutputWriterParaview::writeFile(double currentTime, std::vector<Array2D>& u_gathered, std::vector<Array2D>& v_gathered, std::vector<Array2D>& p_gathered)
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
  dataSet->SetDimensions(size[0] - 1, size[1] - 1, 1);  // we want to have points at each corner of each cell
  
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
  for (int j = 0; j < size[1] - 1; j++)
  {
    for (int i = 0; i < size[0] - 1; i++, index++)
    {
      //! compute process in x direction
      int p_x{i/local_nCells[0]};       //< first guess for process (valid unless on top or right)
      int i_local{i%local_nCells[0]};   //< first guess for index in local array (valid unless on top or right)
      if (p_x==prcs[0])
      {
        //! correct process and index if on top or right
        --p_x;
        i_local+=local_nCells[0];
      }

      //! compute process in y direction
      int p_y{j/local_nCells[1]};       //< first guess for process (valid unless on top or right)
      int j_local{j%local_nCells[1]};   //< first guess for index in local array (valid unless on top or right)
      if (p_y==prcs[1])
      {
        //! correct process and index if on top or right
        --p_y;
        j_local+=local_nCells[1];
      }
      //! compute process number
      int p{p_x+prcs[0]*p_y};

      //! handle inner cell (0,0) index call (load values from process located to the lower left)
      if((j_local==0 && i_local==0) && (p_x>0 && p_y>0))
      {
        const int lower_left_prc{p_x - 1 + prcs[0]*(p_y - 1)};
        u_gathered[p](0,0)=u_gathered[lower_left_prc](local_nCells[0]-2, local_nCells[1]-2);
        v_gathered[p](0,0)=v_gathered[lower_left_prc](local_nCells[0]-2, local_nCells[1]-2);
      }
      
      arrayPressure->SetValue(index, p_gathered[p].interpolate(i_local, j_local, HORIZONTAL_VERTICAL));

      std::array<double,3> velocityVector;

      velocityVector[0] = u_gathered[p].interpolate(i_local, j_local, VERTICAL);
      velocityVector[1] = v_gathered[p].interpolate(i_local, j_local, HORIZONTAL);

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
  
  vtkWriter_->SetDataModeToAscii();     // comment this in to get ascii text files: those can be checked in an editor
  //vtkWriter_->SetDataModeToBinary();      // set file mode to binary files: smaller file sizes

  // finally write out the data
  vtkWriter_->Write();
}