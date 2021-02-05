#include "main.h"

int main(int argc, char *argv[])
{
  //! print detailed results (u, v, p, F, G and RHS) in debug mode
  bool detailed_results{true};

  //! specify and output input file (as argument on Linux and here in Windows)
  std::string file_name{};
  assert(argc == 2);
  file_name = argv[1];
  std::cout << "Settings filename: " << file_name << '\n';

  //! import settings from file
  Settings settings;
  settings.loadFromFile(file_name);

  //! display all settings on console
  settings.printSettings();

  //! setup for execution time measurment
  auto start = std::chrono::high_resolution_clock::now();

  //! run simulation (Navier-Stokes or Lattice-Boltzmann)
  if (settings.navierStokes)
    run_ns(settings, detailed_results);
  else
    run_lbm(settings, detailed_results);

  RELEASE_PRINT('\n'); //< stop overwriting of time status line

  //! end time measurment and output execution time
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed execution time: " << elapsed.count() << " s\n";

  return EXIT_SUCCESS;
}
void run_ns(Settings settings, const bool detailed_results)
{
  //! algorithm setup
  //! initialize discretization
  std::shared_ptr<Discretization> discretization;
  if (settings.useDonorCell)
    discretization = std::make_shared<Donor_cell>(
        settings.nCells, settings.physicalSize, settings.obstExist, settings.obstaclePos,
        settings.re, settings.g, settings.alpha);
  else
    discretization = std::make_shared<Central_differences>(
        settings.nCells, settings.physicalSize, settings.obstExist, settings.obstaclePos,
        settings.re, settings.g);

  //! setup ouput writer
  std::array<double, 2> meshWidth = {discretization->dx(), discretization->dy()};
  OutputWriterParaview OWP{discretization, meshWidth, settings.nCells};
  //OutputWriterText OWT{ discretization };
  int output_counter{0}; //< counter to keep track of output times

  //! initialize pressure solver
  std::shared_ptr<Pressure_solver> pressure_solver;
  //! setup pressure sover implementation
  if (settings.pressureSolver == "SOR")
    pressure_solver = std::make_shared<SOR>(
        settings.epsilon, settings.maximumNumberOfIterations, settings.omega);
  else
  {
    pressure_solver = std::make_shared<CG>(
        settings.epsilon, settings.maximumNumberOfIterations, settings.nCells);
    //! set boundary conditions in p
    pressure_solver->update_boundaries(*discretization,
                                       settings.bcBottom, settings.bcTop,
                                       settings.bcLeft, settings.bcRight,
                                       settings.useDirichletBc);
  }

  //! load boundaries for velocities u and v
  discretization->setup_bound_val_uv(
      settings.bcBottom, settings.bcTop,
      settings.bcLeft, settings.bcRight,
      settings.useDirichletBc);

  //! compute and set boundary values of F and G
  discretization->setup_bound_val_FG();

  //! simulation time setup
  double t{0.0};

  while (t < settings.endTime) // <= iterate through time span
                               //while (t==0.0)
  {

    //! compute and update boundary using boundary conditions of u and v around domain
    discretization->update_bound_val_uv(
        settings.bcBottom, settings.bcTop,
        settings.bcLeft, settings.bcRight,
        settings.useDirichletBc);

    //! time step calculation
    discretization->set_dt(settings.maximumDt, settings.tau); //< set timestep
    //discretization->set_dt(.2);                               //< alternatively set timestep manually
    if ((t + discretization->dt()) > settings.endTime) //< handle last time step
      discretization->set_dt(settings.endTime - t);
    t += discretization->dt(); //< update current time t
    //! inform user about current time step
    DEBUG_PRINT(
        "=========================================================================================\n"
        << "Currently computing time: " << t << "\t\tEnd time: " << settings.endTime << '\n'
        << "Time step has been set to: dt = " << discretization->dt());
    RELEASE_PRINT(
        "Simulation time: " << std::fixed << std::setprecision(2) << t << " / " << settings.endTime << '\r');

    //! compute and update F and G
    discretization->compute_FG();

    discretization->update_bound_val_FG(settings.useDirichletBc);

    //! compute and update RHS
    discretization->compute_RHS();

    //! solve pressure equation and update p
    //pressure_solver->solver(discretization->p_ref(), discretization->RHS());
    pressure_solver->solver(*discretization,
                            settings.bcBottom, settings.bcTop,
                            settings.bcLeft, settings.bcRight,
                            settings.useDirichletBc);

    //! compute and update u and v
    discretization->compute_uv();

    //! compute u and v around boundary
    if (settings.obstExist)
    {
      discretization->compute_bound_val_obstacle();
    }
    //! write results to output every 1/10 s
    //if(output_counter==static_cast<int>(10*t))
    if (output_counter == static_cast<int>(t))
    {
      assert(settings.maximumDt <= .1);
      OWP.writeFile(t);
      //OWT.writeFile(t);
      ++output_counter;
    }

    //! print variables
    if (detailed_results && t == settings.endTime)
    //if(true)
    {
      DEBUG_PRINT(
          "Current result overview:\n"
          << "Velocity u:\n"
          << discretization->u()
          //<< "Velocity v:\n"
          //<< discretization->v()
          << "Pressure p:\n"
          << discretization->p()
          << "F:\n"
          << discretization->F()
          //<< "G:\n"
          //<< discretization->G()
          << "RHS:\n"
          << discretization->RHS());
    }
  }
}

void run_lbm(Settings settings, const bool detailed_results)
{
  //! algorithm setup
  //! initialize lattice_boltzmann
  Lattice_boltzmann lattice_boltzmann{
      settings.nCells, settings.physicalSize,
      settings.re, settings.magicFactor, settings.g};

  //! setup ouput writer
  std::array<double, 2> meshWidth = {lattice_boltzmann.dx(), lattice_boltzmann.dy()};
  OutputWriterParaviewLBM OWP{lattice_boltzmann, meshWidth, settings.nCells};

  // declare f_even and f_odd
  std::array<int, 2> size{settings.nCells[0], settings.nCells[1]};

  std::array<Array2D, 9> f_even{size, size, size, size, size, size, size, size, size};
  std::array<Array2D, 9> f_odd{size, size, size, size, size, size, size, size, size};

  //initialize f_even
  lattice_boltzmann.initialize_f(1.0, f_even);

  // set the obstacle
  if (settings.obstExist)
  {
    lattice_boltzmann.set_obstacle_rect(settings.obstaclePos);
    //lattice_boltzmann.set_obstacle_round({ settings.nCells[0] / 4 , settings.nCells[1] / 2}, settings.nCells[1] / 8);
  }

  //! simulation time setup
  double t{0.0};

  //! set file counter
  int output_counter{1};

  // initialize turn variable
  int turn = 1;

  //compute viscosity
  double viscosity = lattice_boltzmann.compute_viscosity(settings.physicalSize, settings.useDirichletBc,
                                                         settings.bcBottom, settings.bcTop, settings.bcLeft, settings.bcRight);
  //set relaxation parameter
  lattice_boltzmann.set_tau(viscosity);

  //! time step calculation
  lattice_boltzmann.set_dt(); //< set timestep

  while (t < settings.endTime) // <= iterate through time span
  {
    t += lattice_boltzmann.dt(); //< update current time t

    //! inform user about current time step
    DEBUG_PRINT(
        "=========================================================================================\n"
        << "Currently computing time: " << t << "\t\tEnd time: " << settings.endTime << '\n'
        << "Time step has been set to: dt = " << lattice_boltzmann.dt());
    RELEASE_PRINT(
        "Simulation time: " << std::fixed << std::setprecision(4) << t << " / " << settings.endTime << '\r');

    //Two lattice implementation
    if (turn) //input = f_even and output = f_odd
    {
      //do collide - stream for all boundary nodes
      lattice_boltzmann.boundary_treatment(f_even, f_odd, settings.useDirichletBc,
                                           settings.bcBottom, settings.bcTop, settings.bcLeft, settings.bcRight);

      //do collide - stream all inner nodes
      lattice_boltzmann.two_lattice(f_even, f_odd);

      turn = 0;
    }
    else //input = f_odd and output = f_even
    {
      //do collide - stream for all boundary nodes
      lattice_boltzmann.boundary_treatment(f_odd, f_even, settings.useDirichletBc,
                                           settings.bcBottom, settings.bcTop, settings.bcLeft, settings.bcRight);

      //do collide - stream all inner nodes
      lattice_boltzmann.two_lattice(f_odd, f_even);

      turn = 1;
    }

    //! write one output file every seconds
    if (output_counter == static_cast<int>(t))
    {
      ++output_counter;
      OWP.writeFile(t, lattice_boltzmann);
    }

    //! print variables
    if (detailed_results && t == settings.endTime)
    {
      /*
        std::cout << "left boundary:" << '\n';
        for (size_t j = 0; j < size[1]; j++)
        {
          std::cout << "j = " << j << " : " << lattice_boltzmann.u(0,j) << '\n';
        }
        std::cout << "right boundary:" << '\n';
        for (size_t j = 0; j < size[1]; j++)
        {
          std::cout << "j = " << j << " : " << lattice_boltzmann.u(size[0] - 1,j) << '\n';
        }
        */
      DEBUG_PRINT(
          "Current result overview:\n"
          << "Velocity u:\n"
          << lattice_boltzmann.u()
          << "Velocity v:\n"
          << lattice_boltzmann.v()
          << "Pressure p:\n"
          << lattice_boltzmann.p());
    }
  }
}