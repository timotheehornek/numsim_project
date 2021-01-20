#include "main.h"

int main(int argc, char* argv[])
{
   //Lattice Boltzmann
   //! print detailed results (u, v, p) in debug mode
   bool detailed_results{ false };

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


   //! algorithm setup
   //! initialize lattice_boltzmann
   Lattice_boltzmann lattice_boltzmann{
           settings.nCells, settings.physicalSize,
           settings.re, settings.g};

  //! setup ouput writer
  std::array<double, 2> meshWidth = { lattice_boltzmann.dx(), lattice_boltzmann.dy() };
  OutputWriterParaviewLBM OWP{ lattice_boltzmann, meshWidth, settings.nCells };

  // declare f_even and f_odd
  std::array<int, 2> size { settings.nCells[0] , settings.nCells[1] };

  std::array<Array2D, 9> f_even{ size, size, size, size, size, size, size, size, size };
  std::array<Array2D, 9> f_odd{ size, size, size, size, size, size, size, size, size };

  //initialize f_even
  lattice_boltzmann.initialize_f(1.0, f_even);

  // set the obstacle
  lattice_boltzmann.set_obstacle_rect(settings.obstaclePos);
  //lattice_boltzmann.set_obstacle_round({ settings.nCells[0] / 4 , settings.nCells[1] / 2}, settings.nCells[1] / 8);

  //! simulation time setup
  double t{ 0.0 };

  //! set file counter
  int output_counter { 1 };

  // initialize turn variable
  int turn = 1;

  //compute max set velocity
  double max_vel = lattice_boltzmann.get_max_vel(settings.useDirichletBc,
     settings.bcBottom, settings.bcTop, settings.bcLeft, settings.bcRight);

  //set relaxation parameter
  lattice_boltzmann.set_tau(settings.physicalSize[1], max_vel);

  //! time step calculation
  lattice_boltzmann.set_dt();    //< set timestep

  while (t < settings.endTime) // <= iterate through time span
    {
                                       //< alternatively set timestep manually
      if ((t + lattice_boltzmann.dt()) > settings.endTime)       //< handle last time step
        lattice_boltzmann.set_dt(settings.endTime - t);

      t += lattice_boltzmann.dt();                               //< update current time t

      //! inform user about current time step
      DEBUG_PRINT
      (
        "=========================================================================================\n"
         << "Currently computing time: " << t << "\t\tEnd time: " << settings.endTime << '\n'
         << "Time step has been set to: dt = " << lattice_boltzmann.dt()
      );
      RELEASE_PRINT
      (
        "Simulation time: " << std::fixed << std::setprecision(2) << t << " / " << settings.endTime << '\r'
      );

      //Two lattice implementation
      if(turn) //input = f_even and output = f_odd
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

      //! write one output file every second
      if (output_counter == static_cast<int>(t))
      {
        ++output_counter;
        OWP.writeFile(t, lattice_boltzmann);
      }

      //! print variables
      if (detailed_results && t == settings.endTime)
      {
          DEBUG_PRINT(
            "Current result overview:\n"
             << "Velocity u:\n" << lattice_boltzmann.u()
             << "Velocity v:\n" << lattice_boltzmann.v()
             << "Pressure p:\n" << lattice_boltzmann.p()
             /*
             //All even values
             << "f even 0: \n" << f_even[0]
             << "f even 1: \n" << f_even[1]
             << "f even 2: \n" << f_even[2]
             << "f even 3: \n" << f_even[3]
             << "f even 4: \n" << f_even[4]
             << "f even 5: \n" << f_even[5]
             << "f even 6: \n" << f_even[6]
             << "f even 7: \n" << f_even[7]
             << "f even 8: \n" << f_even[8]

            //All odd values
             << "f odd 0: \n" << f_odd[0]
             << "f odd 1: \n" << f_odd[1]
             << "f odd 2: \n" << f_odd[2]
             << "f odd 3: \n" << f_odd[3]
             << "f odd 4: \n" << f_odd[4]
             << "f odd 5: \n" << f_odd[5]
             << "f odd 6: \n" << f_odd[6]
             << "f odd 7: \n" << f_odd[7]
             << "f odd 8: \n" << f_odd[8]
             */
          );
      }
   }
   RELEASE_PRINT('\n'); //< stop overwriting of time status line

   //! end time measurment and output execution time
   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed execution time: " << elapsed.count() << " s\n";

   return EXIT_SUCCESS;
   ///////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
   /*
    //! print detailed results (u, v, p, F, G and RHS) in debug mode
    bool detailed_results{ false };

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


    //! algorithm setup

    //! initialize discretization
    std::shared_ptr<Discretization> discretization;
    if (settings.useDonorCell)
        discretization = std::make_shared< Donor_cell>(
            settings.nCells, settings.physicalSize,
            settings.re, settings.g, settings.alpha);
    else
        discretization = std::make_shared< Central_differences>(
            settings.nCells, settings.physicalSize,
            settings.re, settings.g);

    //! setup ouput writer
    std::array<double, 2> meshWidth = { discretization->dx(), discretization->dy() };
    OutputWriterParaview OWP{ discretization , meshWidth, settings.nCells };
    OutputWriterText OWT{ discretization };

    //! initialize pressure solver
    std::shared_ptr<Pressure_solver> pressure_solver;
    if (settings.pressureSolver == "GaussSeidel")
        pressure_solver = std::make_shared< Gauss_seidel>(
            discretization->dx(), discretization->dy(),
            settings.epsilon, settings.maximumNumberOfIterations);
    else
        pressure_solver = std::make_shared< SOR>(
            discretization->dx(), discretization->dy(),
            settings.epsilon, settings.maximumNumberOfIterations, settings.omega);

    //! load boundaries for velocities u and v
    discretization->setup_bound_val_uv(
        settings.dirichletBcBottom, settings.dirichletBcTop,
        settings.dirichletBcLeft, settings.dirichletBcRight);

    //! compute and set boundary values of F and G
    discretization->compute_bound_val_FG();

    //! simulation time setup
    double t{ 0.0 };

    while (t < settings.endTime) // <= iterate through time span
    {
        //! time step calculation
        discretization->set_dt(settings.maximumDt, settings.tau);    //< set timestep
        //discretization->set_dt(.2);                                  //< alternatively set timestep manually
        if ((t + discretization->dt()) > settings.endTime)       //< handle last time step
            discretization->set_dt(settings.endTime - t);
        t += discretization->dt();                               //< update current time t
        //! inform user about current time step
        DEBUG_PRINT
        (
            "=========================================================================================\n"
            << "Currently computing time: " << t << "\t\tEnd time: " << settings.endTime << '\n'
            << "Time step has been set to: dt = " << discretization->dt()
        );
        RELEASE_PRINT
        (
            "Simulation time: " << std::fixed << std::setprecision(2) << t << " / " << settings.endTime << '\r'
        );


        //! compute and update F and G
        discretization->compute_FG();

        //! compute and update RHS
        discretization->compute_RHS();

        //! solve pressure equation and update p
        pressure_solver->solver(discretization->p_ref(), discretization->RHS());

        //! compute and update u and v
        discretization->compute_uv();



        //! compute and update boundary using boundary conditions of u and v around domain
        discretization->update_bound_val_uv(
            settings.dirichletBcBottom, settings.dirichletBcTop,
            settings.dirichletBcLeft, settings.dirichletBcRight);

        //! write results to output
        OWP.writeFile(t);


        //! print variables
        if (detailed_results&& t == settings.endTime)
        {
            DEBUG_PRINT(
                "Current result overview:\n"
                << "Velocity u:\n" << discretization->u()
                << "Velocity v:\n" << discretization->v()
                << "Pressure p:\n" << discretization->p()
                << "F:\n" << discretization->F()
                << "G:\n" << discretization->G()
                << "RHS:\n" << discretization->RHS()
            );
        }
    }
    RELEASE_PRINT('\n'); //< stop overwriting of time status line

    //! end time measurment and output execution time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed execution time: " << elapsed.count() << " s\n";

    return EXIT_SUCCESS;
    */
}
