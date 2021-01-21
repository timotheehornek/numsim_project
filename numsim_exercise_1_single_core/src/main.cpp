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

    //! algorithm setup

    //! initialize discretization
    std::shared_ptr<Discretization> discretization;
    if (settings.useDonorCell)
        discretization = std::make_shared<Donor_cell>(
            settings.nCells, settings.physicalSize, settings.obstaclePos,
            settings.re, settings.g, settings.alpha);
    else
        discretization = std::make_shared<Central_differences>(
            settings.nCells, settings.physicalSize, settings.obstaclePos,
            settings.re, settings.g);

    //! setup ouput writer
    std::array<double, 2> meshWidth = {discretization->dx(), discretization->dy()};
    OutputWriterParaview OWP{ discretization , meshWidth, settings.nCells };
    //OutputWriterText OWT{ discretization };

    //! initialize pressure solver
	std::array<bool,4> p_0_boundary{
		settings.useDirichletBc[0]==false||settings.bcBottom[1]!=0,
		settings.useDirichletBc[1]==false||settings.bcTop[1]!=0,
		settings.useDirichletBc[2]==false||settings.bcLeft[0]!=0,
		settings.useDirichletBc[3]==false||settings.bcRight[0]!=0,
	}; //< set pressure boundary conditions (p=0 if true; Neumann otherwise)
    std::shared_ptr<Pressure_solver> pressure_solver;
	/*
    if (settings.pressureSolver == "GaussSeidel")
        pressure_solver = std::make_shared<Gauss_seidel>(
            discretization->dx(), discretization->dy(),
            settings.epsilon, settings.maximumNumberOfIterations);
    else
        pressure_solver = std::make_shared<SOR>(
            discretization->dx(), discretization->dy(),
            settings.epsilon, settings.maximumNumberOfIterations, settings.omega);*/
	pressure_solver = std::make_shared<SOR>(
            settings.epsilon, settings.maximumNumberOfIterations, settings.omega, p_0_boundary);

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

        discretization->update_bound_val_FG();

        //! compute and update RHS
        discretization->compute_RHS();

        //! solve pressure equation and update p
        //pressure_solver->solver(discretization->p_ref(), discretization->RHS());
		pressure_solver->solver(*discretization);

        //! compute and update u and v
        discretization->compute_uv();
		
		//! compute u and v around boundary
		discretization->compute_bound_val_obstacle();

        //! write results to output
        OWP.writeFile(t);
        //OWT.writeFile(t);

        //! print variables
        if (detailed_results&&false)
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
                << discretization->RHS()
                );
        }
    }
    RELEASE_PRINT('\n'); //< stop overwriting of time status line

    //! end time measurment and output execution time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed execution time: " << elapsed.count() << " s\n";

    return EXIT_SUCCESS;
}