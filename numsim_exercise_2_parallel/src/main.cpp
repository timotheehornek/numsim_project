#include "main.h"

int main(int argc, char *argv[])
{
    //! initialize MPI
    Message_passer_2D MP{argc, argv};

    //! print detailed results (u, v, p, F, G and RHS) in debug mode
    bool detailed_results{false};

    //! get input file
    std::string file_name{};
    //assert(argc == 2);
    file_name = argv[1];

    //! import settings from file
    Settings settings;
    settings.loadFromFile(file_name);

    //! display all settings and settings file name on console
    if (MP.rank() == 0)
    {
        std::cout << "Settings filename: " << file_name << '\n';
        settings.printSettings();
    }

    //! set processes in each direction
    MP.set_prcs(settings.nCells);

    //! setup for execution time measurment
    auto start = std::chrono::high_resolution_clock::now();

    //! algorithm setup

    //! initialize discretization
    std::shared_ptr<Discretization> discretization;
    if (settings.useDonorCell)
        discretization = std::make_shared<Donor_cell>(
            settings.nCells, settings.physicalSize,
            settings.re, settings.g, MP, settings.alpha);
    else
        discretization = std::make_shared<Central_differences>(
            settings.nCells, settings.physicalSize,
            settings.re, settings.g, MP);

    //! setup ouput writer
    std::vector<Array2D> u_gathered;
    std::vector<Array2D> v_gathered;
    std::vector<Array2D> p_gathered;
    std::vector<int> nCells_local_x;
    std::vector<int> nCells_local_y;
    if (MP.rank() == 0)
    {
        nCells_local_x.resize(MP.size());
        nCells_local_y.resize(MP.size());
    }
    MP.gather(discretization->nCells()[0], nCells_local_x);
    MP.gather(discretization->nCells()[1], nCells_local_y);
    if (MP.rank() == 0)
    {
        u_gathered.reserve(MP.size());
        v_gathered.reserve(MP.size());
        p_gathered.reserve(MP.size());
        for (int i{}; i < MP.size(); ++i)
        {
            u_gathered.push_back({nCells_local_x[i] + 2, nCells_local_y[i] + 2});
            v_gathered.push_back({nCells_local_x[i] + 2, nCells_local_y[i] + 2});
            p_gathered.push_back({nCells_local_x[i] + 2, nCells_local_y[i] + 2});
        }
    }
    int output_counter{1};
    //OutputWriterText OWT{ discretization };

    //! initialize pressure solver
    std::shared_ptr<Pressure_solver> pressure_solver;
    /*
    if (settings.pressureSolver == "GaussSeidel")
        pressure_solver = std::make_shared<Gauss_seidel>(
            discretization->dx(), discretization->dy(),
            settings.epsilon, settings.maximumNumberOfIterations);
    else*/
    pressure_solver = std::make_shared<SOR_parallel>(
        *discretization, MP,
        settings.epsilon, settings.maximumNumberOfIterations, settings.omega);

    //! simulation time setup
    double t{0.0};

    //! compute and update boundary using boundary conditions of u and v around domain
    discretization->communicate_update_bound_val_uv(
        settings.dirichletBcBottom, settings.dirichletBcTop,
        settings.dirichletBcLeft, settings.dirichletBcRight);

    while (t < settings.endTime) // <= iterate through time span
    {
        //! time step calculation
        discretization->set_dt(settings.maximumDt, settings.tau); //< set timestep
        //discretization->set_dt(.2);                        //< alternatively set timestep manually
        if ((t + discretization->dt()) > settings.endTime) //< handle last time step
        {
            const double rem_time{settings.endTime - t};
            if (rem_time<1e-10)
                break;
            discretization->set_dt(rem_time);
        }
            
        t += discretization->dt(); //< update current time t
        if (MP.rank() == 0)
        {
            //! inform user about current time step
            DEBUG_PRINT(
                "=========================================================================================\n"
                << "Currently computing time: " << t << "\t\tEnd time: " << settings.endTime << '\n'
                << "Time step has been set to: dt = " << discretization->dt());
            RELEASE_PRINT(
                //"Simulation time: " << std::fixed << std::setprecision(2) << t << " / " << settings.endTime << '\r');
                "Simulation time: " << std::fixed << std::setprecision(2) << t << " / " << settings.endTime << '\n');
        }

        //! compute and update F and G
        discretization->compute_FG();

        //! communicate F and G required to compute RHS or set boundary values of F and G
        discretization->communicate_update_bound_val_FG();

        //! compute and update RHS
        discretization->compute_RHS();

        //! solve pressure equation and update p
        pressure_solver->solver(discretization->p_ref(), discretization->RHS());

        //! compute and update u and v
        discretization->compute_uv();

        //! compute and update boundary using boundary conditions of u and v around domain
        discretization->communicate_update_bound_val_uv(
            settings.dirichletBcBottom, settings.dirichletBcTop,
            settings.dirichletBcLeft, settings.dirichletBcRight);

        //! write one output file every second
        if (output_counter == static_cast<int>(t))
        {
            ++output_counter;
            //! send results to master (process 0)
            MP.send(discretization->u(), TAG_U, 0);
            MP.send(discretization->v(), TAG_V, 0);
            MP.send(discretization->p(), TAG_P, 0);
            if (MP.rank() == 0)
            {
                for (int i{0}; i < MP.size(); ++i)
                {
                    //! receive and write results
                    MP.receive(u_gathered[i], TAG_U, i);
                    MP.receive(v_gathered[i], TAG_V, i);
                    MP.receive(p_gathered[i], TAG_P, i);
                }
            }
            MP.wait();
            if (MP.rank() == 0)
            {
                static std::array<int, 2> prcs{MP.prcs(0), MP.prcs(1)};
                static std::array<int, 2> nCells_local_0{nCells_local_x[0], nCells_local_y[0]};
                static OutputWriterParaview OWP{discretization, settings.nCells, prcs, nCells_local_0};
                //! write results to output
                OWP.writeFile(t, u_gathered, v_gathered, p_gathered);
                //OWT.writeFile(t);
            }
        }

        //! print variables
        if (detailed_results)
        {
            DEBUG_PRINT(
                "Current result overview:\n"
                << "Velocity u:\n"
                << discretization->u()
                << "Velocity v:\n"
                << discretization->v()
                << "Pressure p:\n"
                << discretization->p()
                << "F:\n"
                << discretization->F()
                << "G:\n"
                << discretization->G()
                << "RHS:\n"
                << discretization->RHS());
        }
    }
    if (MP.rank() == 0)
    {
        RELEASE_PRINT('\n'); //< stop overwriting of time status line
    }

    //! end time measurment and output execution time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed execution time: " << elapsed.count() << " s\n";

    return EXIT_SUCCESS;
}