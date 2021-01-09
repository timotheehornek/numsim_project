#pragma once

#ifndef NDEBUG
#define DEBUG_PRINT(x) do { std::cout << x << '\n'; } while (0)
#define RELEASE_PRINT(x)
#else
#define DEBUG_PRINT(x)
#define RELEASE_PRINT(x) do { std::cout << x ; } while (0)
#endif

#include "output_writer/output_writer_paraview.h"
//#include "output_writer/output_writer_text.h"
#include "array2d/array2d.h"
#include "discretization/central_differences.h"
#include "discretization/donor_cell.h"
//#include "pressure_solver/gauss_seidel.h"
#include "pressure_solver/pressure_solver.h"
#include "settings_parser/settings.h"
#include "pressure_solver/SOR.h"
#include "discretization/staggered_grid.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <memory>