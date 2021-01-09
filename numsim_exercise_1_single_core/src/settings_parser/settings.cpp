#include "settings.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

const bool Settings::extract_bool(std::string &str) const
{
	bool result{};
	//remove white spaces from value string
	str.erase(
		std::remove(str.begin(), str.end(), ' '), str.end());
	str.erase(
		std::remove(str.begin(), str.end(), '\r'), str.end());
	//load value (true or false)
	if (str == "true")
		return true;
	else if (str == "false")
		return false;
	else
	{
		std::cout << "Unable to read bool from input. Return true.\n" << "Input was: "	<< str << ".\n";
		return true;
	}
}

void Settings::loadFromFile(const std::string &filename)
{
	//open file with settings
	std::ifstream settings{filename};
	if (!settings.is_open())
	{
		std::cout << "Unable to open file: " << filename << '\n';
	}
	else
	{
		//read settings line by line
		std::string line{};
		std::string parameter_name{};
		std::string parameter_value_string{};

		//equal sign position; used to find name and value
		size_t eq_sign_pos{};

		while (std::getline(settings, line))
		{
			//check if line is comment or empty or hast no '=' sign
			if ((line[0] != '#') && (line.find_first_of("= \t") != std::string::npos))
			{
				//extract parameter name
				parameter_name = line.substr(0, line.find_first_of('='));
				//remove white spaces from name
				parameter_name.erase(
					std::remove(parameter_name.begin(), parameter_name.end(), ' '),
					parameter_name.end());

				//extract parameter value
				eq_sign_pos = line.find_first_of('=');
				parameter_value_string = line.substr(eq_sign_pos + 1,
													 line.find_first_of("#\t\n") - eq_sign_pos - 1);

				//write parameter values in order of occurence in header file
				if (parameter_name == "nCellsX")
					nCells[0] = std::stoi(parameter_value_string);
				else if (parameter_name == "nCellsY")
					nCells[1] = std::stoi(parameter_value_string);
				else if (parameter_name == "physicalSizeX")
					physicalSize[0] = std::stod(parameter_value_string);
				else if (parameter_name == "physicalSizeY")
					physicalSize[1] = std::stod(parameter_value_string);
				else if (parameter_name == "re")
					re = std::stod(parameter_value_string);
				else if (parameter_name == "endTime")
					endTime = std::stod(parameter_value_string);
				else if (parameter_name == "tau")
					tau = std::stod(parameter_value_string);
				else if (parameter_name == "maximumDt")
					maximumDt = std::stod(parameter_value_string);
				else if (parameter_name == "gX")
					g[0] = std::stod(parameter_value_string);
				else if (parameter_name == "gY")
					g[1] = std::stod(parameter_value_string);
				else if (parameter_name == "useDonorCell")
					useDonorCell = extract_bool(parameter_value_string);
				else if (parameter_name == "alpha")
					alpha = std::stod(parameter_value_string);
				else if (parameter_name == "useDirichletBcBottom")
					useDirichletBc[0] = extract_bool(parameter_value_string);
				else if (parameter_name == "useDirichletBcTop")
					useDirichletBc[1] = extract_bool(parameter_value_string);
				else if (parameter_name == "useDirichletBcLeft")
					useDirichletBc[2] = extract_bool(parameter_value_string);
				else if (parameter_name == "useDirichletBcRight")
					useDirichletBc[3] = extract_bool(parameter_value_string);
				else if (parameter_name == "bcBottomX")
					bcBottom[0] = std::stod(parameter_value_string);
				else if (parameter_name == "bcBottomY")
					bcBottom[1] = std::stod(parameter_value_string);
				else if (parameter_name == "bcTopX")
					bcTop[0] = std::stod(parameter_value_string);
				else if (parameter_name == "bcTopY")
					bcTop[1] = std::stod(parameter_value_string);
				else if (parameter_name == "bcLeftX")
					bcLeft[0] = std::stod(parameter_value_string);
				else if (parameter_name == "bcLeftY")
					bcLeft[1] = std::stod(parameter_value_string);
				else if (parameter_name == "bcRightX")
					bcRight[0] = std::stod(parameter_value_string);
				else if (parameter_name == "bcRightY")
					bcRight[1] = std::stod(parameter_value_string);
				else if (parameter_name == "obstLoLeX")
					obstaclePos[0] = std::stoi(parameter_value_string);
				else if (parameter_name == "obstLoLeY")
					obstaclePos[1] = std::stoi(parameter_value_string);
				else if (parameter_name == "obstUpRiX")
					obstaclePos[2] = std::stoi(parameter_value_string);
				else if (parameter_name == "obstUpRiY")
					obstaclePos[3] = std::stoi(parameter_value_string);
				else if (parameter_name == "pressureSolver")
				{
					//remove white spaces from value string
					parameter_value_string.erase(
						std::remove(parameter_value_string.begin(), parameter_value_string.end(), ' '),
						parameter_value_string.end());
					//load value (SOR or CG)
					if (parameter_value_string == "SOR")
						pressureSolver = "SOR";
					else if (parameter_value_string == "CG")
						pressureSolver = "CG";
					else if (parameter_value_string == "GaussSeidel")
						pressureSolver = "GaussSeidel";
					else
						std::cout << "Unable to read \"pressureSolver\".\n";
				}
				else if (parameter_name == "omega")
					omega = std::stod(parameter_value_string);
				else if (parameter_name == "epsilon")
					epsilon = std::stod(parameter_value_string);
				else if (parameter_name == "maximumNumberOfIterations")
					maximumNumberOfIterations = static_cast<int>(std::stod(parameter_value_string));
			}
		}
		settings.close();
	}
}

void Settings::printSettings()
{
	std::cout << "Settings: " << std::endl
			  << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
			  << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt << std::endl
			  << "  useDirichletBc: bottom: " << useDirichletBc[0] << ", top: " << useDirichletBc[1] << ", left: " << useDirichletBc[2] << ", right: " << useDirichletBc[3] << std::endl
			  << "  dirichletBC: bottom: (" << bcBottom[0] << "," << bcBottom[1] << ")"
			  << ", top: (" << bcTop[0] << "," << bcTop[1] << ")"
			  << ", left: (" << bcLeft[0] << "," << bcLeft[1] << ")"
			  << ", right: (" << bcRight[0] << "," << bcRight[1] << ")" << std::endl
			  << "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl
			  << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}