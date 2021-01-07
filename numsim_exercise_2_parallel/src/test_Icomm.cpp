/**
 *  This file contains test code for the message passer.
 * compile using: 
 * mpic++ -o  test_Icomm test_Icomm.cpp message_passing/message_passer.cpp array2d/array2d.cpp discretization/staggered_grid.cpp
 * run using:
 * mpirun -n 4 test_Icomm
 */

#include "message_passing/message_passer.h"
#include "array2d/array2d.h"
#include "discretization/staggered_grid.h"
#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{
    std::array<int, 2> nCells{2, 2}; //< get from input file, here only for demonstration purposes

    //! initialize message passer
    Message_passer_2D MP(argc, argv, nCells);

    enum Tags
    {
        TAG_U,
        TAG_V,
        TAG_F,
        TAG_G,
        TAG_P,
        TAG_RHS
    };
    assert(MP.size() == 4); //< run with FOUR processes
    enum Prcs
    {
        PRC_MASTER,
        PRC_SLAVE1,
        PRC_SLAVE2,
        PRC_SLAVE3,
    };

    enum Sides
    {
        SIDE_TOP,
        SIDE_LEFT,
        SIDE_RIGHT,
        SIDE_BOTTOM
    };
    double write_max_here;
    double local_max{static_cast<double>(MP.rank())};
    double global_max;
    MP.reduce_max(&global_max,local_max);
    if (MP.rank() == 0)
	{
        std::cout << "global_max= "<<global_max<<'\n';
    }
    if (MP.rank() == 1)
	{
        std::cout << "global_max= "<<global_max<<'\n';
    }
    MP.broadcast(&global_max);

    if (MP.rank() == 1)
	{
        std::cout << "global_max= "<<global_max<<'\n';
    }


    /*
    Staggered_grid sg1{5, 5};
    for (int i{}; i < 5; ++i)
    {
        for (int j{}; j < 5; ++j)
        {
            sg1(i, j) = 10*(j+1)+(i+1);
        }
    }
    Staggered_grid sg2{5, 5};
    for (int i{}; i < 5; ++i)
    {
        for (int j{}; j < 5; ++j)
        {
            sg2(i, j) = 100*(j+1)+10*(i+1);
        }
    }

    std::vector<double> buffer1;
    if (MP.rank() == 0)
    {
        MP.send_right(sg1.get_side(SIDE_BOTTOM), TAG_U);
        //std::cout<<"sg1:\n"<<sg1;
    }
    if (MP.rank() == 1)
    {
        buffer1.resize(3);
        MP.receive_left(buffer1, TAG_U);
    }
    MP.wait();
    if (MP.rank() == 1)
    {
        sg2.write_bound(SIDE_BOTTOM, buffer1);
        std::cout <<"sg2:\n"<< sg2;
    }
    */
}