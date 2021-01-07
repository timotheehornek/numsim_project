/**
 *  This file contains test code for the message passer.
 * compile using: 
 * mpic++ -o  test_main test_main.cpp message_passing/message_passer.cpp array2d/array2d.cpp
 * run using:
 * mpirun -n 4 test_main
 */

#include "message_passing/message_passer.h"
#include <iostream>

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

    //! vector for received data
    std::vector<double> vec_received;
    vec_received.resize(4);

    //! 2d array for received data
    Array2D arr_received{4,4};

    //!###########################################!//
    //! how to send a vector to neighboring cell?
    if (MP.rank() == PRC_MASTER)
    {
        //! some information about initialization
        std::cout << "size = " << MP.size() << '\n'
                  << "prcs_x = " << MP.prcs(0) << "; prcs_y = " << MP.prcs(1) << ";\n";

        //! send a vector to the right and left
        std::vector<double> to_send_right{1.2, 3.4, 5.6, 7.8};
        MP.send_right(to_send_right, TAG_U);

        std::vector<double> to_send_up{-1.2, -3.4, -5.6, -7.8};
        MP.send_top(to_send_up, TAG_U);
    }

    if (MP.rank() == PRC_SLAVE1)
    {
        //! receive vector sent by master
        MP.receive_left(vec_received, TAG_U);

        //! print (not yet received?) vector entries:
        for (int i{0}; i < vec_received.size(); ++i)
        {
            std::cout << "Rank " << MP.rank() << " received vector entry (before wait): " << vec_received[i] << '\n';
        }
    }

    if (MP.rank() == PRC_SLAVE2)
    {
        //! receive vector sent by master
        MP.receive_bottom(vec_received, TAG_U);

        //! print (not yet received?) vector entries:
        for (int i{0}; i < vec_received.size(); ++i)
        {
            std::cout << "Rank " << MP.rank() << " received vector entry (before wait): " << vec_received[i] << '\n';
        }
    }

    //! wait for all message passing requests to finish
    MP.wait();

    if (MP.rank() == PRC_SLAVE1)
    {
        //! print received vector entries:
        for (int i{0}; i < vec_received.size(); ++i)
        {
            std::cout << "Rank " << MP.rank() << " received vector entry: " << vec_received[i] << '\n';
        }
    }

    if (MP.rank() == PRC_SLAVE2)
    {
        //! print received vector entries:
        for (int i{0}; i < vec_received.size(); ++i)
        {
            std::cout << "Rank " << MP.rank() << " received vector entry: " << vec_received[i] << '\n';
        }
    }

    //! wait for all processes to reach this line
    MP.barrier();

    //!###########################################!//
    //! how to send Array2D elements?
    if (MP.rank() == PRC_SLAVE1)
    {
        //! create array and initialize diagonal elements
        Array2D to_send{4, 4};
        for (int i{0}; i < to_send.size()[0];++i)
            to_send(i, i) = i;
        //! send array to master
        MP.send(to_send,TAG_U,PRC_MASTER);
    }

    if (MP.rank() == PRC_MASTER)
    {
        //! receive array
        MP.receive(arr_received,TAG_U,PRC_SLAVE1);
    }

    //! wait for all message passing requests to finish
    MP.wait();

    if (MP.rank() == PRC_MASTER)
    {
        //! receive array
        //MP.receive_2d_array(arr_received,TAG_U,PRC_SLAVE1);
        
        //! print received array:
        std::cout<<arr_received;
    }
    return 0;
}