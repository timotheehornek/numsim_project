#pragma once

#include "../array2d/array2d.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <vector>

class Message_passer_2D
{
private:
    //! some MPI variables
    int m_rank;
    int m_size;

    //! number of processes used for each direction, m_size=m_prcs_x*m_prcs_y
    int m_prcs_x;
    int m_prcs_y;

    //! number of cells in 2D discretization
    //const std::array<int, 2> m_nCells;

    //! variables for requests
    int m_count_reqs{0};
    //MPI_Request m_reqs[48];
    //MPI_Request m_reqs[100];
    std::vector<MPI_Request> m_reqs;

public:
    Message_passer_2D(int argc, char *argv[]);
    ~Message_passer_2D();

    //! getter for size and rank of process
    int size() const;
    int rank() const;

    //! getter for processes in x and y direction
    //! (direction: 0=x, 1=y)
    const int prcs(int dir) const;

    //! set number of processes in two dimensions to minimize communication
    void set_prcs(std::array<int, 2> nCells);

    //! send and receive methods for vectors
    void send_left(const std::vector<double> &vec, const int tag);
    void send_right(const std::vector<double> &vec, const int tag);
    void send_top(const std::vector<double> &vec, const int tag);
    void send_bottom(const std::vector<double> &vec, const int tag);

    void receive_left(std::vector<double> &vec, const int tag);
    void receive_right(std::vector<double> &vec, const int tag);
    void receive_top(std::vector<double> &vec, const int tag);
    void receive_bottom(std::vector<double> &vec, const int tag);

    void send_bottom_right(const double val, const int tag);
    void send_top_left(const double val, const int tag);
    void receive_top_left(double &val, const int tag);
    void receive_bottom_right(double &val, const int tag);

    //! send and receive for 2D arrays and doubles
    void send(const Array2D &arr, const int tag, const int dest);
    void receive(Array2D &arr, const int tag, const int src);
    void send(const double val, const int tag, const int dest);
    void receive(double &val, const int tag, const int src);

    //! reduce to maximum
    void reduce_max(double *result, double local_max, int root = 0) const;

    //! broadcast data
    void broadcast(double *data, int root = 0) const;
    void broadcast(int *data, int root = 0) const;

    //! gather data
    void gather(int val, std::vector<int> &vec, int root = 0);
    void gather(double val, std::vector<double> &vec, int root = 0);

    //! wait for asynchronous communication requests to finish
    void wait();

    //! wait for code in all processes to reach barrier
    void barrier();
};