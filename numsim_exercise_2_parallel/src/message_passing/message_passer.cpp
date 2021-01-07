#include "message_passer.h"

Message_passer_2D::Message_passer_2D(int argc, char *argv[])
{
    //! initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_size);
}

Message_passer_2D::~Message_passer_2D()
{
    MPI_Finalize();
}

//! getter
int Message_passer_2D::size() const
{
    return m_size;
}
int Message_passer_2D::rank() const
{
    return m_rank;
}
const int Message_passer_2D::prcs(int dir) const
{
    assert(0 <= dir && dir <= 1);
    if (dir == 0)
    {
        return m_prcs_x;
    }
    else
    {
        return m_prcs_y;
    }
}

void Message_passer_2D::set_prcs(std::array<int, 2> nCells)
{
    //! compute preliminary optimal rank_y
    double prcs_y_opt{std::pow(m_size * nCells[1] / nCells[0], .5)};

    //! find closest divider of m_size to preliminary optimal rank
    double dist{prcs_y_opt - 1.0};
    double dist_new;
    //! initialize ranks in case divider 1 is the best
    m_prcs_y = 1;
    m_prcs_x = m_size;
    for (int i{2}; i <= m_size; ++i)
    {
        //! check whether i is divider of m_number
        if (m_size % i == 0)
        {
            dist_new = std::abs(i - prcs_y_opt);
            if (dist_new < dist)
            {
                //! update ranks and dist for each dimension
                dist = dist_new;
                m_prcs_y = i;
                m_prcs_x = m_size / m_prcs_y;
            }
        }
    }
}

//! send methods for vectors
void Message_passer_2D::send_left(const std::vector<double> &vec, const int tag)
{
    assert(m_rank % m_prcs_x > 0);

    //m_reqs.emplace_back();
    const int dest{m_rank - 1};
    //MPI_Isend(&vec[0], vec.size(), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &m_reqs.back());

    MPI_Send(&vec[0], vec.size(), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
}

void Message_passer_2D::send_right(const std::vector<double> &vec, const int tag)
{
    assert((m_rank % m_prcs_x) < (m_prcs_x - 1));

    //m_reqs.emplace_back();
    const int dest{m_rank + 1};
    //MPI_Isend(&vec[0], vec.size(), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Send(&vec[0], vec.size(), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
}

void Message_passer_2D::send_top(const std::vector<double> &vec, const int tag)
{
    assert((m_rank + m_prcs_x) < m_size);

    //m_reqs.emplace_back();
    const int dest{m_rank + m_prcs_x};
    //MPI_Isend(&vec[0], vec.size(), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Send(&vec[0], vec.size(), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
}

void Message_passer_2D::send_bottom(const std::vector<double> &vec, const int tag)
{
    assert((m_rank - m_prcs_x) >= 0);

    //m_reqs.emplace_back();
    const int dest{m_rank - m_prcs_x};
    //MPI_Isend(&vec[0], vec.size(), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Send(&vec[0], vec.size(), MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
}

//! receive methods for vectors
void Message_passer_2D::receive_left(std::vector<double> &vec, const int tag)
{
    assert(m_rank % m_prcs_x > 0);

    //m_reqs.emplace_back();
    const int src{m_rank - 1};
    //MPI_Irecv(&vec[0], vec.size(), MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Recv(&vec[0], vec.size(), MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
}

void Message_passer_2D::receive_right(std::vector<double> &vec, const int tag)
{
    assert((m_rank % m_prcs_x) < (m_prcs_x - 1));

    //m_reqs.emplace_back();
    const int src{m_rank + 1};
    //MPI_Irecv(&vec[0], vec.size(), MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Recv(&vec[0], vec.size(), MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
}

void Message_passer_2D::receive_top(std::vector<double> &vec, const int tag)
{
    assert((m_rank + m_prcs_x) < m_size);

    //m_reqs.emplace_back();
    const int src{m_rank + m_prcs_x};
    //MPI_Irecv(&vec[0], vec.size(), MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Recv(&vec[0], vec.size(), MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
}

void Message_passer_2D::receive_bottom(std::vector<double> &vec, const int tag)
{
    assert((m_rank - m_prcs_x) >= 0);

    //m_reqs.emplace_back();
    const int src{m_rank - m_prcs_x};
    //MPI_Irecv(&vec[0], vec.size(), MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Recv(&vec[0], vec.size(), MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
}

//! methods to diagonally send single values
void Message_passer_2D::send_bottom_right(const double val, const int tag)
{
    assert((m_rank - m_prcs_x + 1) >= 0);

    //m_reqs.emplace_back();
    const int dest{m_rank - m_prcs_x + 1};
   // MPI_Isend(&val, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Send(&val, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
}

void Message_passer_2D::send_top_left(const double val, const int tag)
{
    assert((m_rank + (m_prcs_x - 1)) < size());

    //m_reqs.emplace_back();
    const int dest{m_rank + (m_prcs_x - 1)};
    //MPI_Isend(&val, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Send(&val, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
}

void Message_passer_2D::receive_top_left(double &val, const int tag)
{
    assert((m_rank + (m_prcs_x - 1)) < size());

    //m_reqs.emplace_back();
    const int src{m_rank + (m_prcs_x - 1)};
    //MPI_Irecv(&val, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Recv(&val, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
    
}

void Message_passer_2D::receive_bottom_right(double &val, const int tag)
{
    assert((m_rank - m_prcs_x + 1) >= 0);

    //m_reqs.emplace_back();
    const int src{m_rank - m_prcs_x + 1};
    //MPI_Irecv(&val, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &m_reqs.back());
    MPI_Recv(&val, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
}

//! send and receive methods for 2D arrays and doubles
void Message_passer_2D::send(const Array2D &arr, const int tag, const int dest)
{
    m_reqs.emplace_back();
    MPI_Isend(&arr.const_vec_ref()[0], arr.size()[0] * arr.size()[1], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &m_reqs.back());
    //MPI_Send(&arr.const_vec_ref()[0], arr.size()[0] * arr.size()[1], MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
}

void Message_passer_2D::receive(Array2D &arr, const int tag, const int src)
{
    m_reqs.emplace_back();
    MPI_Irecv(&arr.vec_ref()[0], arr.size()[0] * arr.size()[1], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &m_reqs.back());
    //MPI_Recv(&arr.vec_ref()[0], arr.size()[0] * arr.size()[1], MPI_DOUBLE, src, tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
}

void Message_passer_2D::send(const double val, const int tag, const int dest)
{
    m_reqs.emplace_back();
    MPI_Isend(&val, 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &m_reqs.back());
}

void Message_passer_2D::receive(double &val, const int tag, const int src)
{
    m_reqs.emplace_back();
    MPI_Irecv(&val, 1, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &m_reqs.back());
}

void Message_passer_2D::reduce_max(double *result, double local_max, int root) const
{
    MPI_Reduce(&local_max, result, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
}

void Message_passer_2D::broadcast(double *data, int root) const
{
    MPI_Bcast(data, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
}

void Message_passer_2D::broadcast(int *data, int root) const
{
    MPI_Bcast(data, 1, MPI_INT, root, MPI_COMM_WORLD);
}

void Message_passer_2D::gather(int val, std::vector<int> &vec, int root)
{
    MPI_Gather(&val, 1, MPI_INT, &vec[rank()], 1, MPI_INT, root, MPI_COMM_WORLD);
}

void Message_passer_2D::gather(double val, std::vector<double> &vec, int root)
{
    MPI_Gather(&val, 1, MPI_DOUBLE, &vec[rank()], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
}

void Message_passer_2D::wait()
{
    MPI_Waitall(m_reqs.size(), m_reqs.data(), MPI_STATUSES_IGNORE);
}

void Message_passer_2D::barrier()
{
    MPI_Barrier(MPI_COMM_WORLD);
}