
#include "SOR_parallel.h"

SOR_parallel::SOR_parallel(const Discretization &disc, Message_passer_2D &MP, double eps, int max_it, double w)
    : Pressure_solver(disc.dx(), disc.dy(), eps, max_it, MP), m_w{w}, m_disc{disc}, size{disc.nCells()[0] + 2, disc.nCells()[1] + 2}
{
   /*
   p_right_send_buffer.resize(size[1] - 2);
   p_left_send_buffer.resize(size[1] - 2);
   p_top_send_buffer.resize(size[0] - 2);
   p_bottom_send_buffer.resize(size[0] - 2);

   p_right_receive_buffer.resize(size[1] - 2);
   p_left_receive_buffer.resize(size[1] - 2);
   p_top_receive_buffer.resize(size[0] - 2);
   p_bottom_receive_buffer.resize(size[0] - 2);
   */

   p_right_send_buffer.resize(size[1] / 2);
   p_left_send_buffer.resize(size[1] / 2);
   p_top_send_buffer.resize(size[0] / 2);
   p_bottom_send_buffer.resize(size[0] / 2);

   p_right_receive_buffer.resize(size[1] / 2);
   p_left_receive_buffer.resize(size[1] / 2);
   p_top_receive_buffer.resize(size[0] / 2);
   p_bottom_receive_buffer.resize(size[0] / 2);
}

void SOR_parallel::communicate_update_boundary_values(Staggered_grid &p, const bool black)
{
   //! define start values for boundary exchanges (0 or 1)
   const int send_start_top{black ? (size[1] + 1) % 2 : size[1] % 2};
   const int send_start_left{black ? 0 : 1};
   const int send_start_right{black ? (size[0] + 1) % 2 : size[0] % 2};
   const int send_start_bottom{black ? 0 : 1};

   const int recv_start_top{!black ? (size[1] + 1) % 2 : size[1] % 2};
   const int recv_start_left{!black ? 0 : 1};
   const int recv_start_right{!black ? (size[0] + 1) % 2 : size[0] % 2};
   const int recv_start_bottom{!black ? 0 : 1};

   //! communicate with all neighbours / get boundary values of p
   //! right
   if (m_disc.ex_nbrs(SIDE_RIGHT)) // if right neighbour exist, get p values, else set boundary conditions
   {
      //! extract p
      //p_right_send_buffer = p.get_side(SIDE_RIGHT);
      p_right_send_buffer = p.get_side_odd(SIDE_RIGHT, send_start_right);

      //! send p
      m_MP.send_right(p_right_send_buffer, TAG_P);

      //! receive p
      m_MP.receive_right(p_right_receive_buffer, TAG_P);
   }
   else
   {
      //! right boundary
      for (int j{1}; j < size[1] - 1; ++j) // update with homogenous dirichlet Bc
      {
         p(size[0] - 1, j) = p(size[0] - 2, j);
      }
   }

   //! left
   if (m_disc.ex_nbrs(SIDE_LEFT)) // if left neighbour exist, get p values, else set boundary conditions
   {
      //! extract p
      //p_left_send_buffer = p.get_side(SIDE_LEFT);
      p_left_send_buffer = p.get_side_odd(SIDE_LEFT, send_start_left);


      //! send p
      m_MP.send_left(p_left_send_buffer, TAG_P);

      //! receive p
      m_MP.receive_left(p_left_receive_buffer, TAG_P);
   }
   else
   {
      //! left boundary
      for (int j{1}; j < size[1] - 1; ++j) // update with homogenous dirichlet Bc
      {
         p(0, j) = p(1, j);
      }
   }

   //! top
   if (m_disc.ex_nbrs(SIDE_TOP)) // if top neighbour exist, get p values, else set boundary conditions
   {
      //! extract p
      //p_top_send_buffer = p.get_side(SIDE_TOP);
      p_top_send_buffer = p.get_side_odd(SIDE_TOP,send_start_top);

      //! send p
      m_MP.send_top(p_top_send_buffer, TAG_P);

      //! receive p
      m_MP.receive_top(p_top_receive_buffer, TAG_P);
   }
   else
   {
      //! top boundary
      for (int i{1}; i < size[0] - 1; ++i) // update with homogenous dirichlet Bc
      {
         p(i, size[1] - 1) = p(i, size[1] - 2);
      }
   }

   //! bottom
   if (m_disc.ex_nbrs(SIDE_BOTTOM)) // if bottom neighbour exist, get p values, else set boundary conditions
   {
      //! extract p
      //p_bottom_send_buffer = p.get_side(SIDE_BOTTOM);
      p_bottom_send_buffer = p.get_side_odd(SIDE_BOTTOM,send_start_bottom);

      //! send p
      m_MP.send_bottom(p_bottom_send_buffer, TAG_P);

      //! receive p
      m_MP.receive_bottom(p_bottom_receive_buffer, TAG_P);
   }
   else
   {
      //! bottom boundary
      for (int i{1}; i < size[0] - 1; ++i) // update with homogenous dirichlet Bc
      {
         p(i, 0) = p(i, 1);
      }
   }

   //! wait for message passing to complete
   //m_MP.wait();

   //! write boundaries
   //! right
   if (m_disc.ex_nbrs(SIDE_RIGHT))
   {
      //p.write_bound(SIDE_RIGHT, p_right_receive_buffer);
      p.write_bound_odd(SIDE_RIGHT, p_right_receive_buffer, recv_start_right);
   }

   //! left
   if (m_disc.ex_nbrs(SIDE_LEFT))
   {
      //p.write_bound(SIDE_LEFT, p_left_receive_buffer);
      p.write_bound_odd(SIDE_LEFT, p_left_receive_buffer,recv_start_left);
   }

   //! top
   if (m_disc.ex_nbrs(SIDE_TOP))
   {
      //p.write_bound(SIDE_TOP, p_top_receive_buffer);
      p.write_bound_odd(SIDE_TOP, p_top_receive_buffer,recv_start_top);
   }

   //! bottom
   if (m_disc.ex_nbrs(SIDE_BOTTOM))
   {
      //p.write_bound(SIDE_BOTTOM, p_bottom_receive_buffer);
      p.write_bound_odd(SIDE_BOTTOM, p_bottom_receive_buffer,recv_start_bottom);
   }
}

void SOR_parallel::step_black(Staggered_grid &p, const Staggered_grid &RHS)
{
   // run SOR-iteration over all "black" cells
   for (int j{1}; j < size[1] - 1; ++j)
   {
      for (int i{(j + 1) % 2 + 1}; i < size[0] - 1; i = i + 2)
      {
         //! single node iteration step
         p(i, j) =
             (1 - m_w) * p(i, j) + m_w                                                                            //< factor omega
                                       * std::pow(m_dx * m_dy, 2) / (2 * (std::pow(m_dx, 2) + std::pow(m_dy, 2))) //< prefactor
                                       * ((p(i - 1, j) + p(i + 1, j)) / std::pow(m_dx, 2)                         //< horizontal term
                                          + (p(i, j - 1) + p(i, j + 1)) / std::pow(m_dy, 2)                       //< vertical term
                                          - RHS(i, j));                                                           //< right hand side
      }
   }
}

void SOR_parallel::step_red(Staggered_grid &p, const Staggered_grid &RHS)
{
   // run SOR-iteration over all "red" cells
   for (int j{1}; j < size[1] - 1; ++j)
   {
      for (int i{(j % 2) + 1}; i < size[0] - 1; i = i + 2)
      {
         //! single node iteration step
         p(i, j) =
             (1 - m_w) * p(i, j) + m_w                                                                            //< factor omega
                                       * std::pow(m_dx * m_dy, 2) / (2 * (std::pow(m_dx, 2) + std::pow(m_dy, 2))) //< prefactor
                                       * ((p(i - 1, j) + p(i + 1, j)) / std::pow(m_dx, 2)                         //< horizontal term
                                          + (p(i, j - 1) + p(i, j + 1)) / std::pow(m_dy, 2)                       //< vertical term
                                          - RHS(i, j));                                                           //< right hand side
      }
   }
}

void SOR_parallel::run_it_step(Staggered_grid &p, const Staggered_grid &RHS, bool black_red)
{
   //! Assume p is up-to-date
   //! Decide whether to use black-red or red-black scheme
   if (black_red) // choose black-red
   {
      // first black
      step_black(p, RHS);
      // get boundary values from neighbour processes / Bc
      communicate_update_boundary_values(p, true);

      // then red
      step_red(p, RHS);
      // get boundary values from neighbour processes / Bc
      communicate_update_boundary_values(p, false);
   }
   else // choose red-black
   {
      // first red
      step_red(p, RHS);
      // get boundary values from neighbour processes / Bc
      communicate_update_boundary_values(p, false);

      // then red
      step_black(p, RHS);
      // get boundary values from neighbour processes / Bc
      communicate_update_boundary_values(p, true);
   }
}
