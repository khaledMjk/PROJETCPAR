#include <cstddef>
#include <cstdlib>
#include <iostream>
#include "point.hpp"
#include "runge_kutta.hpp"
#include "cartesian_grid_of_speed.hpp"
#include "vortex.hpp"
#include "cloud_of_points.hpp"
#include <mpi.h>
#include <ostream>
#include <tuple>
#include <mpi/mpi.h>
#include <vector>

int main(int argc, char *argv[])
{
    // int rank, nbp;
    MPI_Init(&argc, &argv);
    // MPI_Comm world;
    // MPI_Comm_dup(MPI_COMM_WORLD, &world);

    // MPI_Comm_rank(world, &rank);
    // MPI_Comm_size(world, &nbp);

    // double dt{0.1};
    // std::size_t taille{1};
    // std::size_t nbVortices{1};
    // Numeric::CartesianGridOfSpeed grid;
    // Geometry::CloudOfPoints cloud{taille};

    // cloud[0] = Geometry::Point<double>{1,1};
    // std::cout << "Point before runge-Kutta : " << std::string(cloud[0]) <<std::endl << std::flush;

    // double xleft = 0, ybot = 0, h = 10; 
    // std::size_t nx = 10, ny = 10;

    // grid = Numeric::CartesianGridOfSpeed({nx,ny}, Geometry::Point<double>{xleft,ybot}, h);
    
    // Simulation::Vortices vortices(nbVortices, {grid.getLeftBottomVertex(),
    //                                            grid.getRightTopVertex()});

    // double x = 5,y = 5,force = 5;
    // vortices.setVortex(0,Geometry::Point<double>{x,y}, force);

    // cloud = Numeric::solve_RK4_fixed_vortices(dt, grid, cloud);
    
    // grid.updateVelocityField(vortices);

    // cloud = Numeric::solve_RK4_fixed_vortices(dt, grid, cloud);
    // std::cout << "Point After runge-Kutta : " <<std::string(cloud[0]) << std::endl << std::flush;

    // std::size_t n{10};
    // std::cout << n/3 << std::endl;

    // int c;
    // Get the rank and size in the original communicator
  // Get the rank and size in the original communicator
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int color = (world_rank != 0); // Determine color based on row

    // Split the communicator based on the color and use the
    // original rank for ordering
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &row_comm);

    int row_rank, row_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);

    printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",
        world_rank, world_size, row_rank, row_size);

    MPI_Comm_free(&row_comm);

    int taille_glob = 24;

    if(world_rank == 0)
    {
        MPI_Status stat;
        std::vector<double> buffer(taille_glob/row_size);
        // MPI_Recv(buffer.data(), buffer.size(), MPI_DOUBLE, 1, 101, MPI_COMM_WORLD, &stat);
        for(auto iter : buffer)
        {
            std::cout<< iter << ",";
        }
    }
    else
    {
        if(row_rank !=0)
        {   
            int taille_loc = taille_glob/row_size;
            double buffer = row_rank;
            std::vector<double> buffer_glob(taille_glob);
            MPI_Gather(&buffer, 1, MPI_DOUBLE, buffer_glob.data(), 1, MPI_DOUBLE, 0, row_comm);
            }
        else
        {
            int taille_loc = taille_glob/row_size;
            double buffer_loc = row_rank;
            std::vector<double> buffer_glob(taille_glob/row_size);
            MPI_Gather(&buffer_loc,1, MPI_DOUBLE, buffer_glob.data(), 1, MPI_DOUBLE, 0, row_comm);
            // MPI_Send(buffer_glob.data(), buffer_glob.size(), MPI_DOUBLE, 0, 101, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}