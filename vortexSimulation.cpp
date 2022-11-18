#include <SFML/Window/Keyboard.hpp>
#include <cstddef>
#include <ios>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <chrono>
#include "cartesian_grid_of_speed.hpp"
#include "point.hpp"
#include "vortex.hpp"
#include "cloud_of_points.hpp"
#include "runge_kutta.hpp"
#include "screen.hpp"
#include <mpi/mpi.h>
#include <vector>

std::size_t get_sqr_size(std::size_t nbPoints)
{
    std::size_t sqrtNbPoints = std::size_t(std::sqrt(nbPoints));
    std::size_t nbPointsY    = nbPoints/sqrtNbPoints;
    std::size_t nbPointsX    = sqrtNbPoints + (nbPoints%sqrtNbPoints > 0 ? 1 : 0);
    return nbPointsX*nbPointsY;
}

std::size_t get_sqt_tot_size(std::size_t nbPoints, int nbp)
{
    std::size_t res = 0;
    std::vector<std::size_t> nbPoints_loc(nbp - 1, std::size_t(nbPoints/(nbp - 1)));
    for(std::size_t iter = 0 ; iter < std::size_t(nbPoints%(nbp -1)); iter++)
    {
        nbPoints_loc[iter]++;
    }
    for(auto iter : nbPoints_loc)
    {
        res += get_sqr_size(iter);
    }
    return res;
}


auto readConfigFile( std::ifstream& input, int rank , int nbp)
{
    using point=Simulation::Vortices::point;

    int isMobile;
    std::size_t nbVortices;
    Numeric::CartesianGridOfSpeed cartesianGrid;
    Geometry::CloudOfPoints cloudOfPoints;
    constexpr std::size_t maxBuffer = 8192;
    char buffer[maxBuffer];
    std::string sbuffer;
    std::stringstream ibuffer;
    // Lit la première ligne de commentaire :
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer);// Lecture de la grille cartésienne
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    double xleft, ybot, h;
    std::size_t nx, ny;
    ibuffer >> xleft >> ybot >> nx >> ny >> h;
    cartesianGrid = Numeric::CartesianGridOfSpeed({nx,ny}, point{xleft,ybot}, h);
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit mode de génération des particules
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    int modeGeneration;
    ibuffer >> modeGeneration;
    if (modeGeneration == 0) // Génération sur toute la grille
    {
        std::size_t nbPoints;
        ibuffer >> nbPoints;
        if (rank != 0)
        {
            // On veut sectionner le plan entre les proccessus de 1 à nbp 
            std::size_t nbPoints_loc{nbPoints/(nbp-1)};
            if(rank < int(nbPoints%(nbp -1)))
                nbPoints_loc++;

            Geometry::Point<double> LeftBottom_loc{cartesianGrid.getLeftBottomVertex()};
            Geometry::Point<double> RightTop_loc{cartesianGrid.getRightTopVertex()};

            double longueur = Geometry::Point<double>{xleft,ybot+ny*h}.computeDistance(LeftBottom_loc);

            LeftBottom_loc.y = ybot + (longueur/(nbp - 1)) *(rank - 1);
            RightTop_loc.y = ybot + (longueur/(nbp - 1)) *(rank);

            cloudOfPoints = Geometry::generatePointsIn(nbPoints_loc, {LeftBottom_loc, RightTop_loc});
            // cloudOfPoints = Geometry::generatePointsIn(nbPoints, {cartesianGrid.getLeftBottomVertex(), cartesianGrid.getRightTopVertex()});
        }
        else
        {
            cloudOfPoints = Geometry::CloudOfPoints(get_sqt_tot_size(nbPoints, nbp));
        }
    }
    else
    {
        std::size_t nbPoints;
        double xl, xr, yb, yt;
        ibuffer >> xl >> yb >> xr >> yt >> nbPoints;
        if (rank !=0)
            cloudOfPoints = Geometry::generatePointsIn(nbPoints, {point{xl,yb}, point{xr,yt}});
        else
        {
            std::size_t sqrtNbPoints = std::size_t(std::sqrt(nbPoints));
            std::size_t nbPointsY    = nbPoints/sqrtNbPoints;
            std::size_t nbPointsX    = sqrtNbPoints + (nbPoints%sqrtNbPoints > 0 ? 1 : 0);
            cloudOfPoints = Geometry::CloudOfPoints(nbPointsX*nbPointsY);
        }
    }
    // Lit le nombre de vortex :
    input.getline(buffer, maxBuffer); // Relit un commentaire
    input.getline(buffer, maxBuffer); // Lit le nombre de vortex
    sbuffer = std::string(buffer, maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    try {
        ibuffer >> nbVortices;
    } catch(std::ios_base::failure& err)
    {
        std::cout << "Error " << err.what() << " found" << std::endl;
        std::cout << "Read line : " << sbuffer << std::endl;
        throw err;
    }
    Simulation::Vortices vortices(nbVortices, {cartesianGrid.getLeftBottomVertex(),
                                               cartesianGrid.getRightTopVertex()});
    input.getline(buffer, maxBuffer);// Relit un commentaire
    for (std::size_t iVortex=0; iVortex<nbVortices; ++iVortex)
    {
        input.getline(buffer, maxBuffer);
        double x,y,force;
        std::string sbuffer(buffer, maxBuffer);
        std::stringstream ibuffer(sbuffer);
        ibuffer >> x >> y >> force;
        vortices.setVortex(iVortex, point{x,y}, force);
    }
    input.getline(buffer, maxBuffer);// Relit un commentaire
    input.getline(buffer, maxBuffer);// Lit le mode de déplacement des vortex
    sbuffer = std::string(buffer,maxBuffer);
    ibuffer = std::stringstream(sbuffer);
    ibuffer >> isMobile;
    return std::make_tuple(vortices, isMobile, cartesianGrid, cloudOfPoints);
}



///////////////////////////////////////////////////////////////////////////////////////////////
//  Here is Main                                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////

int main( int nargs, char* argv[] )
{
    MPI_Init(&nargs, &argv);
    int world_rank, world_nbp;
    MPI_Comm glob;

    MPI_Comm_dup(MPI_COMM_WORLD, &glob);
    MPI_Comm_size(glob, &world_nbp);
    MPI_Comm_rank(glob, &world_rank);

    //Here comes the "me trying to do some goofy stuff"
    
    int color = (world_rank != 0);

    MPI_Comm row_comm;
    MPI_Comm_split(glob, color, world_rank, &row_comm);

    int row_rank, row_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);

    ///////////////////////////

    char const* filename;
    if (nargs==1)
    {
        std::cout << "Usage : vortexsimulator <nom fichier configuration>" << std::endl;
        return EXIT_FAILURE;
    }

    filename = argv[1];
    std::ifstream fich(filename);
    auto config = readConfigFile(fich, world_rank, world_nbp);
    fich.close();

    std::size_t resx=800, resy=600;
    if (nargs>3)
    {
        resx = std::stoull(argv[2]);
        resy = std::stoull(argv[3]);
    }

    auto vortices = std::get<0>(config);
    auto isMobile = std::get<1>(config);
    auto grid     = std::get<2>(config);
    auto cloud    = std::get<3>(config);


    if(world_rank == 0)
    {
        std::cout << "######## Vortex simultor ########" << std::endl << std::endl;
        std::cout << "Press P for play animation " << std::endl;
        std::cout << "Press S to stop animation" << std::endl;
        std::cout << "Press right cursor to advance step by step in time" << std::endl;
        std::cout << "Press down cursor to halve the time step" << std::endl;
        std::cout << "Press up cursor to double the time step" << std::endl;
    }

    grid.updateVelocityField(vortices);

    if(world_rank != 0)
    {
        // MPI_Send(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 0, 101, glob);
        if(row_rank != 0)
        {
            MPI_Gather(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0,row_comm);
        }
        else
        {
            Geometry::CloudOfPoints buffer = Geometry::CloudOfPoints(cloud.numberOfPoints()*(row_size));
            MPI_Gather(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, buffer.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 0, row_comm);
            MPI_Send(buffer.data(), buffer.numberOfPoints()*2, MPI_DOUBLE, 0, 101, glob);
        }
    }
    else
    {
        MPI_Status Stat;
        MPI_Recv(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 1, 101, glob, &Stat);    
    }

    if( world_rank == 0)
    {
        // The job that will be sent to 1 to be excecuted
        int job = 0;

        Graphisme::Screen myScreen( {resx,resy}, {grid.getLeftBottomVertex(), grid.getRightTopVertex()} );
        bool animate = false;
        double dt = 0.1;

        while (myScreen.isOpen())
        {
            auto start = std::chrono::system_clock::now();
            bool advance = false;
            // on inspecte tous les évènements de la fenêtre qui ont été émis depuis la précédente itération
            sf::Event event;
            while (myScreen.pollEvent(event))
            {
                // évènement "fermeture demandée" : on ferme la fenêtre
                if (event.type == sf::Event::Closed)
                {
                    myScreen.close();
                    job = 2 ;
                    // // Part 1
                    // MPI_Send(&job, 1, MPI_INT, 1, 101, glob);

                    // Part 3
                    for(int i_dest = 1; i_dest<world_nbp; i_dest++)
                    {
                        MPI_Send(&job, 1, MPI_INT, i_dest, i_dest, glob);
                        MPI_Send(&dt, 1, MPI_DOUBLE, i_dest, i_dest, glob);        
                    }
                }
                if (event.type == sf::Event::Resized)
                {
                    // on met à jour la vue, avec la nouvelle taille de la fenêtre
                    myScreen.resize(event);
                }
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::P)) animate = true;
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::S)) animate = false;
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up)) dt *= 2;
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down)) dt /= 2;
                if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) advance = true;
            }
            if (animate | advance)
            {
                if (isMobile)
                {
                    MPI_Status Stat;
                    job = 1;
                    
                    // Sequential programming
                    // cloud = Numeric::solve_RK4_movable_vortices(dt, grid, vortices, cloud);

                    // Part 2 Sending only to 1
                    // MPI_Send(&job, 1, MPI_INT, 1, 101, glob);
                    // MPI_Send(&dt, 1, MPI_DOUBLE, 1, 101, glob);
                    // MPI_Recv(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 1, 1, glob, &Stat);
                    // MPI_Recv(grid.data(), grid.numberOfPoints()*2, MPI_DOUBLE, 1, 2, glob, &Stat);
                    // MPI_Recv(vortices.data(), vortices.numberOfVortices()*3, MPI_DOUBLE, 1, 3, glob, &Stat);

                    // Part 3 Sending 
                    for(int i_dest = 1; i_dest<world_nbp; i_dest++)
                    {
                        MPI_Send(&job, 1, MPI_INT, i_dest, i_dest, glob);
                        MPI_Send(&dt, 1, MPI_DOUBLE, i_dest, i_dest, glob);        
                    }
                    MPI_Recv(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 1, 1, glob, &Stat);
                    MPI_Recv(grid.data(), grid.numberOfPoints()*2, MPI_DOUBLE, 1, 2, glob, &Stat);
                    MPI_Recv(vortices.data(), vortices.numberOfVortices()*3, MPI_DOUBLE, 1, 3, glob, &Stat);
                }
                else
                {
                    MPI_Status Stat;
                    //compute cloud
                    job = 0;

                    // // Sequential part
                    // cloud = Numeric::solve_RK4_fixed_vortices(dt, grid, cloud);
                    
                    // // Part 2 Sending only to 1
                    // MPI_Send(&job, 1, MPI_INT, 1, 101, glob);
                    // MPI_Send(&dt, 1, MPI_DOUBLE, 1, 101, glob);
                    // MPI_Recv(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 1, 101, glob, &Stat);
                    
                    // Part 3 sending calculation to be done by all procecces
                    for(int i_dest = 1; i_dest<world_nbp; i_dest++)
                    {
                        MPI_Send(&job, 1, MPI_INT, i_dest, i_dest, glob);
                        MPI_Send(&dt, 1, MPI_DOUBLE, i_dest, i_dest, glob);        
                    }
                    MPI_Recv(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 1, 101, glob, &Stat);
                }
            }

            //display stuff
            myScreen.clear(sf::Color::Black);
            std::string strDt = std::string("Time step : ") + std::to_string(dt);
            myScreen.drawText(strDt, Geometry::Point<double>{50, double(myScreen.getGeometry().second-96)});
            myScreen.displayVelocityField(grid, vortices);

            myScreen.displayParticles(grid, vortices, cloud);
            
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> diff = end - start;
            std::string str_fps = std::string("FPS : ") + std::to_string(1./diff.count());
            myScreen.drawText(str_fps, Geometry::Point<double>{300, double(myScreen.getGeometry().second-96)});
            myScreen.display();
        }
    }
    else
    {
        MPI_Status Stat;
        int job = 0;
        double dt = 0.1;
        
        while(job != 2)
        {
            // // Part 2
            // MPI_Recv(&job, 1, MPI_INT, 0, 101, glob, &Stat);
            // MPI_Recv(&dt, 1, MPI_DOUBLE, 0, 101, glob, &Stat); 

            // Part 3
            MPI_Recv(&job, 1, MPI_INT, 0, world_rank, glob, &Stat);
            MPI_Recv(&dt, 1, MPI_DOUBLE, 0, world_rank, glob, &Stat); 

            if(job == 0)
            {
                cloud = Numeric::solve_RK4_fixed_vortices(dt, grid, cloud);
                // // Part 2
                // MPI_Send(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 0, 1, glob);

                // Part 3
                if (row_rank != 0)
                {
                    // Everyone send their result to 1 (or 0 in our local ranking)
                    MPI_Gather(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0,row_comm);
                }
                else
                {
                    // 1 needs to send and to recieve everything
                    Geometry::CloudOfPoints buffer = Geometry::CloudOfPoints(cloud.numberOfPoints()*(row_size));
                    MPI_Gather(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, buffer.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 0, row_comm);
                    MPI_Send(buffer.data(), buffer.numberOfPoints()*2, MPI_DOUBLE, 0, 101, glob);
                }
            }
            else if(job == 1)
            {

                cloud = Numeric::solve_RK4_movable_vortices(dt, grid, vortices, cloud);

                // //Part 2
                // MPI_Send(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 0, 1, glob);
                // MPI_Send(grid.data(), grid.numberOfPoints()*2, MPI_DOUBLE, 0, 2, glob);
                // MPI_Send(vortices.data(), vortices.numberOfVortices()*3, MPI_DOUBLE, 0, 3, glob); 

                //Part 3
                if (row_rank != 0)
                {
                    // Everyone send their result to 1 (or 0 in our local ranking)
                    MPI_Gather(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0,row_comm);
                }
                else
                {
                    // 1 needs to send and to recieve everything
                    Geometry::CloudOfPoints buffer = Geometry::CloudOfPoints(cloud.numberOfPoints()*(row_size));
                    MPI_Gather(cloud.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, buffer.data(), cloud.numberOfPoints()*2, MPI_DOUBLE, 0, row_comm);
                    MPI_Send(buffer.data(), buffer.numberOfPoints()*2, MPI_DOUBLE, 0, 1, glob);
                    MPI_Send(grid.data(), grid.numberOfPoints()*2, MPI_DOUBLE, 0, 2, glob);
                    MPI_Send(vortices.data(), vortices.numberOfVortices()*3, MPI_DOUBLE, 0, 3, glob); 
                }
            }
        }
    }

    //finalize the split
    MPI_Comm_free(&row_comm);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
