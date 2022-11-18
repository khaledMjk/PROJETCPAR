#include "screen.hpp"
#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/PrimitiveType.hpp>
#include <SFML/Graphics/Rect.hpp>
#include <SFML/Graphics/VertexArray.hpp>
#include <iostream>
#include <omp.h>

Graphisme::Screen::Screen( std::pair<std::size_t,std::size_t> const& t_geometry,
                           std::pair<Geometry::Point<double>,Geometry::Point<double>> const& t_domain )
    :   m_window(sf::VideoMode(t_geometry.first, t_geometry.second), "Vortex simulation")
{
    if (!m_font.loadFromFile("Aller_Bd.ttf"))
    {
        throw std::ios_base::failure("Failed to load font file Aller_Bd.ttf");
    }
    auto screenSize = m_window.getSize();
    m_velocityView.setCenter(0.25*screenSize.x,0.5*(screenSize.y-128));
    m_velocityView.setSize  (0.5*screenSize.x, screenSize.y-128);
    m_velocityView.setViewport(sf::FloatRect(0.01f,0.05f,0.48f,0.8f));

    m_particlesView.setCenter(0.25*screenSize.x,0.5*(screenSize.y-128));
    m_particlesView.setSize  (0.5*screenSize.x, screenSize.y-128);
    m_particlesView.setViewport(sf::FloatRect(0.51f,0.05f,0.48f,0.8f));
}

void 
Graphisme::Screen::displayVelocityField( Numeric::CartesianGridOfSpeed const& grid, Simulation::Vortices const& vortices )
{
    using vector=Numeric::CartesianGridOfSpeed::vector;
    m_window.setView(m_velocityView);
    // Affichage moité gauche de l'écran : 
    auto screenSize = m_velocityView.getSize();
    std::size_t width = screenSize.x-1;
    std::size_t height= screenSize.y-1;

    vector domainDimension{grid.getLeftBottomVertex(),grid.getRightTopVertex()};

    auto nbCells = grid.cellGeometry();
    double scalex = width/domainDimension.x;
    double scaley = height/domainDimension.y;
    double hx     = grid.getStep()*scalex;
    double hy     = grid.getStep()*scaley;

    sf::VertexArray drawVertices(sf::Lines, 2*(nbCells.first+nbCells.second+2+nbCells.first*nbCells.second));
    std::size_t indDrawVert = 0;

    for (std::size_t ix=0; ix<=nbCells.first; ++ix)
    {
        drawVertices[indDrawVert] = sf::Vertex(sf::Vector2f(ix*hx, 0));
        drawVertices[indDrawVert].color = sf::Color(64,0,0);
        ++indDrawVert;
        drawVertices[indDrawVert] = sf::Vertex(sf::Vector2f(ix*hx, height));
        drawVertices[indDrawVert].color = sf::Color(64,0,0);
        ++indDrawVert;
        /*sf::Vertex vline[] = {
            sf::Vertex(sf::Vector2f(ix*hx, 0)),
            sf::Vertex(sf::Vector2f(ix*hx, height))
        };
        vline[0].color = sf::Color(64,0,0);
        vline[1].color = sf::Color(64,0,0);*/
        //m_window.draw(vline, 2, sf::Lines);
    }

    for (std::size_t jy=0; jy<=nbCells.second; ++jy)
    {
        drawVertices[indDrawVert] = sf::Vertex(sf::Vector2f(0,jy*hy));
        drawVertices[indDrawVert].color = sf::Color(64,0,0);
        ++indDrawVert;
        drawVertices[indDrawVert] = sf::Vertex(sf::Vector2f(width, jy*hy));
        drawVertices[indDrawVert].color = sf::Color(64,0,0);
        ++indDrawVert;
        /*sf::Vertex hline[] = {
            sf::Vertex(sf::Vector2f(0,jy*hy)),
            sf::Vertex(sf::Vector2f(width, jy*hy))
        };
        hline[0].color = sf::Color(64,0,0);
        hline[1].color = sf::Color(64,0,0);
        m_window.draw(hline, 2, sf::Lines);*/
    } 
    // Affichage des vecteurs :
    for (std::size_t jy=0; jy<nbCells.second; ++jy)
    {
        for (std::size_t ix=0; ix<nbCells.first; ++ix)
        {
            double xdeb = (ix+0.5)*hx;
            double ydeb = (jy+0.5)*hy;
            double xfin = xdeb + scalex*grid.getVelocity(jy, ix).x;
            double yfin = ydeb + scaley*grid.getVelocity(jy, ix).y;
            drawVertices[indDrawVert] = sf::Vertex(sf::Vector2f(xdeb,ydeb));
            drawVertices[indDrawVert].color = sf::Color(0,0,127);
            ++indDrawVert;
            drawVertices[indDrawVert] = sf::Vertex(sf::Vector2f(xfin,yfin));
            drawVertices[indDrawVert].color = sf::Color(255,255,255);
            ++indDrawVert;
            /*sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(xdeb,ydeb)),
                sf::Vertex(sf::Vector2f(xfin,yfin))
            };
            line[0].color = sf::Color(0,0,127);
            line[1].color = sf::Color(255,255,255);*/
            //m_window.draw(line, 2, sf::Lines);
        }
    }
    m_window.draw(drawVertices);
    // Affichage des vortices :
    for (std::size_t iVort=0; iVort<vortices.numberOfVortices(); ++iVort)
    {
        auto c = vortices.getCenter(iVort);
        sf::CircleShape shape{5};
        shape.setPosition(scalex*(c.x-grid.getLeftBottomVertex().x)-5,scaley*(c.y-grid.getLeftBottomVertex().y)-5);
        shape.setFillColor(sf::Color::Red);
        m_window.draw(shape);
    }
    m_window.setView(m_window.getDefaultView());
}
//-----------------------------------------------------------------------------------------------------------
void 
Graphisme::Screen::displayParticles( Numeric::CartesianGridOfSpeed const& grid, Simulation::Vortices const& vortices, Geometry::CloudOfPoints const& points )
{
    using vector=Numeric::CartesianGridOfSpeed::vector;
    m_window.setView(m_particlesView);
    // Affichage moité gauche de l'écran : 
    auto screenSize = m_particlesView.getSize();
    std::size_t width = screenSize.x-1;
    std::size_t height= screenSize.y-1;

    vector domainDimension{grid.getLeftBottomVertex(),grid.getRightTopVertex()};

    double scalex = width/domainDimension.x;
    double scaley = height/domainDimension.y;


    // Affichage des particules en 3/4 transparents :
    sf::VertexArray drawVertices(sf::Points, points.numberOfPoints());
    for (std::size_t iPt=0; iPt<points.numberOfPoints(); ++iPt)
    {
        drawVertices[iPt].position = sf::Vector2f{
            float(scalex*(points[iPt].x-grid.getLeftBottomVertex().x)),
            float(scaley*(points[iPt].y-grid.getLeftBottomVertex().y))};
        drawVertices[iPt].color = sf::Color(255,255,255,128);
    }
    m_window.draw(drawVertices);
    // Affichage des vortices :
    for (std::size_t iVort=0; iVort<vortices.numberOfVortices(); ++iVort)
    {
        auto c = vortices.getCenter(iVort);
        sf::CircleShape shape{5};
        shape.setPosition(scalex*(c.x-grid.getLeftBottomVertex().x)-5,scaley*(c.y-grid.getLeftBottomVertex().y)-5);
        shape.setFillColor(sf::Color::Red);
        m_window.draw(shape);
    }
    m_window.setView(m_window.getDefaultView());
}
//
void 
Graphisme::Screen::drawText( std::string const& t_string, Geometry::Point<double> const& position)
{
    static sf::Text text;
    text.setFont(m_font);
    text.setString(t_string);
    text.setCharacterSize(12);
    text.setFillColor(sf::Color::White);
    text.setPosition(position.x,position.y);
    m_window.draw(text);
}
