#ifndef _NUMERICAL_CARTESIAN_GRID_OF_SPEED_HPP_
#define _NUMERICAL_CARTESIAN_GRID_OF_SPEED_HPP_
#include <cstddef>
#include <vector>
#include <utility>
#include "vector.hpp"
#include "point.hpp"
#include "vortex.hpp"

namespace Numeric
{
    class CartesianGridOfSpeed
    {
    public:
        using vector    = Geometry::Vector<double>;
        using container = std::vector<vector>;
        using point     = Geometry::Point<double>;

        CartesianGridOfSpeed() = default;
        CartesianGridOfSpeed( std::pair<std::size_t,std::size_t> t_dimensions, Geometry::Point<double> m_origin, double t_hStep );
        CartesianGridOfSpeed( CartesianGridOfSpeed const& ) = default;
        CartesianGridOfSpeed( CartesianGridOfSpeed     && ) = default;
        ~CartesianGridOfSpeed() = default;

        point getLeftBottomVertex() const
        {
            return point{m_left, m_bottom};
        }

        point getRightTopVertex() const
        {
            return point{m_left+m_width*m_step, m_bottom+m_height*m_step};
        }

        std::pair<std::size_t,std::size_t> cellGeometry() const
        {
            return {m_width, m_height};
        }

        double getStep() const { return m_step; }

        // Ma modif
        double* data() { return (double*)m_velocityField.data(); }
        double const* data() const { return (double const*)m_velocityField.data(); }
        // Ma modif


        void updateVelocityField( Simulation::Vortices const& t_vortices );


        // Ma modif
        std::size_t numberOfPoints() const { return m_velocityField.size(); }
        // Ma modif

        vector getVelocity( std::size_t iCell, std::size_t jCell ) const 
        {
            return m_velocityField[iCell*m_width+jCell];
        }

        point updatePosition( point const& pt ) const;

        vector computeVelocityFor( point const& p ) const;

        CartesianGridOfSpeed& operator = ( CartesianGridOfSpeed const& ) = default;
        CartesianGridOfSpeed& operator = ( CartesianGridOfSpeed     && ) = default;

    private:
        std::size_t m_width, m_height;
        double      m_left, m_bottom;
        double      m_step;
        container   m_velocityField;
    };
}

#endif