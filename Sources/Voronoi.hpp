/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco, Thomas Le Meur & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_VORONOI_LIGHT
#define DEF_HOA_VORONOI_LIGHT

#include "Defs.hpp"
#include "Math.hpp"

//! @cond

namespace hoa
{
    //! The voronoi class.
    /** The voronoi class discretize the sphere into some points (channels) and some triangles (every channel influence area).
     */
    template <Dimension D> class Voronoi;

    template <> class Voronoi <Hoa3d>
    {

    public:
        //! The voronoi point.
        /** The voronoi point is the ideal place of a channel.
         */
        struct Point
        {
            double x;
            double y;
            double z;
            vector<Point> neightbours;
            vector<Point> bounds;

            //! The point constructor.
            /**	The point constructor allocates and initialize the base classes.
             */
            Point() noexcept :
            x(0), y(0.), z(0.)
            {

            }

            //! The destructor.
            /** The destructor free the memory.
             */
            ~Point() noexcept
            {
                neightbours.clear();
                bounds.clear();
            }

            //! The point constructor.
            /**	The point constructor allocates and initialize the base classes.
            @param x     The x value.
            @param y     The y value.
            @param z     The z value.
            @param index The index of the channel.
             */
            Point(double _x, double _y, double _z) noexcept :
            x(_x), y(_y), z(_z)
            {
                ;
            }

            //! The point constructor by copy.
            /**	The point constructor allocates and initialize the base classes.
            @param other Copy the values of another point into this one.
             */
            Point(Point const& other) noexcept :
            x(other.x), y(other.y), z(other.z)
            {
                ;
            }

            //! Get the distance between this point and another.
            /** Get the distance between this point and another.
            @param other    The other point.
            @return The distance between this point and another.
             */
            double length(Point const& other) const noexcept
            {
                return sqrt((other.x - x) * (other.x - x) + (other.y - y) * (other.y - y) + (other.z - z) * (other.z - z));
            }

            //! Get the distance between this point and the center of the sphere.
            /** Get the distance between this point and the center of the sphere.
            @return The distance between this point and the center of the sphere.
             */
            double length() const noexcept
            {
                return sqrt((x + x) * (x + x) + (y + y) * (y + y) + (z + z) * (z + z));
            }

            //! Get the square of the distance between this point and the center of the sphere.
            /** Get the square of the distance between this point and the center of the sphere.
            @return The square of the distance between this point and the center of the sphere.
             */
            double lenght2() const noexcept
            {
                const double l = length();
                return l * l;
            }

            Point operator-(Point const& other) const noexcept
            {
                return Point(x - other.x, y - other.y, z - other.z);
            }

            Point operator+(Point const& other) const noexcept
            {
                return Point(x + other.x, y + other.y, z + other.z);
            }

            Point operator*(Point const& other) const noexcept
            {
                return Point(x * other.x, y * other.y, z * other.z);
            }

            Point operator*(double val) const noexcept
            {
                return Point(x * val, y * val, z * val);
            }

            Point& operator*=(double val) noexcept
            {
                x *= val; y *= val; z *= val;
                return *this;
            }

            Point operator/(double val) const noexcept
            {
                return Point(x / val, y / val, z / val);
            }

            bool operator==(Point const& other) const noexcept
            {
                return fabs(other.x - x) < HOA_EPSILON && fabs(other.y - y) < HOA_EPSILON && fabs(other.z - z) < HOA_EPSILON;
            }

            bool operator!=(Point const& other) const noexcept
            {
                return fabs(other.x - x) > HOA_EPSILON && fabs(other.y - y) > HOA_EPSILON && fabs(other.z - z) > HOA_EPSILON;
            }

            Point cross(Point const& other) const noexcept
            {
                return Point(other.y * z - other.z * y, other.z * x - other.x * z, other.x * y - other.y * x);
            }

            //! Get the dot product of the point and another.
            /** Get the dot product of the point and another.
            @param  other   The other point.
            @return The dot product of the point and another.
             */
            double dot(Point const& other) const noexcept
            {
                return x * other.x + y * other.y + z * other.z;
            }

            //! Get the dot product of the point and the center of the sphere.
            /** Get the dot product of the point and the center of the sphere.
            @return The dot product of the point and the center of the sphere.
             */
            double dot() const noexcept
            {
                return x * x + y * y + z * z;
            }

            static Point fromPolar(const double r, const double a, const double t) noexcept
            {
                return Point(r * cos(a + HOA_PI2) * cos(t), r * sin(a + HOA_PI2) * cos(t), r * sin(t));
            }

            void normalize() noexcept
            {
                const double l = length();
                if(l != 0.)
                {
                    const double f = (2. / length());
                    x *= f; y *= f; z *= f;
                }
                else
                {
                    x = 0; y = 0; z = 0;
                }
            }

            Point normalized() const noexcept
            {
                Point t = *this;
                t.normalize();
                return t;
            }

            double radius()  const noexcept
            {
                return sqrt(x*x + y*y + z*z);
            }

            double azimuth() const noexcept
            {
                if (x == 0 && y == 0)
                    return 0;
                return atan2(y, x) - HOA_PI2;
            }

            double elevation()  const noexcept
            {
                if(z == 0)
                    return 0;
                return asin(z / sqrt(x*x + y*y + z*z));
            }

            void addNeighbour(Point const& p)
            {
                if(find(neightbours.begin(), neightbours.end(), p) == neightbours.end())
                {
                    neightbours.push_back(p);
                }
            }

            void addBound(Point const& p)
            {
                if(find(bounds.begin(), bounds.end(), p) == bounds.end())
                {
                    bounds.push_back(p);
                }
            }

            void rotateZ(const double _z) noexcept
            {
                const double cosAngle = cos(_z);
                const double sinAngle = sin(_z);
                const double rx = x * cosAngle - y * sinAngle;
                y = x * sinAngle + y * cosAngle;
                x = rx;
            }

            void rotateY(const double _y) noexcept
            {
                const double cosAngle = cos(_y);
                const double sinAngle = sin(_y);
                const double rx = x * cosAngle - z * sinAngle;
                z = x * sinAngle + z * cosAngle;
                x = rx;
            }

            void rotateX(const double _x) noexcept
            {
                const double cosAngle = cos(_x);
                const double sinAngle = sin(_x);
                const double ry = y * cosAngle - z * sinAngle;
                z = y * sinAngle + z * cosAngle;
                y = ry;
            }

            void filterBounds()
            {
                const double el = HOA_PI2 - elevation();
                const double az = azimuth();
                for(ulong i = 0; i < bounds.size(); i++)
                {
                    bounds[i].rotateZ(-az);
                    bounds[i].rotateX(el);
                }
                std::sort(bounds.begin(), bounds.end(), compareAzimuth);
                for(ulong i = 0; i < bounds.size(); i++)
                {
                    bounds[i].rotateX(-el);
                    bounds[i].rotateZ(az);
                }

                ulong size = bounds.size();
                for(ulong i = 0; i < size; i++)
                {
                    const ulong p = i ? i-1 : size-1;
                    const ulong n = (i == size-1) ? 0 : i+1;
                    if(bounds[i].z < 0. && bounds[p].z > 0.)
                    {
                        const double dist = bounds[p].z / (bounds[p].z - bounds[i].z);
                        Point temp((bounds[i].x - bounds[p].x) * dist + bounds[p].x, (bounds[i].y - bounds[p].y) * dist + bounds[p].y, 0.);
                        temp.normalize();
                        bounds.insert(bounds.begin()+int(i), temp);
                        size++;
                    }
                    else if(bounds[i].z < 0. && bounds[n].z > 0.)
                    {
                        const double dist = bounds[n].z / (bounds[n].z - bounds[i].z);
                        Point temp((bounds[i].x - bounds[n].x) * dist + bounds[n].x, (bounds[i].y - bounds[n].y) * dist + bounds[n].y, 0.);
                        temp.normalize();
                        bounds.insert(bounds.begin()+int(i+1), temp);
                        size++;
                    }
                }
                for(ulong i = 0; i < size; i++)
                {
                    if(bounds[i].z < 0.)
                    {
                        bounds.erase(bounds.begin()+int(i));
                        size--; i--;
                    }
                }
            }

            static bool compareAzimuth(Point const& p1, Point const& p2) noexcept
            {
                return p1.azimuth() < p2.azimuth();
            }

            static bool compareElevation(Point const& p1, Point const& p2) noexcept
            {
                return p1.elevation() < p2.elevation();
            }
        };

    private:
        struct Triangle
        {
            Point   a;
            Point   b;
            Point   c;
            Point   p;
            double  r;

            Triangle(Point const& _a, Point const& _b, Point const& _c) noexcept :
            a(_a), b(_b), c(_c), r(0.)
            {
                const Point ac = (c - a);
                const Point ab = (b - a);
                const Point t = ab.cross(ac);
                const double _d = (2. * t.lenght2());
                if(_d > HOA_EPSILON)
                {
                    p = (((t.cross(ab) * ac.lenght2()) + (ac.cross(t) * ab.lenght2())) / _d + a);
                    if(p.length(Point(0., 0., 0)) > HOA_EPSILON)
                    {
                        p.normalize();
                        r = p.length(a);
                    }
                }
            }
        };

        static bool onBottom(Point const& p1) noexcept
        {
            return p1.z < 0.;
        }

        vector<Point>       m_points;
    public:

        Voronoi() noexcept
        {
            ;
        }

        ~Voronoi() noexcept
        {
            clear();
        }

        void add(Point const& p)
        {
            m_points.push_back(p.normalized());
        }

        void clear()
        {
            m_points.clear();
        }

        vector<Point> const& getPoints() const noexcept
        {
            return m_points;
        }

        vector<Point>& getPoints() noexcept
        {
            return m_points;
        }

        vector<Point> const& getBounds(const ulong i) const noexcept
        {
            return m_points[i].bounds;
        }

        vector<Point>& getBounds(const ulong i) noexcept
        {
            return m_points[i].bounds;
        }

        vector<Point> const& getNeightbours(const ulong i) const noexcept
        {
            return m_points[i].neightbours;
        }

        vector<Point>& getNeightbours(const ulong i) noexcept
        {
            return m_points[i].neightbours;
        }

        void compute()
        {
            if(find_if(m_points.begin(), m_points.end(), onBottom) == m_points.end())
            {
                m_points.push_back(Point(0., 0., -1.));
            }
            for(ulong i = 0; i < m_points.size() - 2; i++)
            {
                for(ulong j = i+1; j < m_points.size() - 1; j++)
                {
                    for(ulong k = j+1; k < m_points.size(); k++)
                    {
                        Triangle t(m_points[i], m_points[j], m_points[k]);
                        if(t.r > 0.)
                        {
                            bool valid = true;
                            for(ulong l = 0; l < m_points.size(); l++)
                            {
                                if(l != i && l != j && l != k)
                                {
                                    if(t.p.length(m_points[l]) < t.r - HOA_EPSILON)
                                    {
                                        valid = false;
                                    }
                                }
                            }
                            if(valid)
                            {
                                m_points[i].addNeighbour(m_points[j]);
                                m_points[i].addNeighbour(m_points[k]);
                                m_points[i].addBound(t.p);
                                m_points[j].addNeighbour(m_points[i]);
                                m_points[j].addNeighbour(m_points[k]);
                                m_points[j].addBound(t.p);
                                m_points[k].addNeighbour(m_points[i]);
                                m_points[k].addNeighbour(m_points[j]);
                                m_points[k].addBound(t.p);
                            }
                        }
                    }
                }
            }
            for(ulong i = 0; i < m_points.size(); i++)
            {
                m_points[i].filterBounds();
            }
        }
    };

}
            
//! @endcond

#endif



