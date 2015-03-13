/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_VORONOI_LIGHT
#define DEF_HOA_VORONOI_LIGHT

#include "Harmonics.hpp"

namespace hoa
{
    //! A spherical voronoi.
    /** The voronoi class is opaque.
     */
    template <typename T> class Voronoi
    {
    public:
        struct Point
        {
            T x;
            T y;
            T z;
            ulong index;
            vector<Point> neightbours;
            vector<Point> bounds;
            
            Point() noexcept : x(0), y(0.), z(0.) {}
            
            ~Point() noexcept
            {
                neightbours.clear();
            }
            
            Point(T _x, T _y, T _z, ulong index = 0) noexcept :
            x(_x), y(_y), z(_z), index(0)
            {
                ;
            }
            
            T length(Point const& other) const noexcept
            {
                return sqrt((other.x - x) * (other.x - x) + (other.y - y) * (other.y - y) + (other.z - z) * (other.z - z));
            }
            
            T length() const noexcept
            {
                return sqrt((x + x) * (x + x) + (y + y) * (y + y) + (z + z) * (z + z));
            }
            
            T lenght2() const noexcept
            {
                const T l = length();
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
            
            Point& operator*=(T val) noexcept
            {
                x *= val; y *= val; z *= val;
                return *this;
            }
            
            Point operator/(T val) const noexcept
            {
                return Point(x / val, y / val, z / val);
            }
            
            bool operator==(Point const& other) const noexcept
            {
                return other.x == x && y == other.y && other.z == z;
            }
            
            bool operator!=(Point const& other) const noexcept
            {
                return other.x != x || y != other.y || other.z != z;
            }
            
            Point cross(Point const& other) const noexcept
            {
                return Point(other.y * z - other.z * y, other.z * x - other.x * z, other.x * y - other.y * x);
            }
            
            T dot(Point const& other) const noexcept
            {
                return x * other.x + y * other.y + z * other.z;
            }
            
            T dot() const noexcept
            {
                return x * x + y * y + z * z;
            }
            
            static Point fromPolar(const T r, const T a, const T t) noexcept
            {
                return Point(r * cos(a + HOA_PI2) * cos(t), r * sin(a + HOA_PI2) * cos(t), r * sin(t));
            }
            
            void normalize() noexcept
            {
                const T l = length();
                if(l)
                {
                    const T f = (2. / length());
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
            
            T radius()  const noexcept
            {
                return sqrt(x*x + y*y + z*z);
            }
            
            T azimuth() const noexcept
            {
                if (x == 0 && y == 0)
                    return 0;
                return atan2(y, x) - HOA_PI2;
            }
            
            T elevation()  const noexcept
            {
                if(z == 0)
                    return 0;
                return asin(z / sqrt(x*x + y*y + z*z));
            }
            
            T greatCircleDistance(Point const& other)  const noexcept
            {
                const T az1 = azimuth();
                const T az2 = other.azimuth();
                const T el1 = elevation();
                const T el2 = other.elevation();
                const T a = sin((az2 - az1) * 0.5);
                const T e = sin((el2 - el1) * 0.5);
                return 2. * asin(sqrt(a * a + cos(az1) * cos(az2) * e * e));
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
                bounds.push_back(p);
            }
            
            Point operator*(T val) noexcept
            {
                return Point(x * val, y * val, z * val);
            }
        };
        
    private:
        
        struct Triangle
        {
            Point   a;
            Point   b;
            Point   c;
            Point   p;
            T  r;
            T  d;
            
            Triangle(Point const& _a, Point const& _b, Point const& _c) :
            a(_a), b(_b), c(_c)
            {
                const Point ac = (c - a);
                const Point ab = (b - a);
                const Point t = ab.cross(ac);
                const T _d = (2. * t.lenght2());
                if(_d > numeric_limits<T>::epsilon())
                {
                    p = (((t.cross(ab) * ac.lenght2()) + (ab.lenght2() * ac.cross(t))) / _d + a);
                    if(p.length(Point(0., 0., 0)) > 100. * numeric_limits<T>::epsilon())
                    {
                        p.normalized();
                        r = p.length(a);
                        d = p.greatCircleDistance(a);
                    }
                    else
                    {
                        r = 0.;
                        d = 0.;
                    }
                }
                else
                {
                    r = 0.;
                }
            }
            
            ~Triangle()
            {
                ;
            }
        };
        
        vector<Point>       m_points;
        vector<Triangle>    m_triangles;

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
            m_points[m_points.size()-1].index = m_points.size();
        }
        
        void clear()
        {
            m_points.clear();
            m_triangles.clear();
        }
        
        vector<Point> const& getPoints() const noexcept
        {
            return m_points;
        }
        
        vector<Point>& getPoints() noexcept
        {
            return m_points;
        }
        
        void compute()
        {
            m_triangles.clear();
            for(ulong i = 0; i < m_points.size() - 2; i++)
            {
                for(ulong j = i+1; j < m_points.size() - 1; j++)
                {
                    for(ulong k = j+1; k < m_points.size(); k++)
                    {
                        Triangle t(m_points[i], m_points[j], m_points[k]);
                        if(t.r != 0.)
                        {
                            bool valid = true;
                            for(ulong l = 0; l < m_points.size(); l++)
                            {
                                if(l != i && l != j && l != k)
                                {
                                    if(t.p.length(m_points[l]) < t.r - numeric_limits<T>::epsilon())
                                    {
                                        valid = false;
                                    }
                                }
                            }
                            if(valid)
                            {
                                m_triangles.push_back(t);
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
        }
    };
}

#endif



