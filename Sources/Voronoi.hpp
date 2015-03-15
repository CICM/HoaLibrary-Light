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
#define HOA_EPSILON 1e-6
    
    //! A spherical voronoi.
    /** The voronoi class is opaque.
     */
    template <Dimension D, typename T> class Voronoi;
    
    template <typename T> class Voronoi<Hoa3d, T>
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
            
            Point() noexcept :
            x(0), y(0.), z(0.)
            {
                
            }
            
            ~Point() noexcept
            {
                neightbours.clear();
            }
            
            Point(T _x, T _y, T _z, ulong index = 0) noexcept :
            x(_x), y(_y), z(_z), index(0)
            {
                ;
            }
            
            Point(Point const& other) noexcept :
            x(other.x), y(other.y), z(other.z), index(other.index)
            {
                ;
            }
            
            T distance(Point const& other) const noexcept
            {
                return sqrt((other.x - x) * (other.x - x) + (other.y - y) * (other.y - y) + (other.z - z) * (other.z - z));
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
            
            Point operator*(T val) const noexcept
            {
                return Point(x * val, y * val, z * val);
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
                return fabs(other.x - x) < HOA_EPSILON && fabs(y - other.y) < HOA_EPSILON && fabs(other.z - z) < HOA_EPSILON;
            }
            
            bool operator!=(Point const& other) const noexcept
            {
                return other.x != x || y != other.y || other.z != z;
            }
            
            Point cross(Point const& other) const noexcept
            {
                return Point(other.y * z - other.z * y, other.z * x - other.x * z, other.x * y - other.y * x);
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
            
            void rotateZ(const T _z) noexcept
            {
                const T cosAngle = cos(_z);
                const T sinAngle = sin(_z);
                const T rx = x * cosAngle - y * sinAngle;
                y = x * sinAngle + y * cosAngle;
                x = rx;
            }
            
            void rotateY(const T _y) noexcept
            {
                const T cosAngle = cos(_y);
                const T sinAngle = sin(_y);
                const T rx = x * cosAngle - z * sinAngle;
                z = x * sinAngle + z * cosAngle;
                x = rx;
            }
            
            void rotateX(const T _x) noexcept
            {
                const T cosAngle = cos(_x);
                const T sinAngle = sin(_x);
                const T ry = y * cosAngle - z * sinAngle;
                z = y * sinAngle + z * cosAngle;
                y = ry;
            }
            
            
            void computeView(const bool top = true)
            {
                bool valid = false;
                for(ulong i = 0; i < bounds.size(); i++)
                {
                    if(bounds[i].z > 0.)
                    {
                        valid = true;
                    }
                }
                if(!valid || bounds.size() < 3)
                {
                    bounds.clear();
                }
                else
                {
                    ulong size = bounds.size();
                    for(ulong i = 0; i < size;)
                    {
                        const ulong p = i ? i-1 : size-1;
                        const ulong n = (i == size-1) ? 0 : i+1;
                        if(bounds[i].z < 0. && bounds[p].z >= 0. && bounds[n].z >= 0.)
                        {
                            const T dist = bounds[p].z / (bounds[p].z - bounds[i].z);
                            Point temp1 = (bounds[i] - bounds[p]) * dist + bounds[p];
                            temp1.z = 0.;
                            temp1.normalize();
                            
                            bounds[i] = (bounds[i] - bounds[n]) * dist + bounds[n];
                            bounds[i].z = 0.;
                            bounds[i].normalize();
                            bounds.insert(bounds.begin()+i, temp1);
                            size++;
                            i += 3;
                        }
                        else if(bounds[i].z < 0. && bounds[p].z >= 0.)
                        {
                            const T dist = bounds[p].z / (bounds[p].z - bounds[i].z);
                            Point temp = (bounds[i] - bounds[p]) * dist + bounds[p];
                            temp.z = 0.;
                            temp.normalize();
                            bounds.insert(bounds.begin()+i, temp);
                            size++;
                            i += 2;
                        }
                        else if(bounds[i].z < 0. && bounds[n].z >= 0.)
                        {
                            const T dist = bounds[n].z / (bounds[n].z - bounds[i].z);
                            Point temp = (bounds[i] - bounds[n]) * dist + bounds[n];
                            temp.z = 0.;
                            temp.normalize();
                            bounds.insert(bounds.begin()+n, temp);
                            size++;
                            i += 2;
                        }
                        else
                        {
                            i++;
                        }
                    }
                    size = bounds.size();
                    for(ulong i = 0; i < size;)
                    {
                        const ulong p = i ? i-1 : size-1;
                        const ulong n = (i == size-1) ? 0 : i+1;
                        if(bounds[i].z <= 0. && bounds[p].z <= 0. && bounds[n].z <= 0.)
                        {
                            bounds.erase(bounds.begin()+i);
                            size--;
                        }
                        else
                        {
                            i++;
                        }
                    }
                }
            }
            
            static bool compareAzimuth(Point const& p1, Point const& p2) noexcept
            {
                return p1.azimuth() < p2.azimuth();
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
            
            Triangle(Point const& _a, Point const& _b, Point const& _c) :
            a(_a), b(_b), c(_c)
            {
                const Point ac = (c - a);
                const Point ab = (b - a);
                const Point t = ab.cross(ac);
                const T _d = (2. * t.lenght2());
                if(_d > HOA_EPSILON)
                {
                    p = (((t.cross(ab) * ac.lenght2()) + (ac.cross(t) * ab.lenght2())) / _d + a);
                    if(p.distance(Point(0., 0., 0)) > HOA_EPSILON)
                    {
                        p.normalize();
                        r = p.distance(a);
                    }
                    else
                    {
                        r = 0.;
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
        
    private:
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
            m_points.push_back(p);
            m_points[m_points.size()-1].normalize();
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
                                    if(t.p.distance(m_points[l]) < t.r - HOA_EPSILON)
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
            for(ulong i = 0; i < m_points.size(); i++)
            {
                const T el = HOA_PI2 - m_points[i].elevation();
                const T az = m_points[i].azimuth();
                for(ulong j = 0; j < m_points[i].bounds.size(); j++)
                {
                    m_points[i].bounds[j].rotateZ(-az);
                    m_points[i].bounds[j].rotateX(el);
                }
                sort(m_points[i].bounds.begin(), m_points[i].bounds.end(), Point::compareAzimuth);
                for(ulong j = 0; j < m_points[i].bounds.size(); j++)
                {
                    m_points[i].bounds[j].rotateX(-el);
                    m_points[i].bounds[j].rotateZ(az);
                }
            }
        }
    };

#undef HOA_EPSILON
}

#endif



