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
    template <Dimension D, typename T> class Voronoi;
    
    template <typename T> class Voronoi <Hoa3d, T>
    {
        
    public:
        struct Point
        {
            T x;
            T y;
            T z;
            vector<Point> neightbours;
            vector<Point> bounds;
            
            Point() noexcept :
            x(0), y(0.), z(0.)
            {
                
            }
            
            ~Point() noexcept
            {
                neightbours.clear();
                bounds.clear();
            }
            
            Point(T _x, T _y, T _z, ulong index = 0) noexcept :
            x(_x), y(_y), z(_z)
            {
                ;
            }
            
            Point(Point const& other) noexcept :
            x(other.x), y(other.y), z(other.z)
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
            
            void filterBounds()
            {
                const T el = HOA_PI2 - elevation();
                const T az = azimuth();
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
                        const T dist = bounds[p].z / (bounds[p].z - bounds[i].z);
                        Point temp((bounds[i].x - bounds[p].x) * dist + bounds[p].x, (bounds[i].y - bounds[p].y) * dist + bounds[p].y, 0.);
                        temp.normalize();
                        bounds.insert(bounds.begin()+(i), temp);
                        size++;
                    }
                    else if(bounds[i].z < 0. && bounds[n].z > 0.)
                    {
                        const T dist = bounds[n].z / (bounds[n].z - bounds[i].z);
                        Point temp((bounds[i].x - bounds[n].x) * dist + bounds[n].x, (bounds[i].y - bounds[n].y) * dist + bounds[n].y, 0.);
                        temp.normalize();
                        bounds.insert(bounds.begin()+(i+1), temp);
                        size++;
                    }
                }
                for(ulong i = 0; i < size; i++)
                {
                    if(bounds[i].z < 0.)
                    {
                        bounds.erase(bounds.begin()+i);
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
            T  r;
            
            Triangle(Point const& _a, Point const& _b, Point const& _c) noexcept :
            a(_a), b(_b), c(_c), r(0.)
            {
                const Point ac = (c - a);
                const Point ab = (b - a);
                const Point t = ab.cross(ac);
                const T _d = (2. * t.lenght2());
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

#endif



