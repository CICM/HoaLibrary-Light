/*
// Copyright (c) 2012-2015 Eliott Paris, Julien Colafrancesco & Pierre Guillot, CICM, Universite Paris 8.
// For information on usage and redistribution, and for a DISCLAIMER OF ALL
// WARRANTIES, see the file, "LICENSE.txt," in this distribution.
*/

#ifndef DEF_HOA_VORONOI_LIGHT
#define DEF_HOA_VORONOI_LIGHT

#include "Math.hpp"

namespace hoa
{
#define HOA_EPSILON 1e-15;
    
    template <typename T> class Triangle;
    template <typename T> class Arc;
    
    template <typename T> class Point
    {
    public:
        ulong   i;
        T       x;
        T       y;
        T       z;
        Arc<T>* arc;
        
        Point() noexcept : i(0), x(0), y(0), z(0), arc(NULL)
        {
            ;
        }
        
        Point(const T _x, const T _y, const T _z) noexcept : i(0), x(_x), y(_y), z(_z), arc(NULL)
        {
            ;
        }
        
        static T dot(Point const&a, Point const&b) noexcept
        {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }
        
        static Point subtract(Point const& a, Point const& b) noexcept
        {
            return Point(a.x - b.x, a.y - b.y, a.z - b.z);
        }
        
        static Point cross(Point const& a, Point const& b) noexcept
        {
            return Point(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
        }
        
        static T norm2(Point const& d) noexcept
        {
            return sqrt(dot(d, d));
        }
        
        static T norm(Point const& d) noexcept
        {
            return sqrt(d.norm2());
        }
        
        static Point normalise(Point const& d)  noexcept
        {
            T m = 1 / d.norm();
            return Point(d.x * m, d.y * m, d.z * m);
        }
        
    };
    
    template <typename T> class Edge
    {
    public:
        Triangle<T>* t;
        Point<T>     p;
        Point<T>     neighbor;
        Point<T>     next;
        
        Edge(Triangle<T>* _t, const Point<T>& _p) noexcept :
        t(_t),
        p(_p)
        {
            
        }
        
        static void neighbors(Edge const& a, Edge const&  b)
        {
            (a.neighbor = b).neighbor = a;
        }
    };
    
    template <typename T> class Arc
    {
    public:
        Triangle<T>* t;
        Point<T>     p;
        ulong        index;
        Point<T>*    prev;
        Point<T>*    next;
   
        Arc(Triangle<T>* _t, const Point<T>& _p, const ulong _index) noexcept :
        t(_t),
        p(_p),
        index(_index),
        prev(NULL),
        next(NULL)
        {
            if(p.visible)
            {
                next = p.visible;
            }
            //int zaza;
            
            p.arc = this;
        }
    };
    
    template <typename T> class Triangle
    {
    public:
        ulong       index;
        bool        marked;
        Edge<T>     a;
        Edge<T>     b;
        Edge<T>     c;
        Point<T>    n;
        vector<Arc<T> > arcs;
        
        Triangle(const Point<T> _p1, const Point<T> _p2, const Point<T> _p3, const ulong _index = 0) noexcept :
        index(_index),
        marked(false)
        {
            a = Edge<T>(this, _p1);
            b = Edge<T>(this, _p2);
            c = Edge<T>(this, _p3);
            a.next = _p2;
            b.next = _p3;
            c.next = _p1;
            n = Point<T>::normalise(Point<T>::cross(Point<T>::subtract(_p3, _p1), Point<T>::subtract(_p2, _p1)));
        }
        
        static bool coplanar(Triangle const& t, Point<T> const& p)
        {
            return abs(Point<T>::dot(t.n, p) - Point<T>::dot(t.n, t.a.p)) <= HOA_EPSILON;
        }
        
        static bool visible(Triangle const& t, Point<T> const& p)
        {
            return Point<T>::dot(t.n, p) - Point<T>::dot(t.n, t.a.p) > HOA_EPSILON;
        }
        
        static void addConflict(Triangle const&t, Point<T> const&p, const ulong _index)
        {
            if(visible(t, p))
            {
                t.arcs.push_back(Arc<T>(t, p, _index));
            }
        }


    };
    
    template <Dimension D, typename T> class Voronoi;
    
    template <typename T> class Voronoi<Hoa3d, T>
    {
    private:
        
        void convexhull3d(vector<Point<T> >& points, vector<Triangle<T> >& triangles)
        {
            if(points.size() < 4)
            {
                triangles.clear();
                return;
            }
            
            for(ulong i = 0; i < points.size(); ++i)
            {
                points[i].i = i;
            }
            //int zaza;
            //d3.shuffle(points);
            
            Point<T> a = points[0];
            Point<T> b = points[1];
            Point<T> c = points[2];
            Triangle<T> t(a, b, c);
            
            // Find non-coplanar fourth point.
            ulong i;
            for(i = 3; i < points.size() && coplanar(t, points[i]); ++i)
            {
                ;
            }
            
            if(i == points.size())
            {
                triangles.clear();
                return; // coplanar points
            }
            
            // Create a tetrahedron.
            Point<T> d = points[i];
            points[i] = points[3];
            points[3] = d;
            
            if(visible(t, d))
            {
                Point<T> tmp = b;
                b = c;
                c = tmp;
            }
            
            Triangle<T> ta(a, b, c, 0);
            Triangle<T> tb(d, b, a, 1);
            Triangle<T> tc(c, d, a, 2);
            Triangle<T> td(b, d, c, 3);
            triangles[0] = ta;
            triangles[1] = tb;
            triangles[2] = tc;
            triangles[3] = td;
            
            Edge<T>::neighbors(ta.a, tb.b);
            Edge<T>::neighbors(ta.b, td.c);
            Edge<T>::neighbors(ta.c, tc.c);
            
            Edge<T>::neighbors(tb.a, td.a);
            Edge<T>::neighbors(td.b, tc.a);
            Edge<T>::neighbors(tc.b, tb.c);
            
            // Initialise conflict graph.
            for(ulong i = 4; i < points.size(); ++i)
            {
                Point<T> p = points[i];
                addConflict(ta, p, i);
                addConflict(tb, p, i);
                addConflict(tc, p, i);
                addConflict(td, p, i);
            }
            /*
            for(ulong i = 4; i < points.size(); ++i)
            {
                Point<T> p = points[i];
                Point<T> h = p.visible;
                if(!h)
                {
                    continue;
                }
                
                // Find horizon.
                Point<T> horizon;
                Point<T> a = h;
                do a.t.marked = true; while (a = a.nextF);
                
                a = h; do {
                    var t = a.t;
                    if (horizon = findHorizon(t.a) || findHorizon(t.b) || findHorizon(t.c)) break;
                } while (a = a.nextF);
                
                if (!horizon) continue;
                
                for (var j = 0, m = horizon.length, prev = null, first = null; j < m; ++j) {
                    var e = horizon[j],
                    f1 = e.triangle, f2 = e.neighbor.triangle,
                    t = new Triangle(p, e.neighbor.p, e.p, triangles.length);
                    neighbors(t.b, e);
                    if (prev) neighbors(prev.a, t.c);
                    else first = t;
                    addConflicts(t, f1, f2);
                    triangles.push(prev = t);
                }
                neighbors(prev.a, first.c);
                
                a = h; do {
                    var t = a.t;
                    for (var j = 0, m = t.visible.length; j < m; ++j) t.visible[j].remove();
                    t.visible.length = 0;
                    removeElement(triangles, t.index);
                } while (a = a.nextF);
            }
             */
        }
        
        void delaunay(vector<Point<T> >& points)
        {
            /*
            var p = points.map(cartesian),
            n = points.length,
            triangles = d3.convexhull3d(p);
            
            if (triangles.length) return triangles.forEach(function(t)
            {
                t.coordinates = [points[t.a.p.i], points[t.b.p.i], points[t.c.p.i]];
                t.centre = circumcentre(t);
            }), triangles;
             */
        };
        
    public:
    };
}

#endif


