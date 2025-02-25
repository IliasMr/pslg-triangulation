// header files, libraries and type definitions we use

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include "Custom_Constrained_Delaunay_triangulation_2.h"
#include <CGAL/centroid.h>
#include <CGAL/draw_polygon_2.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <CGAL/Surface_mesh.h>
#include <cstdlib>
#include <ctime>
#include <stack>
#include <set>
#include <map>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Vector_2 Vector;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel> CDT;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Face_handle Face_handle;
typedef Custom_Constrained_Delaunay_triangulation_2<Kernel> MyCDT;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Line_2 Line_2;
typedef CDT::Edge Edge;