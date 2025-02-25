#include "../includes/classes.h"

// the prototypes of the function we use

double calculate_angle(const Vector& v1, const Vector& v2);
int calculate_obtuse_angles(const CDT& cdt,const Polygon_2& region_boundary_Pol );
int in_triangle(const Point_2& A, const Point_2& B, const Point_2& C, const Point_2& D);
MyCDT json_to_cdt(jsonInfo& input, std::vector<Point_2>& reg_bound_points);
Point_2 perpendicularProjection(Point_2 A, Point_2 B, Point_2 C);
bool is_triangle_outside(const Polygon_2& region_boundary, const Point_2& p1, const Point_2& p2, const Point_2& p3);
Point_2* find_polygon_centroid(MyCDT& cdt, Face_handle face);
double circumradius(const CDT::Triangle& triangle);
double longest_edge(const CDT::Triangle& triangle);
double triangle_area(const CDT::Triangle& triangle);
bool has_obtuse_angle(Face_handle face);
void print_output(std::string outputFile, std::vector<Point_2> initial_points, std::vector<Point_2> steiner_points, MyCDT cdt, int obtuse, const jsonInfo& input, const Polygon_2& region_boundary_Pol, bool randomPoint);
void insert_circumcenter(MyCDT& cdt, Point_2 circumcenter, Face_handle original_face);
std::pair<Point_2, int> select_steiner_point(MyCDT& cdt, CDT::Finite_faces_iterator fit, double rand_val, const Polygon_2& region_boundary_Pol, std::vector<double> pheromone, double x, double y);
std::vector<double> compute_heuristic(double ratio);
void UpdatePheromones(std::vector<double>& pheromone, const std::vector<std::tuple<CDT::Finite_faces_iterator, int, Point_2, bool, int>>& selected_points, double lambda, double alpha, double beta, int steiner_count);
CDT::Finite_faces_iterator select_triangle(MyCDT& cdt, double rand_val, const Polygon_2& region_boundary_Pol);
void SaveBestTriangulation(MyCDT &cdt, std::vector<std::tuple<CDT::Finite_faces_iterator, int, Point_2, bool, int>>& triangles_selected, std::vector<Point_2>& steiner_points, int& obtuse_before, const Polygon_2& region_boundary_Pol, std::vector<int>& obtuse_counts);
void ImproveTriangulation(std::vector<std::tuple<CDT::Finite_faces_iterator, int, Point_2, bool, int>>& triangles_selected, CDT::Finite_faces_iterator triangle, std::pair<Point_2, int> selected_steiner, int obtuse_after);
int find_input_type(MyCDT cdt, const Polygon_2& region_boundary_Pol, const jsonInfo& input,const std::vector<Point_2>&  reg_bound_points);
double calculate_convergence_rate(const std::vector<int>& obtuse_counts, int N);
double calculate_energy(int obtuse_count, int steiner_count, double alpha, double beta);
void insert_random_point(MyCDT& cdt, CDT::Finite_faces_iterator fit);


// Functions for the different algorithms
MyCDT sim_Annealing(MyCDT& cdt, double alpha, double beta, int L, const Polygon_2& region_boundary_Pol, std::vector<Point_2>& steiner_points, std::vector<int>& obtuse_counts, bool randomPoint);
MyCDT local_search(MyCDT& cdt, int L, Polygon_2 region_boundary_Pol, std::vector<Point_2>& steiner_points, std::vector<int>& obtuse_counts, bool randomPoint);
MyCDT antColony(MyCDT& cdt, const Polygon_2& region_boundary_Pol, std::vector<Point_2>& steiner_points, double alpha, double beta, double xi, double psi, double lambda, int kappa, int L, std::vector<int>& obtuse_counts, bool randomPoint);