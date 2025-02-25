#include "../includes/functions.h"

int recursion_counter_ant = 0;

// Ant Colony Optimization Method
MyCDT antColony(MyCDT& cdt, const Polygon_2& region_boundary_Pol, 
                std::vector<Point_2>& steiner_points, double alpha, double beta, 
                double xi, double psi, double lambda, int kappa, int L, std::vector<int>& obtuse_counts, bool randomPoint) {

    // Initialize pheromone vector
    std::vector<double> pheromone(4, 1.0); // Initial pheromone values

    // To track selected faces per cycle
    std::vector<std::tuple<CDT::Finite_faces_iterator, int, Point_2, bool, int>> triangles_selected;

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> random_prob(0.0, 1.0);

    int obtuse_before = calculate_obtuse_angles(cdt,region_boundary_Pol);

    // Main loop for cycles
    for (int cycle = 0; cycle < L; ++cycle) {

        triangles_selected.clear();

        // Each ant (agent)
        for (int ant = 0; ant < kappa; ++ant) {

            MyCDT testCdt = cdt;

            // Select triangle
            CDT::Finite_faces_iterator triangle = select_triangle(cdt, random_prob(gen), region_boundary_Pol);

            // Calculate Steiner point
            std::pair<Point_2, int> selected_steiner = select_steiner_point(cdt, triangle, random_prob(gen), region_boundary_Pol, pheromone, xi, psi);

            //Test new configuration
            // if(selected_steiner.second == 1){
            //     insert_circumcenter(testCdt, selected_steiner.first, triangle);
            // } else
            testCdt.insert(selected_steiner.first);
            int obtuse_after = calculate_obtuse_angles(testCdt, region_boundary_Pol);

            if(obtuse_after < obtuse_before)
                ImproveTriangulation(triangles_selected, triangle, selected_steiner, obtuse_after);
        }

        // If no steiner points were added, attempt random point insertion
        if (triangles_selected.size() == 0 && randomPoint) {
            recursion_counter_ant++;
            if (recursion_counter_ant >= 5) {
                randomPoint = false;
            }

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> random_prob(0.0, 1.0);

            // Select a random obtuse triangle
            CDT::Finite_faces_iterator fit = select_triangle(cdt, random_prob(gen), region_boundary_Pol);

            MyCDT testCdt = cdt; // Create a temporary copy of the current CDT
            int obtuse_before_rand = calculate_obtuse_angles(testCdt, region_boundary_Pol);

            // Insert a random point into the test CDT
            insert_random_point(testCdt, fit);

            // Perform ant colony on the test CDT
            std::vector<Point_2> test_steiner_points;
            std::vector<int> test_obtuse_counts;

            testCdt = antColony(testCdt, region_boundary_Pol, test_steiner_points, alpha, beta, xi, psi, lambda, kappa, L, test_obtuse_counts, randomPoint);

            // Recalculate obtuse angles for the ant optimized CDT
            int obt_angles_after_rand = calculate_obtuse_angles(testCdt, region_boundary_Pol);

            // If the optimized CDT is better, update the original CDT
            if (obt_angles_after_rand < obtuse_before_rand) {
                recursion_counter_ant = 0;
                cdt = testCdt;
                steiner_points.insert(steiner_points.end(), test_steiner_points.begin(), test_steiner_points.end());
                obtuse_counts.insert(obtuse_counts.end(), test_obtuse_counts.begin(), test_obtuse_counts.end());
            }
            return cdt;
        }

        SaveBestTriangulation(cdt, triangles_selected, steiner_points, obtuse_before, region_boundary_Pol, obtuse_counts);
        UpdatePheromones(pheromone, triangles_selected, lambda, alpha, beta, steiner_points.size());
    }

    return cdt;
}

void SaveBestTriangulation(MyCDT &cdt, std::vector<std::tuple<CDT::Finite_faces_iterator, int, Point_2, bool, int>>& triangles_selected,
 std::vector<Point_2>& steiner_points, int& obtuse_before, const Polygon_2& region_boundary_Pol, std::vector<int>& obtuse_counts) {
    for(const auto& triangle : triangles_selected){
        if(std::get<3>(triangle) == false){
            continue;
        }
        Point_2 steiner = std::get<2>(triangle);
        int type = std::get<4>(triangle);
        // if(type == 1){
        //     insert_circumcenter(cdt, steiner, std::get<0>(triangle));
        //     steiner_points.push_back(steiner);
        //     continue;
        // }
        cdt.insert(steiner);
        obtuse_counts.push_back(std::get<1>(triangle));
        steiner_points.push_back(steiner);
    }
    obtuse_before = calculate_obtuse_angles(cdt, region_boundary_Pol);
}

void ImproveTriangulation(std::vector<std::tuple<CDT::Finite_faces_iterator, int, Point_2, bool, int>>& triangles_selected, CDT::Finite_faces_iterator triangle, std::pair<Point_2, int> selected_steiner, int obtuse_after) {
    // Find if the selected triangle is already in the list
    auto it = std::find_if(triangles_selected.begin(), triangles_selected.end(),
                           [&triangle](const std::tuple<CDT::Finite_faces_iterator, int, Point_2, bool, int>& temp) {
                               return (std::get<0>(temp) == triangle && std::get<3>(temp) == true);
                           });
    // If the triangle is already in the list, pick the one with the lowest obtuse angles
    if (it != triangles_selected.end()) {
        if (obtuse_after < std::get<1>(*it)) {
            std::get<3>(*it) = false;
            triangles_selected.emplace_back(triangle, obtuse_after, selected_steiner.first, true, selected_steiner.second);
        } else {
            triangles_selected.emplace_back(triangle, obtuse_after, selected_steiner.first, false, selected_steiner.second);
        }
    } else {
        // Else, the triangle is not in the list, add it
        triangles_selected.emplace_back(triangle, obtuse_after, selected_steiner.first, true, selected_steiner.second);
    }
}

void UpdatePheromones(std::vector<double>& pheromone, 
                      const std::vector<std::tuple<CDT::Finite_faces_iterator, int, Point_2, bool, int>>& selected_points, 
                      double lambda, double alpha, double beta, int steiner_count) {
    
    std::vector<double> delta_pheromone(4, 0.0);
    
    // Iterate through all Steiner points
    for (const auto& triangle : selected_points) {
        int type = std::get<4>(triangle);

        int obtuse_count = std::get<1>(triangle);
        delta_pheromone[type] += 1/(1 + alpha * obtuse_count + beta * steiner_count);
    }

    for (int i = 0; i < 4; ++i) {
        pheromone[i] = (1 - lambda) * pheromone[i] + delta_pheromone[i];
    }
}


// Computes heuristic value based on radius-to-height ratio
std::vector<double> compute_heuristic(double ratio) {
    std::vector<double> heuristics(4, 0.0);

    heuristics[0] = std::max(0.0, (3-2*ratio)/3);
    heuristics[1] = ratio/(2+ratio);
    heuristics[2] = std::max(0.0, (ratio-1)/ratio);
    heuristics[3] = 1.0;

    return heuristics;
}

// returns the height_to_radius ratio for the triangle
double height_to_radius_ratio(MyCDT& cdt, CDT::Finite_faces_iterator triangle) {
    return circumradius(cdt.triangle(triangle)) * longest_edge(cdt.triangle(triangle)) / (2.0 * triangle_area(cdt.triangle(triangle)));
}

double calculate_total(const std::vector<double>& pheromone, std::vector<double> heur_value, double x, double y) {
    double case1 = std::pow(pheromone[0], x) * std::pow(heur_value[0], y);
    double case2 = std::pow(pheromone[1], x) * std::pow(heur_value[1], y);
    double case3 = std::pow(pheromone[2], x) * std::pow(heur_value[2], y);
    double case4 = std::pow(pheromone[3], x) * std::pow(heur_value[3], y);

    return case1 + case2 + case3 + case4;
}

// Select Steiner point based on random choice and triangle properties
std::pair<Point_2, int> select_steiner_point(MyCDT& cdt, CDT::Finite_faces_iterator fit, double rand_val, const Polygon_2& region_boundary_Pol, std::vector<double> pheromone, double x, double y) {

    // Psp calculation
    double ratio = height_to_radius_ratio(cdt, fit);

    std::vector<double> heur_value = compute_heuristic(ratio);

    double total = calculate_total(pheromone, heur_value, x, y);

    // Probability of selecting each type of Steiner point option
    double Pmid = std::pow(pheromone[0], x) * std::pow(heur_value[0], y) / total;
    double Pcirc = Pmid + std::pow(pheromone[1], x) * std::pow(heur_value[1], y) / total;
    double Pproj = Pcirc + std::pow(pheromone[2], x) * std::pow(heur_value[2], y) / total;

    Point_2* pol_cen = find_polygon_centroid(cdt, fit);
    if(pol_cen == nullptr){
        heur_value[3] = 0.0;
    }
    double Ppolcen = Pproj + std::pow(pheromone[3], x) * std::pow(heur_value[3], y) / total;
    
    Point_2 p1 = fit->vertex(0)->point();
    Point_2 p2 = fit->vertex(1)->point();
    Point_2 p3 = fit->vertex(2)->point();

    // Calculate vectors for each pair
    Vector v1 = p2 - p1;
    Vector v2 = p3 - p1;
    Vector v3 = p3 - p2;

    // calculate angles 
    double angle1 = calculate_angle(v1, v2);
    double angle2 = calculate_angle(-v1, v3);
    double angle3 = calculate_angle(-v2, -v3);

    int type;

    Point_2 steiner;
    if (rand_val < Pmid) {
        if (angle1 > 90.0) 
            steiner = CGAL::midpoint(p2,p3);
        else if (angle2 > 90.0) 
            steiner = CGAL::midpoint(p1,p3);
        else if (angle3 > 90.0) 
            steiner = CGAL::midpoint(p1,p2);
        type = 0;
    } else if (rand_val < Pcirc) { 
        steiner = CGAL::circumcenter(p1, p2, p3);
        CGAL::Bounded_side res = CGAL::bounded_side_2(region_boundary_Pol.vertices_begin(),
                                                      region_boundary_Pol.vertices_end(), 
                                                      steiner, Kernel());
        if (res == CGAL::ON_UNBOUNDED_SIDE) {
            steiner = CGAL::centroid(p1, p2, p3);
        }
        type = 1;
    } else if (rand_val < Pproj) {
        if (angle1 > 90.0) 
            steiner = perpendicularProjection(p1, p2, p3);
        else if (angle2 > 90.0) 
            steiner = perpendicularProjection(p2, p1, p3);
        else if (angle3 > 90.0) 
            steiner = perpendicularProjection(p3, p2, p1);
        type = 2;
    } else {
        Point_2* polygon_centroid = find_polygon_centroid(cdt, fit);
        steiner = (polygon_centroid != nullptr) ? *polygon_centroid : CGAL::centroid(p1, p2, p3);
        type = 3;
    }
    return {steiner, type};
}

// selects a random triangle and returns an iterator to it 
CDT::Finite_faces_iterator select_triangle(MyCDT& cdt, double rand_val, const Polygon_2& region_boundary_Pol) {
    std::vector<double> probabilities;

    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        if (has_obtuse_angle(fit)) {

            if(is_triangle_outside(region_boundary_Pol, fit->vertex(0)->point(), fit->vertex(1)->point(), fit->vertex(2)->point())){
                 probabilities.push_back(0.0);
            }
            else{
                probabilities.push_back(1.0);
            }
        } else {
            probabilities.push_back(0.0); // Assign probability 0.0 otherwise
        }
    }

    // Create a discrete distribution based on the probabilities
    std::discrete_distribution<> triangle_distribution(probabilities.begin(), probabilities.end());
    
    // Convert rand_val to an integer seed for the random engine
    std::mt19937 generator(static_cast<unsigned int>(rand_val * std::numeric_limits<unsigned int>::max()));
    
    // Select a triangle ID based on the distribution
    int selected_triangle_id = triangle_distribution(generator);

    // Get the iterator to the selected triangle
    CDT::Finite_faces_iterator selected_triangle = cdt.finite_faces_begin();
    std::advance(selected_triangle, selected_triangle_id);

    return selected_triangle;
}

// returns the circumradius of a triangle
double circumradius(const CDT::Triangle& triangle) {
    Point_2 p1 = triangle[0], p2 = triangle[1], p3 = triangle[2];

    double a2 = CGAL::to_double(CGAL::squared_distance(p2, p3));
    double b2 = CGAL::to_double(CGAL::squared_distance(p1, p3));
    double c2 = CGAL::to_double(CGAL::squared_distance(p1, p2));

    double a = CGAL::sqrt(a2), b = CGAL::sqrt(b2), c = CGAL::sqrt(c2);
    double s = (a + b + c) / 2;
    double area = CGAL::sqrt(s * (s - a) * (s - b) * (s - c));

    return (a * b * c) / (4 * area);
}

// returns the longest edge of a triangle
double longest_edge(const CDT::Triangle& triangle) {
    Point_2 p1 = triangle[0], p2 = triangle[1], p3 = triangle[2];

    double length1 = CGAL::to_double(CGAL::squared_distance(p1, p2));
    double length2 = CGAL::to_double(CGAL::squared_distance(p2, p3));
    double length3 = CGAL::to_double(CGAL::squared_distance(p3, p1));

    return CGAL::sqrt(std::max({length1, length2, length3}));
}

// returns the area of a triangle
double triangle_area(const CDT::Triangle& triangle) {
    Point_2 p1 = triangle[0], p2 = triangle[1], p3 = triangle[2];

    auto determinant = CGAL::to_double(p1.x() * (p2.y() - p3.y()) + 
                                       p2.x() * (p3.y() - p1.y()) + 
                                       p3.x() * (p1.y() - p2.y()));

    return std::abs(determinant) / 2.0;
}