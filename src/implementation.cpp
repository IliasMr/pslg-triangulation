#include "../includes/functions.h"

// Returns the converge rate for the saved obtuse counts
double calculate_convergence_rate(const std::vector<int>& obtuse_counts, int N) {
    if (N <= 1) return -1;            // too small N

    double total_p = 0.0;
    for (int n = 0; n < N - 1; ++n) {
        if (obtuse_counts[n] == 0) continue;        
        if (obtuse_counts[n + 1] == 0) break;        
        double p_n = log(static_cast<double>(obtuse_counts[n + 1]) / obtuse_counts[n]) / log(static_cast<double>(n + 1) / n);
        total_p += p_n;
    }
    return total_p / (N - 1);
}

// Returns the Energy
double calculate_energy(int obtuse_count, int steiner_count, double alpha = 4.5, double beta = 1.0) {
    return (1 + alpha * obtuse_count + beta * steiner_count);
}

// returns the degrees of an angle
double calculate_angle(const Vector& v1, const Vector& v2) {
    // dot product of vectors
    double dot_product = CGAL::to_double(v1 * v2);
    
    // magnitudes of the vectors
    double magnitude_v1 = std::sqrt(CGAL::to_double(v1.squared_length()));  //using double for std::sqrt
    double magnitude_v2 = std::sqrt(CGAL::to_double(v2.squared_length()));
    
    // angle in radians
    double angle_radians = std::acos(dot_product / (magnitude_v1 * magnitude_v2));
    
    // convert to degrees 
    double angle_degrees = angle_radians * (180.0 / M_PI);
    return angle_degrees;
}

//returns the number of obtuse angles in a cdt 
int calculate_obtuse_angles(const CDT& cdt, const Polygon_2& region_boundary_Pol){

    int obtCounter = 0;
    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        Face_handle face = fit;

        // Get the vertices of the face
        Point_2 p1 = face->vertex(0)->point();
        Point_2 p2 = face->vertex(1)->point();
        Point_2 p3 = face->vertex(2)->point();

        // go to the next loop if triangle is outside the region boundary
        if (is_triangle_outside(region_boundary_Pol, p1,p2,p3) == true){
            continue;
        }

        // Calculate vectors for each pair of vertices
        Vector v1 = p2 - p1;
        Vector v2 = p3 - p1;
        Vector v3 = p3 - p2;

        // Calculate angles
        double angle1 = calculate_angle(v1, v2);
        double angle2 = calculate_angle(-v1, v3);
        double angle3 = calculate_angle(-v2, -v3);

        // Check if any angle is obtuse (> 90 degrees)
        if (angle1 > 90.0 || angle2 > 90.0 || angle3 > 90.0) {
            obtCounter++;
        }
    }

    return obtCounter;
}



// Returns if point A is inside the triangle B,C,D
int in_triangle(const Point_2& A, const Point_2& B, const Point_2& C, const Point_2& D){

    // Calculate the cross product of the vectors
    Kernel::FT cross1, cross2, cross3;
    cross1 = (B.x()-A.x())*(D.y()-A.y()) - (B.y()-A.y())*(D.x()-A.x());
    cross2 = (C.x()-B.x())*(D.y()-B.y()) - (C.y()-B.y())*(D.x()-B.x());
    cross3 = (A.x()-C.x())*(D.y()-C.y()) - (A.y()-C.y())*(D.x()-C.x());

    if ((cross1 >= 0 && cross2 >= 0 && cross3 >= 0) || (cross1 <= 0 && cross2 <= 0 && cross3 <= 0)){
        return 1;
    } else {
        return 0;
    }
}

// Returns a cdt initialized with the input 
MyCDT json_to_cdt(jsonInfo& input, std::vector<Point_2>& reg_bound_points) {
    // Create a constrained Delaunay triangulation
    MyCDT cdt;

    // Insert the points into the triangulation
    for (const auto& point : input.points){
        cdt.insert(point);
    }

    for(const auto& constr : input.additional_constraints){ 
        cdt.insert_constraint(input.points[constr.first], input.points[constr.second]);
    }

    
    // add region boundary as extra constrains in cdt
    for (size_t i = 0; i < input.region_boundary.size(); ++i) {
        
        // get the indices of the current point and the next point
        int curr_index = input.region_boundary[i];
        int next_index = input.region_boundary[(i + 1) % input.region_boundary.size()];


        cdt.insert_constraint(input.points[curr_index], input.points[next_index]);
        reg_bound_points.push_back(input.points[curr_index]);
    }

    return cdt;
}


// Returns the intersection point between perpendicular projection of point A and the line BC
Point_2 perpendicularProjection(Point_2 A, Point_2 B, Point_2 C) {
    Point_2 projection;
    
    if (B.y() == C.y()) {
        projection = { A.x(), B.y() };  
    }
    else if (B.x() == C.x()) {
        projection = { B.x(), A.y() }; 
    }
    else {
        // slope of BC 
        Kernel::FT m_BC = (C.y() - B.y()) / (C.x() - B.x());
    
        //slope of the perpendicular line
        Kernel::FT m_perpendicular = -1.0 / m_BC;
    
        // find the intersection point
        Kernel::FT x_projection = (m_BC * B.x() - m_perpendicular * A.x() + A.y() - B.y()) / (m_BC - m_perpendicular);
        Kernel::FT y_projection = m_BC * (x_projection - B.x()) + B.y();
    
        projection = { x_projection, y_projection };
    }

    return projection;
}


// Returns true if the triangle is out of bounds 
bool is_triangle_outside(const Polygon_2& region_boundary, const Point_2& p1, const Point_2& p2, const Point_2& p3) {

    // if centroid is out of bounds the triangle will also be out of bounds
    Point_2 centroid = CGAL::centroid(p1,p2,p3);
    CGAL::Bounded_side res = CGAL::bounded_side_2(region_boundary.vertices_begin(), 
    region_boundary.vertices_end(), centroid, Kernel());

    if (res == CGAL::ON_UNBOUNDED_SIDE)
        return true;
    
    return false;
}


// returns true if face is obtuse 
bool has_obtuse_angle(Face_handle face){
    
    Point_2 p1 = face->vertex(0)->point();
    Point_2 p2 = face->vertex(1)->point();
    Point_2 p3 = face->vertex(2)->point();

    // Calculate vectors for each pair
    Vector v1 = p2 - p1;
    Vector v2 = p3 - p1;
    Vector v3 = p3 - p2;

    // calculate angles 
    double angle1 = calculate_angle(v1, v2);
    double angle2 = calculate_angle(-v1, v3);
    double angle3 = calculate_angle(-v2, -v3);

    if (angle1 > 90.0 || angle2 > 90.0 || angle3 > 90.0)
        return true;
    
    return false;
}


//Returns the centroid of the convex polygon formed by adjacent obtuse triangles
//and null if there no such triangles in the cdt 
Point_2* find_polygon_centroid(MyCDT& cdt, Face_handle face) {

    // Check if the input triangle is valid and non-constrained
    for (int i = 0; i < 3; ++i) {
        if (cdt.are_there_incident_constraints(face->vertex(i))) {
            return nullptr;
        }
    }

    // Initialize polygon points with the current face's vertices
    std::vector<Point_2> polygon_points = {
        face->vertex(0)->point(),
        face->vertex(1)->point(),
        face->vertex(2)->point()
    };

    // Process each neighbor
    for (int i = 0; i < 3; ++i) {
        Face_handle neighbor_face = face->neighbor(i);

        // Skip infinite neighbors or non-obtuse neighbors
        if (cdt.is_infinite(neighbor_face) || !has_obtuse_angle(neighbor_face)) {
            continue;
        }

        // Skip neighbors with constraints
        bool has_constraints = false;
        for (int j = 0; j < 3; ++j) {
            if (cdt.are_there_incident_constraints(neighbor_face->vertex(j))) {
                has_constraints = true;
                break;
            }
        }
        if (has_constraints) {
            continue;
        }

        // Add unique points from the neighbor triangle
        for (int j = 0; j < 3; ++j) {
            Point_2 neighbor_point = neighbor_face->vertex(j)->point();
            if (std::find(polygon_points.begin(), polygon_points.end(), neighbor_point) == polygon_points.end()) {
                polygon_points.push_back(neighbor_point);
            }
        }
    }

    // Check if the polygon is convex
    if (!CGAL::is_convex_2(polygon_points.begin(), polygon_points.end())) {
        return nullptr;
    }

    // Calculate and return the centroid of the convex polygon
    Point_2* polygon_centroid = new Point_2(CGAL::centroid(polygon_points.begin(), polygon_points.end()));
    return polygon_centroid;
}



// inserts the cirumcenter in the cdt 
void insert_circumcenter(MyCDT& cdt, Point_2 circumcenter, Face_handle original_face){

    int flag = 0;
    int constraints = 0;
    Face_handle neighbor;
    Point_2 v1,v2,v3;
    Vertex_handle to_remove1, to_remove2;
    
    // Initialize polygon points with original face points 
    std::vector<Point_2> polygon_points = {
        original_face->vertex(0)->point(),
        original_face->vertex(1)->point(),
        original_face->vertex(2)->point()
    };

    //check for incident constraints in the original triangle
    for (int j = 0; j < 3; ++j) {
        if (cdt.are_there_incident_constraints(original_face->vertex(j))) {
            //cdt.insert(CGAL::centroid(cdt.triangle(original_face)));
            return;                                         //immediately return if there are constraints - circumcenter skipped
        }
    }

    //find the neighbor where circumcenter will be 
    for (int i = 0; i < 3; i++) {
        neighbor = original_face->neighbor(i);

        v1 = neighbor->vertex(0)->point();
        v2 = neighbor->vertex(1)->point();
        v3 = neighbor->vertex(2)->point();

        //check for incident constraints in the neighbor
        for (int j = 0; j < 3; ++j) {
            if (cdt.are_there_incident_constraints(neighbor->vertex(j))) {
               // cdt.insert(CGAL::centroid(cdt.triangle(original_face)));
                return;                                         //immediately return if there are constraints - circumcenter skipped
            }
        }

        // check if cirumcenter is inside this neighbor
        CGAL::Triangle_2<Kernel> neighbor_triangle(v1, v2, v3);
        if (neighbor_triangle.bounded_side(circumcenter) == CGAL::ON_BOUNDED_SIDE) {
            flag = 1;
            to_remove1 = original_face->vertex((i + 1) % 3);               // first point of the shared edge
            to_remove2 = original_face->vertex((i + 2) % 3);               // second point of the shared edge
            break;
        }

    }

    // cicrumcenter is not on a neighbor of the triangle 
    if (flag == 0){
        //cdt.insert(CGAL::centroid(cdt.triangle(original_face)));
        return;
    
    }

    // edw exoume eksasfalisei oti to circumcenter einai se geitoniko trigwno 
    // Add the 4th point to the polygon points
    for (int j = 0; j < 3; ++j) {
        Point_2 neighbor_point = neighbor->vertex(j)->point();
        if (std::find(polygon_points.begin(), polygon_points.end(), neighbor_point) == polygon_points.end()) {
            polygon_points.push_back(neighbor_point);
            break;                                                                  // we found the 4th point
        }
    }


    //do not insert if the polygon is not convex
    if (!CGAL::is_convex_2(polygon_points.begin(), polygon_points.end())) {
        //cdt.insert(CGAL::centroid(cdt.triangle(original_face)));
        return;
    }


    std::cout << "polygon sz->" << polygon_points.size()<<std::endl;
    // remove the vertexes from the cdt
    cdt.remove_no_flip(to_remove1);
    cdt.remove_no_flip(to_remove2);

    //CGAL::draw(cdt);

    //add polygons edges as constraints
    for (size_t i = 0; i <=polygon_points.size(); ++i) {
        Point_2 a = polygon_points[i];
        Point_2 b = polygon_points[(i + 1) % polygon_points.size()];

        // Insert the edge as a constraint
        cdt.insert_constraint(a, b);
    }


    //insert circumcenter without fliping
    cdt.insert_no_flip(circumcenter);

    //CGAL::draw(cdt);

    //remove any added constraints
    for (int k=0; k<3; k++){

        //remove added constraints on the original face
        if (cdt.is_constrained(CDT::Edge(original_face, k)))
            //cdt.remove_constraint(original_face, k);
            cdt.remove_constraint_no_flip(original_face, k);

        
        if (cdt.is_constrained(CDT::Edge(neighbor, k)))
            //cdt.remove_constraint(neighbor, k);
            cdt.remove_constraint_no_flip(original_face, k);
    }

}




// returns the type of the given input
int find_input_type(MyCDT cdt, const Polygon_2& region_boundary_Pol, const jsonInfo& input, const std::vector<Point_2>& reg_bound_points){

    // Category A
    if (region_boundary_Pol.is_convex()){  
        if (input.num_constraints == 0)
            return 0;                // 0 -> Type A
    
        // Category Γ - input that containts closed cycles
        Graph constr_graph;                                     //initialize contstraints graph - contains region boundary also   
        for (int k=0; k<input.constr_points.size()-1; k+=2){
            constr_graph.add_edge(input.constr_points[k], input.constr_points[k + 1]);
        }

        //add region boundary points 
        for (size_t k = 0; k < reg_bound_points.size(); ++k) {  
            constr_graph.add_edge(reg_bound_points[k], reg_bound_points[(k + 1) % reg_bound_points.size()]);
        }

        if (constr_graph.has_cycles()) 
            return 2;                       // 2 -> Type Γ
        else 
            return 1;                      // 1 -> Type B


    }

    //Category Δ
    if (input.num_constraints == 0){
        bool only_axis_aligned = true;
        for (auto edge = region_boundary_Pol.edges_begin(); edge != region_boundary_Pol.edges_end(); ++edge) {
            if (!((edge->source().x() == edge->target().x()) || (edge->source().y() == edge->target().y()))) {
                only_axis_aligned = false;
                break;
            }
        }
        if (only_axis_aligned)
            return 3;               //3 -> Type Δ

    }

    //if none of the above categories return type E
    return 4;              

}

void insert_random_point(MyCDT& cdt, CDT::Finite_faces_iterator fit){
    // gaussian distribution initialization
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);

    Point_2 centroid = CGAL::centroid(cdt.triangle(fit));

    // insert the point around the centroid using the gaussian distribution
    Point_2 random_point(centroid.x() + dist(gen), centroid.y() + dist(gen));

    // insert the point in the cdt
    cdt.insert(random_point);
}

// writes output to a json file in json format
void print_output(std::string outputFile, std::vector<Point_2> initial_points, std::vector<Point_2> steiner_points,
 MyCDT cdt, int obtuse, const jsonInfo& input, const Polygon_2& region_boundary_Pol, bool randomPoint) {

    // open output file 
    std::ofstream outFile(outputFile);
    if (!outFile) {
        std::cerr << "Error opening the file.." << std::endl;
        return;
    }
    outFile << "{\n\t\"content_type\": \"CG_SHOP_2025_Solution\",";
    outFile << "\n\t\"instance_uid\": \"" << input.instance_uid << "\",";
    // print Steiner points x
    outFile << "\n\t\"steiner_points_x\": [";
    if(steiner_points.empty())
        outFile << "],";            //close bracket if no steiner points were added

    for (const auto& steiner : steiner_points) {
        outFile << "\"";
        const auto exact_coord = CGAL::exact(steiner.x());
        outFile << exact_coord.get_num() << "/" << exact_coord.get_den();
        if (&steiner != &steiner_points.back()) {
            outFile << "\", ";
        } else {
            outFile << "\"],";
        }
    }
    // print Steiner points y
    outFile << "\n\t\"steiner_points_y\": [";
    if(steiner_points.empty())
        outFile << "],";   
        
    for (const auto& steiner : steiner_points) {
        outFile << "\"";
        const auto exact_coord = CGAL::exact(steiner.y());
        outFile << exact_coord.get_num() << "/" << exact_coord.get_den();
        if (&steiner != &steiner_points.back()) {
            outFile << "\", ";
        } else {
            outFile << "\"],";
        }
    }
    outFile << "\n\t\"edges\": [\n";

    // save edges and print them

    std::vector<Point_2> all_points = initial_points;
    all_points.insert(all_points.end(), steiner_points.begin(), steiner_points.end());

    std::vector<std::pair<int, int>> edges;
    auto boundary_begin = region_boundary_Pol.vertices_begin();
    auto boundary_end = region_boundary_Pol.vertices_end(); 
    for (auto edge_it = cdt.finite_edges_begin(); edge_it != cdt.finite_edges_end(); ++edge_it) {
        CDT::Face_handle face = edge_it->first;
        int index = edge_it->second;
        CDT::Vertex_handle vh1 = face->vertex((index + 1) % 3);
        CDT::Vertex_handle vh2 = face->vertex((index + 2) % 3);

        Point_2 p1 = vh1->point();
        Point_2 p2 = vh2->point();

        // Find indices of vh1->point() and vh2->point() in steiner_points
        int idx1 = std::distance(all_points.begin(), std::find(all_points.begin(), all_points.end(), p1));
        int idx2 = std::distance(all_points.begin(), std::find(all_points.begin(), all_points.end(), p2));
        
        Point_2 midpoint = CGAL::midpoint(p1,p2);
        CGAL::Bounded_side res = CGAL::bounded_side_2(boundary_begin, boundary_end, midpoint, Kernel());

        //skip if edge is out of bounds 
        if (res == CGAL::ON_UNBOUNDED_SIDE){
            continue;
        }
      
        if(idx1 != idx2)
            edges.emplace_back(idx1, idx2);

    }

    for (size_t i = 0; i < edges.size(); ++i) {
        outFile << "\t\t[" << edges[i].first << ", " << edges[i].second << "]";
        if (i != edges.size() - 1) 
            outFile << ",\n";
         else
            outFile << "\n";
        
    }

    outFile << "\t],";
    outFile << "\n\t\"obtuse_count\": " << obtuse << ",";
    outFile << "\n\t\"method\": \"" << input.method << "\",";
    
    //print parameters for the chosen method
    outFile << "\n\t\"parameters\": {";

    if (input.method == "local"){
        outFile << "\n\t\t\"L\": " << input.L;
    }
    else if (input.method == "sa") {
        outFile << "\n\t\t\"alpha\": " << input.alpha << ",";
        outFile << "\n\t\t\"beta\": " << input.beta << ",";
        outFile << "\n\t\t\"L\": " << input.L;
    } else if (input.method == "ant") {
        outFile << "\n\t\t\"alpha\": " << input.alpha << ",";
        outFile << "\n\t\t\"beta\": " << input.beta << ",";
        outFile << "\n\t\t\"xi\": " << input.xi << ",";
        outFile << "\n\t\t\"psi\": " << input.psi << ",";
        outFile << "\n\t\t\"lambda\": " << input.lambda << ",";
        outFile << "\n\t\t\"kappa\": " << input.kappa << ",";
        outFile << "\n\t\t\"L\": " << input.L;
    }

    outFile << "\n\t},";
    outFile << "\n\t\"randomization\": " << (randomPoint ? "true" : "false");
    outFile << "\n}\n";
    outFile.close();

}
