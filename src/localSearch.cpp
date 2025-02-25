#include "../includes/functions.h"

int recursion_counter_ls = 0;

MyCDT local_search(MyCDT& cdt, int L, Polygon_2 region_boundary_Pol, std::vector<Point_2>& steiner_points, std::vector<int>& obtuse_counts, bool randomPoint) {
    int obt_angles_before = calculate_obtuse_angles(cdt, region_boundary_Pol);
    int obt_angles_after = calculate_obtuse_angles(cdt, region_boundary_Pol);
    int obt_angles_before_for = obt_angles_before;

    //L for max iterations
    //If there are no obtuse triangles we break out of the loop
    for(int i = 0; i < L; i++){
        obt_angles_before = calculate_obtuse_angles(cdt, region_boundary_Pol);
        obt_angles_before_for = obt_angles_before;
        for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit){
            // get the vertices of each point
            Point_2 p1 = fit->vertex(0)->point();
            Point_2 p2 = fit->vertex(1)->point();
            Point_2 p3 = fit->vertex(2)->point();

            // skip the triangle if its out of bounds
            if (is_triangle_outside(region_boundary_Pol, p1,p2,p3) == true){
                continue;
            }

            // Calculate vectors for each pair
            Vector v1 = p2 - p1;
            Vector v2 = p3 - p1;
            Vector v3 = p3 - p2;

            // calculate angles 
            double angle1 = calculate_angle(v1, v2);
            double angle2 = calculate_angle(-v1, v3);
            double angle3 = calculate_angle(-v2, -v3);

            // check for obtuse angle in the current face
            if (angle1 > 90.0 || angle2 > 90.0 || angle3 > 90.0) {
                int obt_Centr = 10000; int obt_Mid = 10000; int obt_Proj = 10000; int obt_Circ = 10000; int obt_Pol_Centr = 10000;

                //calculate centroid
                Point_2 centroid = CGAL::centroid(cdt.triangle(fit));

                //calculate midpoint and perpendicular projection
                Point_2 midpoint, projection;                                       
                if (angle1 > 90.0) {
                    midpoint = CGAL::midpoint(p2,p3);
                    projection = perpendicularProjection(p1, p2, p3);
                } else if (angle2 > 90.0) {
                    midpoint = CGAL::midpoint(p1,p3);
                    projection = perpendicularProjection(p2, p1, p3);
                } else if (angle3 > 90.0) {
                    midpoint = CGAL::midpoint(p1,p2);
                    projection = perpendicularProjection(p3, p2, p1);
                }

                //find triangles circumcenter
                Point_2 circumcenter = CGAL::circumcenter(p1,p2,p3);
                Point_2* polygon_centroid = find_polygon_centroid(cdt, fit);
                MyCDT testCdt;
                testCdt = cdt;

                testCdt.insert(centroid);
                obt_Centr = calculate_obtuse_angles(testCdt, region_boundary_Pol);
                
                testCdt = cdt;
                testCdt.insert(midpoint);
                obt_Mid = calculate_obtuse_angles(testCdt, region_boundary_Pol);

                testCdt = cdt;
                testCdt.insert(projection);
                obt_Proj = calculate_obtuse_angles(testCdt, region_boundary_Pol);

                testCdt = cdt;
                testCdt.insert(circumcenter);
                obt_Circ = calculate_obtuse_angles(testCdt, region_boundary_Pol);

                if (polygon_centroid != nullptr){
                    testCdt = cdt;
                    testCdt.insert(*polygon_centroid);
                    obt_Pol_Centr = calculate_obtuse_angles(testCdt, region_boundary_Pol);                    
                } else{
                    obt_Pol_Centr = 10000;
                }

                //find the steiner with the best results
                int obtAngles = std::min({obt_angles_before, obt_Centr, obt_Mid, obt_Proj, obt_Circ, obt_Pol_Centr});

                //insert the point with the best result to the main cdt
                //immediately break and proceed to the next while iteration
                if(obtAngles == obt_angles_before){                        
                    continue;
                } else if (obtAngles == obt_Centr) {
                    cdt.insert(centroid);
                    obtuse_counts.push_back(obtAngles);
                    steiner_points.push_back(centroid);
                    obt_angles_before = obtAngles;
                    if (polygon_centroid != nullptr)              //delete polygon_centroid as we wont use it
                        delete polygon_centroid;
                    break;
                } else if (obtAngles == obt_Mid) {
                    cdt.insert(midpoint);
                    obtuse_counts.push_back(obtAngles);
                    steiner_points.push_back(midpoint);
                    obt_angles_before = obtAngles;
                    if (polygon_centroid != nullptr)
                        delete polygon_centroid;
                    break;

                } else if (obtAngles == obt_Proj) {
                    cdt.insert(projection);
                    obtuse_counts.push_back(obtAngles);
                    steiner_points.push_back(projection);
                    obt_angles_before = obtAngles;
                    if (polygon_centroid != nullptr)
                        delete polygon_centroid;
                    break;

                } else if (obtAngles == obt_Circ){

                    //check if circumcenter is inside the region boundary
                    CGAL::Bounded_side result = CGAL::bounded_side_2(region_boundary_Pol.vertices_begin(),
                    region_boundary_Pol.vertices_end(), circumcenter, Kernel());

                    if (result == CGAL::ON_BOUNDED_SIDE || result == CGAL::ON_BOUNDARY){
                        cdt.insert(circumcenter);
                        obtuse_counts.push_back(obtAngles);
                        steiner_points.push_back(circumcenter);
                        obt_angles_before = obtAngles;
                        if (polygon_centroid != nullptr)
                            delete polygon_centroid;
                        break;
                    }

                } else if (obtAngles == obt_Pol_Centr){
                    cdt.insert(*polygon_centroid);
                    obtuse_counts.push_back(obtAngles);
                    steiner_points.push_back(*polygon_centroid);
                    obt_angles_before = obtAngles;
                    delete polygon_centroid;
                    break;
                }
            }
        }
        obt_angles_after = calculate_obtuse_angles(cdt, region_boundary_Pol);

        // End conditions
        // Break if no obtuse angles were found 
        if (obt_angles_after == 0){
            break;
        }

        // Break if obtuse angles are the same as before, meaning no progress was made (stopping criterion)
        if(obt_angles_after >= obt_angles_before_for && randomPoint){
            recursion_counter_ls++;
            if(recursion_counter_ls >= 5){
                randomPoint = false;
            }

            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<double> dist(0.0, 1.0);

            // Select a random obtuse triangle
            CDT::Finite_faces_iterator fit = select_triangle(cdt, dist(gen), region_boundary_Pol);

            MyCDT testCdt = cdt; // Create a copy of the current CDT
            int obt_angles_before_rand = calculate_obtuse_angles(testCdt, region_boundary_Pol);

            // Insert a random point and update the test CDT
            insert_random_point(testCdt, fit);

            std::vector<Point_2> test_steiner_points;
            std::vector<int> test_obtuse_counts;

            // Perform local search on the test CDT
            testCdt = local_search(testCdt, L, region_boundary_Pol, test_steiner_points, test_obtuse_counts, randomPoint);

            // Recalculate obtuse angles for the optimized CDT
            int obt_angles_after_rand = calculate_obtuse_angles(testCdt, region_boundary_Pol);

            // If the optimized CDT is better, update the original CDT and parameters
            if (obt_angles_after_rand < obt_angles_before_rand) {
                recursion_counter_ls = 0;
                cdt = testCdt;
                steiner_points.insert(steiner_points.end(), test_steiner_points.begin(), test_steiner_points.end());
                obtuse_counts.insert(obtuse_counts.end(), test_obtuse_counts.begin(), test_obtuse_counts.end());
            }
            return cdt;
        }
        
        if(obt_angles_after >= obt_angles_before_for){
            break;
        }
    }
    return cdt;
}
