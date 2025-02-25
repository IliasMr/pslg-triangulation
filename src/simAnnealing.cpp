#include "../includes/functions.h"
#include <chrono>               // for time calculation 


// Simulated Annealing method
MyCDT sim_Annealing(MyCDT& cdt, double alpha, double beta, int L,const Polygon_2& region_boundary_Pol,
 std::vector<Point_2>& steiner_points, std::vector<int>& obtuse_counts, bool randomPoint){

    int obtuse_before = calculate_obtuse_angles(cdt,region_boundary_Pol);
    double Energy = alpha * obtuse_before;                          // no steiner points in the beggining
    double Energy_new, deltaEnergy;
    double T = 1.0;
    int test_obt_angles;
    int steiner_count=0;


    double Tsub = 1/(double)L;
    auto boundary_begin = region_boundary_Pol.vertices_begin();
    auto boundary_end = region_boundary_Pol.vertices_end(); 

    //initialize random values
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr_int(0,4);
    std::uniform_real_distribution<> distr_real(0.1, 1.0);
    // Gaussian (normal) distribution for centroid random steiner
    std::normal_distribution<> gaussian(0.0, 0.15);             

    MyCDT testCdt;
    int itcounter = 0;
    std::vector<double> Pconv; 

    // for calculating run time
    auto start_time = std::chrono::high_resolution_clock::now();
    const auto time_limit = std::chrono::minutes(5);

    int no_improvement = 0;
    while (T >= 0 && itcounter <= L * 2){

        // if elapsed time >5 minutes stop 
        auto current_time = std::chrono::high_resolution_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::minutes>(current_time - start_time);
        if (elapsed_time >= time_limit) {
            break;
        }

        //for each obtuse triangle in the CDT 
        for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit){
            
            // get the vertices of each point
            Point_2 p1 = fit->vertex(0)->point();
            Point_2 p2 = fit->vertex(1)->point();
            Point_2 p3 = fit->vertex(2)->point();

            // if the triangle is outside the region boundary skip it 
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

            int rand_num;
            // check for obtuse angle in the current face
            if (angle1 > 90.0 || angle2 > 90.0 || angle3 > 90.0) {
                            
                //select a random steiner point option    
                //int rand_num = rand() % 5;
                rand_num = distr_int(gen);
                Point_2 steiner;

                if (rand_num == 0){
                    steiner = CGAL::centroid(cdt.triangle(fit));
                }
                else if (rand_num == 1){
                    steiner = CGAL::circumcenter(p1,p2,p3);
                    //check if circumcenter is inside the region boundary pol
                    CGAL::Bounded_side res = CGAL::bounded_side_2(boundary_begin, boundary_end, steiner, Kernel());
                    if (res == CGAL::ON_UNBOUNDED_SIDE){
                        steiner = CGAL::centroid(cdt.triangle(fit));        //use centroid if circumcenter is outside the region boundary
                    }
                }
                else if (rand_num == 2){
                    if (angle1 > 90.0) 
                        steiner = perpendicularProjection(p1, p2, p3);
                    else if (angle2 > 90.0) 
                        steiner = perpendicularProjection(p2, p1, p3);
                    else if (angle3 > 90.0) 
                        steiner = perpendicularProjection(p3, p2, p1);      
                }
                else if (rand_num == 3){
                    if (angle1 > 90.0) 
                        steiner = CGAL::midpoint(p2,p3);
                    else if (angle2 > 90.0) 
                        steiner = CGAL::midpoint(p1,p3);
                    else if (angle3 > 90.0) 
                        steiner = CGAL::midpoint(p1,p2);
                }
                else{ 
                    Point_2* polygon_centroid = find_polygon_centroid(cdt, fit);
                    if (polygon_centroid != nullptr)
                        steiner = *polygon_centroid;
                    else
                        steiner = CGAL::centroid(cdt.triangle(fit));                         //if polygon centroid does not exists we use centroid
                }
                //create temp cdt for calculating new energy
                testCdt = cdt;
                testCdt.insert(steiner);
                test_obt_angles = calculate_obtuse_angles(testCdt, region_boundary_Pol);

                //calculate new Energy and ΔΕ
                Energy_new = alpha * test_obt_angles + beta * (steiner_count+1);
                deltaEnergy = Energy_new - Energy;

                //accept or not the new configuration
                if (deltaEnergy < 0){
                    //accept the new configuration
                    cdt.insert(steiner);
                    obtuse_counts.push_back(calculate_obtuse_angles(cdt, region_boundary_Pol));
                    steiner_points.push_back(steiner);
                    Energy = Energy_new;
                    steiner_count++;
                    // if (rand_num == 0)
                    //     std::cout << "centroid" << std::endl;
                    // else if (rand_num == 1)
                    //     std::cout << "circum" << std::endl;
                    // else if (rand_num == 2)
                    //     std::cout << "proj" << std::endl;
                    // else if (rand_num == 3)
                    //     std::cout << "midpoint" << std::endl;
                    // else if (rand_num == 4)
                    //     std::cout << "pol_centr" << std::endl;

                    no_improvement = 0;
                    break;
                }
                else {
                    double accept_prob = std::exp(-deltaEnergy / T);
                    double rand_prob = distr_real(gen);   
                    //accept the new configuration if random number is lower than e^-ΔΕ/Τ 
                    if (rand_prob <= accept_prob){
                        cdt.insert(steiner);
                        obtuse_counts.push_back(calculate_obtuse_angles(cdt, region_boundary_Pol));
                        steiner_points.push_back(steiner);
                        Energy = Energy_new;
                        steiner_count++;
                    
                        // if (rand_num == 0)
                        //     std::cout << "centroid" << std::endl;
                        // else if (rand_num == 1)
                        //     std::cout << "circum" << std::endl;
                        // else if (rand_num == 2)
                        //     std::cout << "proj" << std::endl;
                        // else if (rand_num == 3)
                        //     std::cout << "midpoint" << std::endl;
                        // else if (rand_num == 4)
                        //     std::cout << "pol_centr" << std::endl;
                        no_improvement = 0;
                        break;
                    }

                    //at this point no steiner point was added 
                    no_improvement++;                                       //increase everytime no steiner points were added in the cdt 

                    //add random "escape" steiner when obtuse_count was not reduced for L+100 iterations
                    if (randomPoint == true){

                        if (no_improvement > 500) {            
                            Point_2 centroid = CGAL::centroid(cdt.triangle(fit));
                            
                            // find offset
                            Vector gaussian_offset(gaussian(gen), gaussian(gen));

                            // add offset to the centroid
                            Point_2 escape_steiner = centroid + gaussian_offset;

                            CGAL::Bounded_side res_centroid = region_boundary_Pol.bounded_side(escape_steiner); 

                            // add the steiner if its inside the region boundary
                            if (res_centroid != CGAL::ON_UNBOUNDED_SIDE){
                                cdt.insert(escape_steiner);
                                steiner_points.push_back(escape_steiner);
                                steiner_count++;
                                obtuse_counts.push_back(calculate_obtuse_angles(cdt, region_boundary_Pol));
                                Energy = alpha * calculate_obtuse_angles(cdt, region_boundary_Pol) + beta * steiner_count;
                                no_improvement = 0;
                                break;
                            }
                        }
                    }
                    //std::cout << "no_improvement" << no_improvement << std::endl;
                }
            }
        }
        //decrease temperature 
        T -= Tsub;
        itcounter++;
       
        //std::cout << "T" << T << std::endl;
        //std::cout << "itcounter ->" <<itcounter <<std::endl;
    }

    return cdt;
}




