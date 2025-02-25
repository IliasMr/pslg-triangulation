#include "../includes/functions.h"


int main(int argc, char** argv) {

    if (argc < 5 || argc > 6) {  
        std::cerr << "Run the program like this:\n\t./opt_triangulation -i /path/to/input.json -o /path/to/output.json [-preselected_params]" << std::endl;
        return 1;  
    }

    std::string inputFile;
    std::string outputFile;
    bool preselectedParams = false;             // true when flag -preselected_params is present

    //arguments check
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-i" && i + 1 < argc) {
            if (std::string(argv[i + 1]) == "-o") {
                std::cerr << "Run the program like this:\n\t./opt_triangulation -i /path/to/input.json -o /path/to/output.json [-preselected_params]" << std::endl;
                return 1;
            }
            inputFile = argv[i + 1];
            ++i;
        } else if (std::string(argv[i]) == "-o" && i + 1 < argc) {
            if (std::string(argv[i + 1]) == "-i") {
                std::cerr << "Run the program like this:\n\t./opt_triangulation -i /path/to/input.json -o /path/to/output.json [-preselected_params]" << std::endl;
                return 1;
            }
            outputFile = argv[i + 1];
            ++i;
        } else if (std::string(argv[i]) == "-preselected_params") {
            preselectedParams = true;
        }
    }

    jsonInfo input(inputFile);
    std::vector<Point_2> reg_bound_points;
    MyCDT cdt = json_to_cdt(input, reg_bound_points);
    //CGAL::draw(cdt);
    
    MyCDT final;
    int obtuse, obtuse_before;
    std::vector<Point_2> steiner_points;                        //to save added steiner points per method
    std::vector<Point_2> steiner_points_total;                        //to save total steiner points added

    // Apply preselected params if any 
    if (preselectedParams) {
        input.alpha = 4.5;
        input.beta = 1.0;
        input.xi = 1.0;
        input.psi = 0.5;
        input.lambda = 0.5;
        input.kappa = 15;
        input.L = 300;
    }

    //initial points for the output
    std::vector<Point_2> initial_points;
    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
        initial_points.push_back(vit->point());
    }
    Polygon_2 region_boundary_Pol(reg_bound_points.begin(), reg_bound_points.end());

    obtuse_before = calculate_obtuse_angles(cdt, region_boundary_Pol);
  
    // find the input category
    int type = find_input_type(cdt, region_boundary_Pol, input, reg_bound_points);
    if (type == 0)
        std::cout << "Input type: A->CONVEX_NO_CONSTRAINTS" << std::endl;
    else if (type == 1)
        std::cout << "Input type: B->CONVEX_NO_CIRCLE_CONSTRAINTS" << std::endl;
    else if (type == 2)
        std::cout << "Input type: C->CONVEX_CIRCLE_CONSTRAINTS" << std::endl;
    else if (type == 3)
        std::cout << "Input type: D->NON_CONVEX_ORTHO" << std::endl;
    else if (type == 4)
        std::cout << "Input type: E->NON_CONVEX_UNSPECIFIED" << std::endl;
 
    bool randomPoint = false;           // change to run with randomization

    std::vector<int> obtuse_counts;
    if(input.method == "local"){
        final = local_search(cdt, input.L, region_boundary_Pol, steiner_points, obtuse_counts, randomPoint);
        obtuse = calculate_obtuse_angles(final, region_boundary_Pol);
        steiner_points_total.insert(steiner_points_total.end(),steiner_points.begin(), steiner_points.end());
        input.alpha = 4.5;          //hardcoded for energy calculation forn local method
        input.beta = 1.0;
    }

    if (input.method == "sa" || input.method == "auto"){
        //check for delaunay - if false we start with cdt returned from local search (hw1)
        MyCDT local = cdt;
        if(input.delaunay == false && input.method != "auto"){
            local = local_search(cdt, input.L, region_boundary_Pol, steiner_points, obtuse_counts, randomPoint);
            steiner_points_total.insert(steiner_points_total.end(),steiner_points.begin(), steiner_points.end());
            steiner_points.clear();
            final = sim_Annealing(local, input.alpha, input.beta, input.L, region_boundary_Pol, steiner_points, obtuse_counts, randomPoint);
            steiner_points_total.insert(steiner_points_total.end(),steiner_points.begin(), steiner_points.end());
        }
        else{
            final = sim_Annealing(cdt, input.alpha, input.beta, input.L, region_boundary_Pol, steiner_points, obtuse_counts, randomPoint);
            steiner_points_total.insert(steiner_points_total.end(),steiner_points.begin(), steiner_points.end());
        }
        obtuse = calculate_obtuse_angles(final, region_boundary_Pol);
    }

    if (input.method == "ant"){
        //check for delaunay - if false we start with cdt returned from local search (hw1)
        MyCDT local = cdt;
        if (input.delaunay == false){
            local = local_search(cdt, input.L, region_boundary_Pol, steiner_points, obtuse_counts, randomPoint);
            
            steiner_points_total.insert(steiner_points_total.end(),steiner_points.begin(), steiner_points.end());
            steiner_points.clear();
            
            final = antColony(local, region_boundary_Pol, steiner_points, input.alpha, input.beta, input.xi, input.psi, input.lambda, input.kappa, input.L, obtuse_counts, randomPoint);
            steiner_points_total.insert(steiner_points_total.end(),steiner_points.begin(), steiner_points.end());

        }
        else {
            final = antColony(cdt, region_boundary_Pol, steiner_points, input.alpha, input.beta, input.xi, input.psi, input.lambda, input.kappa, input.L, obtuse_counts, randomPoint);        
            steiner_points_total.insert(steiner_points_total.end(),steiner_points.begin(), steiner_points.end());
        }
        obtuse = calculate_obtuse_angles(final, region_boundary_Pol);
    }


    // solution info
    std::cout << "Obtuse triangles in cdt before: "<< obtuse_before << std::endl;
    std::cout << "Obtuse triangles in cdt after: "<< obtuse << std::endl;
    std::cout << "Steiner points number: " << steiner_points_total.size() << std::endl;

    //calculate and print conv rate / energy
    double conv_rate = abs(calculate_convergence_rate(obtuse_counts, steiner_points_total.size()));
    double energy = calculate_energy(obtuse, steiner_points_total.size(), input.alpha, input.beta);
    if(obtuse == 0)
        std::cout << "Convergance rate: " << conv_rate << std::endl;
    else
        std::cout << "Energy: " << energy << std::endl;

    print_output(outputFile, initial_points, steiner_points_total, final, obtuse, input, region_boundary_Pol, randomPoint);
    //CGAL::draw(final);

    return 0;
}