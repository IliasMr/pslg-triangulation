#include "../includes/includes.h"

// Class to store the json information from the input
class jsonInfo{
    public:
        int num_points;
        std::string instance_uid;
        std::vector<Point_2> points;
        std::vector<int> region_boundary;
        int num_constraints;
        std::vector<std::pair<int, int>> additional_constraints;
        std::string method;
        bool delaunay;
        std::vector<Point_2> constr_points;                     // points that are vertexes of constrained edges

        //parameters for the different methods 
        double alpha, beta, xi, psi, lambda;
        int kappa, L;



        // Constructor
        jsonInfo(const std::string& filename){
            //load json file
            boost::property_tree::ptree pt;
            boost::property_tree::read_json(filename, pt);

            //extract instance_uid
            instance_uid = pt.get<std::string>("instance_uid");

            //extract num_points
            num_points = pt.get<int>("num_points");

            //extract point_x / point_y
            auto points_x = pt.get_child("points_x");
            auto points_y = pt.get_child("points_y");
            for (auto x_it = points_x.begin(), y_it = points_y.begin(); x_it != points_x.end() && y_it != points_y.end(); ++x_it, ++y_it) {
                int x = x_it->second.get_value<int>();
                int y = y_it->second.get_value<int>();
                points.emplace_back(x, y);
            }
            
            //extract region boundary
            for (auto& item : pt.get_child("region_boundary")) {
                region_boundary.push_back(item.second.get_value<int>());
            }

            //extract num_constraints
            num_constraints = pt.get<int>("num_constraints");

            //extract additional constraints
            for (const auto& item : pt.get_child("additional_constraints")) {
                auto it = item.second.begin();                           // constr_x
                int start = it->second.get_value<int>();  
                ++it;  
                int end = it->second.get_value<int>();                  // constr_y
                additional_constraints.emplace_back(start, end);
                constr_points.push_back(points[start]);
                constr_points.push_back(points[end]);
            }

            method = pt.get<std::string>("method"); 
            delaunay = pt.get<bool>("delaunay");
            
            //read specific parameters for each method
            boost::property_tree::ptree params = pt.get_child("parameters");
            if (method == "ant"){
                alpha = params.get<double>("alpha");
                beta = params.get<double>("beta");
                xi = params.get<double>("xi");
                psi = params.get<double>("psi");
                lambda = params.get<double>("lambda");
                kappa = params.get<int>("kappa");
                L = params.get<int>("L");
            }
            else if (method == "local"){
                L = params.get<int>("L");
            }
            else if (method == "sa" || method == "auto"){
                alpha = params.get<double>("alpha", 3.5);       //default values for auto method with no parameters in the json 
                beta = params.get<double>("beta", 1.0);
                L = params.get<int>("L");
            }



        }
};



// Class to represent the constraints graph and detect any possible cycles
class Graph {
public:
    
    // ads edge p1,p2 in the Graph
    void add_edge(const Point_2& p1, const Point_2& p2) {
        adj_list[p1].insert(p2);
        adj_list[p2].insert(p1);
    }

    // returns true if there is more than 1 cycle in the Graph 
    bool has_cycles() {
        std::set<Point_2> visited_nodes;
        std::set<std::pair<Point_2, Point_2>> visited_edges;
        std::set<std::set<Point_2>> unique_cycles; // unique cycles

        // search all nodes in the adj list
        for (const auto& pair : adj_list) {
            const Point_2& node = pair.first;

            // dfs for every not visited node
            if (visited_nodes.find(node) == visited_nodes.end()) {  
                std::vector<Point_2> current_path;
                dfs(node, Point_2(), visited_nodes, visited_edges, unique_cycles, current_path);
            }
        }
        
        return unique_cycles.size() > 1; 
    }

private:
    std::map<Point_2, std::set<Point_2>> adj_list;

    // dfs search for cycles
    void dfs(const Point_2& node, const Point_2& parent,
             std::set<Point_2>& visited_nodes,
             std::set<std::pair<Point_2, Point_2>>& visited_edges,
             std::set<std::set<Point_2>>& unique_cycles,
             std::vector<Point_2>& current_path) {
        
        visited_nodes.insert(node);             //to hold visited edges
        current_path.push_back(node);           //to hold path 

        for (const auto& neighbor : adj_list[node]) {
            std::pair<Point_2, Point_2> edge = std::make_pair(std::min(node, neighbor), std::max(node, neighbor));  //sorted edge
            
            if (visited_edges.find(edge) != visited_edges.end()) {
                continue;
            }
            
            visited_edges.insert(edge);

            if (visited_nodes.find(neighbor) == visited_nodes.end()) {
                dfs(neighbor, node, visited_nodes, visited_edges, unique_cycles, current_path);
            } 
            else if (neighbor != parent) {
                
                // found new cycle
                auto cycle_start = std::find(current_path.begin(), current_path.end(), neighbor);
                if (cycle_start != current_path.end()) {
                    
                    // save vertices that form the cycle
                    std::set<Point_2> cycle_vertices(cycle_start, current_path.end());
                    unique_cycles.insert(cycle_vertices);
                }
            }
        }
        
        current_path.pop_back();            //remove last path node
    }
};