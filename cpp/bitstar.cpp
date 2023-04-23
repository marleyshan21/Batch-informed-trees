#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <random>
#include <vector>
#include <queue>
#include <bitstar.h>

// #include <opencv2/opencv.hpp>
int main()
{
    Bit_star *Bit_star; 
    return 0;
}


double Bit_star::gt(Node node)
{
    double length = 0.0;

    if (node.x == start.x && node.y == start.y)
    {
        return 0.0;
    }
    //check:  assumed all the nodes in vert are conected?
    if (std::find(this->vert.begin(), this->vert.end(), node) == this->vert.end())
    {
        return std::numeric_limits<double>::infinity();
    }
    double length = 0.0;
    Node *current = &node;
    Node *parent = current->parent;
    while(parent->x != start.x && parent->y != start.y)
    {
        // what is weight here - c_hat between current and parent?
        length = length + c_hat(*current, *parent);
        current = parent;
        parent = current->parent;
    }
    return length;

}


double Bit_star::c_hat(Node node1, Node node2)
{
    return sqrt(pow(node1.x - node2.x, 2) + pow(node1.y - node2.y, 2));
}



std::vector<Node> Bit_star::near(Node node, std::vector<Node> search_list)
{
    std::vector<Node> nearest;

    for (int i = 0; i < search_list.size(); i++)
    {

        if ((c_hat(node, search_list[i]) < Rbit) && (search_list[i].x != node.x && search_list[i].y != node.y))
        {
            nearest.push_back(search_list[i]);
        }
        
    }


    
    return nearest;
}


std::vector<double> Bit_star::sample_unit_ball(int d)
{
    
    // uniform random sample from unit ball
    std::vector<double> u(d);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);



    for (int i = 0; i < d; i++) {
        u[i] = dis(gen);
    }

    double norm = 0.0;
    for (int i = 0; i < d; i++) {
        norm += u[i] * u[i];
    }
    norm = std::sqrt(norm);

    std::uniform_real_distribution<double> dis_r(0.0, 1.0);
    double r = dis_r(gen);

    std::vector<double> x(d);
    for (int i = 0; i < d; i++) {
        x[i] = r * u[i] / norm;
    }

    return x;

}


std::vector<Node> Bit_star::samplePHS()
{
   

    Eigen::Vector2d one_1(1.0, 0.0);
    double cmin = sqrt(pow(goal.x - start.x, 2) + pow(goal.y - start.y, 2));
    Eigen::Vector2d a1{ (goal.x - start.x) / cmin, (goal.y - start.y) / cmin };
    Eigen::Vector2d center { (start.x + goal.x) / 2, (start.y + goal.y) / 2 };

    // SVD
    Eigen::Matrix2d U, S, Vt;
    Eigen::Matrix2d a1_outer_one_1 = a1 * one_1.transpose();
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(a1_outer_one_1, Eigen::ComputeFullU | Eigen::ComputeFullV);

    
    U = svd.matrixU();
    S = svd.singularValues();
    Vt = svd.matrixV().transpose();

    Eigen::Matrix2d Sigma = S.asDiagonal();
    Eigen::Matrix2d lam = Eigen::Matrix2d::Identity(Sigma.rows(), Sigma.cols());
    lam(lam.rows() - 1, lam.cols() - 1) = U.determinant() * Vt.transpose().determinant();
    Eigen::Matrix2d cwe = U * lam * Vt;

    std::vector<double> rn(dim - 1, std::sqrt(std::pow(ci, 2) - std::pow(cmin, 2)) / 2);
    Eigen::MatrixXd r(dim, dim);
    r.setZero();
    r(0, 0) = ci / 2;
    for (int i = 1; i < dim; i++)
    {
        r(i, i) = rn[i - 1];
    }

    Eigen::Vector2d output;
    while(true){

        std::vector<double> xball = sample_unit_ball(dim);
        // convert xball to Eigen::Vector2d
        Eigen::Vector2d xball_eigen(xball[0], xball[1]);
        output = center + (cwe * r) * xball_eigen;

        //  self.intersection - can be checked in another way?
        Eigen::Vector2d int_output{ int(output(0)), int(output(1)) };

        // check for intersection in intersect function
        if(intersection_check(int_output)){
            break;

        }

   }

    Node xrand = Node(output(0), output(1), &start, &goal);
    std::vector<Node> xrand_list;
    xrand_list.push_back(xrand);

    return xrand_list;    
}

void Bit_star::get_f_hat_map(){

    // assume map has been stored in a 2d array
    f_hat_map = Eigen::MatrixXd::Zero(map.rows(), map.cols());

    for(int i = 0; i < map.rows(); i++){
        for(int j = 0; j < map.cols(); j++){
            double f_hat = sqrt(pow(i - goal.x, 2) + pow(j - goal.y, 2)) + sqrt(pow(i - start.x, 2) + pow(j - start.y, 2));
            f_hat_map(i, j) = f_hat;
        }
    }
}





void Bit_star::get_PHS(){

    
    Eigen::Vector2d phs;

    for(int i = 0; i < f_hat_map.rows(); i++){
        for(int j = 0; j < f_hat_map.cols(); j++){
            if(f_hat_map(i, j) < ci){
                phs << i, j;
                xphs.push_back(phs);
            }
        }
    }
    old_ci = ci;

    // TODO: save map's free nodes in a vector somewhere else
    // Assume: free nodes are stored in a vector

    std::sort(xphs.begin(), xphs.end());
    std::sort(free_nodes.begin(), free_nodes.end());

    intersection.clear();
    std::set_intersection(xphs.begin(), xphs.end(),
                      free_nodes.begin(), free_nodes.end(),
                      std::back_inserter(intersection));

}



bool Bit_star::intersection_check(Eigen::Vector2d node){

    if(std::find(intersection.begin(), intersection.end(), node) != intersection.end()){
        return true;
    }
    else{
        return false;
    }
}

void Bit_star::expand_next_vertex(){
    Node vmin = vertex_q.top();

    if vmin unexp_vertex
}



Node Bit_star::sample(){


    if(old_ci != ci){
        get_PHS();
    }

    // TODO: save map size
    // Assume: map size is stored in map_size
    if(xphs.size() <  map_size){

        return samplePHS()[0];


    }
    else{

        // TODO: create the map sample function
        // Assume: map_sample() returns a random node from the map
        return map_sample();        
    }

}

std::vector<Node> Bit_star::prune(){


    std::vector<Eigen::Vector2d> x_reuse;
    std::vector<Node> x_reuse_node;

    // TODO: write a function to save unconnected nodes in a vector
    // Assume: unconnected nodes are stored in a vector

    for(int i = 0; i < unconnected_nodes.size(); i++){
        if(unconnected_nodes[i].f_hat >= ci){
            // TODO: write a function to remove a node from the tree
        }
    }

    // Sort the vector by the gt member variable
    // TODO: write a function to save connected nodes in a vector or how is  the tree stored?
    // Assume: connected nodes are stored in a vector
    std::sort(connected_nodes.begin(), connected_nodes.end(), [](const Node& a, const Node& b) {
        return a.gt < b.gt;
    });

    for(int i = 0; i < connected_nodes.size(); i++){

        if(connected_nodes[i].x != start.x && connected_nodes[i].y != start.y && connected_nodes[i].x != goal.x && connected_nodes[i].y != goal.y){
            if((connected_nodes[i].f_hat > ci) || (connected_nodes[i].gt + connected_nodes[i].f_hat > ci)){

                // TODO: write a function to remove a node from the tree

                // use std::find and check if connected_nodes[i] is in the unexp_vertex vector
                // if it is, remove it from the unexp_vertex vector
                // if it is not, add it to the x_reuse vector
                if(std::find(unexp_vertex.begin(), unexp_vertex.end(), connected_nodes[i]) != unexp_vertex.end()){
                    unexp_vertex.erase(std::remove(unexp_vertex.begin(), unexp_vertex.end(), connected_nodes[i]), unexp_vertex.end());
                }
                
                if(connected_nodes[i].f_hat < ci){
                    x_reuse.push_back(connected_nodes[i]);
                }
        }
        }
    }

    // TODO: write a function to create  node with all aspects and return vector as node
    // Assume: create_node() returns a node with all aspects

    for(int i = 0; i < x_reuse.size(); i++){
        x_reuse_node[i] = create_node(x_reuse[i]);
    }


    return x_reuse_node;
}


std::vector<Eigen::Vector2d> Bit_star::final_solution(){

    std::vector<Node> path;
    double path_length = 0.0;

    Node* current = &goal;
    while(current->x != start.x && current->y != start.y){
        path.push_back(*current);
        path_length += sqrt(pow(current->x - current->parent->x, 2) + pow(current->y - current->parent->y, 2));
        current = current->parent;

    }
    path.push_back(start);

}


int main(){

    Node start(0, 0);
    Node goal(10, 10);

    // TODO: assume map is stored in a 2d array
    // TODO: get the map via the occupancy grid

    Eigen::MatrixXd map = Eigen::MatrixXd::Ones(10, 10);
    Bit_star tree(start, goal, map);

    int iteration = 0;

    while(true){

        iteration++;

        if(tree.edge_q.empty() && tree.vertex_q.empty()){

            tree.x_reuse = tree.prune();

            std::vector<Node> x_sampling;
            // create m samples
            for(int i = 0; i < tree.no_samples; i++){
               
                Node sample = tree.sample();

                if(std::find(x_sampling.begin(), x_sampling.end(), sample) == tree.x_reuse.end()){
                    x_sampling.push_back(sample);
                }

            }

            for(int i = 0; i < tree.x_reuse.size(); i++){
               
               tree.x_new.push_back(tree.x_reuse[i]);

            }

            // can we make this better to remove duplicates?
            for(int i = 0; i < x_sampling.size(); i++){
                // check if x_reuse[i] is in x_new
                if(std::find(tree.x_new.begin(), tree.x_new.end(), x_sampling[i]) == tree.x_new.end()){
                    tree.x_new.push_back(x_sampling[i]);
                    tree.vert.push_back(x_sampling[i]);
                }
            }

            // how to check if the nodes are connected?

            // TODO: write a function to return connected nodes

            // Assume: connected nodes are stored in a vector

            std::vector<Node> connected_nodes = tree.connected_nodes();

            for(int i = 0; i < connected_nodes.size(); i++){
                Node n = connected_nodes[i];
                // tree.qv.put((tree.gt(n) + tree.h_hat(n), n))
                // TODO: push to the priority queue based on a comparison function
                // Assume: assume that the current way is correct
                tree.vertex_q.push(std::make_pair(n.g  + n.h_hat, n));
            }
        }

            while(true){

                if(tree.vertex_q.empty()){
                    break;
                }

                // TODO: write the expand next vertex function
                tree.expand_next_vertex();

                if(tree.edge_q.empty()){
                    continue;
                }

                // TODO: write a functions that returns a partifcualr cost of a node
                if(tree.vertex_q.empty() || (tree.edge_q.top().f_hat <= tree.vertex_q.top().f_hat)){
                    break;

            }
            }

            if(!(tree.edge_q.empty())){

                    //TODO: how to get the nodes in the edge based on the priority queue?
                    Node vmin,  Node xmin = tree.edge_q.top();
                    tree.edge_q.pop();
                    // TODO: write a function to compute the get c_hat
                    if(vmin.g + c_hat(vmin, xmin) + xmin.h_hat < tree.ci){

                        if(vmin.g + c_hat(vmin, xmin) < xmin.g){
                            
                            // CHECK BUG: what should be done here? - unclear in python code
                            double cedge = c_hat(vmin, xmin);
                            if(vmin.g + c_hat(vmin, xmin) + xmin.h_hat < tree.ci){
                                if(vmin.g + c_hat(vmin, xmin) < xmin.g){
                                    // TODO: how to define connected?
                                    if(std::find(tree.connected.begin(), tree.connected.end(), xmin) != tree.connected.end()){
                                   
                                        // TODO: remove edge from edge list
                                        // TODO: add edge with weight

                                    }
                                    else{

                                        tree.vert.push_back(xmin);
                                        // TODO: create edge data structure that takes in the weight too
                                        tree.edge.push_back(std::make_pair(vmin, xmin), cedge);
                                        // TODO: add to vertex queue
                                        tree.unexp_vertex.push_back(xmin);

                                        if(xmin.x == tree.goal.x && xmin.y == tree.goal.y){
                                            tree.vsol.push_back(xmin);
                                        }

                                    }

                                    // CHECK: check if goal has a g.
                                    tree.ci = tree.goal.g;

                                    if(xmin.x == tree.goal.x && xmin.y == tree.goal.y){
                                        std::cout << "Found a solution" << std::endl;
                                        double length = tree.final_solution();
                                        std::cout<< "Path length: " << length << std::endl;

                                    }




                                }


                            }

                        }

                    }
                    else{

                        // create an empty priority queue for the vertex queue
                        std::priority_queue<Node> empty_q;
                        std::swap(tree.vertex_q, empty_q);

                        // TODO: check the datatype of edge queue
                        std::priority_queue<Node> empty_q2;
                        std::swap(tree.edge_q, empty_q2);
                    }

            }
            else
            {

                        // create an empty priority queue for the vertex queue
                        std::priority_queue<Node> empty_q;
                        std::swap(tree.vertex_q, empty_q);

                        // TODO: check the datatype of edge queue
                        std::priority_queue<Node> empty_q2;
                        std::swap(tree.edge_q, empty_q2);



            }
        }



    }


}











