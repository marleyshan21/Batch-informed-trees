#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <random>
#include <vector>
#include <queue>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <set>
#include <opencv2/opencv.hpp>

class Node
{

private:
    void f_hat_cal()
    {
        this->f_hat = this->g_hat + this->h_hat;
    }
    void g_hat_cal()
    {
        this->g_hat = sqrt(pow(this->x - this->start->x, 2) + pow(this->y - this->start->y, 2));
    }
    void h_hat_cal()
    {
        this->h_hat = sqrt(pow(this->x - this->goal->x, 2) + pow(this->y - this->goal->y, 2));
    }
   


public:
    double x, y;
    double f_hat, g_hat, h_hat; // Estimated costs

    double parent_cost;
    // actual costs
    double gt = std::numeric_limits<double>::infinity();
    //inf
    double vertex_weight;           // Actual costs
    
    
    std::vector<Node> children;
    bool is_expanded; // We might use this
    Node *start;
    Node *goal;
    std::vector<Node> parent;
    Node()
    {
        this->x = 0.0;
        this->y = 0.0;
        this->f_hat = 0.0;
        this->g_hat = 0.0;
        this->h_hat = 0.0;
        this->parent = {};
        this->is_expanded = false;
        this->children = {};
        this->parent_cost = 5.0;
    }

    Node(double x, double y)
    {
        this->x = x;
        this->y = y;
        this->f_hat = 0.0;
        this->g_hat = 0.0;
        this->h_hat = 0.0;
        this->parent = {};
        this->children = {};

    }

   


    Node(double x, double y, bool self_calculate)
    {
        this->x = x;
        this->y = y;
        this->start = start;
        // this->gt = 
        if(self_calculate)
        {
            this->g_hat_cal();
            // this->h_hat_cal();
            this->h_hat = 0.0;
            this->f_hat_cal();
            this->vertex_weight = this->gt  + this->h_hat;
        }
    }

    Node(double x, double y, double gt, double g_hat)
    {
        this->x = x;
        this->y = y;
        this->gt = gt;
        this->g_hat = g_hat;       
    }

    Node(double x, double y, Node *start,  bool self_calculate)
    {
        this->x = x;
        this->y = y;
        this->start = start;
        // this->gt = 
        if(self_calculate)
        {
            this->g_hat_cal();
            // this->h_hat_cal();
            this->h_hat = 0.0;
            this->f_hat_cal();
            this->vertex_weight = this->gt  + this->h_hat;
        }
    }
   
    Node(double x, double y, Node *start, Node *goal)
    {
        this->x = x;
        this->y = y;
        this->start = start;
        this->goal = goal;
    }


    Node(double x, double y, Node parent, double gt, double parent_Cost)
    {
        this->x = x;
        this->y = y;
        this->parent.push_back(parent);
        this->gt = gt;
        this->parent_cost = parent_Cost;
    }

    Node(int x, int y, Node *start, Node *goal)
    {
        this->x = static_cast<double>(x);
        this->y = static_cast<double>(y);
        this->start = start;
        this->goal = goal;
    }

    Node(double x, double y, Node *start, Node *goal, bool self_calculate)
    {
        this->start = start;
        this->goal = goal;
        this->x = x;
        this->y = y;
        if(self_calculate)
        {
            this->g_hat_cal();
            this->h_hat_cal();
            this->f_hat_cal();
            this->vertex_weight = this->gt  + this->h_hat;
        }
    }

    Node(double x, double y, Node *start, Node *goal, double gt,  bool self_calculate)
    {
        this->x = x;
        this->y = y;
        this->start = start;
        this->goal = goal;
        this->gt = gt;
        if(self_calculate)
        {
            this->g_hat_cal();
            this->h_hat_cal();
            this->f_hat_cal();
            this->vertex_weight = this->gt  + this->h_hat;
        }


    }
    
   
    Node(double x, double y, double f_hat, double g_hat, double h_hat, double f, double g, double h, Node parent, Node child, Node *start, Node *goal)
    {
        this->x = x;
        this->y = y;
        this->f_hat = f_hat;
        this->g_hat = g_hat;
        this->h_hat = h_hat;
        this->parent.push_back(parent);
        this->children.push_back(child);
        this->start = start;
        this->goal = goal;
    }

   
    bool operator==(const Node& other) const {
        return (this->x == other.x) && (this->y == other.y);
    }
    bool operator<(const Node& other) const {
        // check: how to compare 2 nodes - needed for sets
        return this->g_hat< other.g_hat;
    }
};


class Edge {
    public:
        Node from_node;
        Node to_node;
        double edge_weight;


        Edge(Node from_node, Node to_node, double edge_weight){
            this->from_node = from_node;
            this->to_node = to_node;
            this->edge_weight = edge_weight;
        }

        bool operator<(const Edge& other) const {
            return edge_weight > other.edge_weight;
        }

        bool operator==(const Edge& other) const {
            return ((this->from_node == other.to_node) && (this->to_node == other.from_node) || (this->from_node == other.from_node) && (this->to_node == other.to_node));
        }
};

struct NodeComparator {
    bool operator() (const Node& node1, const Node& node2) {

        // check: node wiht higher cost si pushed to the back

        return node1.gt + node1.h_hat > node2.gt + node2.h_hat;
        }
    };

struct NodeComparatorSort {
    bool operator() (const Node& node1, const Node& node2) {


        return node1.gt + node1.h_hat < node2.gt + node2.h_hat;
        }
    };  
    

class Bit_star
{
public:
    Bit_star(Node start_node, Node goal_node, Eigen::MatrixXd map)
    {

        start = start_node;
        goal = goal_node;
        this->map = map;
        this->dim  = 2;
        this->Rbit = 10;
        this->no_samples = 20;
        this->ci = std::numeric_limits<double>::infinity();
        this->old_ci = std::numeric_limits<double>::infinity();
        this->map_size = map.cols() * map.rows();

        // add node  = vert. and self.V = unconnected_vertex
        this->vert.push_back(start);
        // goal is not connected to the graph and hence unconnected
        this->unconnected_vertex.push_back(goal);
        this->unexp_vertex.push_back(start);
        this->x_new = this->unconnected_vertex;  
    
        // TODO: Read map from file
        // Assuming map is a 2D matrix of 10 x 10 for now

        // For samplePHS
        cmin = sqrt(pow(goal.x - start.x, 2) + pow(goal.y - start.y, 2));
        center = { (start.x + goal.x) / 2, (start.y + goal.y) / 2 };
        a1 = { (goal.x - start.x) / cmin, (goal.y - start.y) / cmin };

        map_width = map.rows();
        map_height = map.cols();
        f_hat_map = Eigen::MatrixXd::Zero(map_width, map_height);

        std::mt19937 gen(0);
        // store free nodes and obstacles
        free_nodes_map();
        f_hat_map_data();
        this->vertex_q.push(start);
        get_PHS();
  
        
    }
    
    // variables
    Node start;
    Node goal;
    Eigen::MatrixXd map;
    int dim;
    double Rbit;
    int no_samples;
    double ci;
    double old_ci;
    double map_size;

    std::vector<Node> vert;
    std::vector<Edge> edges;
    std::vector<Node> x_new;
    std::vector<Node> x_reuse;
    std::vector<Node> unexp_vertex;
    std::vector<Node> unconnected_vertex;
    std::vector<Node> vsol;
     // vertex queue , cost = gt + h_hat of the node.
    std::priority_queue<Node, std::vector<Node>, NodeComparator> vertex_q;
    // edge queue, cost = gt + c_hat + h_hat 
    std::priority_queue<Edge> edge_q;

    double cmin;
    std::vector<Node> free_nodes;
    std::vector<Node> occupied;
    Eigen::MatrixXd f_hat_map;
    int map_width;
    int map_height;
    Eigen::Vector2d a1;
    Eigen::Vector2d center;
    std::vector<Node> xphs;
    std::vector<Node> intersection;

    std::random_device rd;
    std::mt19937 gen;



    // functions
    double a_hat(Node node1, Node node2);
    void get_PHS();
    double gt(Node node);
    double c_hat(Node *node1, Node *node2);
    double c(Node node1, Node node2);
    std::vector<Node> near(Node node, std::vector<Node> search_list);
    std::vector<double> sample_unit_ball(int d);
    void expand_next_vertex();
    Node samplePHS();
    Node sample();
    Node sample_map();
    void prune();
    std::pair<std::vector<Node>, double> final_solution();
    void update_children_gt(Node node);
    std::vector<Node> near(std::vector<Node> search_list, Node node);
    bool nodeEqual(const Node& n1, const Node& n2);
    void f_hat_map_data();
    void free_nodes_map();
    bool intersection_check(Eigen::Vector2d node);  
    void generate_phs();
    void get_f_hat_map();
    
    
    // debug variables
    std::vector<Edge> debug_edges;

};
