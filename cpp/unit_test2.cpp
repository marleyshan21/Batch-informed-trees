#include <bitstar.h>




std::vector<double> Bit_star::sample_unit_ball(int d)
{
    
    // uniform random sample from unit ball
    std::vector<double> u(d);
    // std::random_device rd;
    // std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);



    for (int i = 0; i < d; i++) {
        u[i] = dis(this->gen);
    }

    double norm = 0.0;
    for (int i = 0; i < d; i++) {
        norm += u[i] * u[i];
    }
    norm = std::sqrt(norm);

    std::uniform_real_distribution<double> dis_r(0.0, 1.0);
    double r = dis_r(this->gen);

    std::vector<double> x(d);
    for (int i = 0; i < d; i++) {
        x[i] = r * u[i] / norm;
    }

    return x;

}


bool Bit_star::intersection_check(Eigen::Vector2d node){

    Node node1(node(0), node(1));

    if(std::find(intersection.begin(), intersection.end(), node1) != intersection.end()){
        return true;
    }
    else{
        return false;
    }
}


Node Bit_star::samplePHS(){

    // SVD
    Eigen::Matrix2d U, Vt;
    Eigen::Vector2d one_1(1.0, 0.0);
    Eigen::Matrix2d a1_outer_one_1 = a1 * one_1.transpose();


    Eigen::JacobiSVD<Eigen::Matrix2d> svd(a1_outer_one_1, Eigen::ComputeFullU | Eigen::ComputeFullV);

    U = svd.matrixU();
    Eigen::Vector2d S  = svd.singularValues();
    Vt = svd.matrixV().transpose();

    Eigen::Matrix2d Sigma = S.asDiagonal();
    Eigen::Matrix2d lam = Eigen::Matrix2d::Identity(Sigma.rows(), Sigma.cols());
    lam(lam.rows() - 1, lam.cols() - 1) = U.determinant() * Vt.transpose().determinant();
    Eigen::Matrix2d cwe = U * (lam * Vt);

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
    // Node node = Node(x, y, &start, &goal, true);
    Node xrand = Node(output(0), output(1), &start, &goal, true);
  
    return xrand;


}

Node Bit_star::sample_map(){

    // print sample map started
    // std::cout << "sample map started" << std::endl;
    while(true){

        
        std::uniform_real_distribution<double> dis(0.0, map_width);
        std::uniform_real_distribution<double> dis1(0.0, map_height);
        double x = dis(this->gen);
        double y = dis1(this->gen);

        if (map(int(x), int(y)) > 0){

            Node node = Node(x, y, &start, &goal, true);
            return node;

        }

    }

}


void Bit_star::free_nodes_map(){

    for(int i=0; i<map.rows(); i++){
        for(int j=0; j<map.cols(); j++){
            if(map(i,j) > 0){
                free_nodes.push_back(Node(i,j,&start,&goal));
            }
            else{
                occupied.push_back(Node(i,j));
            }
        }
    }

}


void Bit_star::f_hat_map_data(){

    
    for (int x = 0; x < map.rows(); x++) {
        for (int y = 0; y < map.cols(); y++) {
        double f_hat = std::sqrt(std::pow(x - goal.x, 2) + std::pow(y - goal.y, 2)) 
                       + std::sqrt(std::pow(x - start.x, 2) + std::pow(y - start.y, 2));
        f_hat_map(x, y) = f_hat;
        }
    }        

}

void Bit_star::get_PHS(){

    for(int i=0; i<f_hat_map.rows(); i++){
        for(int j=0; j<f_hat_map.cols(); j++){
            if(f_hat_map(i,j) < ci){
                xphs.push_back(Node(i,j,&start,&goal));
            }
        }
    }

    // old_ci = ci;
   
    std::set<Node> set_xphs(xphs.begin(), xphs.end());
    std::set<Node> set_free_nodes(free_nodes.begin(), free_nodes.end());
    intersection.clear();
    std::set_intersection(set_xphs.begin(), set_xphs.end(),
                        set_free_nodes.begin(), set_free_nodes.end(),
                        std::inserter(intersection, intersection.begin()));
}


Node Bit_star::sample(){


    Node xrand;
    if(old_ci!=ci){
        get_PHS();
    }

    if(xphs.size() < map_size){
        // std::cout << "sample phs" << std::endl;
        xrand = samplePHS();

    }
    else{

        xrand = sample_map();

    }

    return xrand;

}

std::vector<Node> Bit_star::near(std::vector<Node> search_list, Node node){

    std::vector<Node> near_list;
    for (int i = 0; i < search_list.size(); i++)
    {
        if ((c_hat(&search_list[i], &node) <= Rbit) && (search_list[i].x != node.x && search_list[i].y != node.y))
        {
            near_list.push_back(search_list[i]);
        }
    }

    return near_list;

}

void Bit_star::prune(){

    // std::vector<Node> x_reuse;
    std::vector<Node> new_unconnected;
    for(auto node : this->unconnected_vertex){
        if(node.f_hat < ci){
            // x_reuse.push_back(node);
            new_unconnected.push_back(node);
        }
    }
    this->unconnected_vertex = new_unconnected;


    std::vector<Node> sorted_nodes = this->vert;

    std::sort(sorted_nodes.begin(), sorted_nodes.end(), [](const Node& a, const Node& b) {
        return a.g_hat < b.g_hat;
    });


    for(auto node : sorted_nodes){
        if( !(node== this->start) && !(node== goal)){

            if(node.f_hat >= ci || (node.gt + node.h_hat >= ci)){
                // x_reuse.push_back(node);
                // this->remove_node(&node);
                // remove node from vert
                if (std::find(this->vert.begin(), this->vert.end(), node) != this->vert.end())
                {
                    this->vert.erase(std::remove(this->vert.begin(), this->vert.end(), node), this->vert.end());
                }

                // remove node from vsol
                if (std::find(this->vsol.begin(), this->vsol.end(), node) != this->vsol.end())
                {
                    this->vsol.erase(std::remove(this->vsol.begin(), this->vsol.end(), node), this->vsol.end());
                }

                // remove node from unexp_vertex
                if (std::find(this->unexp_vertex.begin(), this->unexp_vertex.end(), node) != this->unexp_vertex.end())
                {
                    this->unexp_vertex.erase(std::remove(this->unexp_vertex.begin(), this->unexp_vertex.end(), node), this->unexp_vertex.end());
                }

                // remove node from edges
                for (int i = 0; i < this->edges.size(); i++)
                {   if(node.parent.size() == 0){
                        break;
                    }
                    //BUG:some values are not having parents and this is stopping then - above line will help?
                    if (this->edges[i].to_node == node.parent[0] || this->edges[i].from_node == node)
                    {
                        this->edges.erase(this->edges.begin() + i);
                    }
                }

                // remove the children of the parent of the node

                if(node.parent.size() > 0){
                    
                    for(auto child : node.parent[0].children){
                        
                        //BUG:some values are not having parents and this is stopping then - above if condition  will help?
                        if(child == node){
                            // remove node from node.parent->children
                            // node.parent->children.erase(std::remove(node.parent->children.begin(), node.parent->children.end(), node), node.parent->children.end());
                            node.parent[0].children.erase(std::remove(node.parent[0].children.begin(), node.parent[0].children.end(), node), node.parent[0].children.end());
                        }
                    }
                }

                if(node.f_hat < ci){

                    x_reuse.push_back(node);

                }


            }
           


        }
           
    }
    this->unconnected_vertex.push_back(this->goal);


    // return x_reuse;

}


std::pair<std::vector<Node>, double> Bit_star::final_solution() {

        std::cout << "final solution" << std::endl;

        std::vector<Node> path;
        double path_length = 0; 
        Node node = this->edges.back().to_node;
        // Node node = this->goal;
        std::vector<Node> parent;
        while (!(node == this->start)) {
            std::cout << "parent cost: " << node.parent_cost << std::endl;


            // double parent_cost = node.parent[0].children[0].parent_cost;
            path_length += node.parent_cost;
            path.push_back(node);
            node = node.parent[0];

        }
        path.push_back(this->start);
        std::reverse(path.begin(), path.end());
        return std::make_pair(path, path_length);
}

void Bit_star::update_children_gt(Node node){

    if(node.children.size() > 0){
        for(auto child : node.children){
            child.gt = child.parent_cost + node.gt;
            this->update_children_gt(child);
        }
    }


}



void Bit_star::expand_next_vertex()
{

    // std::cout << "expand_next_vertex" << std::endl;
    Node vmin = vertex_q.top();
    vertex_q.pop();
    std::vector<Node> near_list;
    // std::cout << "vmin: " << vmin.x << "," << vmin.y << std::endl;

    // p[rint unexp_vertex size
    // std::cout << "unexp_vertex size: " << this->unexp_vertex.size() << std::endl;

    // check if vmin is in unexp_vertex
    if (std::find(this->unexp_vertex.begin(), this->unexp_vertex.end(), vmin) != this->unexp_vertex.end())
    {
        // std::cout << "vmin is in unexp_vertex" << std::endl;
        near_list = near(this->unconnected_vertex, vmin);
    }
    else
    {   
        
        // std::cout << "vmin is not in unexp_vertex" << std::endl;

        std::vector<Node> intersect;

        // // print all valuesin x_new
        // std::cout << "x_new size: " << x_new.size() << std::endl;
        // for(auto node : x_new){
        //     std::cout << node.x << "," << node.y << std::endl;
        // }

        // // print all values in unconnected_vertex
        // std::cout << "unconnected_vertex size: " << unconnected_vertex.size() << std::endl;
        // for(auto node : unconnected_vertex){
        //     std::cout << node.x << "," << node.y << std::endl;
        // }

        

        std::set<Node> set_x_new(x_new.begin(), x_new.end());

        // print set_x_new
        // std::cout << "set_x_new size: " << set_x_new.size() << std::endl;
        // for(auto node : set_x_new){
        //     std::cout << node.x << "," << node.y << std::endl;
        // }

        



        std::set<Node> set_unconnected_vertex(unconnected_vertex.begin(), unconnected_vertex.end());
        
        // print set_unconnected_vertex
        // std::cout << "set_unconnected_vertex size: " << set_unconnected_vertex.size() << std::endl;
        // for(auto node : set_unconnected_vertex){
        //     std::cout << node.x << "," << node.y << std::endl;
        // }
        
        
        std::set_intersection(set_x_new.begin(), set_x_new.end(),
                            set_unconnected_vertex.begin(), set_unconnected_vertex.end(),
                            std::inserter(intersect, intersect.begin()));


        near_list = near(intersect, vmin);

    }

    std::cout << "near_list size: " << near_list.size() << std::endl;

    for(auto near : near_list){
        // std::cout << "near: " << near.x << "," << near.y << std::endl;
        // std::cout << "a_hat(vmin, near): " << a_hat(vmin, near) << std::endl;
        if(a_hat(vmin, near) < ci) {

            std::cout << "test1" << std::endl;
            std::cout << " vim gt: " << vmin.gt << std::endl;
            std::cout << "c(vmin, near): " << c(vmin, near) << std::endl;
            std::cout << "near h_hat: " << near.h_hat << std::endl;
            double cost = vmin.gt + c(vmin, near) + near.h_hat;
            std::cout << "test2" << std::endl;
            Edge edge = Edge(vmin, near, cost);
            // std::cout << "edge cost: " << edge.edge_weight << std::endl;
            // std::cout << "edge: " << edge.from_node.x << "," << edge.from_node.y << " -> " << edge.to_node.x << "," << edge.to_node.y << std::endl;
            this->edge_q.push(edge);
        }
    }

    // std::cout << "edge_q size in expAND: " << this->edge_q.size() << std::endl;

    // check if vmin is in unexp_vertex
    // WHAT? : check what happens if vmin is in unexp_vertex and we havent sampled any new nodes
    if (std::find(this->unexp_vertex.begin(), this->unexp_vertex.end(), vmin) != this->unexp_vertex.end())
    {
        std::vector<Node> v_near = near(this->vert, vmin);
        for(auto near : v_near){

            // check if the edge exist in edges
            if (std::find(this->edges.begin(), this->edges.end(), Edge(vmin, near, c_hat(&vmin, &near))) == this->edges.end())
            {
                if((a_hat(vmin, near) < ci ) && (near.g_hat + c_hat(&vmin, &near) < near.gt)){


                    double cost = vmin.gt + c(vmin, near) + near.h_hat;
                    Edge edge = Edge(vmin, near, cost);
                    this->edge_q.push(edge);
                }
                
            }
            
        }
        // remove vmin from unexp_vertex
        auto new_end = std::remove(this->unexp_vertex.begin(), this->unexp_vertex.end(), vmin);
        this->unexp_vertex.erase(new_end, this->unexp_vertex.end());
    }
}

double Bit_star::c_hat(Node *node1, Node *node2)
{
    return sqrt(pow(node1->x - node2->x, 2) + pow(node1->y - node2->y, 2));
}

double Bit_star::a_hat(Node node1, Node node2)
{
    return node1.g_hat + c_hat(&node1, &node2) + node2.h_hat;
}

double Bit_star::c(Node node1, Node node2)
{
    
    double x1 = node1.x;
    double y1 = node1.y;
    double x2 = node2.x;
    double y2 = node2.y;
    std::cout << "x1, y1: " << x1 << "," << y1 << " x2, y2: " << x2 <<","<< y2 << std::endl;
    // std::cout << sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)) << std::endl;
    int n_divs = int(10 * sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)));
    std::cout << "n_divs: " << n_divs << std::endl;
    if(n_divs == 0){
        return c_hat(&node1, &node2);
    }
    for (double lam = 0; lam <= 1; lam += 1.0/n_divs){

        // std::cout << "lam: " << lam << std::endl;
        int x = int(x1 + lam * (x2 - x1));
        int y = int(y1 + lam * (y2 - y1));
    
        std::cout << "x, y: " << x << "," << y << std::endl;

        if((x==1 && y==1) || (x==1 && y==2) || (x == 3 && y == 1) || (x == 3 && y == 2)){
            std::cout << " OBSTACLE!!!!" << std::endl;
            std::cout << "x1, y1: " << x1 << "," << y1 << " x2, y2: " << x2 <<","<< y2 << std::endl;

        }


        for(auto node : this->occupied){
            if(node.x == x && node.y == y){
                std::cout << "x1, y1: " << x1 << "," << y1 << " x2, y2: " << x2 <<","<< y2 << std::endl;

                return std::numeric_limits<double>::infinity();
            }
        }
    }
    return c_hat(&node1, &node2); 
}
double Bit_star::gt(Node node)
{

    if (node.x == start.x && node.y == start.y)
    {
        return 0.0;
    }
    //check:  assumed all the nodes in vert are conected?
    else if (std::find(this->vert.begin(), this->vert.end(), node) == this->vert.end())
    {
        return std::numeric_limits<double>::infinity();
    }



    return node.parent_cost + node.parent[0].gt;

}

bool Bit_star::nodeEqual(const Node& n1, const Node& n2) {
    if(n1.x == n2.x && n1.y == n2.y){
        return true;
    }
    return false;

}

// TODO: check main function properly
int main()
{
    
    Node* start = new Node(0.0, 0.0, 0.0, 0.0);

    Node* goal = new Node(99.0, 99.0, start, true);
    start->h_hat = sqrt(pow(start->x - goal->x, 2) + pow(start->y -goal->y, 2));

    // print goal's g_hat, h_hat and f_hat
    // std::cout << "goal's g_hat: " << goal->g_hat << std::endl;
    // std::cout << "goal's h_hat: " << goal->h_hat << std::endl;
    // std::cout << "goal's f_hat: " << goal->f_hat << std::endl;

    Eigen::MatrixXd map = Eigen::MatrixXd::Ones(100, 100);
    // map(3,3) = 0;
    // // map(4,4) = 0;
    // // map(4,2) = 0;
    // map(4,3) = 0;

    // map(2,1) = 0;
    // map(2,2) = 0;
    // map(1,2) = 0;

    for(int i = 20; i < 60; i++){
        for(int j = 40; j < 80; j++){
            map(i,j) = 0;
        }
        // std::cout << std::endl;
    }

    // std::Cout << map << std::endl;

    Bit_star *tree = new Bit_star(*start, *goal, map);
    
    int iteration = 0;
    std::cout << "start" << std::endl;

    // check if start is in free_nodes
    if (std::find(tree->free_nodes.begin(), tree->free_nodes.end(), *start) == tree->free_nodes.end())
    {
        std::cout << "start is not in free_nodes" << std::endl;
        return 0;
    }
    // check if goal is in free_nodes
    if (std::find(tree->free_nodes.begin(), tree->free_nodes.end(), *goal) == tree->free_nodes.end())
    {
        std::cout << "goal is not in free_nodes" << std::endl;
        return 0;
    }

    // if staart = goal
    if (tree->nodeEqual(*start, *goal))
    {
        tree->vsol.push_back(*start);
        tree->ci = 0.0;
        std::cout << "start = goal" << std::endl;
        return 0;
    }

    while(tree->ci <= tree->old_ci){
        
        std::cout << "iteration: " << iteration << std::endl;
        iteration++;

        if(tree->edge_q.empty() && tree->vertex_q.empty()){
            
            std::cout << "edge_q and vertex_q are empty" << std::endl;
            tree->prune();
            // std::cout << "x_reuse size: " << tree->x_reuse.size() << std::endl;
            std::vector<Node> x_sampling;

            while(x_sampling.size() < tree->no_samples){
                Node node = tree->sample();
                // check if node is in x_Sampling
                // if (std::find(x_sampling.begin(), x_sampling.end(), node) == x_sampling.end())
                // {
                    x_sampling.push_back(node);
                // }

            }
            // std::cout << "x_sampling size: " << x_sampling.size() << std::endl;
            // print all the nodes in x_sampling
            // for(auto node : x_sampling){
            //     std::cout << "x_sampling: " << node.x << ", " << node.y << std::endl;
            // }

            for(auto node : x_sampling){
                tree->x_new.push_back(node);
            }
            // std::cout << "x_new size: " << tree->x_new.size() << std::endl;

            for(auto node : tree->x_reuse){
                tree->x_new.push_back(node);
            }
            // std::cout << "x_new size: " << tree->x_new.size() << std::endl;
            // remove duplicates
            // sort x_new
            // CHECK:  is this sorting correct?
            std::sort(tree->x_new.begin(), tree->x_new.end(), [](const Node& a, const Node& b) {
                return a.g_hat < b.g_hat;
            });
            auto newEnd = std::unique(tree->x_new.begin(), tree->x_new.end(),
                    [tree](const Node& n1, const Node& n2) {
                        return tree->nodeEqual(n1, n2);
                    });
            tree->x_new.erase(newEnd, tree->x_new.end());
                        

            
            for(auto node : tree->x_new){
                tree->unconnected_vertex.push_back(node);
            }
            std::sort(tree->unconnected_vertex.begin(), tree->unconnected_vertex.end(), [](const Node& a, const Node& b) {
                return a.g_hat < b.g_hat;
            });
            // std::cout<<"unconnected_vertex size: " << tree->unconnected_vertex.size() << std::endl;
            // remove duplicates
            auto newEnd1 = std::unique(tree->unconnected_vertex.begin(), tree->unconnected_vertex.end(),
                    [tree](const Node& n1, const Node& n2) {
                        return tree->nodeEqual(n1, n2);
                    });
            tree->unconnected_vertex.erase(newEnd1, tree->unconnected_vertex.end());

            for(auto node: tree->vert){
                tree->vertex_q.push(node);
            }

            // std::cout<<"unconnected_vertex size: " << tree->unconnected_vertex.size() << std::endl;
            // std::cout<<"vertex_q size: " << tree->vertex_q.size() << std::endl;
            // std::cout<<"edge_q size: " << tree->edge_q.size() << std::endl;


       }

        while(true){

            if(tree->vertex_q.empty()){
                std::cout << "vertex_q is empty" << std::endl;
                break;
            }
            tree->expand_next_vertex();
            std::cout << "vertex_q size: " << tree->vertex_q.size() << std::endl;
            if(tree->edge_q.empty()){
                continue;
            }

            if(tree->vertex_q.empty() || (tree->vertex_q.top().vertex_weight) <= tree->edge_q.top().edge_weight){
                break;
           }

        }

        if(!(tree->edge_q.empty())){
            
            Edge top_edge = tree->edge_q.top();
            

            // print top_edge
            std::cout << "top_edge: " << top_edge.from_node.x << ", " << top_edge.from_node.y << " to " << top_edge.to_node.x << ", " << top_edge.to_node.y << std::endl;
            std::cout << "edge weight: " << top_edge.edge_weight << std::endl;
            std::cout << "top_edge.from_node.gt: " << top_edge.from_node.gt << std::endl;
            // std::cout << "c_hat: " << tree->c_hat(&top_edge.from_node, &top_edge.to_node) << std::endl;
            std::cout << "top_edge.to_node.h_hat: " << top_edge.to_node.h_hat << std::endl;
            std::cout << "tree->ci: " << tree->ci << std::endl;

            if(top_edge.from_node.gt + tree->c_hat(&top_edge.from_node, &top_edge.to_node) + top_edge.to_node.h_hat < tree->ci){
                
                // std::cout << "first if" << std::endl;
                if(top_edge.from_node.gt + tree->c_hat(&top_edge.from_node, &top_edge.to_node) < top_edge.to_node.gt){
                    // std::cout << "second if" << std::endl;
                    double cedge = tree->c(top_edge.from_node, top_edge.to_node);
                    if(top_edge.from_node.gt + cedge + top_edge.to_node.h_hat < tree->ci){
                        // std::cout << "third if" << std::endl;
                        if(top_edge.from_node.gt + cedge < top_edge.to_node.gt){
                            // std::cout << "fourth if" << std::endl;
                        //    check if to_node exists in connected_vertex
                            if (std::find(tree->vert.begin(), tree->vert.end(), top_edge.to_node) != tree->vert.end())
                            {
                                // remove to_node.parent and to_node from edges
                                for(auto edges : tree->edges){
                                    if(!top_edge.to_node.parent.empty()) {
                                        if (edges.from_node.x == top_edge.to_node.parent[0].x &&
                                            edges.from_node.y == top_edge.to_node.parent[0].y &&
                                            edges.to_node.x == top_edge.to_node.x &&
                                            edges.to_node.y == top_edge.to_node.y) {
                                            tree->edges.erase(
                                                    std::remove(tree->edges.begin(), tree->edges.end(), edges),
                                                    tree->edges.end());
                                        }
                                    }
                                }

                                // remove the node from the parent's children
                                 // remove the children of the parent of the node
                                if(!top_edge.to_node.parent.empty()) {
                                    for (auto child: top_edge.to_node.parent[0].children) {
                                        if (child == top_edge.to_node) {
                                            // remove node from node.parent->children
                                            top_edge.to_node.parent[0].children.erase(
                                                    std::remove(top_edge.to_node.parent[0].children.begin(),
                                                                top_edge.to_node.parent[0].children.end(),
                                                                top_edge.to_node),
                                                    top_edge.to_node.parent[0].children.end());
                                        }
                                    }
                                }

                                top_edge.to_node.parent.clear();
                                top_edge.to_node.parent.push_back(top_edge.from_node);
                                //  std::cout << "top_edge.to_node.parent: " << top_edge.to_node.parent[0].x << ", " << top_edge.to_node.parent[0].y << std::endl;
                                top_edge.to_node.parent_cost = cedge;
                                top_edge.to_node.gt = tree->gt(top_edge.to_node);
                                std::cout << "top_edge.to_node.gt: " << top_edge.to_node.gt << std::endl;
                                // add edge between to_node's parent abd to_node
                                double cost = top_edge.to_node.parent[0].gt + tree->c(top_edge.to_node.parent[0], top_edge.to_node) + top_edge.to_node.h_hat;
                                Edge e = Edge(top_edge.from_node, top_edge.to_node, cost);
                                tree->edges.push_back(e);

                                // add to_node to vertex_q
                                tree->vertex_q.push(top_edge.to_node);

                                // add to_node to unexp_vertex
                                tree->unexp_vertex.push_back(top_edge.to_node);

                                // add to_node to to_node.parent->children
                                top_edge.to_node.parent[0].children.push_back(top_edge.to_node);

                                // update to_node's children's gt
                                tree->update_children_gt(top_edge.to_node);

                            }

                            else{


                                //  remove to_node from unconnected_vertex
                                for(auto vert : tree->unconnected_vertex){
                                    if(vert.x == top_edge.to_node.x && vert.y == top_edge.to_node.y){
                                        tree->unconnected_vertex.erase(
                                                std::remove(tree->unconnected_vertex.begin(), tree->unconnected_vertex.end(), vert),
                                                tree->unconnected_vertex.end());
                                    }
                                }

                                // add to_node to vert
                                tree->vert.push_back(top_edge.to_node);
                                
                                // to_nodes parent is from_node
                                top_edge.to_node.parent.clear();
                                top_edge.to_node.parent.push_back(top_edge.from_node);
                                // print top edge to_node's parent
                                // std::cout << "top_edge.to_node.parent: " << top_edge.to_node.parent[0].x << ", " << top_edge.to_node.parent[0].y << std::endl;
                                top_edge.to_node.parent_cost = cedge;
                                top_edge.to_node.gt = tree->gt(top_edge.to_node);
                                std::cout << "top_edge.to_node.gt: " << top_edge.to_node.gt << std::endl;
                                // add edge between to_node's parent abd to_node
                                double cost = top_edge.to_node.parent[0].gt + tree->c(top_edge.to_node.parent[0], top_edge.to_node) + top_edge.to_node.h_hat;
                                Edge e = Edge(top_edge.to_node.parent[0], top_edge.to_node, cost);
                                tree->edges.push_back(e);

                                

                                // add to_node to vertex_q
                                tree->vertex_q.push(top_edge.to_node);

                                // add to_node to unexp_vertex
                                tree->unexp_vertex.push_back(top_edge.to_node);

                                //check if to_node is the goal
                                if(top_edge.to_node.x == goal->x && top_edge.to_node.y == goal->y){
                                    tree->vsol.push_back(top_edge.to_node);
                                }

                                // add to_node to to_node.parent->children
                                top_edge.to_node.parent[0].children.push_back(top_edge.to_node);


                            }


                            
                            
                            if(top_edge.to_node.x == goal->x && top_edge.to_node.y == goal->y){

                                    std::cout << "GOAL Found!" << std::endl;
                                    std::cout << "!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                                    // get results from final_solution()
                                    // std::pair<std::vector<Node>, double> Bit_star::final_solution()
                                    std::pair<std::vector<Node>, double> results = tree->final_solution();
                                    std::vector<Node> path = results.first;
                                    double path_length = results.second;

                                    // path length
                                    std::cout << "path length: " << path_length << std::endl;

                                    // print path
                                    std::cout << "path: " << std::endl;
                                    for(auto node : path){
                                        std::cout << "( " << node.x << ", " << node.y << " ), ";
                                    }
                                    std::cout << std::endl;

                                    tree->ci = top_edge.to_node.gt;
                                    std::cout << "tree->ci: " << tree->ci << std::endl;

                                    std::cout << "GOAL Found!" << std::endl;

                            }   
                        }
                    }
                }
            }
            else{


                // empty the edge queue and vertex queue
                while(!tree->edge_q.empty()){
                    tree->edge_q.pop();
                }
                while(!tree->vertex_q.empty()){
                    tree->vertex_q.pop();
                }

            }

            if(!tree->edge_q.empty()){
                tree->edge_q.pop();
            }
            std::cout << "tree->edge_q.size(): " << tree->edge_q.size() << std::endl;
            // tree->edge_q.pop();

        }
        else{
            
            // empty the edge queue and vertex queue
            while(!tree->edge_q.empty()){
                tree->edge_q.pop();
            }
            while(!tree->vertex_q.empty()){
                tree->vertex_q.pop();
            }
        }

    }

    return 0;

}