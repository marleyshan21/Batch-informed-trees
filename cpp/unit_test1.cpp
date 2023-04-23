#include <bitstar.h>


int main()
{
    
    Node* start = new Node(0.0, 0.0);
    Node* goal = new Node(9.0, 9.0);
    Eigen::MatrixXd map = Eigen::MatrixXd::Ones(10, 10);

    Bit_star *tree = new Bit_star(*start, *goal, map);
    

    return 0;
}