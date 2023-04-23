#! /usr/bin/env python3
#! -*- coding: utf-8 -*-

import numpy as np
from typing import Generic, TypeVar

# Set global variables to avoid circular definitions.
global start_arr
global goal_arr

# Define the generic type T.
T = TypeVar("T")


class Node(Generic[T]):
    """Node class for BIT*. This class is used to represent a node in the tree. It contains the coordinates of the node, the parent node, the edge cost, gt, children set, start, goal, g_hat, h_hat, and f_hat.

    Args:
        Generic (Node): Generic type for the parent node.
    """

    def __init__(
        self,
        coords: tuple,
        parent: T = None,
        gt: float = np.inf,
        par_cost: float = None,
    ) -> None:
        """Initialize the Node class with the coordinates, parent, gt, par_cost, children, start, goal, g_hat, h_hat, and f_hat.

        Args:
            coords (tuple): Coordinates of the node in the form (x, y).
            parent (T, optional): Parent node of the current node. This is an instance of the Node class as the parent node is also a node. Defaults to None.
            gt (float, optional): The cost through the tree to get from the start to the current node. Defaults to np.inf.
            par_cost (float, optional): The cost of the edge from the parent node to the current node. Defaults to None.
        """
        # Extract coordinates from tuple.
        self.x = coords[0]
        self.y = coords[1]
        # Create tuple and numpy array for easy access.
        self.tup = (self.x, self.y)
        self.np_arr = np.array([self.x, self.y])

        # Initialize parent, edge cost (par_cost), and g_t.
        self.parent = parent
        self.par_cost = par_cost
        self.gt = gt
        # Initialize the children set.
        self.children = set()

        # Initialize start and goal.
        global start_arr
        self.start = start_arr
        global goal_arr
        self.goal = goal_arr

        # Generate g_hat and h_hat.
        self.g_hat = self.gen_g_hat()
        self.h_hat = self.gen_h_hat()
        # f_hat is the sum of g_hat and h_hat.
        self.f_hat = self.g_hat + self.h_hat

    def gen_g_hat(self) -> float:
        """Generate the g_hat value for the current node. This is the L2 norm between the current node and the start node.

        Returns:
            g_hat (float): The g_hat value for the current node.
        """
        # Return the L2 norm between the current node and the start node.
        return np.linalg.norm(self.np_arr - self.start)

    def gen_h_hat(self) -> float:
        """Generate the h_hat value for the current node. This is the L2 norm between the current node and the goal node.

        Returns:
            h_hat (float): The h_hat value for the current node.
        """
        # Return the L2 norm between the current node and the goal node.
        return np.linalg.norm(self.np_arr - self.goal)

    def __str__(self) -> str:
        """String representation of the Node class.

        Returns:
            str: String representation of the Node class.
        """
        # Return the string representation of the tuple.
        return str(self.tup)

    def __repr__(self) -> str:
        """String representation of the Node class.

        Returns:
            str: String representation of the Node class.
        """
        # Return the string representation of the tuple.
        return str(self.tup)
