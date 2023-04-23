#! /usr/bin/env python3
#! -*- coding: utf-8 -*-

import random
import numpy as np
from Node import Node


class Map:
    """Map class for BIT*. This class is used to represent the map. It contains the start and goal coordinates, the obstacles, the map, the free and occupied sets, and the f_hat map."""

    def __init__(self, start: Node, goal: Node, occ_grid: np.array) -> None:
        """Initialize the Map class with the start, goal and the occupancy grid.

        Args:
            start (np.array): Start coordinates of the robot in the form np.array([x, y]).
            goal (np.array): Goal coordinates of the robot in the form np.array([x, y]).
            occ_grid (np.array, optional): Occupancy grid of the map which is a 2D numpy array of 0s (occupied) and 1s (free).
        """
        # Set start and goal.
        self.start = start
        self.goal = goal

        self.start_arr = self.start.np_arr
        self.goal_arr = self.goal.np_arr

        # Obstacles set
        self.obstacles = set()
        # Dimensions of the search space.
        self.dim = 2

        # 2D occupancy grid of the map with 0s (occupied) and 1s (free). This is converted from an image.
        self.map = occ_grid
        # Get all the indices of the free cells.
        ind = np.argwhere(self.map > 0)
        # Create a set of tuples of the free cells from the indices for faster lookup.
        self.free = set(list(map(lambda x: tuple(x), ind)))
        # Get all the indices of the occupied cells.
        ind = np.argwhere(self.map == 0)
        # Create a set of tuples of the occupied cells from the indices for faster lookup.
        self.occupied = set(list(map(lambda x: tuple(x), ind)))
        # Get the f_hat map.
        self.get_f_hat_map()

    def sample(self) -> tuple:
        """Sample a random point from the free set. This is used to generate a new node in the tree. We don't use this in the current implementation as it can theoretically be slower than the new_sample function if the obstacle set is large.

        Returns:
            tuple: Random point sampled from the free set.
        """
        # Sample until a point is found in the free set.
        while True:
            # Sample a random point uniformly from the map and in the continuous space.
            x, y = np.random.uniform(0, self.map.shape[0]), np.random.uniform(
                0, self.map.shape[1]
            )
            # Convert the point to an integer tuple and check if it is in the free set.
            if (int(x), int(y)) in self.free:
                # Return the point.
                return (x, y)

    def new_sample(self) -> tuple:
        """Sample a random point from the free set. This is used to generate a new node in the tree. This is a modified version of the sample function which first samples a random point from the free set and then adds a random noise to it to generate a new point and return it if it is in the free set. This could theoretically be faster than the sample function especially if the obstacle set is large.

        Returns:
            tuple: Random point sampled from the free set.
        """
        # Sample until a point is found in the free set.
        while True:
            # Sample from the free set.
            free_node = random.sample(self.free, 1)[0]
            # Add a random noise to the sampled point.
            noise = np.random.uniform(0, 1, self.dim)
            # Add the noise to the sampled point.
            new_node = free_node + noise
            # If the integer tuple of the new node is in the free set, return it.
            if (int(new_node[0]), int(new_node[1])) in self.free:
                return new_node

    def get_f_hat_map(self) -> None:
        """Get the f_hat map which is the heuristic map used by the BIT* algorithm. This is the sum of the L2 norm of the distance from the goal and the start to the current point. This is precomputed for all nodes in the map and stored in a 2D numpy array as a lookup table"""
        # Get the dimensions of the map.
        map_x, map_y = self.map.shape
        # Initialize the f_hat map with zeros of the same dimensions as the map.
        self.f_hat_map = np.zeros((map_x, map_y))
        # For each point in the map, calculate the f_hat value and store it in the f_hat map.
        for x in range(map_x):
            for y in range(map_y):
                # f_hat(x) = g_hat(x) + h_hat(x) for each point x in the map.
                f_hat = np.linalg.norm(
                    np.array([x, y]) - self.goal_arr
                ) + np.linalg.norm(np.array([x, y]) - self.start_arr)
                # Store the f_hat value in the f_hat map.
                self.f_hat_map[x, y] = f_hat
