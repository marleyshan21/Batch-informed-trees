#! /usr/bn/env python3
#! -*- coding: utf-8 -*-

import numpy as np
import json
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import os
import cv2
from typing import List, Tuple
import tqdm


class Visualizer:
    """Visualizer class for visualizing the BIT* algorithm."""

    def __init__(
        self, start: np.array, goal: np.array, occ_map: np.array, output_dir: str
    ) -> None:
        """Initialize the visualizer.

        Args:
            start (np.array): Start point.
            goal (np.array): Goal point.
            occ_map (np.array): Occupancy map of the environment.
            output_dir (str): Output directory to save the plots.
        """
        # Set of all edges in the tree.
        self.edges = set()
        # Final path from start to goal.
        self.final_path = None
        # Initial cost to go from start to goal.
        self.ci = np.inf

        # Simulation number.
        self.sim = 0
        # Map of the environment converted to RGB.
        self.occ_map = cv2.cvtColor(occ_map, cv2.COLOR_BGR2RGB)
        # Start and goal points as numpy arrays.
        self.start = start
        self.goal = goal

        # All the new edges, removed edges, final paths, costs, and final edges over all the simulations.
        (
            self.all_new_edges,
            self.all_rem_edges,
            self.all_final_paths,
            self.all_cis,
            self.all_final_edge_list,
        ) = (
            [],
            [],
            [],
            [],
            [],
        )
        # A dictionary of all the lines in the plot.
        self.lines = {}

        # Initialize the plot.
        self.fig, self.ax = plt.subplots(figsize=(20, 20))
        # The output directory to save the plots.
        self.output_dir = output_dir

    def read_json(self, folder: str, max_iter: int = np.inf) -> None:
        """Reads the json files from the log directory.

        Args:
            folder (str): Name of the log directory.
            max_iter (int, optional): If nothing is given it reads all the simulations logs if a number is given 0 <= max_iter <= len(simulations) then that amount of simulations are read. Defaults to np.inf meaning all files are read.
        """
        # Sort the files in the folder.
        files = sorted(os.listdir(folder))
        # If max_iter is not given then read all the files.
        max_iter = min(max_iter, len(files))
        for i in tqdm.tqdm(range(max_iter), desc="Reading JSON files", total=max_iter):
            # Open one file at a time and read the data.
            with open(os.path.join(folder, files[i]), "r") as f:
                # Load the data from the json file.
                data = json.load(f)
                # Get all new edges and append it to the list.
                self.all_new_edges.append(data["new_edges"])
                # Get all remove edges and append it to the list.
                self.all_rem_edges.append(data["rem_edges"])
                # Get all final paths and append it to the list.
                self.all_final_paths.append(data["final_path"])
                # Get all the costs and append it to the list.
                self.all_cis.append(np.array(data["ci"]))
                # Get all the final edges and append it to the list.
                self.all_final_edge_list.append(data["final_edge_list"])

    def draw_final_path(self, path: List[Tuple[float, float]]) -> None:
        """Draw the final path from start to goal.

        Args:
            path (list): Path from start to goal.
        """
        # If the path is empty then return.
        if len(path) == 0:
            return
        # Convert to a numpy array for easy plotting.
        path = np.array(path)
        # Split the path into x and y coordinates.
        x, y = path[:, 0], path[:, 1]
        # Plot the path.
        self.ax.plot(
            y,
            x,
            color="darkorchid",
            lw=4,
            label="Final Path",
        )

    def draw_ellipse(self, ci: float, colour: str = "dimgrey") -> None:
        """Draws the Prolate Hyperspheroid (PHS) for the given ci.

        Args:
            ci (float): Cost to go from start to goal.
            colour (str, optional): The color of the ellipse while plotting. Defaults to "dimgrey".
        """
        # Return if ci is infinity.
        if ci == np.inf:
            return
        # Calculate the minimum ci which is the L2 norm from the start to the goal.
        cmin = np.linalg.norm(self.goal - self.start)
        # Get the center of the ellipse.
        center = (self.start + self.goal) / 2.0
        # Get the first radius of the Prolate Hyperspheroid (PHS).
        r1 = ci
        # Get the second radius of the Prolate Hyperspheroid (PHS).
        r2 = np.sqrt(ci**2 - cmin**2)
        # Get the angle of the Prolate Hyperspheroid (PHS) with respect to the x-axis.
        theta = np.arctan2(self.goal[0] - self.start[0], self.goal[1] - self.start[1])
        # Convert the angle from radians to degrees.
        theta = np.degrees(theta)
        # Use matplotlib patches to draw the ellipse.
        patch = Ellipse(
            (center[1], center[0]),
            r1,
            r2,
            theta,
            color=colour,
            fill=False,
            lw=5,
            ls="--",
            label="Prolate Hyperspheroid (PHS)",
        )
        # Add the ellipse to the plot.
        self.ax.add_patch(patch)

    def draw_edge(self, edge: List[List[List[float]]]) -> None:
        """Draws the edge between two points.

        Args:
            edge (List[List[List[float]]]): List of tuples of the form [[[x1, y1], [x2, y2]], [[x2, y2], [x3, y3]], ... ].
        """
        # Convert the edge to a tuple for easy indexing.
        edge_tup = tuple(map(tuple, edge))

        # Plot the edge between the two points.
        l = self.ax.plot(
            [edge[0][1], edge[1][1]],
            [edge[0][0], edge[1][0]],
            color="sandybrown",
            lw=2,
            marker="x",
            markersize=4,
            markerfacecolor="darkcyan",
            markeredgecolor="darkcyan",
            label="Branches",
        )
        # Add the line to the dictionary.
        self.lines[edge_tup] = l

    def draw_tree(self, sim: int) -> None:
        """This is a slow drawing function. It draws the tree from the start to the goal as the algorithm progresses through a simulation.

        Args:
            sim (int): Simulation number.
        """
        # Get the start and goal.
        start = self.start
        goal = self.goal
        # Calls redraw_map to redraw the map.
        self.redraw_map(sim)

        # Get the figure and axes.
        fig = self.fig
        ax = self.ax

        # Loop through all the edges and draw them.
        for i in range(len(self.all_new_edges[sim])):
            # Get the new edge, removed edge, final path and cost to go with each step in the simulation.
            new_edge = self.all_new_edges[sim][i]
            rem_edge = self.all_rem_edges[sim][i]
            path = self.all_final_paths[sim][i]
            ci = self.all_cis[sim][i]

            # Remove the edges that are no longer in the tree.
            if rem_edge:
                for rem_e in rem_edge:
                    rem_e_tup = tuple(map(tuple, rem_e))
                    try:
                        ax.lines.remove(self.lines[rem_e_tup][0])
                        self.edges.remove(rem_e_tup)
                    except:
                        continue

            # Add the new edges to the tree.
            if new_edge is not None:
                new_e_tup = tuple(map(tuple, new_edge))
                self.edges.add(new_e_tup)
                self.draw_edge(new_edge)

            # Plot the final path if it exists.
            if path is None:
                if self.final_path is not None:
                    self.draw_final_path(self.final_path)
            else:
                self.final_path = path
                self.draw_final_path(path)

            # Plot the PHS for the current cost to go.
            self.draw_ellipse(ci, colour="dimgrey")
            # Plot the start and goal.
            ax.plot(
                start[1],
                start[0],
                color="red",
                marker="*",
                markersize=20,
            )
            ax.plot(
                goal[1],
                goal[0],
                color="blue",
                marker="*",
                markersize=20,
            )
            # Remove the legend to avoid duplicates and add the legend again.
            handles, labels = self.ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            plt.legend(
                by_label.values(),
                by_label.keys(),
                bbox_to_anchor=(1.05, 1.0),
                loc="upper left",
            )

            # Show the plot and pause for a short time without blocking.
            plt.show(block=False)
            plt.pause(0.0001)

    def draw_fast(self, sim: int) -> None:
        """This is a fast drawing function. It draws the final tree, final path and PHS for a given simulation.

        Args:
            sim (int): Simulation number.
        """
        # Get the final path, edges, and cost to go for the given simulation.
        self.edges = self.all_final_edge_list[sim]
        path = self.all_final_paths[sim][-1]
        ci = self.all_cis[sim][0]

        # Redraw the map.
        self.redraw_map(sim)

        # If the final path exists, draw it.
        if path is not None:
            self.final_path = path
            self.draw_final_path(path)

        # Draw the PHS.
        self.draw_ellipse(ci, colour="dimgrey")
        # Draw the start and goal.
        self.ax.plot(
            self.start[1],
            self.start[0],
            color="red",
            marker="*",
            markersize=20,
        )
        self.ax.plot(
            self.goal[1],
            self.goal[0],
            color="blue",
            marker="*",
            markersize=20,
        )
        # Remove the legend to avoid duplicates and add the legend again.
        handles, labels = self.ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(
            by_label.values(),
            by_label.keys(),
            bbox_to_anchor=(1.05, 1.0),
            loc="upper left",
        )
        # Save the plot at the given output directory and show the plot.
        plt.savefig(f"{self.output_dir}/Bitstar_Simulation_{sim:02d}.png")
        # Show the plot and pause for a short time without blocking.
        plt.show(block=False)
        plt.pause(1)

    def draw(self, sim: int, fast: bool = False) -> None:
        """This function calls either the fast or slow drawing function depending on the value of fast flag.

        Args:
            sim (int): Simulation number.
            fast (bool, optional): If True, the fast drawing function is called. Defaults to False.
        """
        # If the fast flag is True, call the fast drawing function.
        if fast:
            self.draw_fast(sim)
        # Else, call the slow drawing function.
        else:
            self.draw_tree(sim)

    def redraw_map(self, sim: int) -> None:
        # Clear the current plot.
        plt.close()
        self.fig, self.ax = plt.subplots(figsize=(20, 20))
        # Get the occupancy map for this simulation.
        im = self.ax.imshow(self.occ_map, cmap=plt.cm.gray, extent=[0, 100, 100, 0])

        # Loop through all the edges and draw them.
        for e in self.edges:
            self.draw_edge(e)

        # Plot the start and goal.
        self.ax.plot(
            self.start[1],
            self.start[0],
            color="red",
            marker="*",
            markersize=20,
            label="Start",
        )
        self.ax.plot(
            self.goal[1],
            self.goal[0],
            color="blue",
            marker="*",
            markersize=20,
            label="Goal",
        )
        # Set the title and axis labels.
        self.ax.set_title(f"BIT* - Simulation {sim}", fontsize=30)
        self.ax.set_xlabel(r"X$\rightarrow$", fontsize=10)
        self.ax.set_ylabel(r"Y$\rightarrow$", fontsize=10)

        # Remove the legend to avoid duplicates and add the legend again.
        handles, labels = self.ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(
            by_label.values(),
            by_label.keys(),
            bbox_to_anchor=(1.05, 1.0),
            loc="upper left",
        )
        # Set the axis limits.
        self.ax.set_xlim(-10, 110)
        self.ax.set_ylim(-10, 110)
