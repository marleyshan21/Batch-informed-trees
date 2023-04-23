import argparse
from BIT_Star import *
import BIT_Star as BIT_Star
import Node as node
from Node import *
from Map import *
from Visualize import *
import sys
import shutil
from datetime import datetime
from PIL import Image
from PrintColours import *


def main(
    map_name: str,
    vis: bool,
    start: list,
    goal: list,
    rbit: float,
    samples: int,
    dim: int,
    seed: int,
    stop_time: int,
    fast: bool,
) -> None:
    """The main function for the BIT* algorithm this function will run the algorithm after parsing the command line arguments. Call the bitstar class/function to run the algorithm and if the vis or fast flag is set then call the visualizer to show the path. If the fast flag is set then the fast visualizer is called else the normal visualizer is called. This also creates a text file that contains the path length and time taken for each seed. All the output including plots is saved in the Output folder. The logs are saved in the Logs folder and deleted after the algorithm is run.

    Args:
        map_name (str): The name of the map to run the algorithm on.
        vis (bool): If set to true then the visualizer will be called to show the path after the algorithm is run.
        start (list): The start coordinates in the form ['x', 'y'].
        goal (list): The goal coordinates in the form ['x', 'y'].
        rbit (float): The radius of the ball (rbit*) to be used in the algorithm.
        samples (int): The number of samples to be used in the algorithm.
        dim (int): The dimension of the search space.
        seed (int): The number of seeds to be used in the algorithm.
        stop_time (int): The time in seconds after which the algorithm will stop.
        fast (bool): If set to true then the fast visualizer will be called to show the path after the algorithm is run.
    """
    # Get the working directory of this file.
    pwd = os.path.abspath(os.path.dirname(__file__))

    # Empty list to store the time taken and path length for each seed.
    time_taken_all = []
    path_lengths = []

    # The text file to store the path length and time taken for each seed.
    text_path = f"{pwd}/../Output/path_lengths_and_times_{map_name}.txt"
    # Make the directory if it does not exist.
    os.makedirs(os.path.dirname(text_path), exist_ok=True)

    # Run the algorithm for each seed.
    for seed in range(seed):
        # Set the seed for the random number generator.
        random.seed(seed)
        np.random.seed(seed)

        # Get the start and goal coordinates.
        start = []
        goal = []
        for i in range(opt.dim):
            # Convert the strings to floats.
            start.append(float(opt.start[i]))
            goal.append(float(opt.goal[i]))

        # Convert the lists to numpy arrays.
        start = np.array(start)
        goal = np.array(goal)

        # Create the log directory.
        log_dir = f"{pwd}/../Logs/{map_name}"
        os.makedirs(log_dir, exist_ok=True)

        # Get the occupancy grid map.
        map_path = f"{pwd}/../gridmaps/{map_name}.png"
        # Open the image and convert it to a numpy array.
        occ_map = np.array(Image.open(map_path))

        # Set the start and goal coordinates in the node class. This is done so that the node class can access the start and goal coordinates.
        node.start_arr = start
        node.goal_arr = goal

        # Create the start and goal nodes. The start node has a cost of 0.
        start_node = Node(tuple(start), gt=0)
        goal_node = Node(tuple(goal))

        # Create the map object and pass the occupancy grid map, start node, and goal node.
        map_obj = Map(start=start_node, goal=goal_node, occ_grid=occ_map)
        planner = None

        # Decide which planner to use based on the vis and fast flags. If the vis or fast flag is set then a log directory must be passed to the planner.
        if vis or fast:
            planner = bitstar(
                start=start_node,
                goal=goal_node,
                occ_map=map_obj,
                no_samples=samples,
                rbit=rbit,
                dim=dim,
                log_dir=log_dir,
                stop_time=stop_time,
            )
        # Else the planner is called without the log directory so that the logs are not saved.
        else:
            planner = bitstar(
                start=start_node,
                goal=goal_node,
                occ_map=map_obj,
                no_samples=samples,
                rbit=rbit,
                dim=dim,
                stop_time=stop_time,
            )
        # Make the plan.
        path, path_length, time_taken = planner.make_plan()

        print(
            f"{CGREEN2}Seed: {seed}\t\tFinal CI: {planner.ci}\tOld CI: {planner.old_ci}\tFinal Path Length: {path_length}\nPath:{CEND} {path}\n{CGREEN2}Time Taken per iteration:{CEND} {time_taken}\n{CEND}"
        )
        # Append the time taken and path length to the lists.
        time_taken_all.append(time_taken)
        path_lengths.append(path_length)
        # Convert the lists to strings and write them to the text file.
        time_taken_str = ", ".join([str(t) for t in time_taken])
        with open(text_path, "a") as f:
            f.write(
                f"Seed: {seed}\nPath Length: {path_length}\nTime Taken: {time_taken_str}\n"
            )

        if vis or fast:
            # Create the output directory. The directory name is the map name and the current date and time.
            output_dir = f"{pwd}/../Output/{map_name} - {str(datetime.now())}/"
            os.makedirs(output_dir, exist_ok=True)

            # Invert the occupancy grid map so that the free space is white and the occupied space is black. Weird bug in matplotlib which requires us to do this or use cv2.CvtColor function.
            inv_map = np.where((occ_map == 0) | (occ_map == 1), occ_map ^ 1, occ_map)

            # Create the visualizer object and pass the start and goal coordinates, the occupancy grid map, and the output directory.
            visualizer = Visualizer(start, goal, inv_map, output_dir)

            print(
                f"{CGREEN2}{CBOLD}{len(os.listdir(log_dir))} files in Log Directory:{CEND} {log_dir}/{sorted(os.listdir(log_dir))}\n"
            )
            # Read the log files.
            visualizer.read_json(log_dir, max_iter=np.inf)

            for i in range(len(os.listdir(log_dir))):
                # For each simulation draw with the fast visualizer if the fast flag is set.
                visualizer.draw(i, fast)
            # After all the simulations are drawn, set the title of the plot and show it is done drawing.
            visualizer.ax.set_title("BIT* - Final Path", fontsize=30)
            # Wait for the user to close the plot.
            plt.show()

        # Delete the log directory. This is done so the results of this experiment does not affect the results of the next experiment.
        print(
            f"{CRED2}================== Deleting log directory: {log_dir} =================={CEND}"
        )
        # Remove the log directory.
        shutil.rmtree(log_dir)


def parse_opt() -> argparse.Namespace:
    """Parse the command line arguments and return the arguments.

    Returns:
        argparse.Namespace: The arguments passed from the command line.
    """
    # Create the parser.
    parser = argparse.ArgumentParser()
    # Adding the other arguments.
    parser.add_argument(
        "--map_name",
        type=str,
        default="Default",
        help="Name of the map file. The map must be in the gridmaps folder. The map must be a png file. The map must be a black and white image where black is the obstacle space and white is the free space. Name of the map must be in the form 'map_name.png'. Do not enter the file extension in the map_name argument. The included maps are Default (or empty), Enclosure, Maze, Random, Symmetry, and Wall_gap. If none, will Default 100x100 empty grid will be used. Eg: --map_name Default",
    )
    parser.add_argument(
        "--vis",
        action="store_true",
        help="Whether or not to save and visualize outputs",
    )
    parser.add_argument("--start", nargs="+", help="Start coordinates. Eg: --start 0 0")
    parser.add_argument("--goal", nargs="+", help="Goal coordinates. Eg: --goal 99 99")
    parser.add_argument(
        "--rbit", type=float, default=10, help="Maximum Edge length. Eg: --rbit 10"
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=50,
        help="Number of new samples per iteration. Eg: --samples 50",
    )
    parser.add_argument(
        "--dim", type=int, default=2, help="Dimensions of working space. Eg: --dim 2"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed for reproducibility. Eg: --seed 1",
    )
    parser.add_argument(
        "--stop_time",
        type=int,
        default=60,
        help="When to stop the algorithm. Eg: --stop_time 60",
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Whether or not to only plot the final edge list of each iteration. Note: when this flag is set the vis is also considered set. Eg: --fast",
    )
    # Parse the arguments.
    opt = parser.parse_args()
    # Return the arguments.
    if opt.fast:
        opt.vis = True
    return opt


if __name__ == "__main__":
    # If no arguments are passed, then print the help message.
    if len(sys.argv) == 1:
        sys.argv.append("--help")

    # Parse the arguments.
    opt = parse_opt()
    print(f"{CYELLOW2}{opt}{CEND}")

    # Make sure the start and goal coordinates are of the correct dimension.
    assert len(opt.start) == opt.dim
    assert len(opt.goal) == opt.dim

    # Start the experiment.
    main(**vars(opt))
