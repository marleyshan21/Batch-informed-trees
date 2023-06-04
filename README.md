# Batch Informed Trees

[**Project Report**](Batch_informed_trees.pdf) | [**Project Slides**](https://1drv.ms/p/s!AiYkRRrhfnuCkVzbbknyJE-WsNeL?e=4lao3E) | [**Turtlebot Demo - (Wall Gap Scenario)**](https://youtu.be/JurQ1YLwIVY)

## BIT* Algorithm - Introduction


 Path Planning in Robotics has always relied on simple approximations to identify solutions. This is due to the difficulty to find one due to the high dimensional nature of the problem. 
 
 Generally, we can divide the approximations into 2 types: 
 Search-based and Sampling-based. 
 
 Heuristics are used by search-based planners like A* to effectively search across graphs, but their efficiency is limited by the resolution of the selected approximation. To approximate the problem, sample-based planners such as RRT* employ random sampling. Here, the resolution can be increased until we find a suitable solution. These random samples approximate the region in all directions at the same time, making the search ineffective. 
 
 A recent approach called Batch Informed Trees (BIT*) combines the strengths of both Search-based sampling-based planners. Heuristics and Sampling is used by BIT* to alternate between searching and approximating. In this work, we have used the pseudo-code from the paper and coded the algorithm from scratch, and tested its performance in R2 space for different motion planning scenarios.

## Installation

```bash
git clone https://github.com/marleyshan21/Batch-informed-trees.git
cd BIT-star
pip install -r requirements.txt
```

## Usage

<!-- <details> -->
<summary> Python </summary>

To use our python implementation of BIT-star, we provide a run.py file with options to specify all the arguments passed to the algorithm. A full list of options can be seen by running 
```bash
cd python
python3 run.py --help
```

### Example Usage
To run BIT-star on a default map and only obtain the final path once the stop time has been reached, simply run:
```bash
python3 run.py --map_name Default --start 0 0 --goal 99 99 --seed 1 --stop_time 20
```

To run BIT-star on a more complex map (Maze) and obtain visualizations once the stop time ahs been reached, run:

```bash
python3 run.py --map_name Maze --start 0 0 --goal 99 99 --seed 1 --stop_time 60 --vis --fast
```

We also provide options to change the rbit (maximum edge length), and number of samples when running. You can also visualize every edge addition and removal by disabling the --fast option.
<!-- </details> -->



## Examples
- Empty Scenario
<img src="https://github.com/marleyshan21/Batch-informed-trees/blob/master/Output/default_gif.gif"  alt="Empty Scenario">

- Enclosure Scenario
<img src="https://github.com/marleyshan21/Batch-informed-trees/blob/master/Output/enc_gif.gif"  alt="Enclosure Scenario">

- Maze Scenario
<img src="https://github.com/marleyshan21/Batch-informed-trees/blob/master/Output/maze_gif.gif"  alt="Maze Scenario">

- Wall Scenario
<img src="https://github.com/marleyshan21/Batch-informed-trees/blob/master/Output/wall_gif.gif"  alt="Wall Scenario">


## Collaboration

Done in collaboration with - Saharajit Anantharamakrishnan, Drake Moore, Thejaswini Gopiparthi, and Francis Jacob Kalliath