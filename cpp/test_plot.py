import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# create a 5x5 matrix of ones
matrix = np.ones((100, 100))

# create a figure and axes
fig, ax = plt.subplots()

# create a rectangle patch with position (2,2), width=1, height=1
rect = patches.Rectangle((1,2), 2, 1, linewidth=1, edgecolor='black', facecolor='black')
rect2 = patches.Rectangle((2,1), 1, 1, linewidth=1, edgecolor='black', facecolor='black')

# rect = patches.Rectangle((40,20), 40, 40, linewidth=1, edgecolor='black', facecolor='black')

# add the rectangle patch to the plot
ax.add_patch(rect)
ax.add_patch(rect2)

# set the x and y limits of the plot to include the rectangle patch
ax.set_xlim([0, 5])
ax.set_ylim([0, 5])


coords = [( 0, 0 ), ( 0.899529, 0.0831445 ), ( 2.15803, 0.514926 ), ( 3.28236, 1.11079 ), ( 3.5259, 1.963 ), ( 3.94908, 3.87296 ), ( 4, 4 )]
# extract x-coordinates and y-coordinates from the list of tuples
y_coords, x_coords = zip(*coords)
# x_coords, y_coords = zip(*coords)
# plot the points
ax.invert_yaxis()
ax.plot(x_coords, y_coords, '-ro')

# show the plot
plt.show()