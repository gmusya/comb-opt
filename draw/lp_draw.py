import cities
import argparse
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

coords = []
xs = []
ys = []
for key in cities.cities:
    coords.append((cities.cities[key][1], cities.cities[key][0]))
    
for x, y in coords:
    xs.append(x)
    ys.append(y)
    
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input')
parser.add_argument('-o', '--output')

args = parser.parse_args()

f = open(args.input, "r")
vals = int(f.readline())
weights = []
for i in range(vals):
    weights.append([float(x) for x in f.readline().split()])

for i in range(vals):
    for j in range(vals):
        if weights[i][j] <= 0:
            continue
        if weights[i][j] == 1:
            plt.plot([coords[i][0], coords[j][0]], [coords[i][1], coords[j][1]], c='black')
        elif weights[i][j] == 0.5:
            plt.plot([coords[i][0], coords[j][0]], [coords[i][1], coords[j][1]], c='red')
        else:
            plt.plot([coords[i][0], coords[j][0]], [coords[i][1], coords[j][1]], c='blue')
# plt.plot([1, 2], [3, 4])
# plt.plot([3, 4], [1, 2])

# ax.add_patch(polygon1)

# name = args.input[args.input.rfind('/') + 1 + len('output-'):-4]
# plt.title(name + ', length = ' + str(total_length))

# print(path)
plt.scatter(x=xs, y=ys, s=10, c='red')

# plt.ylabel('latitude')
# plt.xlabel('longitude')

# plt.ylim(40,70)
# plt.xlim(0,200)


plt.savefig(args.output)
