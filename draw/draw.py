import cities
import argparse
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

coords = []
for key in cities.cities:
    coords.append((cities.cities[key][1], cities.cities[key][0]))
    
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input')
parser.add_argument('-o', '--output')

args = parser.parse_args()

f = open(args.input, "r")
vals = int(f.readline())
order = [int(x) for x in f.readline().split()]
total_length = int(f.readline())

path = []
xs = []
ys = []
for v in order:
    path.append(coords[v])
    xs.append(coords[v][0])
    ys.append(coords[v][1])

polygon1 = Polygon(path, fill=False, color='black')

fig, ax = plt.subplots(1,1)

ax.add_patch(polygon1)

name = args.input[args.input.rfind('/') + 1 + len('output-'):-4]
plt.title(name + ', length = ' + str(total_length))

print(path)
plt.scatter(x=xs, y=ys, s=10, c='red')

plt.ylabel('latitude')
plt.xlabel('longitude')

plt.ylim(40,70)
plt.xlim(0,200)


plt.savefig(args.output)
