# Getting started
Required packages to compile the code, how to compile the code, and run an experiment.

### Prerequisites
Required packages are: git cmake libpq-dev libpqxx-dev

Installed by 
```sh
sudo apt-get update && sudo apt-get install git cmake lipq-dev libpqxx-dev
```

### Building
The code can be build by moving to the build directory and calling the CMakeList.txt
```sh
cd build/; cmake -DCMAKE_BUILD_TYPE=Release .. && make
```


When using vscode as IDE you can alternatively use the shortcut ctrl + shift + b, as defined in *.vscode/tasks.json*
Build code by ctrl + shift + b (in vscode) or $ cd build/; cmake -DCMAKE_BUILD_TYPE=Release .. && make

Note that the general configuration is to run in Release (optimized), for debugging change *RELEASE* to *DEBUG* and it is suggested to remove the O3 optimizer from the CMakeList.txt (string(APPEND CMAKE_CXX_FLAGS " -O3")) as otherwise variables may still be optimized out when using break points. ()

### Running any binary
Executing the compiled binaries can be done by running them from the build directory:
```sh
./experiments
or
./example
```
Running the code in *main/experiments.cpp* or *main/example.cpp* respectively, the first used for running actual experiments on many queries and collecting statistics, the latter is an example to display how to initialize SpatialSketch, perform insertions and queries.

Again if vscode is used as IDE, a shortcut to run experiments is done via the F5 key, as defined in *launch.json*

# Running Experiments
Running the experiments binary will run it with the default configurations as defined at the top of *experiments.cpp*, alternatively these can be passed as parameters when calling the binary, which is especially useful when running multiple experiments with different parameters.

*Relevant parameters and how to pass them coming soon*

### Statistics logging (to be updated)
Every experiment that starts up succesfully will output a Statistics.csv with the current date and time. The content currently consists of data related to the experiment, grid size, polygon area, polygon name, number of vertices. But also statistics derived during query execution, number of rectangles (a polygon is partitioned into), partition time, (synopsis) query time, number of (synopsis) subqueries, and the query answer. Thus note that the total query time for one polygon is the partition time + query time.

Additional notes, a query position sample will be a unique record in the statistics csv. For the partition/query samples a query is done multiple time and the averaged time is put in the record.


#### Datasets
Currently experiments.cpp expects datasets consisting of four columns <timestamp, ip, longitude, latitude> of types <long, long, int, int>, where the longitude and latitude fall in the range [0, N-1].
To generate synthetic datasets complying with this specification, which are also the datasets used in the experimental evaluation of the paper, the *Synthetic.ipynb* jupyter notebook in the *datasets* can be used. 

To generate the real-world dataset, also used in the paper, first the IP to geolocation dataset is required of which the free version is available at: https://dev.maxmind.com/geoip/geolite2-free-geolocation-data
This dataset will be joined with Caida network traces. In the paper the 2019 passive dataset was used, for which access can be requested via https://www.caida.org/catalog/datasets/passive_dataset/, although any kind of network trace dataset could suffice.

Combining these into a useful dataset requires the following steps:
1. The Caida traces have to be pre-processed from .pcap.gz to .csv. *datasets/ProcessCaida.py* does this for you by unzipping the file, splitting it into smaller .pcap chunks, of which the timestamp, source, and target ip are extracted and converted to .csv by use of tshark.
2. The unique IP addresses need to be computed, and joined with the GeoLite dataset in order to create a mapping from source IP address to geolocation.
3. Each pre-processed Caida trace has to be pulled through the mapping from source IP to geolocation and combined in one large dataset.
Steps 2. and 3. are covered in the *datasets/RealWorld.ipynb* notebook.

The current code uses pre-processed shape file to determine what region to query. Examples of these are found in *datasets/ExampleRegions* and in general contain some information on the shape, their vertices, and by how much they can be offset to stay within bounds. The polygon created by the vertices outlines the shape, therefore aligns with half integers. For example, if an inclusive square (1,1), (2,2) is queried, the respective vertices will be (0.5,0.5),(2.5,2.5) and the *Partitioner.cpp* will return partitions covering the square in integer space. To generate new shapes, either for simple shapes such as rectangles, or real-world shapes such as lakes, US states, or other geometries, the *PolygonPreprocessing.ipynb* can be used. 

# MARQ
MARQ implementation details (to be made more neat):
Reservoir sampling adapted from: http://hadjieleftheriou.com/sketches/index.html which falls under the GNU Lesser General Public License, therefore, our version of the MARQ implementation is made available under the same license.
Max heap: https://www.geeksforgeeks.org/binary-heap/ (free to share and adapt as described in https://www.geeksforgeeks.org/legal/copyright-information/)
MARQ implementation started from https://github.com/r4n4sh/sub_cube_MD_queries (no information on license), who took their range tree implementation from https://github.com/Lucaweihs/range-tree/tree/master (free to use under MIT license). which again is adapted in this repository.
Count-min and random numbers from ?

# Points to be added to this readme and repository in the upcoming week
- [ ] Experiment parameter passing
- [ ] Explanation on datasets and how to create real-world dataset with jupyter notebook as guide.
- [ ] Explanation of MARQ implementation, which code is ours and which is used from existing repositories.
- [ ] Experiments.cpp will be simplified for easier use, non-relevant implementation will be removed or moved to seperate files.