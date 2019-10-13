# RDyn: graph benchmark handling community dynamics

[![Build Status](https://travis-ci.org/GiulioRossetti/RDyn.svg?branch=master)](https://travis-ci.org/GiulioRossetti/RDyn)
[![Coverage Status](https://coveralls.io/repos/github/GiulioRossetti/RDyn/badge.svg?branch=master)](https://coveralls.io/github/GiulioRossetti/RDyn?branch=master)
[![pyversions](https://img.shields.io/pypi/pyversions/rdyn.svg)](https://badge.fury.io/py/RDyn)
[![PyPI version](https://badge.fury.io/py/rdyn.svg)](https://badge.fury.io/py/RDyn)
[![Updates](https://pyup.io/repos/github/GiulioRossetti/RDyn/shield.svg)](https://pyup.io/repos/github/GiulioRossetti/RDyn/)
[![BCH compliance](https://bettercodehub.com/edge/badge/GiulioRossetti/RDyn?branch=master)](https://bettercodehub.com/)


Graph models provide an understanding of the dynamics of network formation and evolution; as a direct consequence, synthesizing graphs having controlled topology and planted partitions has been often identified as a strategy to describe benchmarks able to assess the performances of community discovery algorithm. However, one relevant aspect of real-world networks has been ignored by benchmarks proposed so far: community dynamics. As time goes by network communities rise, fall and may interact with each other generating merges and splits. Indeed, during the last decade dynamic community discovery has become a very active research field: in order to provide a coherent environment to test novel algorithms aimed at identifying mutable network partitions we introduce 
RDyn, an approach able to generates dynamic networks along with time-dependent ground-truth partitions having tunable quality.

## Citation
If you use our algorithm please cite the following works:

> Giulio Rossetti
> ["RDyn: graph benchmark handling community dynamics."](https://academic.oup.com/comnet/article-abstract/doi/10.1093/comnet/cnx016/3925036/text-RD-small-text-YN-graph-benchmark-handling) 
> Journal of Complex Networks 2017. 
> doi:10.1093/comnet/cnx016

## Installation

In order to install the package just download (or clone) the current project and copy the demon folder in the root of your application.

Alternatively use pip:
```bash
sudo pip install rdyn
```

# Generator Characteristics
For the main specifics of the RDyn algorithm please refer to the original publication.
The current release is shipped with *only* **conductance** as community quality function.

In order to speed up the computation a ``simplified`` flag has been introduced: setting it to ``True`` community quality is checked incrementally only on communities selected for merge/split operations.

**NB:** when set the ``simplified=True`` an approximation of the original process is executed: some network characteristics can be not preserved (e.g. degree distribution, interactions time to live).

# Execution

The algorithm can be used as standalone program as well as integrated in python scripts.

## Standalone

RDyn can be executed as standalone script with the following parameters:

**Positional arguments**

Name  |  Type | Description | Default 
-------------  | ------------- |------------- | -------------
nodes  | Integer | Number of nodes | 1000
iterations |Integer | Number of iterations| 1000
simplified | Boolean |Simplified execution | True

**Optional arguments**

Flag | Extended Name  |  Type | Description | Default 
-------------  | ------------- |------------- | ------------- | -------------
-d | --avg_degree | Integer | Average node degree | 15
-s | --sigma | Float | Percentage of node's edges within a community | .7
-l | --lbd | Float | Lambda community size distribution | 1
-a | --alpha | Float |Alpha degree distribution | 3
-p | --prob_action | Float |Probability of node action | 1
-r | --prob_renewal | Float |Probability of edge renewal | .8
-q | --quality_threshold | Float | Conductance quality threshold | .3
-n | --new_nodes | Float |Probability of node appearance | 0
-j | --delete_nodes | Float |Probability of node vanishing | 0
-e | --max_events | Integer |Max number of community events for stable iteration | 1

All parameters have a default value.
In order to generate a dynamic graph of 1000 nodes for 1000 iterations applying the simplified version of the algorithm just use:

```bash
python rdyn 1000 1000 True
```

## As python library

RDyn can be also executed within a python program.
In order to generate a dynamic network with the default parameter values just use the following snippet

```python
from rdyn import RDyn
rdb = RDyn()
rdb.execute(simplified=True)
```

To custumize the execution specify the usual parameters during object instantiation.

# Output

RDyn generates a folder named ``results`` and, for each specific model configuration a subfolder having the following naming convention:
 - ``nodes_iterations_avgDegree_sigma_renewal_qualityThreshold_maxEvts`` 
 - Example (default parameters): ``results/1000_1000_15_0.7_0.8_0.3_1``

Within such folder the following files are generated:
 - ``graph-*.txt``: Edgelist representation of the generated graph. One file for stable iteration.
 - ``communities-*.txt``: community description. One file for stable iteration.
 - ``events.txt``: summary of merge\split action per stable iteration.
 - ``interaction.txt``: dynamic graph description as edge stream.
 
The syntax of each class of output files is the following:

**Communities**

A community per line descibed as:
```bash
community_id	[node1, node2, node3, ..., nodeN]
```

**Events**

An block of events per stable iteration descibed as:

```bash
iteration_id:
 	Event1
 	Event2
 	...
 	EventN
```

Where the available events are:
 - ``START``, used for the first stable iteration
 - ``SPLIT id_origina_community [id_new_com1, id_new_com2]``
 - ``MERGE [id_old_com1, id_old_com2]`` the new com will inherit ``id_old_com_1``
 
**Interactions**
 
One interaction per line with the syntax:

``iteration_id	interaction_seq	operation	node1	node2``

Where:
 - ``iteration_id`` identify the iteration in which the interaction occurs
 - ``interaction_seq`` describe an absolute ordering among all the interactions
 - ``operation`` define if the interaction produces a new edge ``+`` or destroy an existing one ``-``
 - ``node1`` and ``node2`` are interaction endpoints
  
Example:
```bash
123	5361	+	385	390
123	5362	-	385	379
```

# Dependencies

RDyn is written in python and requires the following package to run:
- python>=2.7.11
- networkx==1.11
- numpy==1.11.1
- tqdm
- six
- future
