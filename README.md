# DynaMo: Dynamic Community Detection by Incrementally Maximizing Modularity

- If you find this dataset and code useful in your research, please consider citing:

`@ARTICLE{zhuang2019dynamo,
  author={Zhuang, Di and Chang, J. Morris and Li, Mingchen},
  journal={IEEE Transactions on Knowledge and Data Engineering}, 
  title={DynaMo: Dynamic Community Detection by Incrementally Maximizing Modularity}, 
  year={2021},
  volume={33},
  number={5},
  pages={1934-1945},
  doi={10.1109/TKDE.2019.2951419}}`

- A preprint version of our paper: [https://arxiv.org/abs/1709.08350](https://arxiv.org/abs/1709.08350).

- This is a sample code of our DynaMo paper, including all the necessary implementations for our experiments on the synthetic dynamic networks: 1) synthetic_exp_1.java: using the ground truth communities as the initial community structure; 2) synthetic_exp_2.java: using the result of static algorithm (i.e., Louvain) as the initial community structure.  

- The necessary implementations for our experiments on the real-world dynamic networks is in real_world_network_exp.java. The real-world dynamic network datasets will be added soon.  

#### Requirements
- We use [RDyn](https://github.com/GiulioRossetti/RDyn/blob/master/README.md) to generate the synthetic dynamic networks.

- We use [xmeasures](https://github.com/eXascaleInfolab/xmeasures) to evaluate the performance of the community detection results.

#### Quick Start (only tested using Eclipse on Ubuntu, JDK8)  

- Download or clone the whole repository.  
- Import our code as a Maven project into Eclipse.  
- Go to the folder ./xmeasures-master and build xmeasures using ```make release```.
- If you are using Python 2.7, nothing need to be changed. If you are using Python 3.XX, changing ``python RDyn-master/rdyn`` in synthetic_exp_1.java and synthetic_exp_2.java to ``python3 RDyn-master/rdyn``.
- Run synthetic_exp_1.java and synthetic_exp_2.java for the synthetic dynamic network experiments.
- Run real_world_network_exp.java for [the real-world dynamic network](https://www.kaggle.com/nogrady/realworld-dynamic-networks-dynamo) experiments.  

Contact info: Di Zhuang - ‬zhuangdi1990@gmail.com
