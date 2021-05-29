# Developing best practices for single-cell analysis: data integration

## Intoduction

Single-cell sequencing is paving the way for precision medicine. It is the next srep towards making precision medicine more accurate. However, the analysis of single-cell data is incredibly complex with numerous distinct approaches resulting in more than 500 Python and R libraries existing today.

## Project goals

The goal of this project is to tackle the complexity of data analysis by identifying the best approaches. The single-cell transcriptomics analysis has multiple steps, but we have focused on data integration &mdash; a crucial step when working with clinical data coming from patients.

## Objectives

1. Single-cell sequencing data search and preprocessing
2. Selection of metric to compare performance of different algorithms
3. Simulation of artificial scRNA-seq datasets to compare the performances of the algorithms on real data and similar simulated
4. Comparative analysis of algorithms aimed at batch correction in single-cell sequencing data
5. Create a tool based on existing algorithms for data integration

## Materials and methods

1. We found 3 public datasets approximately 500, 15000, 80000 cells and 2, 4, 2 batches respectively. It may be found [here](https://hub.docker.com/r/jinmiaochenlab/batch-effect-removal-benchmarking). Also you need docker to be installed on your computer.

2. Silhouette-score with cosine distance was chosen as a metric. If we denote through a &mdash; the average distance from this object to objects from the same cluster, through b &mdash; the average distance from this object to objects from the nearest cluster (different from the one in which the object itself lies). In our case labels of clusters are *cell type* or *batch*. Then the silhouette of this object is called the value:

![equation](https://latex.codecogs.com/gif.latex?s%20%3D%20%7B%5Cfrac%20%7Bb-a%7D%7B%5Cmax%5C%28a%2Cb%5C%29%7D%7D)
    
3. Generated 3 artificial datasets 500, 15000, 80000 cells and 2, 5, 2 batches respectively. We used [SymSim](https://github.com/YosefLab/SymSim) R library, the code can be found [here]().

4. For the downstream work we have chosen 5 different approaches: [Combat](https://github.com/zhangyuqing/ComBat-seq) and [Regress out](https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.regress_out.html) which use linear models for batch correction and [MNN](https://github.com/MarioniLab/MNN2017), [BBKNN](https://github.com/Teichlab/bbknn) and [Scanorama](https://github.com/brianhie/scanorama) which in turn are looking for mutual nearest neighbors in others batches for every cell.

5. Based on these 5 libraries, we made a program that for each specific dataset can choose the best approach and apply it to correct the batch effect.

## Results

We analyzed 5 algorithms on 6 datasets (3 real and 3 simulated). We provided batch correction of expression matrix, then calculated silhouettes for each cell in data and finally run MU-test to compare silhouettes of baseline (data before correction) and corrected data. We did it in two ways:
1. Calculate silhouettes based on *cell type* labels and then run MU-test with *alternative* argument equal "less", because in this situation correction should increase our metric (cell types must form separate clusters)
2. Calculate silhouettes based on *batch* labels and then run MU-test with *alternative* argument equal "greater", because in this situation correction should decrease our metric (different batches should mix)

Figure 1 shows the results obtained. Blue color means that the silhouette according to the batch is significantly reduced after correction, and yellow, that the silhouette according to the cell type is significantly increased.

<p align="center" width="100%">
    <img width="60%" src="visualization/significance.png"> 
</p>
<p>
    <em>Fig. 1. Significance of differences between silhouettes.</em>
</p>

Also in Figure 2 you can see the time it took for each algorithm to work on the corresponding data. It is expected that MNN based on the search for mutual nearest neighbors worked most often longer, but at the same time its modifications (BBKNN and Scanorama), which preliminarily reduce the dimension of the data using canonical correlation analysis worked much faster.

<p align="center" width="100%">
    <img src="visualization/time_performance.png"> 
</p>
<p>
    <em>Fig. 2. Time performance of different algorithms on different datasets.</em>
</p>

Figure 3 depicts how clusters change after correction: cell types form separate clusters, and batches are stirred.

<p align="center" width="100%">
    <img src="visualization/before_after.png">
</p>
<p>
    <em>Fig. 3. UMAP before and after correction with Scanorama, labeled by cell type and batch.</em>
</p>

As a result, based on the existing algorithms for correcting the batch effect, we have created a program that can apply the algorithm chosen by the user or run all algorithms and select a list of the most successful ones (See **Usage** for more information about downloading and running).

## Usage

You can download and run program on our small example using following commands in your terminal:

```bash
git clone https://github.com/immunomind/bi2021spring.git
git checkout dev
cd bi2021spring/
pip install -r requirements.txt
cd source/
python data_integrator.py --adata ../example_data/example.h5ad \
                          --celltype celltype batch \
                          --batch batch \
                          --do_filter True \
                          --algtorun all \
                          --out example_result
```

Where **adata** is expression matrix with genes as columns and cells as rows and also with annotated cell and batch types (you can see more information in --help); **celltype** and **batch** are names of corresponding columns; **do_filter** tells if program apply casual filters (the gene is expressed in more than 100 cells and more than 20 genes have nonzero expression in a cell); **algtorun** means which algorithms will be used for correction or *all* and at the end there will be a list of the best ones; **out** is the directory with results.

We ran this program on Ubuntu 18.04 with 8 threads and 64 Gb RAM. With **algtorun** option equals *all* it has taken approximately 1 m on the small dataset, 10 m on the medium one and 75 m on the big. If your data is more than 100000 cells then we do not recommend using MNN due to too long working time.

## References

1. Haghverdi, L., Lun, A., Morgan, M. et al. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. Nat Biotechnol 36, 421–427 (2018).
2. Hie, B., Bryson, B. & Berger, B. Efficient integration of heterogeneous single-cell transcriptomes using Scanorama. Nat Biotechnol 37, 685–691 (2019).
3. Krzysztof Polański, Matthew D Young, Zhichao Miao, Kerstin B Meyer, Sarah A Teichmann, Jong-Eun Park, BBKNN: fast batch alignment of single cell transcriptomes, Bioinformatics 36, 3, (2020).
4. Tran, H.T.N., Ang, K.S., Chevrier, M. et al. A benchmark of batch-effect correction methods for single-cell RNA sequencing data. Genome Biol 21, 12 (2020).
5. Zhang, X., Xu, C. & Yosef, N. Simulating multiple faceted variability in single cell RNA sequencing. Nat Commun 10, 2611 (2019).
6. Yuqing Zhang, Giovanni Parmigiani, W Evan Johnson, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics  2, 3, (2020). 