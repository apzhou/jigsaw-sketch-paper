# Jigsaw-Sketch

Jigsaw-Sketch is a new sketch for finding the top-k items. It is submitted to the International Conference on Very Large Databases ([VLDB](https://vldb.org/2022/)). This is the repository for the source codes.

## Introduction
Finding top-k frequent items is one of the most fundamental tasks in data stream processing. With limited memory, it is challenging to track top-k items from the single-pass and unbounded data streams. Furthermore, different from frequency estimation, top-k finding requires storing the IDs of top-k items, enlarging the difficulties of this task. Existing solutions employ sketches to fit in the limited memory and focus on tracking the hot items. However, they pay much less attention to improving the ID storage scheme. Moreover, existing ID storage schemes cannot simultaneously achieve high accuracy and high speed. In this work, we propose Jigsaw-Sketch based on a novel item ID storage scheme called jigsaw storage scheme, which can achieve both high accuracy and high speed with limited memory. Extensive experimental results show that Jigsaw-Sketch can achieve a 99.6% precision on finding top-1000 items with small memory, and is always the most accurate and fastest algorithm compared to the state-of-the-art.

## About this repo
Jigsaw-Sketch and the baselines are implemented on CPU, the codes are compiled using gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04). We employ the default optimization options when compiling codes. For each algorithm, the compiling command is as follows.
```shell
g++ main.cpp MurmurHash3.cpp -O -m64
```

