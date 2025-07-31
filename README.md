# pe-lcd
Pure Expansion based Local Community Detection (PE-LCD)
A community detection algorithm for identifying overlapping community structures in unweighted networks.
To be able to use this code, first compile it using the C++ compiler on a Linux/Unix machine as follows:
g++ -std=c++11 pe-lcd.cpp -o pe-lcd
Then an executable file named pe-lcd will be created. This can be executed, only for unweighted graphs, as follows:
./pe-lcd network_file -rh [option] 

The -rh option lets you set a parameter for the algorithm. If you don’t specify a value, it uses the default: 0.30. 
If you use this code for research purpose, please cite the following article:

@article{kumar2024pure,
  title={Pure expansion-based local community detection},
  author={Kumar, Abhinav and Kumar, Pawan and Dohare, Ravins},
  journal={International Journal of Data Science and Analytics},
  pages={1--19},
  year={2024},
  publisher={Springer}
}
