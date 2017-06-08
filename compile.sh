cd src
javac -d ../class ParseNetwork.java
g++ AdjList2.cpp -o ../bin/AdjList2
g++ SortBlastp.cpp -o ../bin/SortBlastp
g++ BlastpTopHits.cpp -o ../bin/BlastpTopHits
javac -d ../class CombineAdjLists.java
javac -d ../class GraphletAlign.java
g++ AlignmentClustering.cpp -o ../bin/AlignmentClustering