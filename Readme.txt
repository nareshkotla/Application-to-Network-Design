An Application to Network Design
--------------------------------
The theme of this project is to implement the basic network design model and experiment with it.


 As input, it receives the number of nodes (N), the trac demand values (bij) between pairs of nodes, and the unit cost values for
the potential links (aij).

 As output, the program generates a network topology, with capacities assigned to the links using the shortest path based fast solution method.



Platform/Compiler: Windows, javac
------------------


Files:
------
FloydWarshall.java - The program on top of which code has been added based on the requirement (see references below). 



Compilation:
------------
- Download the files specified in the References.
- Modify the file FloydWarshall.java with the code specified in the Report.
- Navigate to the folder through the command prompt.
- Type the below command for compiling all the java files in the folder.

javac *.java


Run:
----
- Type the below command

java FloydWarshall 

Output:
-------
For each value of k 
- Prints TrafficDemand, Cost matrices of the network.
- Prints All Pairs Shortest Paths with costs for each edge along the path.
- Prints Network Cost and Network Density values.

References: 
-----------
Thomas H. Cormen Charles E. Leiserson Ronald L. Rivest Clifford Stein, Section 25.2, The Floyd-Warshall algorithm Introduction to Algorithms (Third ed.)
http://algs4.cs.princeton.edu/44sp