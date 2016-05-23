import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.Stack;

public class FloydWarshall {
    private boolean hasNegativeCycle;  // is there a negative cycle?
    private static double[][] distTo;  // distTo[v][w] = length of shortest v->w path
    private static DirectedEdge[][] edgeTo;  // edgeTo[v][w] = last edge on shortest v->w path
    private static int V = 35; //Number of Nodes
    private static AdjMatrixEdgeWeightedDigraph G;
    static FloydWarshall spt;
    private static int[][] trafficDemand;
    private static int[][] cost;
	private static int[][] finalCost;
	private static int[][] capacity;
	private static Map<Integer, Double> networkCostMap = new HashMap<Integer, Double>();
	private static Map<Integer, Double> densityMap = new HashMap<Integer, Double>();
	
    /**
     * Computes a shortest paths tree from each vertex to to every other vertex in
     * the edge-weighted digraph <tt>G</tt>. If no such shortest path exists for
     * some pair of vertices, it computes a negative cycle.
     * @param G the edge-weighted digraph
     */
    public FloydWarshall(int k) {
    	initialize(k);
        int V = G.V();
        distTo = new double[V][V];
        edgeTo = new DirectedEdge[V][V];

        // initialize distances to infinity
        for (int v = 0; v < V; v++) {
            for (int w = 0; w < V; w++) {
                distTo[v][w] = Double.POSITIVE_INFINITY;
            }
        }

        // initialize distances using edge-weighted digraph's
        for (int v = 0; v < G.V(); v++) {
            for (DirectedEdge e : G.adj(v)) {
                distTo[e.from()][e.to()] = e.weight();
                edgeTo[e.from()][e.to()] = e;
            }
            // in case of self-loops
            if (distTo[v][v] >= 0.0) {
                distTo[v][v] = 0.0;
                edgeTo[v][v] = null;
            }
        }

        // Floyd-Warshall updates
        for (int i = 0; i < V; i++) {
            // compute shortest paths using only 0, 1, ..., i as intermediate vertices
            for (int v = 0; v < V; v++) {
                if (edgeTo[v][i] == null) continue;  // optimization
                for (int w = 0; w < V; w++) {
                    if (distTo[v][w] > distTo[v][i] + distTo[i][w]) {
                        distTo[v][w] = distTo[v][i] + distTo[i][w];
                        edgeTo[v][w] = edgeTo[i][w];
                    }
                }
                // check for negative cycle
                if (distTo[v][v] < 0.0) {
                    hasNegativeCycle = true;
                    return;
                }
            }
        }
    }

    /**
     * Is there a negative cycle?
     * @return <tt>true</tt> if there is a negative cycle, and <tt>false</tt> otherwise
     */
    public boolean hasNegativeCycle() {
        return hasNegativeCycle;
    }

    /**
     * Returns a negative cycle, or <tt>null</tt> if there is no such cycle.
     * @return a negative cycle as an iterable of edges,
     * or <tt>null</tt> if there is no such cycle
     */
    public Iterable<DirectedEdge> negativeCycle() {
        for (int v = 0; v < distTo.length; v++) {
            // negative cycle in v's predecessor graph
            if (distTo[v][v] < 0.0) {
                int V = edgeTo.length;
                EdgeWeightedDigraph spt = new EdgeWeightedDigraph(V);
                for (int w = 0; w < V; w++)
                    if (edgeTo[v][w] != null)
                        spt.addEdge(edgeTo[v][w]);
                EdgeWeightedDirectedCycle finder = new EdgeWeightedDirectedCycle(spt);
                assert finder.hasCycle();
                return finder.cycle();
            }
        }
        return null;
    }

    /**
     * Is there a path from the vertex <tt>s</tt> to vertex <tt>t</tt>?
     * @param s the source vertex
     * @param t the destination vertex
     * @return <tt>true</tt> if there is a path from vertex <tt>s</tt>
     * to vertex <tt>t</tt>, and <tt>false</tt> otherwise
     */
    public boolean hasPath(int s, int t) {
        return distTo[s][t] < Double.POSITIVE_INFINITY;
    }

    /**
     * Is there an edge from the vertex s to vertext ?
     * Mainly used for checking self edges
     * @param s the source vertex
     * @param t the destination vertex
     * @return true if there is an edge from vertex s to vertex t
     * and false otherwise
     */
    public boolean hasEdge(int s, int t) {
        return edgeTo[s][t] != null;
    }
    
    /**
     * Returns the length of a shortest path from vertex <tt>s</tt> to vertex <tt>t</tt>.
     * @param s the source vertex
     * @param t the destination vertex
     * @return the length of a shortest path from vertex <tt>s</tt> to vertex <tt>t</tt>;
     * <tt>Double.POSITIVE_INFINITY</tt> if no such path
     * @throws UnsupportedOperationException if there is a negative cost cycle
     */
    public double dist(int s, int t) {
        if (hasNegativeCycle())
            throw new UnsupportedOperationException("Negative cost cycle exists");
        return distTo[s][t];
    }

    /**
     * Returns a shortest path from vertex <tt>s</tt> to vertex <tt>t</tt>.
     * @param s the source vertex
     * @param t the destination vertex
     * @return a shortest path from vertex <tt>s</tt> to vertex <tt>t</tt>
     * as an iterable of edges, and <tt>null</tt> if no such path
     * @throws UnsupportedOperationException if there is a negative cost cycle
     */
    public Iterable<DirectedEdge> path(int s, int t) {
        if (hasNegativeCycle())
            throw new UnsupportedOperationException("Negative cost cycle exists");
        if (!hasPath(s, t)) return null;
        Stack<DirectedEdge> path = new Stack<DirectedEdge>();
        for (DirectedEdge e = edgeTo[s][t]; e != null; e = edgeTo[s][e.from()]) {
            path.push(e);
        }
        return path;
    }

    // check optimality conditions
    private boolean check(EdgeWeightedDigraph G, int s) {

        // no negative cycle
        if (!hasNegativeCycle()) {
            for (int v = 0; v < G.V(); v++) {
                for (DirectedEdge e : G.adj(v)) {
                    int w = e.to();
                    for (int i = 0; i < G.V(); i++) {
                        if (distTo[i][w] > distTo[i][v] + e.weight()) {
                            System.err.println("edge " + e + " is eligible");
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }


    /**
     *Initializes trafficDemand and cost matrices
     *and adds edge to the graph
     *@param k can range from 3 to 16*
     */
    public static void initialize(int k){
    	finalCost = new int[V][V];
    	capacity = new int[V][V];
    	trafficDemand = new int[V][V];
    	cost = new int[V][V];
    	int randomNumber;
    	Random random = new Random();
    	HashSet<Integer> hashSet = new HashSet<Integer>();
    	
    	G = new AdjMatrixEdgeWeightedDigraph(V);
    	
    	 for (int i = 0; i < V; i++) {
         	hashSet.clear();
         	//Assigning traffic demand values for all edges
     		for(int v = 0; v < V; v++){
     			trafficDemand[i][v] = random.nextInt(4);    			
     		}

     		//Generating k random indices
         	while(hashSet.size() < k){
         		randomNumber = random.nextInt(V);
         		if(randomNumber != i){
         			hashSet.add(randomNumber);	
         		}
     		}

         	//Assigning cost values for all edges
     		for(int v=0; v < V; v++){
     			if(hashSet.contains(v)){
     				cost[i][v] = 1;
     			}
     			else{
     				cost[i][v] = 250;
     			}
     		}
         	
     		//Adding edges to the graph
         	for(int v=0; v < V; v++){
         		double weight = cost[i][v];
         		G.addEdge(new DirectedEdge(i, v, weight));
         	}
         }

//         StdOut.println(G);
    }
    
    /**
     *Prints the trafficDemand matrix
     *@param input array with double values*
     */
    public void printGraph(double[][] input){
        for (int v = 0; v < G.V(); v++) {
            StdOut.printf("%5d ", v);
        }
        
        StdOut.println();
        
        //Print Traffic Demand Matrix
        for(int v=0; v<G.V(); v++){
        	StdOut.printf("%3d: ", v);
        	for(int w=0; w<G.V(); w++){
        		StdOut.print(input[v][w] + "     ");
        	}
        	StdOut.println();
        }
    }
    
    /**
     *Prints the trafficDemand matrix
     *@param input array with integer values*
     */
    public void printGraph(int[][] input){
        for (int v = 0; v < G.V(); v++) {
            StdOut.printf("%5d ", v);
        }
        
        StdOut.println();
        
        //Print Traffic Demand Matrix
        for(int v=0; v<G.V(); v++){
        	StdOut.printf("%3d: ", v);
        	for(int w=0; w<G.V(); w++){
        		StdOut.print(input[v][w] + "     ");
        	}
        	StdOut.println();
        }
    }
    
    /**
     * Computes the capacity and cost of each of the link
     * @returns networkCost cost of the network for the
     * value of k when the method is called
     */
    public double computeLinkCapacityAndCost(){
        int pathTraffic = 0;
        double networkCost = 0;
        int fromNode = 0, toNode = 0;
        
        //Prepare Capacity Matrix
        for(int v=0; v<G.V(); v++){
        	for(int w=0; w<G.V(); w++){
        		pathTraffic = trafficDemand[v][w];
        		for(DirectedEdge edge : spt.path(v, w)){
        			if(edge == null){
        				break;
        			}
        			fromNode = edge.from();
        			toNode = edge.to();
        			finalCost[fromNode][toNode] += (pathTraffic*cost[fromNode][toNode]);
        			capacity[fromNode][toNode] += pathTraffic;
        		}
        		networkCost += finalCost[fromNode][toNode];
        	}
        }
        return networkCost;
    }
    
    /**
     *Computes density of the network
     *@returns count/(V*(V-1)) density of the network for the
     *value of the k when the method is called 
     */
    public double computeDensity(){
    	int count = 0;
    	
    	//Count Number of Non-Zero capacity edges
    	for(int v=0; v<G.V(); v++){
        	for(int w=0; w<G.V(); w++){
        		if(capacity[v][w]!=0)
        			count++;
        	}
        }
    	
    	return (double) count/(V*(V-1));
    }
    
    /**
     * Unit tests the <tt>FloydWarshall</tt> data type.
     */
    public static void main(String[] args) {
    	for(int k=3; k<=16; k++){
    		StdOut.println("***************");
    		StdOut.println("*****K="+k+"****");
    		StdOut.println("***************");
//    		run Floyd-Warshall algorithm
            spt = new FloydWarshall(k);
            StdOut.println("\nTraffic Demand Matrix: ");
            spt.printGraph(trafficDemand);
            
            StdOut.println("\nInitial Cost Matrix: ");
            spt.printGraph(cost);

         // Print edges matrix
            StdOut.println("\nEdges Matrix: ");
            StdOut.printf("  ");
            for (int v = 0; v < G.V(); v++) {
                StdOut.printf("%6d ", v);
            }
            StdOut.println();
            for(int v=0; v<G.V(); v++){
            	StdOut.printf("%d: ", v);
            	for(int w=0; w<G.V(); w++){
            		if(edgeTo[v][w] == null){
            			StdOut.printf("  NIL");
            			continue;
            		}
            		StdOut.printf("%2d", edgeTo[v][w].from());
            	}
        		StdOut.println();
            }
            
            StdOut.println("\nAll Pairs SPs: ");
            // print negative cycle
            if (spt.hasNegativeCycle()) {
                StdOut.println("Negative cost cycle:");
                for (DirectedEdge e : spt.negativeCycle())
                    StdOut.println(e);
                StdOut.println();
            }

            // Print all-pairs shortest paths
            else {
                for (int v = 0; v < G.V(); v++) {
                    for (int w = 0; w < G.V(); w++) {
                        if (spt.hasPath(v, w)) {
                            StdOut.printf("%d to %d (%5.2f)  ", v, w, spt.dist(v, w));
                            for (DirectedEdge e : spt.path(v, w))
                                StdOut.print(e + "  ");
                            StdOut.println();
                        }
                        else {
                            StdOut.printf("%d to %d no path\n", v, w);
                        }
                    }
                }
            }
            
            networkCostMap.put(k, spt.computeLinkCapacityAndCost());
            
            densityMap.put(k, spt.computeDensity());
            
    	}
    	
    	//Printing network traffic map
        StdOut.println("\nNetwork Traffic Map: ");
        for(Map.Entry<Integer, Double> entry : networkCostMap.entrySet()){
        	StdOut.printf(entry.toString() + "\n");
        }
        
        //Printing network density map
        StdOut.println("\nDensity Map: ");
        for(Map.Entry<Integer, Double> entry : densityMap.entrySet()){
        	StdOut.printf(entry.toString() + "\n");
        }
    }
   }