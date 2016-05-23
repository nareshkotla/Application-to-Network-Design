import java.util.HashSet;
import java.util.Random;


public class Test {
	static int nodeCount = 5;
	static int[][] trafficDemand = new int[nodeCount][nodeCount];
	static int[][] cost = new int[nodeCount][nodeCount];
	static int k = 3; //Need to run for k=3,...16
	static double[][] distance = new double[nodeCount][nodeCount];
	static double[][] edge = new double[nodeCount][nodeCount];
	
	public Test(){
		Random random = new Random();
    	HashSet<Integer> hashSet = new HashSet<Integer>();
    	
    	//Assigning traffic demand values for all edges
    	for(int i=0; i < nodeCount; i++){
    		for(int j=0; j < nodeCount; j++){
    			trafficDemand[i][j] = random.nextInt(4);    			
    		}
    	}

    	//Assigning cost values for all edges
    	while(hashSet.size()<k){
			hashSet.add(random.nextInt(k));
		}

    	for(int i=0; i < nodeCount; i++){
    		for(int j=0; j < nodeCount; j++){
    			if(hashSet.contains(j)){
    				cost[i][j] = 1;
    			}
    			else{
    				cost[i][j] = 250;
    			}
    		}
    	}
    	
    	 // initialize distances to infinity
        for (int v = 0; v < nodeCount; v++) {
            for (int w = 0; w < nodeCount; w++) {
                 distance[v][w] = Double.POSITIVE_INFINITY;
            }
        }
        
        // initialize edges to NULL/0
        for (int v = 0; v < nodeCount; v++) {
            for (int w = 0; w < nodeCount; w++) {
                 edge[v][w] = 0;
            }
        }
	}
	
	private void computeAllPairsSP() {
        for(int i=0; i < nodeCount; i++){
        	for(int v=0; v<nodeCount; v++){
        		for (int w = 0; w < nodeCount; w++) {
                    if (cost[v][w] > cost[v][i] + cost[i][w]) {
                        cost[v][w] = cost[v][i] + cost[i][w];
                        edge[v][w] = edge[i][w];
                    }
                }
        	}
        }
	}

	public static void main(String args[]){
    	Test obj =new Test();
    	obj.computeAllPairsSP();
    	
    	for (int v = 0; v < nodeCount; v++) {
            StdOut.printf("%5d ", v);
        }
    	
    	StdOut.println();
    	for (int v = 0; v < nodeCount; v++) {
            StdOut.printf("%d: ", v);
            for (int w = 0; w < nodeCount; w++) {
                System.out.printf("%5d", cost[v][w]);
            }
            System.out.println();
        }
    	
    	// print all-pairs shortest paths
        for (int v = 0; v < nodeCount; v++) {
            for (int w = 0; w < nodeCount; w++) {
            	System.out.printf("%d to %d (%5.2f)  ", v, w, cost[v][w]);
            }
            System.out.println();
        }
	}
}
