import java.util.Iterator;
import java.util.NoSuchElementException;

public class AdjMatrixEdgeWeightedDigraph {
    private int V;
    private int E;
    private DirectedEdge[][] adj;
    
    /**
     * Initializes an empty edge-weighted digraph with <tt>V</tt> vertices and 0 edges.
     * param V the number of vertices
     * @throws java.lang.IllegalArgumentException if <tt>V</tt> < 0
     */
    public AdjMatrixEdgeWeightedDigraph(int V) {
        if (V < 0) throw new RuntimeException("Number of vertices must be nonnegative");
        this.V = V;
        this.E = 0;
        this.adj = new DirectedEdge[V][V];
    }

    /**
     * Initializes a random edge-weighted digraph with <tt>V</tt> vertices and <em>E</em> edges.
     * param V the number of vertices
     * param E the number of edges
     * @throws java.lang.IllegalArgumentException if <tt>V</tt> < 0
     * @throws java.lang.IllegalArgumentException if <tt>E</tt> < 0
     */
    public AdjMatrixEdgeWeightedDigraph(int V, int E) {
        this(V);
        if (E < 0) throw new RuntimeException("Number of edges must be nonnegative");
        if (E > V*V) throw new RuntimeException("Too many edges");

        // can be inefficient
        while (this.E != E) {
            int v = (int) (V * Math.random());
            int w = (int) (V * Math.random());
            double weight = Math.round(100 * Math.random()) / 100.0;
            addEdge(new DirectedEdge(v, w, weight));
        }
    }

    /**
     * Returns the number of vertices in the edge-weighted digraph.
     * @return the number of vertices in the edge-weighted digraph
     */
    public int V() {
        return V;
    }

    /**
     * Returns the number of edges in the edge-weighted digraph.
     * @return the number of edges in the edge-weighted digraph
     */
    public int E() {
        return E;
    }

    /**
     * Adds the directed edge <tt>e</tt> to the edge-weighted digraph (if there
     * is not already an edge with the same endpoints).
     * @param e the edge
     */
    public void addEdge(DirectedEdge e) {
        int v = e.from();
        int w = e.to();
        if (adj[v][w] == null) {
            E++;
            adj[v][w] = e;
        }
    }

    /**
     * Returns the directed edges incident from vertex <tt>v</tt>.
     * @return the directed edges incident from vertex <tt>v</tt> as an Iterable
     * @param v the vertex
     * @throws java.lang.IndexOutOfBoundsException unless 0 <= v < V
     */
    public Iterable<DirectedEdge> adj(int v) {
        return new AdjIterator(v);
    }

    // support iteration over graph vertices
    private class AdjIterator implements Iterator<DirectedEdge>, Iterable<DirectedEdge> {
        private int v, w = 0;
        public AdjIterator(int v) { this.v = v; }

        public Iterator<DirectedEdge> iterator() { return this; }

        public boolean hasNext() {
            while (w < V) {
                if (adj[v][w] != null) return true;
                w++;
            }
            return false;
        }

        public DirectedEdge next() {
            if (hasNext()) { return adj[v][w++];                 }
            else           { throw new NoSuchElementException(); }
        }

        public void remove()  { throw new UnsupportedOperationException();  }
    }

    /**
     * Returns a string representation of the edge-weighted digraph. This method takes
     * time proportional to <em>V</em><sup>2</sup>.
     * @return the number of vertices <em>V</em>, followed by the number of edges <em>E</em>,
     *   followed by the <em>V</em> adjacency lists of edges
     */
    public String toString() {
        String NEWLINE = System.getProperty("line.separator");
        StringBuilder s = new StringBuilder();
        s.append(V + " " + E + NEWLINE);
        for (int v = 0; v < V; v++) {
            s.append(v + ": ");
            for (DirectedEdge e : adj(v)) {
                s.append(e + "  ");
            }
            s.append(NEWLINE);
        }
        return s.toString();
    }

    /**
     * Unit tests the <tt>AdjMatrixEdgeWeightedDigraph</tt> data type.
     */
    public static void main(String[] args) {
        int V = Integer.parseInt(args[0]);
        int E = Integer.parseInt(args[1]);
        AdjMatrixEdgeWeightedDigraph G = new AdjMatrixEdgeWeightedDigraph(V, E);
        StdOut.println(G);
    }

}