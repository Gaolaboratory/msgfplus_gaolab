package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.AminoAcidSet;
import edu.ucsd.msjava.msutil.Annotation;
import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.msutil.Peptide;

import java.util.ArrayList;

public abstract class DeNovoGraph<T extends Matter> {
    protected T source;
    protected T pmNode;
    protected ArrayList<T> sinkNodes;
    protected ArrayList<T> intermediateNodes;

    public T getSource() {
        return source;
    }

    public T getPMNode() {
        return pmNode;
    }

    public ArrayList<T> getSinkList() {
        return sinkNodes;
    }

    public ArrayList<T> getIntermediateNodeList() {
        return intermediateNodes;
    }

    public abstract boolean isReverse();

    public abstract int getScore(Peptide pep);

    public abstract int getScore(Annotation annotation);

    public abstract int getNodeScore(T node);

    //	public abstract int getEdgeScore(T curNode, T prevNode);
    public abstract ArrayList<Edge<T>> getEdges(T curNode);

    /**
     * Optional dense numbering of nodes, used by GeneratingFunction to keep its tables in arrays.
     * Returns an index in [0, getNodeIndexSize()) or -1 if the node has no index.
     * A node's index must not change during the graph's lifetime.
     */
    public int getNodeIndex(T node) {
        return -1;
    }

    /**
     * Size of the dense node numbering; 0 if the graph does not provide one.
     */
    public int getNodeIndexSize() {
        return 0;
    }

    public abstract T getComplementNode(T node);

    public abstract AminoAcidSet getAASet();

    public static class Edge<T extends Matter> {
        private T prevNode;
        private float probability;
        private int index;
        private float mass;

        // scores
        private int cleavageScore;
        private int errorScore;

        // Dense index of prevNode in the owning graph (see DeNovoGraph.getNodeIndex), or
        // NO_PREV_NODE_INDEX if not set by the graph
        public static final int NO_PREV_NODE_INDEX = Integer.MIN_VALUE;
        private int prevNodeIndex = NO_PREV_NODE_INDEX;

        public int getPrevNodeIndex() {
            return prevNodeIndex;
        }

        public void setPrevNodeIndex(int prevNodeIndex) {
            this.prevNodeIndex = prevNodeIndex;
        }

        public Edge(T prevNode, float probability, int index, float mass) {
            this.prevNode = prevNode;
            this.probability = probability;
            this.index = index;
            this.mass = mass;
        }

        public T getPrevNode() {
            return prevNode;
        }

        public void setCleavageScore(int cleavageScore) {
            this.cleavageScore = cleavageScore;
        }

        public void setErrorScore(int errorScore) {
            this.errorScore = errorScore;
        }

        public void setEdgeMass(float mass) {
            this.mass = mass;
        }

        public int getEdgeScore() {
            return cleavageScore + errorScore;
        }

        public int getErrorScore() {
            return errorScore;
        }

        public float getEdgeProbability() {
            return probability;
        }

        public int getEdgeIndex() {
            return index;
        }

        public float getEdgeMass() {
            return mass;
        }
    }
}
