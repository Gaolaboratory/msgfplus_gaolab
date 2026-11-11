package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.Annotation;
import edu.ucsd.msjava.msutil.Enzyme;
import edu.ucsd.msjava.msutil.Matter;
import edu.ucsd.msjava.suffixarray.SuffixArray;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;


public class GeneratingFunction<T extends Matter> implements GF<T> {
    private final DeNovoGraph<T> graph;

    private boolean backtrack = true;
    private boolean calcNumber = true;
    private boolean calcProb = true;
    private Enzyme enzyme = Enzyme.TRYPSIN;

    //	private int numScoreBinsPerNode = 1000;
    private int gfTableCapacity;

    private ScoreDist distribution = null;
    private BacktrackTable<T> backtrackTable = null;

    private class GFTable extends LinkedHashMap<T, ScoreDist> {

        /**
         *
         */
        private static final long serialVersionUID = 1L;
        private final int capacity;

        public GFTable(int capacity) {
            super(capacity + 1, 1.1f, false);
            this.capacity = capacity;
        }

        @Override
        protected boolean removeEldestEntry(Map.Entry<T, ScoreDist> eldest) {
            return size() > capacity;
        }
    }

    private HashMap<T, ScoreDist> fwdTable;
    private HashMap<T, Integer> minScoreTable = null;

    // Array versions of fwdTable and minScoreTable for nodes that have a dense index
    // (DeNovoGraph.getNodeIndex); nodes without one keep using the maps
    private ScoreDist[] fwdArr = null;
    private int[] minScoreArr = null;
    private boolean[] hasMinScoreArr = null;

    private boolean useNodeIndex() {
        // fwdTable evicts entries only if its capacity is below the number of nodes;
        // the arrays never evict, so use them only when no eviction can happen
        return graph.getNodeIndexSize() > 0
                && gfTableCapacity >= 1 + graph.intermediateNodes.size() + graph.sinkNodes.size();
    }

    private ScoreDist getFwd(T node) {
        if (fwdArr != null) {
            int i = graph.getNodeIndex(node);
            if (i >= 0)
                return fwdArr[i];
        }
        return fwdTable.get(node);
    }

    private void putFwd(T node, ScoreDist dist) {
        fwdTable.put(node, dist);
        if (fwdArr != null) {
            int i = graph.getNodeIndex(node);
            if (i >= 0)
                fwdArr[i] = dist;
        }
    }

    // Index of an edge's previous node: cached on the edge by the graph when available
    private int prevNodeIndex(DeNovoGraph.Edge<T> edge) {
        int i = edge.getPrevNodeIndex();
        if (i == DeNovoGraph.Edge.NO_PREV_NODE_INDEX)
            i = graph.getNodeIndex(edge.getPrevNode());
        return i;
    }

    // Same as getFwd(edge.getPrevNode())
    private ScoreDist getFwdOfPrev(DeNovoGraph.Edge<T> edge) {
        if (fwdArr != null) {
            int i = prevNodeIndex(edge);
            if (i >= 0)
                return fwdArr[i];
        }
        return fwdTable.get(edge.getPrevNode());
    }

    // Reusable buffers for setCurNode (replaces a per-node ArrayList of edges)
    @SuppressWarnings("unchecked")
    private DeNovoGraph.Edge<T>[] edgeBuffer = new DeNovoGraph.Edge[32];
    private ScoreDist[] prevDistBuffer = new ScoreDist[32];

    private boolean hasMinScore(T node) {
        if (minScoreArr != null) {
            int i = graph.getNodeIndex(node);
            if (i >= 0)
                return hasMinScoreArr[i];
        }
        return minScoreTable.get(node) != null;
    }

    // only valid if hasMinScore(node)
    private int getMinScoreOf(T node) {
        if (minScoreArr != null) {
            int i = graph.getNodeIndex(node);
            if (i >= 0)
                return minScoreArr[i];
        }
        return minScoreTable.get(node);
    }

    private void putMinScore(T node, int score) {
        if (minScoreArr != null) {
            int i = graph.getNodeIndex(node);
            if (i >= 0) {
                minScoreArr[i] = score;
                hasMinScoreArr[i] = true;
                return;
            }
        }
        minScoreTable.put(node, score);
    }

    private boolean isGFComputed = false;

//	private HashMap<T, Integer> srmScore = null;

    public GeneratingFunction(DeNovoGraph<T> graph) {
        this.graph = graph;
        this.gfTableCapacity = 1 + graph.intermediateNodes.size() + graph.sinkNodes.size();
    }

    // Builder
    public GeneratingFunction<T> doNotBacktrack() {
        this.backtrack = false;
        return this;
    }

    public GeneratingFunction<T> doNotCalcNumber() {
        this.calcNumber = false;
        return this;
    }

    public GeneratingFunction<T> doNotCalcProb() {
        this.calcProb = false;
        return this;
    }

    public GeneratingFunction<T> enzyme(Enzyme enzyme) {
        this.enzyme = enzyme;
        return this;
    }

    //	public GeneratingFunction<T> numScoreBinsPerNode(int numBins)	{ this.numScoreBinsPerNode = numBins; return this; }
    public GeneratingFunction<T> gfTableCapacity(int gfTableCapacity) {
        this.gfTableCapacity = gfTableCapacity;
        return this;
    }

    public boolean backtrack() {
        return backtrack;
    }

    public boolean calcNumber() {
        return calcNumber;
    }

    public boolean calcProb() {
        return calcProb;
    }

    public Enzyme getEnzyme() {
        return enzyme;
    }

    //	public int getNumScoreBinsPerNode()	{ return numScoreBinsPerNode; }
    public boolean isGFComputed() {
        return this.isGFComputed;
    }

    public DeNovoGraph<T> getGraph() {
        return graph;
    }

    protected HashMap<T, ScoreDist> getFwdTable() {
        return fwdTable;
    }

    protected BacktrackTable<T> getBacktrackTable() {
        return backtrackTable;
    }

    public int getScore(Annotation annotation) {
        return graph.getScore(annotation);
    }

    public int getEnergy(Annotation annotation) {
        return getMaxScore() - getScore(annotation);
    }

    public double getSpectralProbability(Annotation annotation) {
        int score = getScore(annotation);
        return getSpectralProbability(score);
    }

    // score: inclusive
    public double getSpectralProbability(int score) {
        if (!this.distribution.isProbSet())
            return 100;
        return distribution.getSpectralProbability(score);
    }

    public double getNumEqualBetterPeptides(Annotation annotation) {
        int score = getScore(annotation);
        return getNumEqualOrBetterPeptides(score);
    }

    public double getNumEqualOrBetterPeptides(int score) {
        if (!this.distribution.isNumSet())
            return -1.;
        return distribution.getNumEqualOrBetterPeptides(score);
    }

    public double getDictionarySize(float specProb) {
        return getNumEqualOrBetterPeptides(getThresholdScore(specProb));
    }

    // returns t where totalProb(t) > specProb && totalProb(t+1) <= specProb
    public static int getThresholdScore(float specProb, ScoreDist distribution) {
        if (!distribution.isProbSet())
            return -1;
        float totalProb = 0;

        for (int t = distribution.getMaxScore() - 1; t >= distribution.getMinScore(); t--) {
            totalProb += distribution.getProbability(t);
            if (totalProb > specProb)
                return t;
        }
        return -1;
    }

    // returns t where totalProb(t) > specProb && totalProb(t+1) <= specProb
    public int getThresholdScore(float specProb) {
        return getThresholdScore(specProb, distribution);
    }

    public ScoreDist getScoreDist() {
        return distribution;
    }

    /**
     * Generate reconstructions with score "score" and have match with "sa" and put it in "reconstructions".
     *
     * @param score           the score of reconstructions to be generated
     * @param reconstructions a container where reconstructions will be stored
     * @param sa              suffix array that will filter reconstructions
     * @return
     */
    private void generateReconstructions(int score, ArrayList<String> reconstructions, SuffixArray sa) {
        if (backtrackTable == null)
            return;
        if (enzyme == null) {
            for (T sink : graph.getSinkList())
                backtrackTable.getReconstructions(sink, score, "", reconstructions, sa);
        } else {
            //TODO: add prefix info?
            for (T sink : graph.getSinkList())
                backtrackTable.getReconstructions(sink, score - graph.getAASet().getNeighboringAACleavageCredit(), "R.", reconstructions, sa);
            for (T sink : graph.getSinkList())
                backtrackTable.getReconstructions(sink, score - graph.getAASet().getNeighboringAACleavagePenalty(), "L.", reconstructions, sa);
        }
    }

    public String getOneReconstruction(int score) {
        if (backtrackTable == null)
            return null;
        return backtrackTable.getOneReconstruction(graph.getPMNode(), score, "");
    }

    public ArrayList<String> getReconstructions(int score) {
        ArrayList<String> reconstructions = new ArrayList<String>();
        generateReconstructions(score, reconstructions, null);
        return reconstructions;
    }

    public ArrayList<String> getReconstructionsEqualOrAboveScore(int score) {
        ArrayList<String> reconstructions = new ArrayList<String>();
        for (int t = this.getMaxScore() - 1; t >= score; t--)
            generateReconstructions(t, reconstructions, null);
        return reconstructions;
    }

    public ArrayList<String> getDictionary(float specProbThreshold) {
        assert (calcProb);
        int threshold = getThresholdScore(specProbThreshold);
        return getReconstructionsEqualOrAboveScore(threshold + 1);
    }

    public ArrayList<String> getReconstructions(float specProbThreshold, float numRecsThreshold, boolean isNumInclusive, SuffixArray sa) {
        assert (calcProb && calcNumber);
        ArrayList<String> recs = new ArrayList<String>();
        int threshold = getThresholdScore(specProbThreshold);
        float numRecs = 0;
        for (int t = getMaxScore() - 1; t > threshold; t--) {
            numRecs += distribution.getNumberRecs(t);
            if (!isNumInclusive) {
                if (numRecs <= numRecsThreshold)
                    generateReconstructions(t, recs, sa);
                else
                    break;
            } else {
                generateReconstructions(t, recs, sa);
                if (numRecs >= numRecsThreshold)
                    break;
            }
        }
        return recs;
    }

    public int getMinScore() {
        return this.distribution.getMinScore();
    }

    public int getMaxScore() {
        return this.distribution.getMaxScore();
    }

    public void setUpScoreThreshold(int score) {
        minScoreTable = new HashMap<T, Integer>();
        if (useNodeIndex()) {
            minScoreArr = new int[graph.getNodeIndexSize()];
            hasMinScoreArr = new boolean[graph.getNodeIndexSize()];
        } else {
            minScoreArr = null;
            hasMinScoreArr = null;
        }
        if (enzyme != null)
            score -= graph.getAASet().getNeighboringAACleavageCredit();

        for (T sink : graph.getSinkList()) {
            putMinScore(sink, score);
            for (DeNovoGraph.Edge<T> edge : graph.getEdges(sink)) {
                int newPrevMinScore = score - edge.getEdgeScore();
                lowerPrevMinScore(edge, newPrevMinScore);
            }
        }

        ArrayList<T> intermediateNodeList = graph.getIntermediateNodeList();

        for (int i = intermediateNodeList.size() - 1; i >= 0; i--) {
            T curNode = intermediateNodeList.get(i);
            if (!hasMinScore(curNode))
                continue;
            int curScore = getMinScoreOf(curNode);
            int curNodeScore = graph.getNodeScore(curNode);
            ArrayList<DeNovoGraph.Edge<T>> edges = graph.getEdges(curNode);
            for (int e = 0, numEdges = edges.size(); e < numEdges; e++) {
                DeNovoGraph.Edge<T> edge = edges.get(e);
                int newPrevMinScore = curScore - (curNodeScore + edge.getEdgeScore());
                lowerPrevMinScore(edge, newPrevMinScore);
            }
        }
    }

    // minScore(prevNode) = min(minScore(prevNode), newPrevMinScore), creating the entry if absent
    private void lowerPrevMinScore(DeNovoGraph.Edge<T> edge, int newPrevMinScore) {
        if (minScoreArr != null) {
            int i = prevNodeIndex(edge);
            if (i >= 0) {
                if (!hasMinScoreArr[i] || minScoreArr[i] > newPrevMinScore) {
                    minScoreArr[i] = newPrevMinScore;
                    hasMinScoreArr[i] = true;
                }
                return;
            }
        }
        T prevNode = edge.getPrevNode();
        Integer prevMinScore = minScoreTable.get(prevNode);
        if (prevMinScore == null || prevMinScore > newPrevMinScore)
            minScoreTable.put(prevNode, newPrevMinScore);
    }

    public boolean computeGeneratingFunction() {
        ScoreDistFactory factory = new ScoreDistFactory(calcNumber, calcProb);
        // initialization of the source
        ScoreDist sourceDist = factory.getInstance(0, 1);
        if (calcNumber)
            sourceDist.setNumber(0, 1);
        if (calcProb)
            sourceDist.setProb(0, 1);
        fwdTable = new GFTable(gfTableCapacity);
        fwdArr = useNodeIndex() ? new ScoreDist[graph.getNodeIndexSize()] : null;
        putFwd(graph.getSource(), sourceDist);
        if (backtrack) {
            backtrackTable = new BacktrackTable<T>(graph);
            BacktrackPointer sourcePointer = new BacktrackPointer(0, 1, 0);
            sourcePointer.setBacktrack(0, 0);
            backtrackTable.put(graph.getSource(), sourcePointer);
        }

        // dynamic programming, source node (i=0) is excluded
        ArrayList<T> intermediateNodeList = graph.getIntermediateNodeList();

        for (int i = 1; i < intermediateNodeList.size(); i++) {
            T curNode = intermediateNodeList.get(i);
            setCurNode(curNode, factory);
        }

        // process dest node
        int minScore = Integer.MAX_VALUE;
        int maxScore = Integer.MIN_VALUE;

        for (T curNode : graph.getSinkList()) {
            setCurNode(curNode, factory);
            ScoreDist curDist = getFwd(curNode);
            if (curDist == null)    // curNode is not connected from the source
                continue;
            if (curDist.getMinScore() < minScore)
                minScore = curDist.getMinScore();
            if (curDist.getMaxScore() > maxScore)
                maxScore = curDist.getMaxScore();
        }

        if (maxScore <= minScore)
            return false;

        if (minScore < -10000 || maxScore > 10000) {
            System.err.println("Error! MinScore: " + minScore + ", MaxScore: " + maxScore + " ");
            System.exit(-1);
        }

        // merge distributions of dest nodes
        ScoreDist mergedDist = factory.getInstance(minScore, maxScore);
        for (T sinkNode : graph.getSinkList()) {
            if (calcNumber)
                mergedDist.addNumDist(fwdTable.get(sinkNode), 0);
            if (calcProb)
                mergedDist.addProbDist(fwdTable.get(sinkNode), 0, 1);
        }

        // process neighboring amino acid
        ScoreDist finalDist;
        if (enzyme != null && enzyme.getResidues() != null) {
            int neighboringAACleavageCredit = graph.getAASet().getNeighboringAACleavageCredit();
            int neighboringAACleavagePenalty = graph.getAASet().getNeighboringAACleavagePenalty();
            finalDist = factory.getInstance(mergedDist.getMinScore() + neighboringAACleavagePenalty, mergedDist.getMaxScore() + neighboringAACleavageCredit);
            if (calcNumber) {
                finalDist.addNumDist(mergedDist, neighboringAACleavageCredit, enzyme.getResidues().length);
                finalDist.addNumDist(mergedDist, neighboringAACleavagePenalty, graph.getAASet().size() - enzyme.getResidues().length);
            }
            if (calcProb) {
                finalDist.addProbDist(mergedDist, neighboringAACleavageCredit, graph.getAASet().getProbCleavageSites());
                finalDist.addProbDist(mergedDist, neighboringAACleavagePenalty, 1 - graph.getAASet().getProbCleavageSites());
            }
        } else {
            finalDist = mergedDist;
        }

        this.distribution = finalDist;
        isGFComputed = true;
        return true;
    }

    // scoreThreshold : inclusive
    public HashMap<T, Float> getDestProfile(int scoreThreshold) {
        assert (calcNumber);
        HashMap<T, Float> destProf = new HashMap<T, Float>();
        for (T sinkNode : graph.getSinkList()) {
            float num = 0;
            ScoreDist dist = fwdTable.get(sinkNode);
            for (int t = dist.getMaxScore() - 1; t >= dist.getMinScore() && t >= scoreThreshold; t--)
                num += dist.getNumberRecs(t);
            if (num > 0)
                destProf.put(sinkNode, num);
        }
        return destProf;
    }

    private void setCurNode(T curNode, ScoreDistFactory scoreDistFactory) {
        int curNodeScore = graph.getNodeScore(curNode);
        int curMaxScore = Integer.MIN_VALUE;
        int curMinScore;
        if (minScoreTable == null)
            curMinScore = Integer.MAX_VALUE;
        else {
            if (!hasMinScore(curNode))
                return;
            curMinScore = getMinScoreOf(curNode);
        }

        // determine minScore and maxScore
        // (edges with a predecessor distribution, and those distributions, are kept in edgeBuffer / prevDistBuffer)
        ArrayList<DeNovoGraph.Edge<T>> allEdges = graph.getEdges(curNode);
        int numAllEdges = allEdges.size();
        if (edgeBuffer.length < numAllEdges) {
            edgeBuffer = java.util.Arrays.copyOf(edgeBuffer, Math.max(numAllEdges, edgeBuffer.length * 2));
            prevDistBuffer = java.util.Arrays.copyOf(prevDistBuffer, edgeBuffer.length);
        }
        int numEdges = 0;
        for (int e = 0; e < numAllEdges; e++) {
            DeNovoGraph.Edge<T> edge = allEdges.get(e);
            ScoreDist prevDist = getFwdOfPrev(edge);
            if (prevDist != null) {
                int edgeScore = edge.getEdgeScore();
                int combinedScore = curNodeScore + edgeScore;
                if (prevDist.getMaxScore() + combinedScore > curMaxScore)
                    curMaxScore = prevDist.getMaxScore() + combinedScore;
                if (minScoreTable == null) {
                    if (prevDist.getMinScore() + combinedScore < curMinScore)
                        curMinScore = prevDist.getMinScore() + combinedScore;
                }
                edgeBuffer[numEdges] = edge;
                prevDistBuffer[numEdges] = prevDist;
                numEdges++;
            }
        }
        if (curMinScore >= curMaxScore)
            return;

        if (curMinScore < -10000) {
            System.err.println("Warning, MinScore is abnormally low; "
                    + "MinScore: " + curMinScore + ", MaxScore: " + curMaxScore + ", "
                    + "CurNode: " + curNode.getNominalMass() + ", CurNodeScore: " + curNodeScore);
            // Could abort processing
            // System.exit(-1);
            // Instead, skip this node
            return;
        }

        if (curMaxScore > 10000) {
            System.err.println("Warning, MaxScore is abnormally high; "
                    + "MinScore: " + curMinScore + ", MaxScore: " + curMaxScore + ", "
                    + "CurNode: " + curNode.getNominalMass() + ", CurNodeScore: " + curNodeScore);
            // Could abort processing
            // System.exit(-1);
            // Instead, skip this node
            return;
        }

        ScoreDist curDist = scoreDistFactory.getInstance(curMinScore, curMaxScore);
        BacktrackPointer backPointer = null;
        if (backtrack)
            backPointer = new BacktrackPointer(curMinScore, curMaxScore, curNodeScore);
        for (int e = 0; e < numEdges; e++) {
            DeNovoGraph.Edge<T> edge = edgeBuffer[e];
            T prevNode = edge.getPrevNode();
            ScoreDist prevDist = prevDistBuffer[e];    // = getFwd(prevNode); fwd tables do not change within this loop
            if (prevDist != null) {
                int edgeScore = edge.getEdgeScore();
                int combinedScore = curNodeScore + edgeScore;

                if (calcNumber)
                    curDist.addNumDist(prevDist, combinedScore, 1);
                if (calcProb)
                    curDist.addProbDist(prevDist, combinedScore, edge.getEdgeProbability());
                if (backtrack) {
                    BacktrackPointer prevPointer = backtrackTable.get(prevNode);
                    backPointer.addBacktrackPointers(prevPointer, edge.getEdgeIndex(), edgeScore);
                }
            }
        }
        if (calcProb) {
            if (curDist.getProbability(curDist.maxScore - 1) == 0)    // to avoid underflow
            {
                assert (false) : "Underflow! " + curNode.getNominalMass() + " " + curDist.getProbability(curDist.maxScore - 1);
                curDist.setProb(curDist.maxScore - 1, Float.MIN_VALUE);
            }
        }
        putFwd(curNode, curDist);
        if (backtrack)
            backtrackTable.put(curNode, backPointer);
    }
}
