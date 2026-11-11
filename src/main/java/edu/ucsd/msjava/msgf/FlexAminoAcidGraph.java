package edu.ucsd.msjava.msgf;

import edu.ucsd.msjava.msutil.*;
import edu.ucsd.msjava.msutil.Modification.Location;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class FlexAminoAcidGraph extends DeNovoGraph<NominalMass> {
    public static final int MODIFIED_EDGE_PENALTY = 0;
    private ScoredSpectrum<NominalMass> scoredSpec;
    private Enzyme enzyme;
    private boolean direction;    // true: forward (e.g. Lys-C), false: reverse (e.g. Trypsin)
    private AminoAcidSet aaSet;
    private boolean useProtNTerm;
    private boolean useProtCTerm;

    private HashMap<NominalMass, ArrayList<DeNovoGraph.Edge<NominalMass>>> edgeMap;
    private HashMap<NominalMass, Integer> nodeScore;

    // Array copies of edgeMap and nodeScore indexed by (nominal mass - nodeIndexOffset),
    // built once the graph is complete (see buildNodeIndex)
    private int nodeIndexOffset;
    private int nodeIndexSize = 0;
    private ArrayList<DeNovoGraph.Edge<NominalMass>>[] edgeArr;
    private int[] nodeScoreArr;
    private boolean[] hasNodeScoreArr;

    private static AtomicInteger negativeCompNodeMassWarnCount;
    private static AtomicInteger negativeNodeMassWarnCount;

    private static AtomicInteger nullNodeCountGetNodeScore;
    private static AtomicInteger exceptionCountGetNodeScore;

    public FlexAminoAcidGraph(
            AminoAcidSet aaSet,
            int peptideMass,
            Enzyme enzyme,
            ScoredSpectrum<NominalMass> scoredSpec
    ) {
        this(aaSet, peptideMass, enzyme, scoredSpec, false, false);
    }

    public FlexAminoAcidGraph(
            AminoAcidSet aaSet,
            int peptideMass,
            Enzyme enzyme,
            ScoredSpectrum<NominalMass> scoredSpec,
            boolean useProteinNTerm,
            boolean useProteinCTerm
    ) {
        this.enzyme = enzyme;
        this.direction = scoredSpec.getMainIonDirection();
        this.scoredSpec = scoredSpec;
        this.aaSet = aaSet;
        this.useProtNTerm = useProteinNTerm;
        this.useProtCTerm = useProteinCTerm;

        super.source = new NominalMass(0);

        super.pmNode = new NominalMass(peptideMass);

        if (negativeNodeMassWarnCount == null) {
            negativeNodeMassWarnCount = new AtomicInteger();
        }

        if (negativeCompNodeMassWarnCount == null) {
            negativeCompNodeMassWarnCount = new AtomicInteger();
        }

        if (nullNodeCountGetNodeScore == null) {
            nullNodeCountGetNodeScore = new AtomicInteger();
        }

        if (exceptionCountGetNodeScore == null) {
            exceptionCountGetNodeScore = new AtomicInteger();
        }

        edgeMap = new HashMap<NominalMass, ArrayList<DeNovoGraph.Edge<NominalMass>>>();
        initBuildArrays();
        putBuildEdges(source, new ArrayList<DeNovoGraph.Edge<NominalMass>>());
        setForwardEdgesFromSource();
        setForwardEdgesFromIntermediateNodes();
        // the forward edges were collected in arrays indexed by nominal mass; copy them to edgeMap
        for (int i = 0; i < buildEdges.length; i++) {
            if (buildEdges[i] != null)
                edgeMap.put(buildNodes[i], buildEdges[i]);
        }
        buildEdges = null;
        buildNodes = null;
        super.intermediateNodes = new ArrayList<NominalMass>(edgeMap.keySet());
        Collections.sort(super.intermediateNodes);
        this.setBackwardEdgesFromSink();

        super.sinkNodes = new ArrayList<NominalMass>();
        sinkNodes.add(pmNode);

        computeNodeScores();
        buildNodeIndex();
    }

    // Construction-time arrays (indexed by nominal mass - buildOffset) holding the forward edges
    // and the node object of each node created so far; equivalent to edgeMap during construction
    private int buildOffset;
    private ArrayList<DeNovoGraph.Edge<NominalMass>>[] buildEdges;
    private NominalMass[] buildNodes;

    @SuppressWarnings("unchecked")
    private void initBuildArrays() {
        // Forward edges start at the source (0) or at an intermediate node 1..pm-1, so every node
        // created is >= min(0, smallest amino acid nominal mass) and < pm
        int minAANominalMass = 0;
        for (Location location : new Location[]{Location.N_Term, Location.Protein_N_Term, Location.C_Term, Location.Protein_C_Term, Location.Anywhere}) {
            for (AminoAcid aa : aaSet.getAAList(location))
                minAANominalMass = Math.min(minAANominalMass, aa.getNominalMass());
        }
        buildOffset = minAANominalMass;
        int size = pmNode.getNominalMass() - buildOffset + 1;
        if (size < 1)
            size = 1;
        buildEdges = new ArrayList[size];
        buildNodes = new NominalMass[size];
    }

    private ArrayList<DeNovoGraph.Edge<NominalMass>> getBuildEdges(int nominalMass) {
        int i = nominalMass - buildOffset;
        if (i < 0 || i >= buildEdges.length)
            return getBuildEdgesOutOfRange(nominalMass);
        return buildEdges[i];
    }

    private void putBuildEdges(NominalMass node, ArrayList<DeNovoGraph.Edge<NominalMass>> edges) {
        int i = node.getNominalMass() - buildOffset;
        if (i < 0 || i >= buildEdges.length) {
            putBuildEdgesOutOfRange(node, edges);
            return;
        }
        buildEdges[i] = edges;
        buildNodes[i] = node;
    }

    // Not expected to be used (see initBuildArrays); kept so behavior matches the map-based version
    private HashMap<NominalMass, ArrayList<DeNovoGraph.Edge<NominalMass>>> outOfRangeEdges = null;

    private ArrayList<DeNovoGraph.Edge<NominalMass>> getBuildEdgesOutOfRange(int nominalMass) {
        return outOfRangeEdges == null ? null : outOfRangeEdges.get(new NominalMass(nominalMass));
    }

    private void putBuildEdgesOutOfRange(NominalMass node, ArrayList<DeNovoGraph.Edge<NominalMass>> edges) {
        if (outOfRangeEdges == null)
            outOfRangeEdges = new HashMap<NominalMass, ArrayList<DeNovoGraph.Edge<NominalMass>>>();
        outOfRangeEdges.put(node, edges);
        edgeMap.put(node, edges);
    }

    @SuppressWarnings("unchecked")
    private void buildNodeIndex() {
        int min = pmNode.getNominalMass();
        int max = pmNode.getNominalMass();
        for (NominalMass n : edgeMap.keySet()) {
            min = Math.min(min, n.getNominalMass());
            max = Math.max(max, n.getNominalMass());
        }
        for (NominalMass n : nodeScore.keySet()) {
            min = Math.min(min, n.getNominalMass());
            max = Math.max(max, n.getNominalMass());
        }
        long size = (long) max - min + 1;
        if (size <= 0 || size > (1 << 22))
            return;    // keep using the maps
        int n = (int) size;
        ArrayList<DeNovoGraph.Edge<NominalMass>>[] edges = new ArrayList[n];
        for (java.util.Map.Entry<NominalMass, ArrayList<DeNovoGraph.Edge<NominalMass>>> e : edgeMap.entrySet())
            edges[e.getKey().getNominalMass() - min] = e.getValue();
        int[] scores = new int[n];
        boolean[] hasScore = new boolean[n];
        for (java.util.Map.Entry<NominalMass, Integer> e : nodeScore.entrySet()) {
            if (e.getValue() == null)
                continue;
            int i = e.getKey().getNominalMass() - min;
            scores[i] = e.getValue();
            hasScore[i] = true;
        }
        nodeIndexOffset = min;
        edgeArr = edges;
        nodeScoreArr = scores;
        hasNodeScoreArr = hasScore;
        nodeIndexSize = n;

        // cache the index of each edge's previous node on the edge (used by GeneratingFunction)
        for (ArrayList<DeNovoGraph.Edge<NominalMass>> list : edgeMap.values()) {
            if (list == null)
                continue;
            for (DeNovoGraph.Edge<NominalMass> edge : list)
                edge.setPrevNodeIndex(getNodeIndex(edge.getPrevNode()));
        }
    }

    @Override
    public int getNodeIndex(NominalMass node) {
        int i = node.getNominalMass() - nodeIndexOffset;
        return (i >= 0 && i < nodeIndexSize) ? i : -1;
    }

    @Override
    public int getNodeIndexSize() {
        return nodeIndexSize;
    }

    @Override
    public NominalMass getComplementNode(NominalMass node) {
        return new NominalMass(pmNode.getNominalMass() - node.getNominalMass());
    }

    @Override
    public ArrayList<DeNovoGraph.Edge<NominalMass>> getEdges(NominalMass curNode) {
        if (curNode != null) {
            int i = getNodeIndex(curNode);
            if (i >= 0)
                return edgeArr[i];
        }
        return edgeMap.get(curNode);
    }

    @Override
    public int getNodeScore(NominalMass node) {

        if (node == null) {
            int errorCount = nullNodeCountGetNodeScore.addAndGet(1);
            if (notifyError(errorCount)) {
                System.out.println("Note: null node encountered in getNodeScore");
            }
            return 0;
        }

        int i = getNodeIndex(node);
        if (i >= 0 && hasNodeScoreArr[i])
            return nodeScoreArr[i];

        try {
            return nodeScore.get(node);
        } catch (Exception ex) {
            int errorCount = exceptionCountGetNodeScore.addAndGet(1);
            if (notifyError(errorCount)) {
                System.out.println("Note: Exception in getNodeScore retrieving node at nominal mass " +
                        node.getNominalMass() + ": " + ex.getMessage());
            }
            return 0;
        }

    }

    @Override
    public int getScore(Peptide pep) {
        int score = 0;

        NominalMass prevNode = source;
        int nominalMass = 0;
        for (int i = 0; i < pep.size() - 1; i++) {
            AminoAcid aa;
            if (direction == true)
                aa = pep.get(i);
            else
                aa = pep.get(pep.size() - 1 - i);

            nominalMass += aa.getNominalMass();
            NominalMass curNode = new NominalMass(nominalMass);
            int nodeScore = getNodeScore(curNode);
            int edgeScore = scoredSpec.getEdgeScore(curNode, prevNode, aa.getMass());
            if (prevNode == source && direction == false && enzyme != null) {
                if (enzyme.isCleavable(aa))
                    edgeScore += aaSet.getPeptideCleavageCredit();
                else
                    edgeScore += aaSet.getPeptideCleavagePenalty();
            }
            prevNode = curNode;
            score += nodeScore + edgeScore;
        }
        if (direction == true && enzyme != null) {
            if (enzyme.isCleavable(pep.get(pep.size() - 1)))
                score += aaSet.getPeptideCleavageCredit();
            else
                score += aaSet.getPeptideCleavagePenalty();
        }
        if (direction == true)
            nominalMass += pep.get(pep.size() - 1).getNominalMass();
        else
            nominalMass += pep.get(0).getNominalMass();

        if (nominalMass != pmNode.getNominalMass())
            return Integer.MIN_VALUE;
        else
            return score;
    }

    @Override
    public int getScore(Annotation annotation) {
        int score = getScore(annotation.getPeptide());
        if (enzyme != null) {
            AminoAcid neighboringAA;
            if (enzyme.isCTerm())
                neighboringAA = annotation.getPrevAA();
            else
                neighboringAA = annotation.getNextAA();
            if (neighboringAA == null || enzyme.isCleavable(neighboringAA))
                score += aaSet.getNeighboringAACleavageCredit();
            else
                score += aaSet.getNeighboringAACleavagePenalty();
        }
        return score;
    }

    @Override
    public boolean isReverse() {
        return !direction;
    }

    @Override
    public AminoAcidSet getAASet() {
        return aaSet;
    }

    private void computeNodeScores() {
        nodeScore = new HashMap<NominalMass, Integer>();
        nodeScore.put(source, 0);

        boolean warnNegativeNodeMass = false;
        boolean warnNegativeCompNodeMass = false;

        for (int i = 1; i < intermediateNodes.size(); i++) {
            NominalMass node = intermediateNodes.get(i);
            NominalMass compNode = this.getComplementNode(node);
            if (node.getNominalMass() < 0 && !warnNegativeNodeMass) {
                warnNegativeNodeMass = true;
                // Mass of the node is negative
                // This can happen if we have a negative dynamic mod at the C-terminus, for example Lys-Loss
                int warnCount = negativeNodeMassWarnCount.addAndGet(1);
                if (notifyError(warnCount)) {
                    System.out.println("Note: negative node mass in computeNodeScores " +
                            "(count = " + Integer.toString(warnCount) + ")");
                }
            }
            if (compNode.getNominalMass() < 0 && !warnNegativeCompNodeMass) {
                warnNegativeCompNodeMass = true;
                int warnCount = negativeCompNodeMassWarnCount.addAndGet(1);
                if (notifyError(warnCount)) {
                    System.out.println("Note: negative compnode mass in computeNodeScores " +
                            "(count = " + Integer.toString(warnCount) + ")");
                }
            }
            int score;
            if (isReverse())
                score = scoredSpec.getNodeScore(compNode, node);
            else
                score = scoredSpec.getNodeScore(node, compNode);
            nodeScore.put(node, score);
        }
        for (NominalMass node : this.sinkNodes)
            nodeScore.put(node, 0);
    }

    private boolean notifyError(int errorCount) {
        if (errorCount < 5 || errorCount == 100 || errorCount == 1000 || errorCount % 10000 == 0) {
            return true;
        } else {
            return false;
        }
    }

    private void setForwardEdgesFromSource() {
        Location location;
        if (direction) {
            if (!this.useProtNTerm)
                location = Location.N_Term;
            else
                location = Location.Protein_N_Term;
        } else {
            if (!this.useProtCTerm)
                location = Location.C_Term;
            else
                location = Location.Protein_C_Term;
        }

        ArrayList<AminoAcid> aaList = aaSet.getAAList(location);
        makeForwardEdges(source, aaList, aaIndices(aaList), enzyme != null && direction == enzyme.isNTerm());
    }

    private int[] aaIndices(ArrayList<AminoAcid> aaList) {
        int[] indices = new int[aaList.size()];
        for (int k = 0; k < indices.length; k++)
            indices[k] = aaSet.getIndex(aaList.get(k));
        return indices;
    }

    private void setForwardEdgesFromIntermediateNodes() {
        ArrayList<AminoAcid> aaList = aaSet.getAAList(Location.Anywhere);
        int[] aaIndices = aaIndices(aaList);
        for (int i = 1; i < pmNode.getNominalMass(); i++) {
            int b = i - buildOffset;
            if (b >= 0 && b < buildEdges.length) {
                // only nodes that exist have outgoing edges; reuse the node object
                if (buildEdges[b] != null)
                    makeForwardEdges(buildNodes[b], aaList, aaIndices, false);
            } else {
                makeForwardEdges(new NominalMass(i), aaList, aaIndices, false);
            }
        }
    }

    private void setBackwardEdgesFromSink() {
        Location location;
        if (direction) {
            if (!this.useProtCTerm)
                location = Location.C_Term;
            else
                location = Location.Protein_C_Term;
        } else {
            if (!this.useProtNTerm)
                location = Location.N_Term;
            else
                location = Location.Protein_N_Term;
        }

        ArrayList<AminoAcid> aaList = aaSet.getAAList(location);

//		if(enzymaticCleavageOnly && direction != enzyme.isNTerm())
//			aaList = aaSet.getEnzymeAAList();
//		else
//			aaList = aaSet.getAAList(location);

        int peptideNominalMass = pmNode.getNominalMass();
        ArrayList<DeNovoGraph.Edge<NominalMass>> edges = new ArrayList<DeNovoGraph.Edge<NominalMass>>();
        for (AminoAcid aa : aaList) {
            NominalMass prevNode = new NominalMass(peptideNominalMass - aa.getNominalMass());
            if (edgeMap.containsKey(prevNode)) {
                DeNovoGraph.Edge<NominalMass> edge = new DeNovoGraph.Edge<NominalMass>(prevNode, aa.getProbability(), aaSet.getIndex(aa), aa.getMass());
                edges.add(edge);
                if (enzyme != null && direction != enzyme.isNTerm()) {
                    if (enzyme.isCleavable(aa))
                        edge.setCleavageScore(aaSet.getPeptideCleavageCredit());
                    else
                        edge.setCleavageScore(aaSet.getPeptideCleavagePenalty());
                }
                if (aa.isModified())
                    edge.setErrorScore(MODIFIED_EDGE_PENALTY);
            }
        }
        edgeMap.put(pmNode, edges);
    }

    private void makeForwardEdges(NominalMass curNode, ArrayList<AminoAcid> aaList, int[] aaIndices, boolean addCleavageScore) {
        if (getBuildEdges(curNode.getNominalMass()) == null)
            return;
        int curNominalMass = curNode.getNominalMass();
        for (int k = 0, numAA = aaList.size(); k < numAA; k++) {
            AminoAcid aa = aaList.get(k);
            int nextNodeNominalMass = curNominalMass + aa.getNominalMass();
            if (nextNodeNominalMass >= pmNode.getNominalMass())
                continue;
            ArrayList<DeNovoGraph.Edge<NominalMass>> edges = getBuildEdges(nextNodeNominalMass);
            NominalMass nextNode;
            if (edges == null) {
                nextNode = new NominalMass(nextNodeNominalMass);
                edges = new ArrayList<DeNovoGraph.Edge<NominalMass>>();
                putBuildEdges(nextNode, edges);
            } else {
                int b = nextNodeNominalMass - buildOffset;
                nextNode = (b >= 0 && b < buildNodes.length) ? buildNodes[b] : new NominalMass(nextNodeNominalMass);
            }

            DeNovoGraph.Edge<NominalMass> edge = new DeNovoGraph.Edge<NominalMass>(
                    curNode,
                    aa.getProbability(),
                    aaIndices[k],
                    aa.getMass());
//			if(curNode.getNominalMass() == 57 && nextNode.getNominalMass() == 114)
//				System.out.println("Debug");
            int errorScore = scoredSpec.getEdgeScore(nextNode, curNode, aa.getMass());
//			if(aa.isModified())
//				errorScore += MODIFIED_EDGE_PENALTY;
            if (errorScore < -100 || errorScore > 100) {
                System.err.println("Warning, invalid ErrorScore: " + errorScore);

                // Could abort the search
                // System.exit(-1);

                // Instead, use a score of -4
                errorScore = -4;
            }
            edge.setErrorScore(errorScore);
            if (addCleavageScore) {
                if (enzyme.isCleavable(aa))
                    edge.setCleavageScore(aaSet.getPeptideCleavageCredit());
                else
                    edge.setCleavageScore(aaSet.getPeptideCleavagePenalty());
            }

            edges.add(edge);
        }
    }

//	private void computeEdgeScores()
//	{
//		for(NominalMass curNode : intermediateNodes)
//		{
//			ArrayList<DeNovoGraph.Edge<NominalMass>> edges = edgeMap.get(curNode);
//			for(DeNovoGraph.Edge<NominalMass> edge : edges)
//			{
//				NominalMass prevNode = edge.getPrevNode();
//				int errorScore = scoredSpec.getEdgeScore(curNode, prevNode, edge.getEdgeMass());
//				assert(errorScore == edge.getErrorScore());
//				edge.setErrorScore(errorScore);
//			}
//			edgeMap.put(curNode, edges);
//		}
//	}	
}
