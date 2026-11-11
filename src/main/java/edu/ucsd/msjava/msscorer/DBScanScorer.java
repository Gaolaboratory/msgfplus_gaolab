package edu.ucsd.msjava.msscorer;

import edu.ucsd.msjava.msgf.NominalMass;

// Fast scorer for DB search, consider edges
public class DBScanScorer extends FastScorer {

    private float[] nodeMass = null;
    private NewRankScorer scorer = null;
    private Partition partition;
    private float probPeak;
    private boolean isNodeMassPRM;    // prefix: true, suffix: false

    // Edge-score lookups precomputed per spectrum, so the search loop does not touch the
    // scorer's shared parameter maps (see NewRankScorer.getIonExistenceScore / getErrorScore)
    private final float[] ionExistenceScore = new float[4];
    private float[] errorScoreTable;

    public DBScanScorer(NewScoredSpectrum<NominalMass> scoredSpec, int peptideMass) {
        super(scoredSpec, peptideMass);
        this.scorer = scoredSpec.getScorer();

        nodeMass = new float[peptideMass];

        for (int i = 0; i < nodeMass.length; i++)
            nodeMass[i] = -1;

        isNodeMassPRM = scoredSpec.getMainIonDirection();
        // assign node mass
        nodeMass[0] = 0;
        for (int nominalMass = 1; nominalMass < nodeMass.length; nominalMass++) {
            nodeMass[nominalMass] = scoredSpec.getNodeMass(new NominalMass(nominalMass));
        }

        partition = scoredSpec.getPartition();
        probPeak = scoredSpec.getProbPeak();

        for (int i = 0; i < ionExistenceScore.length; i++)
            ionExistenceScore[i] = scorer.getIonExistenceScore(partition, i, probPeak);
        errorScoreTable = scorer.getErrorScoreTable(partition);
    }

    // fromIndex: inclusive, toIndex: exclusive
    @Override
    public int getScore(double[] prefixMassArr, int[] nominalPrefixMassArr, int fromIndex, int toIndex, int numMods) {
        int nodeScore = super.getScore(prefixMassArr, nominalPrefixMassArr, fromIndex, toIndex, numMods);
        int edgeScore = 0;
        if (!isNodeMassPRM)    // reverse
        {
            int nominalPeptideMass = nominalPrefixMassArr[toIndex - 1];
            for (int i = toIndex - 2; i >= fromIndex; i--)
                edgeScore += getEdgeScoreInt(nominalPeptideMass - nominalPrefixMassArr[i], nominalPeptideMass - nominalPrefixMassArr[i + 1], (float) (prefixMassArr[i + 1] - prefixMassArr[i]));
        } else                    // forward
        {
            for (int i = fromIndex; i <= toIndex - 2; i++)
                edgeScore += getEdgeScoreInt(nominalPrefixMassArr[i], nominalPrefixMassArr[i - 1], (float) (prefixMassArr[i] - prefixMassArr[i - 1]));
        }
        return nodeScore + edgeScore;
    }

    @Override
    public int getEdgeScore(NominalMass curNode, NominalMass prevNode, float theoMass) {
        return getEdgeScoreInt(curNode.getNominalMass(), prevNode.getNominalMass(), theoMass);
    }

    private int getEdgeScoreInt(int curNominalMass, int prevNominalMass, float theoMass) {
        // Debug
//		if(curNominalMass == 114 && prevNominalMass == 57)
//			System.out.println("Debug");
        if (curNominalMass >= nodeMass.length || prevNominalMass >= nodeMass.length || curNominalMass < 0 || prevNominalMass < 0)
            return 0;
        int ionExistenceIndex = 0;
        float curMass = nodeMass[curNominalMass];
        if (curMass >= 0)
            ionExistenceIndex += 1;
        float prevMass = nodeMass[prevNominalMass];
        if (prevMass >= 0)
            ionExistenceIndex += 2;

        float edgeScore = ionExistenceScore[ionExistenceIndex];
        if (ionExistenceIndex == 3) {
            float error = curMass - prevMass - theoMass;
            if (errorScoreTable != null)
                edgeScore += errorScoreTable[scorer.getErrorIndex(error)];
            else
                edgeScore += scorer.getErrorScore(partition, error);
        }
        return Math.round(edgeScore);
    }
}
