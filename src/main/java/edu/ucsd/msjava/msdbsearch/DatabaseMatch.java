package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.msutil.ActivationMethod;

import java.util.SortedSet;
import java.util.TreeSet;

public class DatabaseMatch extends Match {
    private int index;
    private byte length;

    // optional
    private boolean isProteinNTerm;
    private boolean isProteinCTerm;
    private boolean isNTermMetCleaved = false;

    private Float psmQValue = null;
    private Float pepQValue = null;

    // for degenerate peptides: other suffix-array indices of the same peptide, in insertion order
    // (a plain int array; getIndices() returns them sorted and de-duplicated together with index)
    private int[] extraIndices;
    private int numExtraIndices;

    public DatabaseMatch(
            int index,
            byte length,
            int score,
            float peptideMass,
            int nominalPeptideMass,
            int charge,
            String pepSeq,
            ActivationMethod[] actMethodArr
    ) {
        super(score, peptideMass, nominalPeptideMass, charge, pepSeq, actMethodArr);
        this.index = index;
        this.length = length;
        isProteinNTerm = false;
        isProteinCTerm = false;
    }

    public DatabaseMatch setProteinNTerm(boolean isProteinNTerm) {
        this.isProteinNTerm = isProteinNTerm;
        return this;
    }

    public DatabaseMatch setProteinCTerm(boolean isProteinCTerm) {
        this.isProteinCTerm = isProteinCTerm;
        return this;
    }

    public DatabaseMatch setNTermMetCleaved(boolean isNTermMetCleaved) {
        this.isNTermMetCleaved = isNTermMetCleaved;
        return this;
    }

    public boolean isNTermMetCleaved() {
        return this.isNTermMetCleaved;
    }

    public void setPSMQValue(float psmQValue) {
        this.psmQValue = psmQValue;
    }

    public Float getPSMQValue() {
        return this.psmQValue;
    }

    public void setPepQValue(Float pepQValue) {
        this.pepQValue = pepQValue;
    }

    public Float getPepQValue() {
        return this.pepQValue;
    }

    public void addIndex(int index) {
        if (extraIndices == null)
            extraIndices = new int[4];
        else if (numExtraIndices == extraIndices.length)
            extraIndices = java.util.Arrays.copyOf(extraIndices, numExtraIndices * 2);
        extraIndices[numExtraIndices++] = index;
    }

    public SortedSet<Integer> getIndices() {
        SortedSet<Integer> indices = new TreeSet<Integer>();
        indices.add(index);
        for (int i = 0; i < numExtraIndices; i++)
            indices.add(extraIndices[i]);
        return indices;
    }

    public int getIndex() {
        return index;
    }

    public int getLength() {
        return length;
    }

    public boolean isProteinNTerm() {
        return isProteinNTerm;
    }

    public boolean isProteinCTerm() {
        return isProteinCTerm;
    }

    public int hashCode() {
        return index * length;
    }

    public boolean equals(Object obj) {
        if (obj instanceof DatabaseMatch) {
            DatabaseMatch other = (DatabaseMatch) obj;
            if (index == other.index && length == other.length)
                return true;
        }
        return false;
    }
}
