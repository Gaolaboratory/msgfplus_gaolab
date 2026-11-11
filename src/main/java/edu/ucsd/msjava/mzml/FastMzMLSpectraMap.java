package edu.ucsd.msjava.mzml;

import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.msutil.SpectrumAccessorBySpecIndex;

import java.util.ArrayList;

/**
 * SpectrumAccessorBySpecIndex backed by {@link FastMzMLReader}; behaves like {@link MzMLSpectraMap}.
 */
public class FastMzMLSpectraMap implements SpectrumAccessorBySpecIndex {
    private final FastMzMLReader reader;

    public FastMzMLSpectraMap(FastMzMLReader reader) {
        this.reader = reader;
    }

    @Override
    public Spectrum getSpectrumBySpecIndex(int specIndex) {
        return reader.getSpectrumBySpecIndex(specIndex);
    }

    @Override
    public Spectrum getSpectrumById(String specId) {
        return reader.getSpectrumById(specId);
    }

    @Override
    public String getID(int specIndex) {
        return reader.getID(specIndex);
    }

    @Override
    public Float getPrecursorMz(int specIndex) {
        return reader.getPrecursorMz(specIndex);
    }

    @Override
    public String getTitle(int specIndex) {
        return null;
    }

    @Override
    public ArrayList<Integer> getSpecIndexList() {
        return reader.getSpecIndexList();
    }
}
