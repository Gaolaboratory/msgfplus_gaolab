package edu.ucsd.msjava.mzml;

import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.CvParamInfo;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Spectrum;
import edu.ucsd.msjava.mzid.Constants;
import uk.ac.ebi.jmzidml.model.mzidml.CvParam;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * Lightweight mzML reader used in place of jmzML for reading spectra.
 * <p>
 * The file is scanned once to record the byte range, id and index of every &lt;spectrum&gt; element.
 * Individual spectra are then parsed on demand with StAX and converted exactly as
 * {@link SpectrumConverter#getSpectrumFromJMzMLSpec} converts a jmzML spectrum.
 * Spectra using features this reader does not handle (referenceable parameter groups, numpress or other
 * compression, unknown binary precision, missing attributes) are read with jmzML instead, so the results
 * are the same as before; jmzML is only initialized if such a spectrum is encountered.
 * <p>
 * Thread-safe: spectra can be read concurrently.
 */
public class FastMzMLReader {

    /**
     * Set the system property msgfplus.fastMzML=false to always use jmzML.
     */
    public static boolean isEnabled() {
        return !"false".equalsIgnoreCase(System.getProperty("msgfplus.fastMzML"));
    }

    private final File specFile;
    private final Charset charset;
    private final FileChannel channel;

    // Spectrum elements in file order
    private final long[] startOffset;
    private final int[] length;
    private final String[] ids;
    private final int[] indexAttr;    // value of the index attribute, -1 if missing
    private final HashMap<String, Integer> idToPos;
    private final HashMap<Integer, Integer> indexToPos;

    // Same defaults as MzMLAdapter
    private final int minMSLevel = 2;        // inclusive
    private final int maxMSLevel = Integer.MAX_VALUE;

    private final ThreadLocal<XMLInputFactory> xmlInputFactory = ThreadLocal.withInitial(() -> {
        XMLInputFactory f = XMLInputFactory.newInstance();
        f.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, false);
        f.setProperty(XMLInputFactory.SUPPORT_DTD, false);
        f.setProperty(XMLInputFactory.IS_COALESCING, true);
        return f;
    });

    // jmzML, created only if needed
    private volatile MzMLAdapter fallbackAdapter = null;
    private volatile MzMLSpectraMap fallbackMap = null;
    private volatile boolean fallbackNoticePrinted = false;

    private CvParam spectrumIDFormatCvParam = null;

    public FastMzMLReader(File specFile) throws IOException {
        this.specFile = specFile;
        this.channel = FileChannel.open(specFile.toPath(), StandardOpenOption.READ);
        this.charset = detectCharset();

        SpectrumIndexBuilder builder = new SpectrumIndexBuilder();
        builder.scan();
        int n = builder.starts.size();
        startOffset = new long[n];
        length = new int[n];
        ids = new String[n];
        indexAttr = new int[n];
        idToPos = new HashMap<String, Integer>(n * 2);
        indexToPos = new HashMap<Integer, Integer>(n * 2);
        for (int i = 0; i < n; i++) {
            startOffset[i] = builder.starts.get(i);
            long len = builder.ends.get(i) - startOffset[i];
            if (len > Integer.MAX_VALUE)
                throw new IOException("Spectrum element too large in " + specFile);
            length[i] = (int) len;
            ids[i] = builder.idList.get(i);
            indexAttr[i] = builder.indexList.get(i);
            if (ids[i] != null)
                idToPos.putIfAbsent(ids[i], i);
            if (indexAttr[i] >= 0)
                indexToPos.putIfAbsent(indexAttr[i], i);
        }
    }

    public int getNumSpectra() {
        return startOffset.length;
    }

    public int getMinMSLevel() {
        return minMSLevel;
    }

    public int getMaxMSLevel() {
        return maxMSLevel;
    }

    /**
     * @param specIndex 1-based spectrum index (index attribute + 1)
     */
    public String getID(int specIndex) {
        Integer pos = indexToPos.get(specIndex - 1);
        return pos == null ? null : ids[pos];
    }

    public ArrayList<Integer> getSpecIndexList() {
        ArrayList<Integer> list = new ArrayList<Integer>(new TreeSet<Integer>(indexToPos.keySet()));
        return list;
    }

    /**
     * Spectrum by 1-based index, or null if not found or outside the MS level range
     */
    public Spectrum getSpectrumBySpecIndex(int specIndex) {
        String specId = getID(specIndex);
        if (specId != null)
            return getSpectrumById(specId);
        return null;
    }

    /**
     * Spectrum by native ID, or null if not found or outside the MS level range
     */
    public Spectrum getSpectrumById(String specId) {
        Integer pos = idToPos.get(specId);
        if (pos == null)
            return null;
        Spectrum spec = readSpectrumAt(pos);
        if (spec == null)
            return null;
        if (spec.getMSLevel() < minMSLevel || spec.getMSLevel() > maxMSLevel)
            return null;
        return spec;
    }

    public Float getPrecursorMz(int specIndex) {
        String specId = getID(specIndex);
        if (specId == null)
            return null;
        Integer pos = idToPos.get(specId);
        if (pos == null)
            return null;
        try {
            Elem spectrum = parseElement(pos);
            return precursorMzOf(spectrum);
        } catch (Exception e) {
            return getFallbackMap().getPrecursorMz(specIndex);
        }
    }

    /**
     * Iterator over spectra in file order within the MS level range (like MzMLSpectraIterator).
     */
    public Iterator<Spectrum> iterator() {
        return new FastIterator();
    }

    /**
     * Parses and converts the spectrum at the given position (file order); falls back to jmzML if needed.
     * Does not apply the MS level filter.
     */
    Spectrum readSpectrumAt(int pos) {
        try {
            Elem spectrum = parseElement(pos);
            return convert(spectrum);
        } catch (UnsupportedFeature | XMLStreamException | IOException | RuntimeException e) {
            return readWithJmzML(pos, e);
        }
    }

    /**
     * Fast path only, without the jmzML fallback (for testing)
     */
    Spectrum readSpectrumAtWithoutFallback(int pos) throws Exception {
        return convert(parseElement(pos));
    }

    String getIdAt(int pos) {
        return ids[pos];
    }

    private Spectrum readWithJmzML(int pos, Exception reason) {
        if (!fallbackNoticePrinted) {
            fallbackNoticePrinted = true;
            System.out.println("Note: reading some spectra with jmzML (" + reason.getClass().getSimpleName()
                    + (reason.getMessage() != null ? ": " + reason.getMessage() : "") + ")");
        }
        String id = ids[pos];
        if (id == null)
            return null;
        MzMLSpectraMap map = getFallbackMap();
        uk.ac.ebi.jmzml.model.mzml.Spectrum jmzSpec;
        try {
            jmzSpec = map.getJmzMLSpectrumById(id);
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        if (jmzSpec == null)
            return null;
        return SpectrumConverter.getSpectrumFromJMzMLSpec(jmzSpec);
    }

    private MzMLSpectraMap getFallbackMap() {
        if (fallbackMap == null) {
            synchronized (this) {
                if (fallbackMap == null) {
                    fallbackAdapter = new MzMLAdapter(specFile);
                    fallbackMap = new MzMLSpectraMap(fallbackAdapter);
                }
            }
        }
        return fallbackMap;
    }

    // ---------------------------------------------------------------------------------------------
    // Spectrum ID format (replaces MzMLAdapter.getSpectrumIDFormatCvParam without building a jmzML index)

    public synchronized CvParam getSpectrumIDFormatCvParam() {
        if (spectrumIDFormatCvParam != null)
            return spectrumIDFormatCvParam;

        try (InputStream in = new BufferedInputStream(new FileInputStream(specFile), 1 << 16)) {
            XMLStreamReader r = xmlInputFactory.get().createXMLStreamReader(in, charset.name());
            // path of local names from the root
            ArrayList<String> path = new ArrayList<String>();
            try {
                while (r.hasNext()) {
                    int event = r.next();
                    if (event == XMLStreamConstants.START_ELEMENT) {
                        String name = localName(r.getLocalName());
                        path.add(name);
                        if (name.equals("run"))
                            break;    // sourceFileList is in the header
                        if (name.equals("cvParam") && path.size() >= 4
                                && path.get(path.size() - 2).equals("sourceFile")
                                && path.get(path.size() - 3).equals("sourceFileList")
                                && path.get(path.size() - 4).equals("fileDescription")) {
                            String acc = r.getAttributeValue(null, "accession");
                            String name2 = r.getAttributeValue(null, "name");
                            if (acc != null) {
                                long accNum = Long.parseLong(acc.substring(acc.lastIndexOf(':') + 1));
                                if (accNum >= 1000768 && accNum <= 1000777
                                        || accNum == 1000823 || accNum == 1000824 || accNum == 1000929
                                        || accNum == 1001508 || accNum == 1001526 || accNum == 1001528
                                        || accNum == 1001531 || accNum == 1001532 || accNum == 1001559
                                        || accNum == 1001562 || accNum == 1002818 || accNum == 1001480
                                        || accNum == 1002303 || accNum == 1002532 || accNum == 1002898
                                ) {
                                    spectrumIDFormatCvParam = Constants.makeCvParam(acc, name2);
                                    return spectrumIDFormatCvParam;
                                }
                            }
                        }
                    } else if (event == XMLStreamConstants.END_ELEMENT) {
                        path.remove(path.size() - 1);
                    }
                }
            } finally {
                r.close();
            }
        } catch (Exception e) {
            // fall through to jmzML below
        }
        // Not found (or not parseable): let jmzML decide, including its error handling
        if (fallbackAdapter == null)
            getFallbackMap();
        spectrumIDFormatCvParam = fallbackAdapter.getSpectrumIDFormatCvParam();
        return spectrumIDFormatCvParam;
    }

    private static String localName(String name) {
        int colon = name.indexOf(':');
        return colon < 0 ? name : name.substring(colon + 1);
    }

    // ---------------------------------------------------------------------------------------------
    // Index building: one pass over the bytes of the file

    private Charset detectCharset() throws IOException {
        ByteBuffer buf = ByteBuffer.allocate(512);
        channel.read(buf, 0);
        String head = new String(buf.array(), 0, buf.position(), StandardCharsets.ISO_8859_1);
        int decl = head.indexOf("<?xml");
        if (decl >= 0) {
            int end = head.indexOf("?>", decl);
            String xmlDecl = end > 0 ? head.substring(decl, end) : head.substring(decl);
            int enc = xmlDecl.indexOf("encoding");
            if (enc >= 0) {
                int q = enc + 8;
                while (q < xmlDecl.length() && xmlDecl.charAt(q) != '"' && xmlDecl.charAt(q) != '\'')
                    q++;
                if (q < xmlDecl.length()) {
                    char quote = xmlDecl.charAt(q);
                    int q2 = xmlDecl.indexOf(quote, q + 1);
                    if (q2 > q) {
                        try {
                            return Charset.forName(xmlDecl.substring(q + 1, q2));
                        } catch (Exception e) {
                            // use default
                        }
                    }
                }
            }
        }
        return StandardCharsets.UTF_8;
    }

    private class SpectrumIndexBuilder {
        final ArrayList<Long> starts = new ArrayList<Long>();
        final ArrayList<Long> ends = new ArrayList<Long>();
        final ArrayList<String> idList = new ArrayList<String>();
        final ArrayList<Integer> indexList = new ArrayList<Integer>();

        private final byte[] buf = new byte[1 << 22];
        private int bufLen = 0;
        private int bufPos = 0;
        private long bufStart = 0;    // file offset of buf[0]
        private InputStream in;

        private int read() throws IOException {
            if (bufPos >= bufLen) {
                bufStart += bufLen;
                bufLen = in.read(buf, 0, buf.length);
                bufPos = 0;
                if (bufLen <= 0) {
                    bufLen = 0;
                    return -1;
                }
            }
            return buf[bufPos++] & 0xFF;
        }

        private long position() {
            return bufStart + bufPos;
        }

        void scan() throws IOException {
            in = new FileInputStream(specFile);
            try {
                boolean inSpectrum = false;
                ByteArrayOutputStream tag = new ByteArrayOutputStream(256);
                while (true) {
                    // find '<'
                    int c;
                    while (bufPos < bufLen && buf[bufPos] != '<')
                        bufPos++;
                    c = read();
                    if (c < 0)
                        break;
                    if (c != '<')
                        continue;
                    long tagStart = position() - 1;
                    // tag name
                    StringBuilder name = new StringBuilder(16);
                    while (true) {
                        c = read();
                        if (c < 0 || c == '>' || c == ' ' || c == '\t' || c == '\n' || c == '\r' || (c == '/' && name.length() > 0))
                            break;
                        name.append((char) c);
                        if (name.length() > 64)
                            break;
                    }
                    if (c < 0)
                        break;
                    String tagName = localName(name.toString());
                    if (tagName.equals("spectrum") && name.charAt(0) != '/') {
                        // read the rest of the start tag, respecting quotes
                        tag.reset();
                        int quote = 0;
                        boolean selfClosing = false;
                        int prev = c;
                        if (c != '>') {
                            while (true) {
                                c = read();
                                if (c < 0)
                                    break;
                                if (quote != 0) {
                                    if (c == quote)
                                        quote = 0;
                                } else if (c == '"' || c == '\'') {
                                    quote = c;
                                } else if (c == '>') {
                                    selfClosing = prev == '/';
                                    break;
                                }
                                tag.write(c);
                                prev = c;
                            }
                        }
                        if (c < 0)
                            break;
                        Map<String, String> attrs = parseAttributes(tag.toByteArray());
                        starts.add(tagStart);
                        idList.add(attrs.get("id"));
                        String indexStr = attrs.get("index");
                        int index = -1;
                        if (indexStr != null) {
                            try {
                                index = Integer.parseInt(indexStr.trim());
                            } catch (NumberFormatException e) {
                                index = -1;
                            }
                        }
                        indexList.add(index);
                        if (selfClosing)
                            ends.add(position());
                        else
                            inSpectrum = true;
                    } else if (name.length() > 0 && name.charAt(0) == '/' && localName(name.substring(1)).equals("spectrum") && inSpectrum) {
                        // end tag: skip to '>'
                        while (c >= 0 && c != '>')
                            c = read();
                        if (c < 0)
                            break;
                        ends.add(position());
                        inSpectrum = false;
                    }
                }
                if (inSpectrum)
                    throw new IOException("Unterminated spectrum element in " + specFile);
            } finally {
                in.close();
            }
        }

        private Map<String, String> parseAttributes(byte[] bytes) {
            String s = new String(bytes, charset);
            HashMap<String, String> attrs = new HashMap<String, String>();
            int i = 0;
            int n = s.length();
            while (i < n) {
                while (i < n && Character.isWhitespace(s.charAt(i)))
                    i++;
                int nameStart = i;
                while (i < n && s.charAt(i) != '=' && !Character.isWhitespace(s.charAt(i)))
                    i++;
                String attrName = s.substring(nameStart, i);
                while (i < n && Character.isWhitespace(s.charAt(i)))
                    i++;
                if (i >= n || s.charAt(i) != '=')
                    break;
                i++;
                while (i < n && Character.isWhitespace(s.charAt(i)))
                    i++;
                if (i >= n)
                    break;
                char quote = s.charAt(i);
                if (quote != '"' && quote != '\'')
                    break;
                int valueEnd = s.indexOf(quote, i + 1);
                if (valueEnd < 0)
                    break;
                attrs.put(localName(attrName), decodeEntities(s.substring(i + 1, valueEnd)));
                i = valueEnd + 1;
            }
            return attrs;
        }
    }

    // XML attribute value normalization: entities and character references
    // (whitespace characters are normalized to spaces as an XML parser would)
    static String decodeEntities(String v) {
        if (v.indexOf('&') < 0 && v.indexOf('\t') < 0 && v.indexOf('\n') < 0 && v.indexOf('\r') < 0)
            return v;
        StringBuilder sb = new StringBuilder(v.length());
        for (int i = 0; i < v.length(); i++) {
            char c = v.charAt(i);
            if (c == '\t' || c == '\n' || c == '\r') {
                if (c == '\r' && i + 1 < v.length() && v.charAt(i + 1) == '\n')
                    i++;
                sb.append(' ');
            } else if (c == '&') {
                int semi = v.indexOf(';', i);
                if (semi < 0) {
                    sb.append(c);
                    continue;
                }
                String ent = v.substring(i + 1, semi);
                if (ent.equals("amp")) sb.append('&');
                else if (ent.equals("lt")) sb.append('<');
                else if (ent.equals("gt")) sb.append('>');
                else if (ent.equals("quot")) sb.append('"');
                else if (ent.equals("apos")) sb.append('\'');
                else if (ent.startsWith("#x") || ent.startsWith("#X"))
                    sb.appendCodePoint(Integer.parseInt(ent.substring(2), 16));
                else if (ent.startsWith("#"))
                    sb.appendCodePoint(Integer.parseInt(ent.substring(1)));
                else {
                    sb.append(c);
                    continue;
                }
                i = semi;
            } else {
                sb.append(c);
            }
        }
        return sb.toString();
    }

    // ---------------------------------------------------------------------------------------------
    // Parsing one spectrum element into a small tree

    static class UnsupportedFeature extends Exception {
        private static final long serialVersionUID = 1L;

        UnsupportedFeature(String msg) {
            super(msg);
        }
    }

    static class Elem {
        final String name;
        final String[] attrNames;
        final String[] attrValues;
        final ArrayList<Elem> children = new ArrayList<Elem>(4);
        StringBuilder text = null;

        Elem(String name, String[] attrNames, String[] attrValues) {
            this.name = name;
            this.attrNames = attrNames;
            this.attrValues = attrValues;
        }

        String attr(String n) {
            for (int i = 0; i < attrNames.length; i++)
                if (attrNames[i].equals(n))
                    return attrValues[i];
            return null;
        }

        Elem child(String n) {
            for (Elem c : children)
                if (c.name.equals(n))
                    return c;
            return null;
        }

        List<Elem> children(String n) {
            ArrayList<Elem> list = new ArrayList<Elem>(children.size());
            for (Elem c : children)
                if (c.name.equals(n))
                    list.add(c);
            return list;
        }
    }

    private Elem parseElement(int pos) throws IOException, XMLStreamException, UnsupportedFeature {
        byte[] bytes = new byte[length[pos]];
        ByteBuffer bb = ByteBuffer.wrap(bytes);
        long filePos = startOffset[pos];
        while (bb.hasRemaining()) {
            int n = channel.read(bb, filePos + bb.position());
            if (n < 0)
                throw new EOFException();
        }
        XMLStreamReader r = xmlInputFactory.get().createXMLStreamReader(new ByteArrayInputStream(bytes), charset.name());
        try {
            ArrayDeque<Elem> stack = new ArrayDeque<Elem>();
            Elem root = null;
            while (r.hasNext()) {
                int event = r.next();
                if (event == XMLStreamConstants.START_ELEMENT) {
                    String name = localName(r.getLocalName());
                    if (name.equals("referenceableParamGroupRef"))
                        throw new UnsupportedFeature("referenceableParamGroupRef in a spectrum");
                    int na = r.getAttributeCount();
                    String[] an = new String[na];
                    String[] av = new String[na];
                    for (int i = 0; i < na; i++) {
                        an[i] = localName(r.getAttributeLocalName(i));
                        av[i] = r.getAttributeValue(i);
                    }
                    Elem e = new Elem(name, an, av);
                    if (name.equals("binary"))
                        e.text = new StringBuilder();
                    if (stack.isEmpty())
                        root = e;
                    else
                        stack.peek().children.add(e);
                    stack.push(e);
                } else if (event == XMLStreamConstants.END_ELEMENT) {
                    stack.pop();
                    if (stack.isEmpty())
                        break;
                } else if (event == XMLStreamConstants.CHARACTERS || event == XMLStreamConstants.CDATA) {
                    Elem top = stack.peek();
                    if (top != null && top.text != null)
                        top.text.append(r.getTextCharacters(), r.getTextStart(), r.getTextLength());
                }
            }
            if (root == null || !root.name.equals("spectrum"))
                throw new UnsupportedFeature("not a spectrum element");
            return root;
        } finally {
            r.close();
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Conversion: mirrors SpectrumConverter.getSpectrumFromJMzMLSpec step by step

    private static List<Elem> cvParams(Elem e) {
        return e.children("cvParam");
    }

    private static String acc(Elem cvParam) {
        return cvParam.attr("accession");
    }

    private static int requiredInt(Elem e, String attr) throws UnsupportedFeature {
        String v = e.attr(attr);
        if (v == null)
            throw new UnsupportedFeature("missing attribute " + attr + " on " + e.name);
        try {
            return Integer.parseInt(v.trim());
        } catch (NumberFormatException ex) {
            throw new UnsupportedFeature("bad attribute " + attr + " on " + e.name);
        }
    }

    private Spectrum convert(Elem jmzMLSpec) throws UnsupportedFeature {
        Spectrum spec = new Spectrum();

        // ID
        String id = jmzMLSpec.attr("id");
        if (id == null)
            throw new UnsupportedFeature("spectrum without id");
        spec.setID(id);

        // scan number
        String[] idToken = id.split("\\s+");
        if (idToken.length > 0 && idToken[idToken.length - 1].matches("scan=\\d+")) {
            int scanNum = Integer.parseInt(idToken[idToken.length - 1].substring(5));
            spec.setScanNum(scanNum);
        }

        // MS Level
        Elem msLevelParam = null;

        Boolean isCentroided = false;
        Spectrum.Polarity scanPolarity = Spectrum.Polarity.POSITIVE;

        for (Elem cvParam : cvParams(jmzMLSpec)) {
            String a = acc(cvParam);
            if (a == null)
                throw new UnsupportedFeature("cvParam without accession");
            if (a.equals("MS:1000511"))    // MS level
            {
                msLevelParam = cvParam;
            } else if (a.equals("MS:1000127"))    // centroid spectrum
            {
                isCentroided = true;
            } else if (a.equals("MS:1000128"))    // profile spectrum
            {
                isCentroided = false;
            } else if (a.equals("MS:1000129"))    // negative mode scan
            {
                scanPolarity = Spectrum.Polarity.NEGATIVE;
            } else if (a.equals("MS:1000130"))    // positive mode scan
            {
                scanPolarity = Spectrum.Polarity.POSITIVE;
            }
        }

        spec.setIsCentroided(isCentroided);
        spec.setScanPolarity(scanPolarity);

        float precursorMz = -1;
        float scanStartTime = -1;
        boolean scanStartTimeIsSeconds = true;

        // Scan list to get monoisotopic m/z
        Elem scanList = jmzMLSpec.child("scanList");
        if (scanList != null) {
            List<Elem> scans = scanList.children("scan");
            if (requiredInt(scanList, "count") > 0 && scans.size() > 0) {
                Elem scan = scans.get(0);
                for (Elem cvParam : cvParams(scan)) {
                    String a = acc(cvParam);
                    if (a == null)
                        throw new UnsupportedFeature("cvParam without accession");
                    if (a.equals("MS:1000016")) {
                        String value = cvParam.attr("value");
                        String unitAccession = cvParam.attr("unitAccession");
                        if (value == null || unitAccession == null)
                            throw new UnsupportedFeature("scan start time without value or unit");
                        scanStartTime = Float.parseFloat(value);
                        if (unitAccession.equals("UO:0000031")) {
                            // in minutes
                            scanStartTimeIsSeconds = false;
                        } else if (unitAccession.equals("UO:0000010")) {
                            // in seconds
                            scanStartTimeIsSeconds = true;
                        }
                    }
                    // is_a: MS:1002892 ! ion mobility attribute cvParams
                    else if (a.equals("MS:1001581")      // FAIMS compensation voltage
                            || a.equals("MS:1002476")    // ion mobility drift time
                            || a.equals("MS:1002815")) { // inverse reduced ion mobility
                        CvParamInfo cvParamInfo;
                        String unitAccession = cvParam.attr("unitAccession");
                        if (unitAccession != null && !unitAccession.isEmpty()) {
                            cvParamInfo = new CvParamInfo(a, cvParam.attr("name"), cvParam.attr("value"), unitAccession, cvParam.attr("unitName"));
                        } else {
                            cvParamInfo = new CvParamInfo(a, cvParam.attr("name"), cvParam.attr("value"));
                        }

                        spec.addAddlCvParam(cvParamInfo);
                    }
                }
                for (Elem param : scan.children("userParam")) {
                    String name = param.attr("name");
                    if (name == null)
                        throw new UnsupportedFeature("userParam without name");
                    if (name.equals("[Thermo Trailer Extra]Monoisotopic M/Z:")) {
                        String value = param.attr("value");
                        if (value == null)
                            throw new UnsupportedFeature("userParam without value");
                        precursorMz = Float.parseFloat(value);
                        break;
                    }
                }
            }
        }
        spec.setRt(scanStartTime);
        spec.setRtIsSeconds(scanStartTimeIsSeconds);

        int msLevel;
        if (msLevelParam != null) {
            String value = msLevelParam.attr("value");
            if (value == null)
                throw new UnsupportedFeature("ms level without value");
            msLevel = Integer.parseInt(value);
        } else {
            msLevel = 0;
        }
        spec.setMsLevel(msLevel);

        // Precursor
        Elem precursorList = jmzMLSpec.child("precursorList");
        if (precursorList != null) {
            List<Elem> precursors = precursorList.children("precursor");
            if (requiredInt(precursorList, "count") > 0) {
                if (precursors.isEmpty())
                    throw new UnsupportedFeature("empty precursorList");
            }
            if (requiredInt(precursorList, "count") > 0 && precursors.get(0).child("selectedIonList") != null) {
                Elem precursor = precursors.get(0);    // consider only the first precursor

                Elem isolationWindowParams = precursor.child("isolationWindow");
                if (isolationWindowParams != null) {
                    Float isolationWindowTargetMz = null;
                    for (Elem param : cvParams(isolationWindowParams)) {
                        String a = acc(param);
                        if (a == null)
                            throw new UnsupportedFeature("cvParam without accession");
                        if (a.equals("MS:1000827"))    // selected ion m/z
                        {
                            isolationWindowTargetMz = Float.parseFloat(param.attr("value"));    // assume that unit is m/z (MS:1000040)
                        }
                    }
                    spec.setIsolationWindowTargetMz(isolationWindowTargetMz);
                }

                // precursor mz, charge
                int precursorCharge = 0;
                float precursorIntensity = 0;

                List<Elem> selectedIons = precursor.child("selectedIonList").children("selectedIon");
                if (selectedIons.isEmpty())
                    throw new UnsupportedFeature("empty selectedIonList");
                Elem paramGroup = selectedIons.get(0);
                for (Elem param : cvParams(paramGroup)) {
                    String a = acc(param);
                    if (a == null)
                        throw new UnsupportedFeature("cvParam without accession");
                    if (precursorMz < 0.01 && a.equals("MS:1000744"))    // selected ion m/z
                    {
                        precursorMz = Float.parseFloat(param.attr("value"));    // assume that unit is m/z (MS:1000040)
                    } else if (a.equals("MS:1000041"))    // charge state
                    {
                        precursorCharge = Integer.parseInt(param.attr("value"));
                    } else if (a.equals("MS:1000042"))    // peak intensity
                    {
                        precursorIntensity = Float.parseFloat(param.attr("value"));
                    }
                }

                spec.setPrecursor(new Peak(precursorMz, precursorIntensity, precursorCharge));

                // activation method
                Elem actMethodParams = precursor.child("activation");
                if (actMethodParams == null)
                    throw new UnsupportedFeature("precursor without activation");
                boolean isETD = false;
                for (Elem param : cvParams(actMethodParams)) {
                    String a = acc(param);
                    if (a == null)
                        throw new UnsupportedFeature("cvParam without accession");
                    ActivationMethod am = ActivationMethod.getByCV(a);
                    if (am != null) {
                        if (am == ActivationMethod.ETD) {
                            isETD = true;
                            break;
                        }
                        if (spec.getActivationMethod() == null)
                            spec.setActivationMethod(am);
                    }
                }
                if (isETD)
                    spec.setActivationMethod(ActivationMethod.ETD);
            }
        }

        // Peak list
        Elem mzArray = null, intenArray = null;

        Elem binaryDataArrayList = jmzMLSpec.child("binaryDataArrayList");
        if (binaryDataArrayList != null) {
            for (Elem array : binaryDataArrayList.children("binaryDataArray")) {
                if (requiredInt(array, "encodedLength") == 0)
                    continue;
                // check the cvParams
                for (Elem param : cvParams(array)) {
                    String a = acc(param);
                    if (a == null)
                        throw new UnsupportedFeature("cvParam without accession");
                    if (a.equals("MS:1000514")) {
                        mzArray = array;
                        break;
                    }
                    if (a.equals("MS:1000515")) {
                        intenArray = array;
                        break;
                    }
                }
                if (mzArray != null && intenArray != null)
                    break;
            }
        }

        if (mzArray != null && intenArray != null) {
            float[] mzNumbers = decodeBinaryDataArray(mzArray);
            float[] intenNumbers = decodeBinaryDataArray(intenArray);

            if (mzNumbers.length != intenNumbers.length)
                throw new UnsupportedFeature("different sizes for m/z and intensity arrays");    // jmzML path reports this

            for (int i = 0; i < mzNumbers.length; i++)
                spec.add(new Peak(mzNumbers[i], intenNumbers[i], 1));
        }

        // SpecIndex
        spec.setSpecIndex(requiredInt(jmzMLSpec, "index") + 1);    // 1-based spectrum index

        // sort peaks by increasing order of m/z
        Collections.sort(spec);

        spec.determineIsCentroided();

        return spec;
    }

    private Float precursorMzOf(Elem jmzMLSpec) throws UnsupportedFeature {
        // mirrors SpectrumConverter.getPrecursorMzFromJMzMLSpec
        Elem precursorList = jmzMLSpec.child("precursorList");
        if (precursorList != null && requiredInt(precursorList, "count") > 0) {
            List<Elem> precursors = precursorList.children("precursor");
            if (precursors.isEmpty())
                throw new UnsupportedFeature("empty precursorList");
            Elem selectedIonList = precursors.get(0).child("selectedIonList");
            if (selectedIonList != null) {
                List<Elem> selectedIons = selectedIonList.children("selectedIon");
                if (selectedIons.isEmpty())
                    throw new UnsupportedFeature("empty selectedIonList");
                for (Elem param : cvParams(selectedIons.get(0))) {
                    String a = acc(param);
                    if (a == null)
                        throw new UnsupportedFeature("cvParam without accession");
                    if (a.equals("MS:1000744"))
                        return Float.parseFloat(param.attr("value"));
                }
            }
        }
        return null;
    }

    // Compression and encoding accessions other than zlib / none: let jmzML handle them
    private static final Set<String> UNSUPPORTED_ARRAY_ACCESSIONS = new HashSet<String>(Arrays.asList(
            "MS:1002312", "MS:1002313", "MS:1002314",    // MS-Numpress
            "MS:1002746", "MS:1002747", "MS:1002748",    // MS-Numpress followed by zlib
            "MS:1003088", "MS:1003089", "MS:1003090",    // truncation / zstd variants
            "MS:1003091", "MS:1003092", "MS:1003093",
            "MS:1000520"                                 // 16-bit float
    ));

    private static float[] decodeBinaryDataArray(Elem array) throws UnsupportedFeature {
        boolean zlib = false;
        int precision = 0;    // 1: 32-bit float, 2: 64-bit float, 3: 32-bit int, 4: 64-bit int
        for (Elem param : cvParams(array)) {
            String a = acc(param);
            if (a == null)
                throw new UnsupportedFeature("cvParam without accession");
            if (UNSUPPORTED_ARRAY_ACCESSIONS.contains(a))
                throw new UnsupportedFeature("binary array encoding " + a);
            switch (a) {
                case "MS:1000574":
                    zlib = true;
                    break;
                case "MS:1000576":
                    break;
                case "MS:1000521":
                    precision = setPrecision(precision, 1);
                    break;
                case "MS:1000523":
                    precision = setPrecision(precision, 2);
                    break;
                case "MS:1000519":
                    precision = setPrecision(precision, 3);
                    break;
                case "MS:1000522":
                    precision = setPrecision(precision, 4);
                    break;
                default:
                    break;
            }
        }
        if (precision == 0)
            throw new UnsupportedFeature("binary array without precision");

        Elem binary = array.child("binary");
        if (binary == null || binary.text == null)
            throw new UnsupportedFeature("binary array without data");
        byte[] raw;
        try {
            raw = Base64.getMimeDecoder().decode(binary.text.toString());
        } catch (IllegalArgumentException e) {
            throw new UnsupportedFeature("invalid base64");
        }
        if (zlib)
            raw = inflate(raw);

        int bytesPerValue = (precision == 1 || precision == 3) ? 4 : 8;
        if (raw.length % bytesPerValue != 0)
            throw new UnsupportedFeature("binary array length is not a multiple of the value size");
        int n = raw.length / bytesPerValue;
        ByteBuffer bb = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
        float[] values = new float[n];
        switch (precision) {
            case 1:
                for (int i = 0; i < n; i++)
                    values[i] = bb.getFloat(i * 4);
                break;
            case 2:
                for (int i = 0; i < n; i++)
                    values[i] = (float) bb.getDouble(i * 8);
                break;
            case 3:
                for (int i = 0; i < n; i++)
                    values[i] = (float) bb.getInt(i * 4);
                break;
            default:
                for (int i = 0; i < n; i++)
                    values[i] = (float) bb.getLong(i * 8);
                break;
        }
        return values;
    }

    private static int setPrecision(int current, int value) throws UnsupportedFeature {
        if (current != 0 && current != value)
            throw new UnsupportedFeature("conflicting binary precision");
        return value;
    }

    private static byte[] inflate(byte[] compressed) throws UnsupportedFeature {
        Inflater inflater = new Inflater();
        try {
            inflater.setInput(compressed);
            ByteArrayOutputStream out = new ByteArrayOutputStream(Math.max(64, compressed.length * 3));
            byte[] buf = new byte[1 << 16];
            while (!inflater.finished()) {
                int n = inflater.inflate(buf);
                if (n == 0) {
                    if (inflater.needsInput() || inflater.needsDictionary())
                        break;
                }
                out.write(buf, 0, n);
            }
            return out.toByteArray();
        } catch (DataFormatException e) {
            throw new UnsupportedFeature("invalid zlib data");
        } finally {
            inflater.end();
        }
    }

    // ---------------------------------------------------------------------------------------------

    private class FastIterator implements Iterator<Spectrum> {
        // Spectra are parsed ahead in parallel batches; results are returned in file order
        private static final int BATCH = 512;
        private final int numThreads = Math.max(1, Math.min(8, Runtime.getRuntime().availableProcessors()));
        private final java.util.concurrent.ExecutorService pool;
        private int nextPos = 0;
        private ArrayDeque<java.util.concurrent.Future<Spectrum[]>> pending = new ArrayDeque<java.util.concurrent.Future<Spectrum[]>>();
        private Spectrum[] batch = null;
        private int batchIndex = 0;
        private Spectrum current;

        FastIterator() {
            pool = java.util.concurrent.Executors.newFixedThreadPool(numThreads, r -> {
                Thread t = new Thread(r, "mzML-reader");
                t.setDaemon(true);
                return t;
            });
            fill();
            current = advance();
        }

        private void fill() {
            while (pending.size() < numThreads * 2 && nextPos < startOffset.length) {
                final int from = nextPos;
                final int to = Math.min(startOffset.length, from + BATCH);
                nextPos = to;
                pending.add(pool.submit(() -> {
                    Spectrum[] out = new Spectrum[to - from];
                    for (int p = from; p < to; p++)
                        out[p - from] = readSpectrumAt(p);
                    return out;
                }));
            }
        }

        private Spectrum advance() {
            while (true) {
                if (batch != null && batchIndex < batch.length) {
                    Spectrum spec = batch[batchIndex++];
                    if (spec == null)
                        continue;
                    if (spec.getMSLevel() < minMSLevel || spec.getMSLevel() > maxMSLevel)
                        continue;
                    return spec;
                }
                if (pending.isEmpty()) {
                    pool.shutdown();
                    return null;
                }
                try {
                    batch = pending.poll().get();
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    pool.shutdownNow();
                    return null;
                } catch (java.util.concurrent.ExecutionException e) {
                    pool.shutdownNow();
                    throw new RuntimeException(e.getCause());
                }
                batchIndex = 0;
                fill();
            }
        }

        public boolean hasNext() {
            return current != null;
        }

        private long negativePolarityWarningCount = 0;

        public Spectrum next() {
            if (current == null)
                throw new NoSuchElementException();
            Spectrum spec = current;
            current = advance();
            if (spec.getScanPolarity() == Spectrum.Polarity.NEGATIVE)
                warnNegativePolarity(spec);
            return spec;
        }

        // same messages as MzMLSpectraIterator
        private void warnNegativePolarity(Spectrum currentSpectrum) {
            negativePolarityWarningCount++;
            if (negativePolarityWarningCount > edu.ucsd.msjava.parser.SpectrumParser.MAX_NEGATIVE_POLARITY_WARNINGS)
                return;
            if (negativePolarityWarningCount == 1) {
                System.out.println("Warning: negative polarity spectrum found; you likely need to use a negative charge carrier");
            }
            System.out.println("Negative polarity spectrum found, scan " + Long.toString(currentSpectrum.getScanNum()));
            if (negativePolarityWarningCount == edu.ucsd.msjava.parser.SpectrumParser.MAX_NEGATIVE_POLARITY_WARNINGS) {
                System.out.println("Additional warnings regarding negative polarity will not be shown");
            }
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
}
