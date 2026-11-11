# Why MS-GF+ is slow on the E4 baseline search

Measured 2026-11-11 on `harrier`, MS-GF+ v2024.03.26, OpenJDK 21.0.13, 16 threads.

## Finding

The search threads spend most of their time waiting for one lock, not computing.
In 77 thread dumps (1,103 worker-thread samples) 73% of samples were `BLOCKED`, and 80% had
`java.util.Hashtable.get` as the innermost frame, called from
`NewRankScorer.getIonExistenceScore`. `Hashtable.get` is `synchronized`, and the table is shared by
all worker threads, so every edge score takes a global lock.

- Lock site: `src/main/java/edu/ucsd/msjava/msscorer/NewRankScorer.java:135`
  (`this.ionExistenceTable.get(part)`; the field is declared at line 45 and built at line 381).
- Caller: `src/main/java/edu/ucsd/msjava/msscorer/DBScanScorer.java:71` (`getEdgeScoreInt`), reached
  from both `DBScanner.dbSearch` and `DBScanner.computeSpecEValue`.
- 762 of 801 "waiting to lock" entries name the same `Hashtable` object.
- `/usr/bin/time` on the benchmark run agrees: 504% CPU on 16 threads and 64 million voluntary
  context switches in 22 minutes.

The other tables in the same class (`rankDistTable`, `ionErrDistTable`, `noiseErrDistTable`, lines
30–44) are also `Hashtable`s read on the scoring path; `getMissingIonScore` (line 112) shows up in
4% of samples. They are read-only after the parameter file is loaded.

Fixed: see [Changes and results](#changes-and-results). Removing the lock cut the search stage from
1,171 s to about 230 s; the other changes there bring the whole run to 171 s.

## Where the time goes

Wall time by stage, baseline cell, HCD run (102,726 spectra, human Swiss-Prot + reversed decoys):

| Stage | Benchmark run, rep 0 (s) | Diagnostic run (s) |
|---|---|---|
| Load / build database index | 17 | 14 |
| Read mzML (single thread) | 106 | 112 |
| Search (48 tasks on 16 threads) | 1,171 | about 1,380 |
| Write mzIdentML | 41 | 50 |
| Total | 1,318 | 1,555 |

The diagnostic run is slower because it prints stage messages and takes a thread dump every 20 s;
the host is shared, so wall times vary between runs.

Inside the search stage, summed over the 48 tasks of the diagnostic run (`-verbose 1` messages):

| Sub-stage | Thread-seconds | Share |
|---|---|---|
| Preprocessing spectra | 1,099 | 5% |
| Database search | 3,860 | 19% |
| Computing spectral E-values | 15,617 | 76% |

Thread samples by stage give the same picture: `DBScanner.computeSpecEValue` 71%,
`DBScanner.dbSearch` 18%, `ScoredSpectraMap.preProcessIndividualSpectra` 5%. These thread-seconds
include time spent blocked on the lock, so they are not CPU cost.

For scale, the same search takes 23 s in Sage, 25 s in MSFragger and 52 s in Comet on the same
machine and thread count.

## Changes and results

Measured 2026-11-11 on `harrier` with the patched source, same commands and 16 pinned CPUs. The host
was heavily loaded (load average 75–100 on 96 cores), so wall times are approximate; user CPU time
from `/usr/bin/time -v` is given as a second measure.

| Search | Benchmark rep 0, wall (s) | Patched, wall (s) | Benchmark user CPU (s) | Patched user CPU (s) |
|---|---|---|---|---|
| HCD, 10 ppm (`baseline`) | 1,318 | 171 | 6,355 | 1,402 |
| Velos, 10 ppm (`frag_1.0005`) | 511 | 106 | 2,646 | 958 |
| Velos, 3 Da (`x_3da_frag1`) | 668 | 356 | 4,788 | 2,636 |

Stage times for the patched HCD run: load database 5 s, read mzML 15 s, search about 125 s,
write mzIdentML 27 s.

Outputs are unchanged. For all three searches the patched `.mzid` has the same lines as the
benchmark `rep0/out.mzid` once dates, versions and file paths are removed and the lines are sorted
(99,977, 61,972 and 63,449 results); only the order of results with tied scores differs. The new
mzML reader was also compared spectrum by spectrum with jmzML on all four files in `data/runs`
(511,806 spectra): every field and peak matched and none needed the jmzML fallback.

What changed, in order of effect (paths under `src/main/java/edu/ucsd/msjava/`):

1. Lock contention (`msscorer/NewRankScorer.java`, `msscorer/DBScanScorer.java`): the parameter
   tables are `HashMap`s; `DBScanScorer` computes the four ion-existence scores once per spectrum and
   reads error scores from a per-partition array cached in a `ConcurrentHashMap`.
2. mzML reading (`mzml/FastMzMLReader.java`, `mzml/FastMzMLSpectraMap.java`,
   `msutil/SpectraAccessor.java`): replaces jmzML's index (79–288 s per file) with a byte scan
   (4–9 s) and parses spectra with StAX, converting them exactly as `SpectrumConverter` does.
   Spectra with referenceable parameter groups, numpress or other unexpected content are read with
   jmzML instead. `-Dmsgfplus.fastMzML=false` turns the reader off. `mzid/MZIdentMLGen.java` reads
   the spectra it needs for output in parallel.
3. Database scan (`msdbsearch/DBScanner.java` and related): the suffix array is memory-mapped once
   instead of streamed by each of the 48 tasks; the precursor-mass window is a binary search over a
   sorted array instead of `TreeMap.subMap(...).size()`; the per-suffix loop skips peptide lengths
   with no matches; `StringBuffer`, `TreeSet<Integer>` and boxed-key `HashMap` lookups on the hot
   path are replaced.
4. Spectral E-values (`msgf/GeneratingFunction.java`, `msgf/FlexAminoAcidGraph.java`,
   `msgf/DeNovoGraph.java`): graph construction and the generating-function tables use arrays
   indexed by nominal mass instead of `HashMap<NominalMass, …>`.

What is left: in the 3 Da search 64% of CPU is the spectral E-value computation, mostly the
score-distribution convolution in `ScoreDist.addProbDist`, which runs once per candidate peptide
mass. In the 10 ppm searches the database scan and E-values are about 45% and 35% of CPU. Further
gains need algorithmic changes (for example sharing graph construction between neighbouring
masses).

`TestMSGFPlus.TestMzML` fails before and after the change because it reads from
`C:\DMS_WorkDir1`; the other unit tests pass (`mvn -DskipTests=false test`).

## Commands

All paths are under `/home/ygao/work/mass_complexity` unless absolute.

Benchmark run (E4 harness; three repeats under `/usr/bin/time -v`):

```
.venv/bin/python analysis/scripts/e4_bench.py --engine msgfplus --cell baseline \
  --mzml data/runs/20161213_NGHF_DBJ_SA_Exp3A_HeLa_1ug_60min_7500_02.mzML --dataset PXD006932 \
  --fasta data/fasta/human_sp_2026_04.fasta --threads 16 --repeats 3 --cpus 1,3,...,31
```

The MS-GF+ command it issues (from the log):

```
java -Xmx200G -jar tools/msgfplus/MSGFPlus.jar -s <run>.mzML -d human_sp_2026_04.fasta -o out.mzid \
  -mod params.mods -t 10ppm -ti 0,0 -thread 16 -tda 1 -m 0 -inst 1 -e 1 -ntt 2 \
  -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 4 -n 1 -maxMissedCleavages 2 -addFeatures 0
```

Diagnostic run (same command plus `-verbose 1`, pinned to CPUs 0,2,…,30, `kill -3` every 20 s so the
JVM writes thread dumps into the log; there is no `jstack` on this host):

```
data/work/msgf_diag/run.sh
```

Summaries:

```
.venv/bin/python analysis/scripts/msgf_threaddump_summary.py data/work/msgf_diag/log.txt
grep -E 'elapsed time' data/work/msgf_diag/log.txt
```

## Files

| What | Path |
|---|---|
| Diagnostic log with 77 thread dumps and stage messages | `data/work/msgf_diag/log.txt` |
| Diagnostic run script | `data/work/msgf_diag/run.sh` |
| Thread-dump summary script | `analysis/scripts/msgf_threaddump_summary.py` |
| Benchmark logs (stage times, `/usr/bin/time -v`) | `data/work/e4/msgfplus/20161213_NGHF_DBJ_SA_Exp3A_HeLa_1ug_60min_7500_02/baseline/rep{0,1,2}/log.txt` |
| Ion-trap run logs (511 s and 623 s) | `data/work/e4/msgfplus/20121012_Velos1_PJ_Arne_Hela_Ibaq/{frag_1.0005,x_3da_frag1}/rep*/log.txt` |
| Modification file | `analysis/params/msgfplus/baseline.mods` |
| Benchmark harness | `analysis/scripts/e4_bench.py` |
| MS-GF+ jar | `tools/msgfplus/MSGFPlus.jar` |
| MS-GF+ source (commit `f0bc79e`) | `/home/ygao/work/msgfplus/src/main/java/edu/ucsd/msjava/` |
