## Patch Note

### v0.3.1 (2023-09-07)

#### New Feature

- Add multithreading support

#### Changes

- More conserved memory estimation

#### Bug fix

- Now rapid job queuing will no longer cause thread pool to exit prematurely

### v0.2.0 (2023-09-04)

#### New Feature

- Restriction of memory usage (-M) by limiting the size of read chunk.
- More general parsing ability that accounts for sciRNAseq3 _et al._.
  - Now read names and read tags are both supported to be used as CBC or UMI.

#### Bug fix

- Not providing output folder will no longer crash the app with segmentaion fault.
- MAPQ and primary/secondary mapping information are considered in deduplication of reads.
- Fix memory leaks.
