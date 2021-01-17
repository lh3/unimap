## Getting Started
```sh
# compile and install
git clone https://github.com/lh3/unimap
cd unimap && make

# simple use cases (-b24 below saves memory for this toy example)
./unimap -b24 -c test/MT-human.fa test/MT-orang.fa > out.paf

# prebuild index; recommended as unimap indexing is slower than minimap2 indexing
./unimap -b24 -d MT-human.umi test/MT-human.fa
./unimap -c MT-human.umi test/MT-orang.fa > out.paf

# use presets (no test data)
unimap -cxasm5 --cs -t8 ref.fa asm.fa         # if asm.fa is near identical to ref.fa
unimap -cxhifi --cs -t8 ref.fa hifi-reads.fa  # HiFi reads to reference alignment
unimap -cxont  --cs -t8 ref.fa ont-reads.fa   # Nanopore reads to reference alignment
```

## Introduction

Unimap is a fork of [minimap2][mm2] optimized for assembly-to-reference
alignment. It integrates the [minigraph][mg] chaining algorithm and can align
through long INDELs (up to 100kb by default) much faster than minimap2. Unimap
is a better fit for resolving segmental duplications and is recommended over
minimap2 for alignment between high-quality assemblies.

Unimap does not replace minimap2 for other types of alignment. It drops the
support of multi-part index and short-read mapping. Its long-read alignment is
different from minimap2 but is not necessarily better. Unimap is more of a
specialized minimap2 at the moment.

## Notes

* With the default `asm5` preset, unimap may align a highly diverged region
  as a long insertions followed by a long deletion. [Truvari][tru] may identify
  two false positive calls in this case, but these arguably are not errors.

* Unimap takes ~5 minutes to index a human genome, slower than minimap2. It is
  recommended to save the index for faster startup.

* The default `ont` preset has been tuned for more recent Nanopore reads at
  ~95% accuracy.

[mm2]: https://github.com/lh3/minimap2
[mg]: https://github.com/lh3/minigraph
[tru]: https://github.com/spiralgenetics/truvari
