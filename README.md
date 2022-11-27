# Quintet Rooting via Unrooting

Quintet Rooting via Unrooting (QRU) is an attempt to root species trees
using a method based on, but different from [Quintet Rooting](https://github.com/ytabatabaee/Quintet-Rooting) (QR).
We also use a lot of code from QR, so that repository is included as a submodule here.
Additionally, we use the "ASTRAL-modified" code from the [FASTRAL](https://github.com/PayamDiba/FASTRAL) repository,
which is included as a subdirectory here.
To get QRU into your system properly, please use
```
git clone --recurse-submodules https://github.com/bakhil/qru
```
Then run the `make.sh` in `ASTRAL-modified` to create an ASTRAL executable.
You'll need to have a java compiler installed for running `make.sh`.
You also need `dendropy`, `table-five`, and `numpy` for running QRU (these can be installed using `pip`).

You can also find a few scripts on the `scripts` branch of this repository.
You might find them useful for iterating over the [QR-paper](https://github.com/ytabatabaee/QR-paper) datasets,
computing error metrics over the datasets, and plotting.
You'll also find a PDF report describing this work in more detail on that branch.
(This was done in Fall 2022 as a course project for [CS581 at UIUC](http://tandy.cs.illinois.edu/CS581-Fa2022-short.html).)

QRU uses the same input-output format as QR.
There just one additional optional argument for specifying a path
to an ASTRAL executable (ignore if you just want to use QRU as is).
Please see instructions from the [QR repository](https://github.com/ytabatabaee/Quintet-Rooting)
for details on the other arguments and an example (same should work with QRU as well,
just copy the `example` directory to the correct location).
Also, the `-cfs` option of QR does not work for QRU.

### Rooting an unrooted species tree
**Input:** A file containing an unrooted species tree (with at least 5 taxa) and a file containing a set of unrooted gene trees, both in newick format (with or without branch lengths).

**Output:** A file containing the rooted species tree in newick format 
```
$ python3 quintet_rooting.py -t <species-topology.tre> -g <input-genes.tre> -o <output-tree.tre>
```
**Arguments**
- **Required**
```
 -t,  --speciestree        input unrooted species tree in newick format
 -g,  --genetrees          input gene trees in newick format
 -o,  --output             output file containing a rooted species tree
```
- **Optional**
```
 -h,  --help               show this help message and exit
 -sm, --samplingmode       TC for triplet cover, LE for linear encoding, EXH for exhaustive
 -c,  --cost               cost function (STAR for QR-STAR)
 -mult, --multiplicity     multiplicity (number of quintets mapped to each edge) in QR-LE
 -norm, --normalized       using normalization for unresolved gene trees or missing taxa
 -coef, --coef             shape coefficient in QR-STAR
 -abratio, --abratio       ratio of invariants to inequalities in QR-STAR
 -rs,  --seed              random seed
 -apath, --astralpath      path to the ASTRAL executable, relative to qru.py
```

