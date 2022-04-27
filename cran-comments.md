rspca v0.1.0
=======
***

## Test environments
* windows-x86_64-devel (r-devel)
* macos-highsierra-release-cran (r-release)
* fedora-clang-devel (r-devel)

## R CMD check results
* There is no ERRORs or WARNINGs.
* There is one NOTE found in both Windows (Server 2022, R-devel, 64 bit) and Fedora (Linux, R-devel, clang, gfortran) regarding new submission and the spelling of the word *transcriptome*: 

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'X. Steve Niu <xin2001@med.cornell.edu>'
New submission

Possibly misspelled words in DESCRIPTION:
  omics (10:70)
  transcriptome (10:118)
```
This is a new release and the word *transcriptome* is spelt correctly as in [Transcriptome_Wikipedia](https://en.wikipedia.org/wiki/Transcriptome), the same is *omics* as in [Omics_Wikipedia](https://en.wikipedia.org/wiki/omics).

* There is another NOTE that was found only on Windows (Server 2022, R-devel, 64 bit):

```
> checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
```
As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this was due to a bug/crash in MiKTeX and can be ignored.

***

Thanks!

Steve X. Niu
