# D0jetTreeEmbedding

The original repository for this jet framework is based on https://github.com/joelmazer/star-jetframework, then modified by Diptanil Roy, then taken over by me.

Macro for analisis is `RunHeavyIonOverlayAngularities.C`, to run it for testing before submitting, one needs to change the parameters `doTEST` in [RunHeavyIonOverlayAngularities.C#L34](./RunHeavyIonOverlayAngularities.C#L34) and execute it with

```bash
root -l RunHeavyIonOverlayAngularities.C
```

The main source file for StHIOverlayAngularities [RunHeavyIonOverlayAngularities.C#L281](RunHeavyIonOverlayAngularities.C#L281) is in [StRoot/StMyAnalysisMaker/StHIOverlayAngularities.cxx](./StRoot/StMyAnalysisMaker/StHIOverlayAngularities.cxx)

Some important links in code :

* First, find original PYTHIA jets (MC jets): [line 621](./StRoot/StMyAnalysisMaker/StHIOverlayAngularities.cxx#L621)

* Then, find reconstructed PYTHIA jets (MC Reco jets): [line 788](./StRoot/StMyAnalysisMaker/StHIOverlayAngularities.cxx#L621)

* The last part is to embed those MC Reco jets into Heavy-Ion environments and try to find $D^0$ jets there (Reco jets): [line 861](./StRoot/StMyAnalysisMaker/StHIOverlayAngularities.cxx#L861)

* Background subtraction is done here : [line 1135](./StRoot/StMyAnalysisMaker/StHIOverlayAngularities.cxx#L1135)

* The important function for finding the exact $D^0$ jet is `DoesItHaveAGoodD0Jet()` [line 1748](./StRoot/StMyAnalysisMaker/StHIOverlayAngularities.cxx#L1748)

---
To run the whole chain, just use

```bash
star-submit send2014.xml
```

or try automatization with merging using [run.csh](./run.csh) with atomatic merging and waiting for the jobs using [wait_for_jobs.sh](./wait_for_jobs.sh)

---

### FASTJET install 

This framework requires the installation of FastJet + FastJet contrib.

One can use my installations by creating links:

```bash
ln -s /gpfs01/star/pwg/prozorov/install/fastjet-install/include/siscone/ siscone
ln -s /gpfs01/star/pwg/prozorov/install/fastjet-install/include/fastjet/ fastjet
```
or install it manually. FastJet can be found at: http://fastjet.fr/all-releases.html -  download FastJet 3.0 +

FastJet contrib can be found at: http://fastjet.hepforge.org/contrib/

RCAS & PDSF both require 32-bit build to fit into other software framework.  To force this install on the 64-bit systems, please follow the below prescription to avoid headaches.

FASTJET installation example:

``` bash
  tar zxvf fastjet-3.4.2.tar.gz
  cd fastjet-3.4.2/
  ./configure --prefix=$PWD/../fastjet-install CXXFLAGS="-m32 -fPIC -fno-inline" CFLAGS="-m32 -fPIC -fno-inline" LDFLAGS="-m32"

  make
  make check
  make install
```

FastJet contrib installation example:
``` bash
  tar zxvf fjcontrib-1.046.tar.gz
  cd fjcontrib-1.046
  ./configure --fastjet-config=$PWD/../fastjet-3.4.2/fastjet-config --prefix=$PWD/../fastjet-install CXXFLAGS="-m32 -fPIC -fno-inline" CFLAGS="-m32 -fPIC -fno-inline" LDFLAGS="-m32"

  make
  make install
  make fragile-shared
  make fragile-shared-install
```

unset flags after:
``` bash
  unset CXXFLAGS CFLAGS LDFLAGS
```

You will additionally need to set FastJet environmental variables in your startup file (mine is .bashrc, and the below is an example):

``` bash
export FASTJET='/path/to/your/FastJet/fastjet-install'
export FASTJET_CONTRIB='/path/to/your/FastJet/fjcontrib-1.046'
```

and soft link the include paths to your working directory which contains your StRoot directory and where you run your readPicoDst.C file from:

``` bash
ln -s /path/to/your/FastJet/include/siscone siscone
ln -s /path/to/your/FastJet/include/fastjet fastjet
```