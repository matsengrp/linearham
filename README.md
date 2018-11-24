# linearham

```
     __      _.._
  .-'__`-._.'.--.'.__.,
 /--'  '-._.'    '-._./
/__.--._.--._.'``-.__/
'._.-'-._.-._.-''-..'
```

[![wercker status](https://app.wercker.com/status/284280f33f13e936de0d544a332121af/m/master "wercker status")](https://app.wercker.com/project/byKey/284280f33f13e936de0d544a332121af)

Developer documentation: http://matsengrp.github.io/linearham (or build locally by running `doxygen Doxyfile`)

## dependencies

On Debian/Ubuntu, execute all the `RUN` instructions in the `Dockerfile`.

## usage

Running `scons --run-partis --fasta-path=data/liao_dataset.fasta --all-clonal-seqs` runs `partis` assuming all the sequences in the Liao (2013) dataset are clonal.
Remove that flag if you want to run `partis` in `partition` mode.

Running `scons --run-linearham --template-path=templates/revbayes_template.rev --mcmc-iter=25 --mcmc-thin=1 --tune-iter=0` runs the `linearham` phylogenetic inference pipeline.

Running `_build/test/test` runs the `linearham` tests.

Running `_build/linearham/linearham` runs the `linearham` binary built from the source files in `src/`.

For more information on the `scons` arguments, run `scons --help`.
