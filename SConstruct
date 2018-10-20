#!/usr/bin/env python
# -*- coding: utf-8 -*-
import nestly
from nestly import scons as nestly_scons
import os
from SCons.Script import Environment
import SCons.Script as Script


#### Set up command line arguments/options

# partis arguments

Script.AddOption("--build-partis",
        dest="build_partis",
        action="store_true",
        default=False,
        help="Should we build partis?")

Script.AddOption("--run-partis",
        dest="run_partis",
        action="store_true",
        default=False,
        help="Should we run partis?")

Script.AddOption("--fasta-path",
        dest="fasta_path",
        type="str",
        default=None,
        help="The repertoire FASTA file path.")

# linearham arguments

Script.AddOption("--build-linearham",
        dest="build_linearham",
        action="store_true",
        default=False,
        help="Should we build linearham?")

Script.AddOption("--run-linearham",
        dest="run_linearham",
        action="store_true",
        default=False,
        help="Should we run linearham?")

Script.AddOption("--cluster-ind",
        dest="cluster_ind",
        type="str",
        default="0",
        help="An index specifying the clonal family of interest.")

Script.AddOption("--template-path",
        dest="template_path",
        type="str",
        default=None,
        help="The Rev template path.")

Script.AddOption("--mcmc-iter",
        dest="mcmc_iter",
        type="str",
        default="10000",
        help="How many RevBayes MCMC iterations should we use?")

Script.AddOption("--mcmc-thin",
        dest="mcmc_thin",
        type="str",
        default="10",
        help="What RevBayes MCMC thinning frequency should we use?")

Script.AddOption("--tune-iter",
        dest="tune_iter",
        type="str",
        default="5000",
        help="How many RevBayes tuning iterations should we use?")

Script.AddOption("--tune-thin",
        dest="tune_thin",
        type="str",
        default="100",
        help="What RevBayes tuning thinning frequency should we use?")

Script.AddOption("--num-rates",
        dest="num_rates",
        type="str",
        default="4",
        help="The number of gamma rate categories.")

Script.AddOption("--burnin-frac",
        dest="burnin_frac",
        type="str",
        default="0.1",
        help="What fraction of MCMC burnin should we use?")

Script.AddOption("--subsamp-frac",
        dest="subsamp_frac",
        type="str",
        default="0.05",
        help="What bootstrap sampling fraction should we use?")

Script.AddOption("--num-cores",
        dest="num_cores",
        type=int,
        default=1,
        help="The number of cores to use for ASR sampling.")

Script.AddOption("--seed",
        dest="seed",
        type="str",
        default="0",
        help="The RNG seed.")

# partis/linearham arguments

Script.AddOption("--outdir",
        dest="outdir",
        default="output",
        help="The output directory.")


#### Process command line arguments/options

def process_multiarg(arg_str, type, delim):
    return [type(arg) for arg in arg_str.split(delim)]

def get_options(env):
    return dict(
        # partis arguments
        build_partis = env.GetOption("build_partis"),
        run_partis = env.GetOption("run_partis"),
        fasta_path = env.GetOption("fasta_path"),

        # linearham arguments
        build_linearham = env.GetOption("build_linearham"),
        run_linearham = env.GetOption("run_linearham"),
        cluster_ind = process_multiarg(env.GetOption("cluster_ind"), int, ","),
        template_path = env.GetOption("template_path"),
        mcmc_iter = process_multiarg(env.GetOption("mcmc_iter"), int, ","),
        mcmc_thin = process_multiarg(env.GetOption("mcmc_thin"), int, ","),
        tune_iter = process_multiarg(env.GetOption("tune_iter"), int, ","),
        tune_thin = process_multiarg(env.GetOption("tune_thin"), int, ","),
        num_rates = process_multiarg(env.GetOption("num_rates"), int, ","),
        burnin_frac = process_multiarg(env.GetOption("burnin_frac"), float, ","),
        subsamp_frac = process_multiarg(env.GetOption("subsamp_frac"), float, ","),
        num_cores = env.GetOption("num_cores"),
        seed = process_multiarg(env.GetOption("seed"), int, ","),

        # partis/linearham arguments
        outdir = env.GetOption("outdir")
    )

env = Environment()
options = get_options(env)


#### Set up the nesting structure

nest = nestly_scons.SConsWrap(nestly.Nest(), dest_dir=options["outdir"], alias_environment=env)

def default_label(d):
    return d.get("id")


#### Install partis (if necessary)

if options["build_partis"]:

    partis_env = env.Clone()
    partis_env.Command(["lib/partis/packages/samtools/samtools",
                        "lib/partis/packages/ig-sw/src/ig_align/ig-sw",
                        "lib/partis/packages/ham/bcrham"], "",
                       "cd lib/partis && ./bin/build.sh")


#### Run partis (if necessary)

if options["run_partis"]:

    @nest.add_target()
    def partis_output(outdir, c):
        partis_output = env.Command(
            [os.path.join(outdir, filename) for filename in
                ["partis_output.yaml", "partis_output.log"]],
            options["fasta_path"],
            "lib/partis/bin/partis partition --linearham" \
                + " --infname $SOURCE" \
                + " --parameter-dir " + os.path.join(outdir, "hmm_param_dir/") \
                + " --outfname ${TARGETS[0]}" \
                + " > ${TARGETS[1]}")
        env.Depends(partis_output, c["partis_build"])
        return partis_output


#### Install linearham (if necessary)

if options["build_linearham"]:

    yamlcpp_env = env.Clone()
    yamlcpp_env.VariantDir("_build/yaml-cpp", "lib/yaml-cpp", duplicate=0)
    yamlcpp_env.Command("_build/yaml-cpp/libyaml-cpp.a", "",
                        "cp -r lib/yaml-cpp _build/ && " + \
                        "cd _build/yaml-cpp && cmake . && make")

    libptpll_env = env.Clone()
    libptpll_env.VariantDir("_build/libptpll", "lib/libptpll", duplicate=0)
    libptpll_env.Command("_build/libptpll/libptpll_static.a", "",
                         "cp -r lib/libptpll _build/ && cd _build/libptpll && " + \
                         "make && cp -t . _build/src/libptpll_static.a " + \
                         "_build/lib/_prefix/lib/libpll_algorithm.a " + \
                         "_build/lib/_prefix/lib/libpll_optimize.a " + \
                         "_build/lib/_prefix/lib/libpll_tree.a " + \
                         "_build/lib/_prefix/lib/libpll_util.a " + \
                         "_build/lib/_prefix/lib/libpll.a " + \
                         "_build/lib/lesplace/src/liblesplace-static.a")

    linearham_env = env.Clone()
    linearham_env.Append(CPPPATH=["lib/eigen", "lib/yaml-cpp/include",
                                  "lib/fast-cpp-csv-parser", "lib/libptpll/src",
                                  "lib/tclap/include",
                                  "_build/libptpll/_build/lib/_prefix/include"])
    linearham_env.Append(CCFLAGS=["-pthread", "-std=c++11", "-g"])
    linearham_env.Append(LIBPATH=["_build/yaml-cpp", "_build/libptpll"])
    linearham_env.Append(LIBS=["ptpll_static", "pll_algorithm", "pll_optimize",
                               "pll_tree", "pll_util", "pll", "lesplace-static",
                               "gsl", "blas", "yaml-cpp", "pthread"])
    linearham_env.Append(LINKFLAGS=["-g"])
    linearham_env.VariantDir("_build/linearham", "src")
    linearham_env.StaticLibrary(
        target="_build/linearham/liblinearham.a",
        source=Glob("_build/linearham/*.cpp",
                    exclude=["_build/linearham/linearham.cpp"])
    )
    linearham_env.Append(CPPPATH=["src"])
    linearham_env.Append(LIBPATH=["_build/linearham"])
    linearham_env.Prepend(LIBS=["linearham"])
    linearham_env.Program(target="_build/linearham/linearham",
                          source="_build/linearham/linearham.cpp")

    test_env = linearham_env.Clone()
    test_env.VariantDir("_build/test", "test")
    test_env.Program(target="_build/test/test", source="_build/test/test.cpp")


#### Run linearham (if necessary)

if options["run_linearham"]:

    @nest.add_target()
    def partis_yaml_file(outdir, c):
        return os.path.join(outdir, "partis_output.yaml")

    @nest.add_target()
    def hmm_param_dir(outdir, c):
        return os.path.join(outdir, "hmm_param_dir/")

    @nest.add_nest(label_func=default_label)
    def partition(c):
        return [{"id": "partition" + str(i), "index": i} for i in options["cluster_ind"]]

    @nest.add_target()
    def partition_fasta_file(outdir, c):
        return os.path.join(outdir, "input_seqs.fasta")

    @nest.add_nest(label_func=default_label)
    def revbayes_settings(c):
        return [{"id": "mcmc-iter" + str(mcmc_iter) + "_mcmc-thin" + str(mcmc_thin) + \
                       "_tune-iter" + str(tune_iter) + "_tune-thin" + str(tune_thin) + \
                       "_seed" + str(seed),
                 "mcmc_iter": mcmc_iter, "mcmc_thin": mcmc_thin,
                 "tune_iter": tune_iter, "tune_thin": tune_thin, "seed": seed}
                for mcmc_iter in options["mcmc_iter"]
                for mcmc_thin in options["mcmc_thin"]
                for tune_iter in options["tune_iter"]
                for tune_thin in options["tune_thin"]
                for seed in options["seed"]]

    @nest.add_target()
    def revbayes_rev_file(outdir, c):
        revbayes_rev_file = env.Command(
            os.path.join(outdir, "revbayes_run.rev"), options["template_path"],
            "scripts/generate_revbayes_rev_file.py $SOURCE" \
                + " --fasta-path " + str(c["partition_fasta_file"]) \
                + " --mcmc-iter " + str(c["revbayes_settings"]["mcmc_iter"]) \
                + " --mcmc-thin " + str(c["revbayes_settings"]["mcmc_thin"]) \
                + " --tune-iter " + str(c["revbayes_settings"]["tune_iter"]) \
                + " --tune-thin " + str(c["revbayes_settings"]["tune_thin"]) \
                + " --seed " + str(c["revbayes_settings"]["seed"]) \
                + " --output-path $TARGET")
        env.Depends(revbayes_rev_file, "scripts/generate_revbayes_rev_file.py")
        return revbayes_rev_file
