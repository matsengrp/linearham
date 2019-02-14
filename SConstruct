#!/usr/bin/env python
# -*- coding: utf-8 -*-
import nestly
from nestly import scons as nestly_scons
import os
from SCons.Script import Environment
import SCons.Script as Script


#### Set up command line arguments/options

# partis arguments

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

Script.AddOption("--all-clonal-seqs",
        dest="all_clonal_seqs",
        action="store_true",
        default=False,
        help="Should we assume all the sequences in the FASTA file are clonal?")

# linearham arguments

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

Script.AddOption("--seed-seq",
        dest="seed_seq",
        type="str",
        default=None,
        help="The name of the seed sequence.")

Script.AddOption("--asr-pfilters",
        dest="asr_pfilters",
        type="str",
        default="0.1",
        help="The ancestral sequence posterior probability threshold.")

# partis/linearham arguments

Script.AddOption("--build-partis-linearham",
        dest="build_partis_linearham",
        action="store_true",
        default=False,
        help="Should we build partis and linearham?")

Script.AddOption("--hmm-param-dir",
        dest="hmm_param_dir",
        type="str",
        default=None,
        help="An optional directory of partis HMM germline parameter files.")

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
        run_partis = env.GetOption("run_partis"),
        fasta_path = env.GetOption("fasta_path"),
        all_clonal_seqs = env.GetOption("all_clonal_seqs"),

        # linearham arguments
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
        seed_seq = process_multiarg(env.GetOption("seed_seq"), str, ",") if env.GetOption("seed_seq") is not None else None,
        asr_pfilters = process_multiarg(env.GetOption("asr_pfilters"), float, ","),

        # partis/linearham arguments
        build_partis_linearham = env.GetOption("build_partis_linearham"),
        hmm_param_dir = env.GetOption("hmm_param_dir"),
        outdir = env.GetOption("outdir")
    )

env = Environment(ENV = os.environ)
options = get_options(env)
env.SConsignFile(os.path.join(options["outdir"], ".sconsign"))


#### Set up the nesting structure

nest = nestly_scons.SConsWrap(nestly.Nest(), dest_dir=options["outdir"], alias_environment=env)

def default_label(d):
    return d.get("id")


#### Install partis and linearham (if necessary)

if options["build_partis_linearham"]:

    @nest.add_target()
    def partis_build(outdir, c):
        return env.Command("lib/partis/packages/ham/bcrham", "",
                           "cd lib/partis && ./bin/build.sh")

    @nest.add_target()
    def linearham_build(outdir, c):
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
        linearham_env.Append(CPPPATH=["lib/eigen", "lib/fast-cpp-csv-parser",
                                      "lib/libptpll/src", "lib/tclap/include",
                                      "_build/libptpll/_build/lib/_prefix/include"])
        linearham_env.Append(CCFLAGS=["-pthread", "-std=c++11", "-g"])
        linearham_env.Append(LIBPATH=["_build/libptpll"])
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
        linearham_bin = linearham_env.Program(target="_build/linearham/linearham",
                                              source="_build/linearham/linearham.cpp")

        test_env = linearham_env.Clone()
        test_env.VariantDir("_build/test", "test")
        test_env.Program(target="_build/test/test", source="_build/test/test.cpp")
        return linearham_bin


#### Run partis (if necessary)

if options["run_partis"]:

    @nest.add_target()
    def partis_output(outdir, c):
        partis_mode = "annotate" if options["all_clonal_seqs"] else "partition"
        partis_hmm_param_dir = options["hmm_param_dir"].rstrip("/") if options["hmm_param_dir"] is not None else os.path.join(outdir, "hmm_param_dir")
        partis_output = env.Command(
            [os.path.join(outdir, filename) for filename in
                ["partis_run.yaml", "partis_run.stdout.log"]],
            options["fasta_path"],
            "lib/partis/bin/partis " + partis_mode + " --linearham" \
                + (" --all-seqs-simultaneous" if options["all_clonal_seqs"] else "") \
                +  " --infname $SOURCE" \
                +  " --parameter-dir " + partis_hmm_param_dir \
                +  " --outfname ${TARGETS[0]}" \
                +  " > ${TARGETS[1]}")
        env.Depends(partis_output, "lib/partis/packages/ham/bcrham")
        return partis_output


#### Run linearham (if necessary)

if options["run_linearham"]:

    @nest.add_target()
    def partis_yaml_file(outdir, c):
        return os.path.join(outdir, "partis_run.yaml")

    @nest.add_target()
    def hmm_param_dir(outdir, c):
        linearham_hmm_param_dir = options["hmm_param_dir"].rstrip("/") if options["hmm_param_dir"] is not None else os.path.join(outdir, "hmm_param_dir")
        return linearham_hmm_param_dir + "/hmm/hmms"

    @nest.add_nest(label_func=default_label)
    def partition(c):
        return [{"id": "partition" + str(i), "index": i} for i in options["cluster_ind"]]

    @nest.add_target()
    def partition_fasta_file(outdir, c):
        return os.path.join(outdir, "input_seqs.fasta")

    @nest.add_nest(label_func=default_label)
    def revbayes_setting(c):
        return [{"id": "mcmciter" + str(mcmc_iter) + "_mcmcthin" + str(mcmc_thin) + \
                       "_tuneiter" + str(tune_iter) + "_tunethin" + str(tune_thin) + \
                       "_numrates" + str(num_rates) + "_seed" + str(seed),
                 "mcmc_iter": mcmc_iter, "mcmc_thin": mcmc_thin,
                 "tune_iter": tune_iter, "tune_thin": tune_thin,
                 "num_rates": num_rates, "seed": seed}
                for mcmc_iter in options["mcmc_iter"]
                for mcmc_thin in options["mcmc_thin"]
                for tune_iter in options["tune_iter"]
                for tune_thin in options["tune_thin"]
                for num_rates in options["num_rates"]
                for seed in options["seed"]]

    @nest.add_target()
    def revbayes_rev_file(outdir, c):
        revbayes_rev_file = env.Command(
            os.path.join(outdir, "revbayes_run.rev"), options["template_path"],
            "scripts/generate_revbayes_rev_file.py $SOURCE" \
                + " --fasta-path " + c["partition_fasta_file"] \
                + " --mcmc-iter " + str(c["revbayes_setting"]["mcmc_iter"]) \
                + " --mcmc-thin " + str(c["revbayes_setting"]["mcmc_thin"]) \
                + " --tune-iter " + str(c["revbayes_setting"]["tune_iter"]) \
                + " --tune-thin " + str(c["revbayes_setting"]["tune_thin"]) \
                + " --num-rates " + str(c["revbayes_setting"]["num_rates"]) \
                + " --seed " + str(c["revbayes_setting"]["seed"]) \
                + " --output-path $TARGET")
        env.Depends(revbayes_rev_file, "scripts/generate_revbayes_rev_file.py")
        return revbayes_rev_file

    @nest.add_target()
    def revbayes_output(outdir, c):
        revbayes_output = env.Command(
            [os.path.join(outdir, filename) for filename in
                ["revbayes_run.trees", "revbayes_run.log", "revbayes_run.stdout.log"]],
            c["revbayes_rev_file"],
            "lib/revbayes/projects/cmake/rb $SOURCE > ${TARGETS[2]}")
        env.Depends(revbayes_output, "lib/revbayes/projects/cmake/rb")
        return revbayes_output

    @nest.add_target()
    def linearham_intermediate_output(outdir, c):
        linearham_intermediate_output = env.Command(
            os.path.join(outdir, "lh_revbayes_run.trees"), c["revbayes_output"][0],
            "_build/linearham/linearham --pipeline" \
                + " --yaml-path " + c["partis_yaml_file"] \
                + " --cluster-ind " + str(c["partition"]["index"]) \
                + " --hmm-param-dir " + c["hmm_param_dir"] \
                + " --seed " + str(c["revbayes_setting"]["seed"]) \
                + " --num-rates " + str(c["revbayes_setting"]["num_rates"]) \
                + " --input-path $SOURCE" \
                + " --output-path $TARGET")
        env.Depends(linearham_intermediate_output, "_build/linearham/linearham")
        return linearham_intermediate_output

    @nest.add_nest(label_func=default_label)
    def linearham_setting(c):
        return [{"id": "burninfrac" + str(burnin_frac) + "_subsampfrac" + str(subsamp_frac),
                 "burnin_frac": burnin_frac, "subsamp_frac": subsamp_frac}
                for burnin_frac in options["burnin_frac"]
                for subsamp_frac in options["subsamp_frac"]]

    @nest.add_target()
    def linearham_final_output(outdir, c):
        linearham_final_output = env.Command(
            os.path.join(outdir, "linearham_run.trees"),
            c["linearham_intermediate_output"],
            "Rscript --slave --vanilla scripts/run_bootstrap_asr.R" \
                + " $SOURCE" \
                + " " + c["partition_fasta_file"] \
                + " " + str(c["linearham_setting"]["burnin_frac"]) \
                + " " + str(c["linearham_setting"]["subsamp_frac"]) \
                + " " + str(options["num_cores"]) \
                + " " + str(c["revbayes_setting"]["seed"]) \
                + " $TARGET")
        env.Depends(linearham_final_output, "scripts/run_bootstrap_asr.R")
        return linearham_final_output

    if options["seed_seq"] is not None:

        @nest.add_nest(label_func=default_label)
        def seed_seq(c):
            return [{"id": "seedseq" + seed_seq, "name": seed_seq}
                    for seed_seq in options["seed_seq"]]

        @nest.add_target()
        def linearham_visualization(outdir, c):
            outbase = os.path.join(outdir, c["seed_seq"]["name"])
            linearham_visualization = env.Command(
                [outbase + ".pfilter" + str(pfilter) + ".aa_lineage_graph.png" for pfilter in options["asr_pfilters"]],
                c["linearham_final_output"],
                "scripts/trees_to_counted_ancestors.py $SOURCE" \
                    + " --seed-seq " + c["seed_seq"]["name"] \
                    + " --pfilters " + " ".join(str(pfilter) for pfilter in options["asr_pfilters"]) \
                    + " --output-base " + outbase)
            env.Depends(linearham_visualization, "scripts/trees_to_counted_ancestors.py")
            return linearham_visualization
