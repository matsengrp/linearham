env = Environment()

## dependencies ##

yamlcpp_env = env.Clone()
yamlcpp_env.VariantDir('_build/yaml-cpp', 'lib/yaml-cpp', duplicate=0)
yamlcpp_env.Command('_build/yaml-cpp/libyaml-cpp.a', '',
                    'cp -r lib/yaml-cpp _build/ && ' + \
                    'cd _build/yaml-cpp && cmake . && make')

libptpll_env = env.Clone()
libptpll_env.VariantDir("_build/libptpll", "lib/libptpll", duplicate=0)
libptpll_env.Command('_build/libptpll/libptpll_static.a', '',
                     'cp -r lib/libptpll _build/ && cd _build/libptpll && ' + \
                     'make && cp -t . _build/src/libptpll_static.a ' + \
                     '_build/lib/_prefix/lib/libpll_algorithm.a ' + \
                     '_build/lib/_prefix/lib/libpll_optimize.a ' + \
                     '_build/lib/_prefix/lib/libpll_tree.a ' + \
                     '_build/lib/_prefix/lib/libpll_util.a ' + \
                     '_build/lib/_prefix/lib/libpll.a ' + \
                     '_build/lib/lesplace/src/liblesplace-static.a')

## linearham ##

linearham_env = env.Clone()
linearham_env.Append(CPPPATH=['lib/eigen', 'lib/yaml-cpp/include',
                              'lib/fast-cpp-csv-parser', 'lib/libptpll/src',
                              'lib/tclap/include',
                              '_build/libptpll/_build/lib/_prefix/include'])
linearham_env.Append(CCFLAGS=['-pthread', '-std=c++11', '-g'])
linearham_env.Append(LIBPATH=['_build/yaml-cpp', '_build/libptpll'])
linearham_env.Append(LIBS=['ptpll_static', 'pll_algorithm', 'pll_optimize',
                           'pll_tree', 'pll_util', 'pll', 'lesplace-static',
                           'gsl', 'blas', 'yaml-cpp', 'pthread'])
linearham_env.Append(LINKFLAGS=['-g'])
linearham_env.VariantDir('_build/linearham', 'src')
linearham_env.StaticLibrary(
    target='_build/linearham/liblinearham.a',
    source=Glob('_build/linearham/*.cpp',
                exclude=['_build/linearham/linearham.cpp'])
)
linearham_env.Append(CPPPATH=['src'])
linearham_env.Append(LIBPATH=['_build/linearham'])
linearham_env.Prepend(LIBS=['linearham'])
linearham_env.Program(target='_build/linearham/linearham',
                      source='_build/linearham/linearham.cpp')

test_env = linearham_env.Clone()
test_env.VariantDir('_build/test', 'test')
test_env.Program(target='_build/test/test', source='_build/test/test.cpp')
