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

common_env = env.Clone()
common_env.Append(CPPPATH=['lib/eigen', 'lib/yaml-cpp/include',
                           'lib/fast-cpp-csv-parser', 'lib/libptpll/src',
                           '_build/libptpll/_build/lib/_prefix/include'])
common_env.Append(CCFLAGS=['-pthread', '-std=c++11', '-g'])
common_env.Append(LIBPATH=['_build/yaml-cpp', '_build/libptpll'])
common_env.Append(LIBS=['ptpll_static', 'pll_algorithm', 'pll_optimize',
                        'pll_tree', 'pll_util', 'pll', 'lesplace-static', 'gsl',
                        'blas', 'yaml-cpp', 'pthread'])
common_env.Append(LINKFLAGS=['-g'])

# Doubles compilation time.
#common_env.Append(CCFLAGS=['-O3', '-msse2'])

###### FIX BELOW!!!

linearham_env = common_env.Clone()
linearham_env.VariantDir('_build/linearham', 'src')
linearham_env.StaticLibrary(target='_build/linearham/linearham',
                            source=Glob('_build/linearham/*.cpp'))

test_env = common_env.Clone()
test_env.VariantDir('_build/test', 'test')
test_env.Append(CPPPATH=['src'])
test_env.Append(LIBPATH=['_build/linearham'])
test_env.Prepend(LIBS=['linearham'])
test_env.Program(target='_build/test/test',
                 source=Glob('_build/test/test.cpp'))
