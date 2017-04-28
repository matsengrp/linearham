env = Environment()

## dependencies ##

yaml_cpp_env = env.Clone()
yaml_cpp_env.VariantDir('_build/yaml-cpp', 'lib/yaml-cpp', duplicate=0)
yaml_cpp_env.Command('_build/yaml-cpp/libyaml-cpp.a', '', 'cd _build/yaml-cpp && cmake ../../lib/yaml-cpp && make')

pt_env = env.Clone()
pt_env.Command('_build/libptpll/libptpll_static.a', '', 'cp -r lib/libptpll _build && cd _build/libptpll && make && mv _build/src/libptpll_static.a .')

brent_env = env.Clone()
brent_env.VariantDir('_build/brent', 'lib/brent', duplicate=0)
brent_env.StaticLibrary(target='_build/brent/brent',
                        source='_build/brent/brent.cpp')

## linearham ##

common_env = env.Clone()
common_env.Append(CPPPATH=['lib/eigen', 'lib/yaml-cpp/include', 'lib/fast-cpp-csv-parser', 'lib/libptpll/src', 'lib/brent'])
common_env.Append(CCFLAGS=['-pthread', '-std=c++11', '-g'])
common_env.Append(LIBPATH=['_build/yaml-cpp', '_build/libptpll', '_build/brent'])
common_env.Append(LIBS=['ptpll_static', File('/usr/local/lib/libpll.a'), 'yaml-cpp', 'pthread', 'brent'])
common_env.Append(LINKFLAGS=['-g'])

# Doubles compilation time.
#common_env.Append(CCFLAGS=['-O3', '-msse2'])

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
                 source=Glob('_build/test/*.cpp'))
