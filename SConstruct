env = Environment()

## dependencies ##

yaml_cpp_env = env.Clone()
yaml_cpp_env.VariantDir('_build/yaml-cpp', 'lib/yaml-cpp', duplicate=0)
yaml_cpp_env.Command('_build/yaml-cpp/libyaml-cpp.a', '', 'cd _build/yaml-cpp && cmake ../../lib/yaml-cpp && make')

libpll_env = env.Clone()
libpll_env.Command('_build/libpll/libpll.a', '', 'cp -r lib/libpll _build && cd _build/libpll && ./autogen.sh && ./configure && make && mv src/.libs/libpll.a .')

## linearham ##

common_env = env.Clone()
common_env.Append(CPPPATH=['lib/eigen', 'lib/yaml-cpp/include', 'lib/fast-cpp-csv-parser', 'lib/libpll/src'])
common_env.Append(CCFLAGS=['-pthread', '-std=c++11', '-g'])
common_env.Append(LIBPATH=['_build/yaml-cpp', '_build/libpll'])
common_env.Append(LIBS=['pll', 'yaml-cpp', 'pthread'])
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
                 source=Glob('_build/test/*.cpp'))
