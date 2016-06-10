env = Environment()

## dependencies ##

yaml_cpp_env = env.Clone()
yaml_cpp_env.VariantDir('_build/yaml-cpp', 'lib/yaml-cpp', duplicate=0)
yaml_cpp_env.Command('_build/yaml-cpp/libyaml-cpp.a', '', 'cd _build/yaml-cpp && cmake ../../lib/yaml-cpp && make')

## linearham ##

common_env = env.Clone()
common_env.Append(CPPPATH=['lib/eigen', 'lib/yaml-cpp/include'])
common_env.Append(CCFLAGS=['-std=c++11', '-g'])
common_env.Append(LIBPATH=['_build/yaml-cpp'])
common_env.Append(LIBS=['yaml-cpp'])
common_env.Append(LINKFLAGS=['-g'])

# Doubles compilation time.
#common_env.Append(CCFLAGS=['-O3', '-msse2'])

linearham_env = common_env.Clone()
linearham_env.VariantDir('_build/linearham', 'src')
linearham_env.Library(target='_build/linearham/linearham',
                      source=Glob('_build/linearham/*.cpp'))

test_env = common_env.Clone()
test_env.VariantDir('_build/test', 'test')
test_env.Append(CPPPATH=['src'])
test_env.Program(target='_build/test/test',
                 LIBPATH=['_build/linearham', '_build/yaml-cpp'],
                 LIBS=['linearham', 'yaml-cpp'],
                 source=Glob('_build/test/*.cpp'))
