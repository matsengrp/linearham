env = Environment(
    CPPPATH = ['eigen', 'src', 'yaml-cpp/include'],
    CCFLAGS='-std=c++11 -g',
    LIBPATH='.',
    LIBS='yaml-cpp',
    LINKFLAGS='-g')

# Doubles compilation time.
#env.Append(CCFLAGS='-O3 -msse2')

env.Library(target='yaml-cpp', source=Glob('yaml-cpp/src/*.cpp'))

env.Program(
    target='test',
    source=[Glob('src/*.cpp')])
