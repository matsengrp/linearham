env = Environment(
    CPPPATH = ['eigen', 'src'],
    CCFLAGS='-std=c++11 -g',
    LINKFLAGS='-g')

# Doubles compilation time.
#env.Append(CCFLAGS='-O3 -msse2')

env.Program(
    target='test',
    source=[Glob('src/*.cpp')])

