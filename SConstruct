env = Environment(
    CPPPATH = ['eigen', 'src'])

# Doubles compilation time.
#env.Append(CCFLAGS='-O3 -msse2')

env.Program(
    target='test',
    source=Split('src/linalg.cpp src/test.cpp'))

