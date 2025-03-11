mkdir -p bin
for f in cpp/*.cpp ; do
  BASE=$( basename ${f} .cpp )
  g++ ${f} -o bin/${BASE}_cpp
done

