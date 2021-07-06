cc = g++
targets = FiltersSpatial FiltersTransform ProcessImage

all: stacklim $(targets)

stacklim:
	ulimit -s unlimited

clean:
	rm $(targets)