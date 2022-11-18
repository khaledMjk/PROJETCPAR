# Configuration pour MSYS 2
# Peut être modifié à votre convenance...
CXX = mpic++
CXXFLAGS = -std=c++17 -fopenmp
ifdef DEBUG
CXXFLAGS += -g -O0 -Wall -fbounds-check -pedantic -D_GLIBCXX_DEBUG -fopenmp
CXXFLAGS2 = CXXFLAGS
else
CXXFLAGS2 = ${CXXFLAGS} -O2 -march=native -Wall -fopenmp
CXXFLAGS += -O3 -march=native -Wall
endif
# LIB=-lsfml-graphics -lsfml-window -lsfml-main -lsfml-system
LIB=-lsfml-graphics -lsfml-window -lsfml-system
