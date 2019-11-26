cmake_minimum_required(VERSION 3.0)

project(dyet)

# Crée des variables avec les fichiers à compiler
set(SRCS
    base.cpp
    genome.cpp
    parameters.cpp
    system.cpp
    element.cpp
    importsystem.cpp
    population.cpp
    tools.cpp
    exportsystem.cpp
    indivgroup.cpp
    recombination.cpp
    gamete.cpp
    individual.cpp
    simul.cpp
    generator.cpp
    main.cpp
    species.cpp
    )

set(HEADERS
    argstream.h
    base.h
    config.h
    dyet.h
    element.h
    exception.h
    exportsystem.h
    gamete.h
    generator.h
    genome.h
    importsystem.h
    indivgroup.h
    individual.h
    parameters.h
    population.h
    recombination.h
    simul.h
    species.h
    system.h
    tools.h
    gamete.inl
    individual.inl
    individu.inl
    population.inl
    recombination.inl
    species.inl
    system.inl
    )

add_executable(dyet_core ${SRCS} ${HEADERS})

find_package(PkgConfig REQUIRED)
find_package(LibXml2 REQUIRED)
#~ find_package(LibXml++ REQUIRED)
pkg_check_modules(LibXml++ REQUIRED libxml++-2.6)

target_include_directories(dyet_core PRIVATE ${LibXml2_INCLUDE_DIRS})
target_include_directories(dyet_core PRIVATE ${LibXml++_INCLUDE_DIRS})

target_link_libraries(dyet_core ${LibXml2_LIBRARIES} ${LibXml++_LIBRARIES})

