/***************************************************************************
                          exportsystem.h  -  description
                             -------------------
    begin                : Wed Feb 12 2003
    copyright            : (C) 2003 by Arnaud Le Rouzic
    email                : lerouzic@cnrs-gif.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 
#ifndef EXPORTSYSTEM_H
#define EXPORTSYSTEM_H 

#include <vector>
#include <fstream>
#include <string>

#include "element.h"


// forward declarations, since exportsystem.h is included in element.h
class System; 
class Species;
class Population;
class Individual;
class Gamete;
class Element;
template <typename T> struct ElementFeature;

class ExportSystem
{
  public:
  ExportSystem(void);
  ExportSystem(const std::string);
  ExportSystem(const System *);
  ExportSystem(const std::string, const System *);

  void Save(void);
  // Impression de l'en-tête du fichier xml
  static std::string DTD(void);
  
  private:
  void Print_header(void);
  void Print_genome(void);
  void Print_memory(void);
  void Print_systeme(const System *);
  void Print_espece(const Species *);
  void Print_population(const Population *);
  void Print_individu(const Individual *, const std::string &);
  void Print_gamete(const Gamete *);
  void Print_ref_element(const Ptr2Element);
  void Print_element(const Ptr2Element);
  template <typename T> void Print_ElementFeature(const ElementFeature<T> &, const std::string &);
  void Print_ccl(void);

  bool Ouv_fic(const std::string);

  void UpdateMemory(void);

  void _e_ouvfic(const std::string &) const;
  void _e_systnull(void) const;

  const System * sys;
  std::vector<Ptr2Element> memory;
  std::ofstream ficsx;
};

#endif // EXPORTSYSTEM_H
