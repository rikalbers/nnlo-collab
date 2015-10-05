//FJSTARTHEADER
// $Id: JadePlugin.cc 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2007-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER

// fastjet stuff
#include "../fjcore.hh"
#include "JadePlugin.hh"
#include <iostream>
#include <sstream>  // needed for internal io
#include <iomanip>

// other stuff
#include <vector>
#include <sstream>
#include <limits>




namespace fj = fjcore;
using namespace std;

//----------------------------------------------------------------------
// forward declaration for printing out info about a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream &, const fj::PseudoJet &);

namespace fjcore {



//----------------------------------------------------------------------
/// class to help run a JADE algorithm distance measured as in the E-scheme:
// y = (pa + pb)^2/S
class JadeEBriefJet {
public:
  void init(const PseudoJet & jet) {
    double norm;
    if (jet.modp2() !=  0.0) {
      norm = 1.0/sqrt(jet.modp2());
      nx = jet.px() * norm;
      ny = jet.py() * norm;
      nz = jet.pz() * norm;
    } else {
      nx = 0.0;
      ny = 0.0;
      nz = 0.0;
    };
    pabs = jet.modp();
    E    = jet.E();
    m2   = jet.m2();
  }

  double distance(const JadeEBriefJet * jet) const {
    double dij = m2 + jet->m2 + 2.0*E*jet->E - 2.0*pabs*jet->pabs*(nx*jet->nx + ny*jet->ny + nz*jet->nz);
    return dij;
  }

  double beam_distance() const {
    return numeric_limits<double>::max();
  }

private:
  double m2, E, pabs, nx, ny, nz;
};

//----------------------------------------------------------------------
/// class to help run a JADE algorithm distance is measured as in E0:
// E0 : y = (Ea Eb (1 - cos theta(a,b)))/S
class JadeE0BriefJet {
public:
  void init(const PseudoJet & jet) {
    double norm;
    if (jet.modp2() !=  0.0) {
      norm = 1.0/sqrt(jet.modp2());
      nx = jet.px() * norm;
      ny = jet.py() * norm;
      nz = jet.pz() * norm;
    } else {
      nx = 0.0;
      ny = 0.0;
      nz = 0.0;
    };
    rt2E = sqrt(2.0)*jet.E();
  }

  double distance(const JadeE0BriefJet * jet) const {
    double dij = 1 - nx*jet->nx
                   - ny*jet->ny
                   - nz*jet->nz;
    if (nx == 0.0 && ny == 0.0 && nz == 0.0) dij = 0.0;
    dij *= rt2E*jet->rt2E;
    return dij;
  }

  double beam_distance() const {
    return numeric_limits<double>::max();
  }

private:
  double rt2E, nx, ny, nz;
};


//----------------------------------------------------------------------
string JadeE0Plugin::description () const {
  ostringstream desc;
  desc << "e+e- JADE algorithm plugin E0-scheme";
  return desc.str();
}

//----------------------------------------------------------------------
void JadeE0Plugin::run_clustering(ClusterSequence & cs) const {
  int njets = cs.jets().size();
  NNH<JadeE0BriefJet> nnh(cs.jets());

  // if testing against Hoeth's implementation, need to rescale the
  // dij by Q^2.
  //double Q2 = cs.Q2(); 

  while (njets > 0) {
    int i, j, k;
    double dij = nnh.dij_min(i, j);

    if (j >= 0) {
      cs.plugin_record_ij_recombination(i, j, dij, k);
      nnh.merge_jets(i, j, cs.jets()[k], k);
    } else {
      double diB = cs.jets()[i].E()*cs.jets()[i].E(); // get new diB
      cs.plugin_record_iB_recombination(i, diB);
      nnh.remove_jet(i);
    }
    njets--;
  }
}

//----------------------------------------------------------------------
string JadeEPlugin::description () const {
  ostringstream desc;
  desc << "e+e- JADE algorithm plugin E-scheme";
  return desc.str();
}

//----------------------------------------------------------------------
void JadeEPlugin::run_clustering(ClusterSequence & cs) const {
  int njets = cs.jets().size();
  NNH<JadeEBriefJet> nnh(cs.jets());

  // if testing against Hoeth's implementation, need to rescale the
  // dij by Q^2.
  //double Q2 = cs.Q2(); 

  while (njets > 0) {
    int i, j, k;
    double dij = nnh.dij_min(i, j);

    if (j >= 0) {
      cs.plugin_record_ij_recombination(i, j, dij, k);
      nnh.merge_jets(i, j, cs.jets()[k], k);
    } else {
      double diB = cs.jets()[i].E()*cs.jets()[i].E(); // get new diB
      cs.plugin_record_iB_recombination(i, diB);
      nnh.remove_jet(i);
    }
    njets--;
  }
}

}

//----------------------------------------------------------------------
// does the actual work for printing out a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream & ostr, const fj::PseudoJet & jet) {
  ostr << "pt, y, phi ="
       << " " << setw(10) << jet.perp()
       << " " << setw(6) <<  jet.rap()
       << " " << setw(6) <<  jet.phi()
       << ", mass = " << setw(10) << jet.m()
       << ", btag = " << jet.user_index();
  return ostr;
}
