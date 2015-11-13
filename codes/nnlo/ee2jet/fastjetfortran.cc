//STARTHEADER
// $Id: fastjetfortran.cc 1570 2009-05-25 10:45:18Z salam $
//
// Copyright (c) 2005-2007, Matteo Cacciari, Gavin Salam and Gregory Soyez
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
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet; if not, write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//----------------------------------------------------------------------
//ENDHEADER

#include <iostream>
#include "../fjcore.hh"
#include "JadePlugin.hh"

namespace fj = fjcore;
using namespace std;

namespace fjcore {

typedef JetDefinition::DefaultRecombiner DefRecomb;

class E0Recombiner : public  DefRecomb {
public:
  E0Recombiner(RecombinationScheme recomb_scheme = E_scheme) : 
    DefRecomb(recomb_scheme) {};

  virtual std::string description() const {
    return DefRecomb::description()+"";}

  virtual void recombine(const PseudoJet & pa, const PseudoJet & pb, 
                         PseudoJet & pab) const {
    double Eab  = pa.e() + pb.e();
    double norm = sqrt(pow(pa.px() + pb.px(),2) +
	               pow(pa.py() + pb.py(),2) +
	               pow(pa.pz() + pb.pz(),2))/Eab;
    if (norm != 0.0) {
      PseudoJet comb_jet((pa.px() + pb.px())/norm,
	                 (pa.py() + pb.py())/norm,
	                 (pa.pz() + pb.pz())/norm,
			 Eab);
      pab = comb_jet;
    }
    else {
      PseudoJet comb_jet(0.0,0.0,0.0,Eab);
      pab = comb_jet;
    };
    pab.set_user_index(max(pa.user_index(),0) + max(pb.user_index(),0));
  }
};

}


extern "C" {   

// To make the sequence globally reachable:
fj::ClusterSequence * cluster;
fj::JetDefinition jet_def_my;
fj::JetDefinition::Plugin * plugin;

// f77 interface to the pp generalised-kt (sequential recombination)
// algorithms, as defined in arXiv.org:0802.1189, which includes
// kt, Cambridge/Aachen and anti-kt as special cases.
//
// Corresponds to the following Fortran subroutine
// interface structure:
//
//   SUBROUTINE FASTJETPPSEQREC(P,NPART,R,PALG,F77JETS,NJETS)
//   DOUBLE PRECISION P(4,*), R, PALG, F, F77JETS(4,*)
//   INTEGER          NPART, NJETS
// 
// where on input
//
//   P        the input particle 4-momenta
//   NPART    the number of input momenta
//   R        the radius parameter
//   PALG     the power for the generalised kt alg 
//            (1.0=kt, 0.0=C/A,  -1.0 = anti-kt)
//
// and on output 
//
//   F77JETS  the output jet momenta (whose second dim should be >= NPART)
//            sorted in order of decreasing p_t.
//   NJETS    the number of output jets 
//
// For the values of PALG that correspond to "standard" cases (1.0=kt,
// 0.0=C/A, -1.0 = anti-kt) this routine actually calls the direct
// implementation of those algorithms, whereas for other values of
// PALG it calls the generalised kt implementation.
//
void fastjetppgenkt_(const double * p, const int & npart,                   
                     const double & R, const double & palg, const double & ptmin,
                     double * f77jets, int & njets, int * f77jetvec) {

    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<fj::PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      fj::PseudoJet psjet(mom);
      input_particles.push_back(psjet);    
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
    }
    
    // prepare jet def and run fastjet
    fj::JetDefinition jet_def;
    if (palg == 1.0) {
      jet_def = fj::JetDefinition(fj::kt_algorithm, R);
    }  else if (palg == 0.0) {
      jet_def = fj::JetDefinition(fj::cambridge_algorithm, R);
    }  else if (palg == -1.0) {
      jet_def = fj::JetDefinition(fj::antikt_algorithm, R);
    } else {
      jet_def = fj::JetDefinition(fj::genkt_algorithm, R, palg);
    }

    
    // perform clustering
    fj::ClusterSequence cs(input_particles, jet_def);
    // extract jets (pt-ordered)
    vector<fj::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
    njets = jets.size();

    // find particles inside i-th jet
    vector<fj::PseudoJet> *constit;
    constit=new vector<fj::PseudoJet>[njets];
    for(int i=0; i<njets; i++) {
      constit[i] = cs.constituents(jets[i]); 
      //cout<<"jet "<<i<<endl;
      //cout<<"mult "<<constit[i].size()<<endl;
      for(int j=0; j<constit[i].size(); j++) {
	*(f77jetvec + constit[i][j].user_index()-1) = i+1;
      }
    }


    // transfer jets -> f77jets[4*ijet+0..3]
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
        *f77jets = jets[i][j];
        f77jets++;
      } 
    }

    // clean up
    delete [] constit;
    
   }

// scheme = 0 E0 scheme E0 distance, E0 recombination
// scheme = 1 E scheme E distance, E recombination
// scheme = 2 E0tilde scheme E0 distance, E recombination
void fastjeteegenjade_(const double * p, const int & npart, const double & ycut,
                       const double & scheme,
                       double * f77jets, int & njets, int * f77jetvec) {

    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<fj::PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      fj::PseudoJet psjet(mom);
      input_particles.push_back(psjet);    
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
    }

    //fj::JadePlugin jade;
    fj::JetDefinition::Plugin * pluginE = new fj::JadeEPlugin();
    fj::JetDefinition::Plugin * pluginE0 = new fj::JadeE0Plugin();
    fj::JetDefinition::Recombiner *recomb = new fj::E0Recombiner();
    fj::JetDefinition jet_def;
    // E0 scheme:
    if (scheme == 0.0) {
      jet_def = fj::JetDefinition(pluginE0);
      jet_def.set_recombiner(recomb);
    }
    // E scheme:
    else if (scheme == 1.0) {
      jet_def = fj::JetDefinition(pluginE);
    }
    // E0 scheme w/ E recombination scheme:
    else if (scheme == 2.0) {
      jet_def = fj::JetDefinition(pluginE0);
    };
    fj::ClusterSequence cs(input_particles, jet_def);

    vector<fj::PseudoJet> jets = cs.exclusive_jets_ycut(ycut);

    njets = jets.size();

    // find particles inside i-th jet
    vector<fj::PseudoJet> *constit;
    constit=new vector<fj::PseudoJet>[njets];
    for(int i=0; i<njets; i++) {
      constit[i] = cs.constituents(jets[i]); 
      //cout<<"jet "<<i<<endl;
      //cout<<"mult "<<constit[i].size()<<endl;
      for(int j=0; j<constit[i].size(); j++) {
	*(f77jetvec + constit[i][j].user_index()-1) = i+1;
      }
    }

    // transfer jets -> f77jets[4*ijet+0..3]
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
        *f77jets = jets[i][j];
        f77jets++;
      } 
    }

    // clean up
    delete [] constit;
    delete pluginE;
    delete pluginE0;
    delete recomb;
    
  }

// Creates a cluster sequence corresponding to the JADE algo, that is scheme = 2
void fastjeteegenjadecluster_(const double * p, const int & npart) {

    // transfer p[4*ipart+0..3] -> input_particles[i]
    vector<fj::PseudoJet> input_particles;   
    for (int i=0; i<npart; i++) {
      valarray<double> mom(4); // mom[0..3]
      for (int j=0;j<=3; j++) {
         mom[j] = *(p++);
      }
      fj::PseudoJet psjet(mom);
      input_particles.push_back(psjet);    
      // label input_particles entries
      input_particles[i].set_user_index(i+1);
    }

    plugin = new fj::JadeE0Plugin();

    jet_def_my = fj::JetDefinition(plugin);

    cluster = new fj::ClusterSequence(input_particles, jet_def_my);

  }

void fastjetgetjets_(const double & ycut, 
                     double * f77jets, int & njets, int * f77jetvec) {

    vector<fj::PseudoJet> jets = cluster->exclusive_jets_ycut(ycut);

    njets = jets.size();

    // find particles inside i-th jet
    vector<fj::PseudoJet> *constit;
    constit=new vector<fj::PseudoJet>[njets];
    for(int i=0; i<njets; i++) {
      constit[i] = cluster->constituents(jets[i]); 
      //cout<<"jet "<<i<<endl;
      //cout<<"mult "<<constit[i].size()<<endl;
      for(int j=0; j<constit[i].size(); j++) {
	*(f77jetvec + constit[i][j].user_index()-1) = i+1;
      }
    }

    // transfer jets -> f77jets[4*ijet+0..3]
    for (int i=0; i<njets; i++) {
      for (int j=0;j<=3; j++) {
        *f77jets = jets[i][j];
        f77jets++;
      } 
    }

    delete [] constit;

  }

void cleancluster_() {
    delete cluster;
    delete plugin;
  }
}

